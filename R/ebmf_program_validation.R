#' @title Validate EBMF Programs Against Filtering and Evidence Metrics
#' @description One-stop summary of EBMF programs for filtering and reporting.
#' Folds the graph-reliability metrics (`mean_internal_similarity`,
#' `connectedness`) together with trait-subsampling **stability**, and the
#' descriptive **additional scores** (factor strength, aggregate posterior
#' evidence and factor distinctiveness) into a single per-program table with
#' pass/fail flags and an overall `status`.
#'
#' Filters applied (a program must pass all three to get `status == "valid"`):
#' \itemize{
#'   \item **Size / internal coherence** — at least `min_module_size`
#'   EBMF-supported SNPs (`n_snps`), mean internal SNP similarity
#'   `>= min_mean_internal`, and pair connectedness `>= min_connectedness` on
#'   the SNP similarity graph. The mean (and median) pairwise similarity and
#'   the connectedness score come from the existing SNP-SNP similarity
#'   framework; `min_module_size` ensures the metrics are computed on enough
#'   filtered SNPs to be meaningful. `n_snps_filtered` reports the subset of
#'   EBMF-supported SNPs that also passes the lFSR/magnitude filter.
#'   \item **Trait-subsampling stability** — the top-loading SNPs must be
#'   recovered with replication `>= stability_threshold` when `1 - frac_traits`
#'   of the traits are held out and flashier is refit. Because this refits
#'   flashier `n_rep` times, it is only run for programs that pass the
#'   internal-coherence checks; programs failing coherence are already rejected
#'   and are not re-checked.
#' }
#'
#' Additional scores (reported, but not used to gate `status`):
#' \itemize{
#'   \item `factor_strength` / `factor_strength_per_snp` — the proportion of
#'   the observed matrix signal attributable to the factor, computed from its
#'   fitted reconstruction on the lFSR/magnitude-filtered (assigned) SNPs;
#'   `factor_strength_per_snp = factor_strength / sqrt(n_snps_filtered)`
#'   normalises it by the square root of the number of filtered SNPs.
#'   \item `posterior_evidence` — aggregate (`median` or `mean`) of
#'   `-log10(lfsr)` across a program's assigned SNPs. Distinguishes programs
#'   dominated by borderline `lfsr ~ 0.05` memberships from programs whose
#'   memberships carry much stronger posterior sign evidence.
#'   \item `max_abs_factor_corr` / `redundant` — the maximum absolute
#'   correlation of the program's SNP loading vector (`F_pm`) with any other
#'   program; `redundant` flags programs whose maximum correlation exceeds
#'   `redundancy_threshold`.
#' }
#' @param clustering_result Result of `run_univariate_clustering()` with
#'   `cluster_type = "ebmf"`.
#' @param s_matrix SNP-by-SNP similarity matrix used for the coherence metrics.
#'   Defaults to `clustering_result$s_matrix`.
#' @param edge_threshold Absolute similarity used for connectedness. Defaults to
#'   `0.2`.
#' @param min_module_size Minimum assigned SNPs for the coherence metrics to be
#'   meaningful and for `size_pass`. Defaults to `5`.
#' @param min_mean_internal Minimum mean internal SNP similarity for
#'   `internal_pass`. Defaults to `0.3`.
#' @param min_connectedness Minimum pair-connectedness for `internal_pass`.
#'   Defaults to `0.5`.
#' @param n_null Number of permutation nulls for membership calibration.
#' @param alpha_membership Membership FDR level passed to
#'   `calibrate_ebmf_programs()` for `core` cell labelling.
#' @param n_candidate_tier Per program, how many non-core cells to label as the
#'   uncertified "candidate" tier.
#' @param n_rep Number of trait-subsample replicates for stability. Set to `0`
#'   to disable the stability check entirely.
#' @param frac_traits Fraction of trait rows sampled per stability replicate.
#'   Defaults to `0.8` (hold out 20% of traits).
#' @param top_n Size of each program's top-loading SNP set used for stability
#'   matching.
#' @param stability_threshold Minimum replication for `stability_pass`.
#' @param posterior_stat Aggregate statistic for the posterior evidence score:
#'   `"median"` or `"mean"` of `-log10(lfsr)`.
#' @param posterior_evidence_cap Upper bound for `-log10(lfsr)` when computing
#'   the posterior evidence score. flashier can return lFSR values that underflow
#'   towards 0, so `-log10(lfsr)` can reach hundreds or `Inf`; flooring lFSR at
#'   `10^-posterior_evidence_cap` keeps the score finite and readable while
#'   preserving the ordering. Defaults to `50`.
#' @param redundancy_threshold Maximum absolute factor-loading correlation with
#'   another program before `redundant` is flagged.
#' @param seed RNG seed.
#' @param cores Number of cores for the parallel trait-subsample refits (passed
#'   to `stability_ebmf_programs()`). Defaults to `1` (serial).
#' @param verbose Print progress messages.
#' @return A list with:
#'   \itemize{
#'     \item programs: the folded per-program table (filter metrics, additional
#'       scores, pass/fail flags, and `status`)
#'     \item memberships: calibrated SNP x program posterior from
#'       `calibrate_ebmf_programs()`
#'     \item assigned: gated assigned memberships per program
#'     \item factor_correlation: K x K correlation matrix of SNP loading vectors
#'     \item null_summary: permutation-null summary
#'     \item settings: settings used
#'   }
#' @export
summarise_ebmf_programs <- function(clustering_result,
                                    s_matrix = NULL,
                                    edge_threshold = 0.2,
                                    min_module_size = 5L,
                                    min_mean_internal = 0.3,
                                    min_connectedness = 0.5,
                                    n_null = 10,
                                    alpha_membership = 0.05,
                                    n_candidate_tier = 25L,
                                    n_rep = 10,
                                    frac_traits = 0.8,
                                    top_n = 25L,
                                    stability_threshold = 0.7,
                                    posterior_stat = c("median", "mean"),
                                    posterior_evidence_cap = 50,
                                    redundancy_threshold = 0.9,
                                    seed = 1,
                                    cores = 1,
                                    verbose = TRUE) {
  params <- clustering_result$parameters
  if (is.null(params) || !identical(params$cluster_type, "ebmf")) {
    stop("clustering_result must come from cluster_type = 'ebmf'")
  }
  posterior_stat <- match.arg(posterior_stat)

  if (is.null(s_matrix)) {
    s_matrix <- clustering_result$s_matrix
  }
  if (is.null(s_matrix)) {
    stop(
      "s_matrix is required (supply it or ensure clustering_result$s_matrix ",
      "is set)"
    )
  }

  lfsr_threshold <- params$ebmf_lfsr_threshold
  mag_threshold <- params$ebmf_magnitude_threshold

  cal <- calibrate_ebmf_programs(
    clustering_result,
    n_null = n_null,
    alpha_membership = alpha_membership,
    n_candidate_tier = n_candidate_tier,
    seed = seed,
    verbose = verbose
  )

  posterior <- cal$memberships
  ebmf_memberships <- posterior |>
    dplyr::filter(is.finite(abs_loading), abs_loading > 0)
  assigned <- ebmf_memberships |>
    dplyr::filter(
      !is.na(lfsr),
      lfsr < lfsr_threshold,
      is.na(mag_threshold) | abs_loading > mag_threshold
    )

  programs <- cal$programs
  programs$program <- as.integer(programs$program)

  n_snps_by_program <- ebmf_memberships |>
    dplyr::group_by(program) |>
    dplyr::summarise(
      n_snps = dplyr::n_distinct(snp_id),
      .groups = "drop"
    )
  n_snps_filtered_by_program <- assigned |>
    dplyr::group_by(program) |>
    dplyr::summarise(
      n_snps_filtered = dplyr::n_distinct(snp_id),
      .groups = "drop"
    )

  coherence_rows <- lapply(programs$program, function(pg) {
    snps <- assigned$snp_id[assigned$program == pg]
    return(.program_internal_coherence(s_matrix, snps, edge_threshold))
  })
  coherence <- dplyr::bind_rows(coherence_rows)
  coherence$program <- programs$program
  coherence_columns <- c(
    "n_snps", "mean_internal_similarity", "median_internal_similarity",
    "mean_external_similarity", "separation", "connectedness",
    "n_internal_pairs"
  )
  for (column in coherence_columns) {
    if (!column %in% names(coherence)) {
      coherence[[column]] <- numeric(0)
    }
  }
  coherence$n_snps <- NULL

  out <- programs |>
    dplyr::left_join(coherence, by = "program") |>
    dplyr::left_join(n_snps_by_program, by = "program") |>
    dplyr::left_join(n_snps_filtered_by_program, by = "program") |>
    dplyr::mutate(
      n_snps = as.integer(tidyr::replace_na(n_snps, 0L)),
      n_snps_filtered = as.integer(tidyr::replace_na(n_snps_filtered, 0L)),
      size_pass = n_snps_filtered >= min_module_size,
      internal_pass = size_pass &
        is.finite(mean_internal_similarity) &
        mean_internal_similarity >= min_mean_internal &
        is.finite(connectedness) &
        connectedness >= min_connectedness
    )

  fit <- clustering_result$cluster_details$flash_fit
  obs_input <- clustering_result$ebmf_input
  if (is.null(obs_input)) {
    obs_input <- clustering_result$x_star
  }
  total_signal <- sum(abs(obs_input), na.rm = TRUE)
  F_pm <- fit$F_pm
  traits_signal <- if (!is.null(fit$L_pm)) {
    colSums(abs(fit$L_pm), na.rm = TRUE)
  } else {
    NULL
  }
  out$factor_strength <- vapply(out$program, function(pg) {
    snps <- assigned$snp_id[assigned$program == pg]
    if (is.null(F_pm) || is.null(traits_signal) || length(snps) == 0) {
      return(NA_real_)
    }
    snp_signal <- sum(abs(F_pm[snps, pg]), na.rm = TRUE)
    raw <- traits_signal[pg] * snp_signal
    if (is.finite(total_signal) && total_signal > 0) {
      raw / total_signal
    } else {
      NA_real_
    }
  }, numeric(1))
  out$factor_strength_per_snp <- out$factor_strength / sqrt(out$n_snps_filtered)

  out$replication <- rep(NA_real_, nrow(out))
  out$sd_replication <- rep(NA_real_, nrow(out))
  out$stability_checked <- rep(FALSE, nrow(out))
  out$stability_pass <- rep(TRUE, nrow(out))
  need_stability <- any(out$internal_pass %in% TRUE) && n_rep > 0L
  if (need_stability) {
    st <- stability_ebmf_programs(
      clustering_result,
      n_rep = n_rep,
      frac_traits = frac_traits,
      top_n = top_n,
      seed = seed,
      cores = cores,
      verbose = verbose
    )
    if (nrow(st) > 0) {
      st <- st |>
        dplyr::filter(program %in% out$program[out$internal_pass %in% TRUE]) |>
        dplyr::select(program, replication, sd_replication)
      out <- out |>
        dplyr::left_join(st, by = "program", suffix = c("", ".stb")) |>
        dplyr::mutate(
          replication = dplyr::coalesce(replication, replication.stb),
          sd_replication = dplyr::coalesce(sd_replication, sd_replication.stb),
          stability_checked = !is.na(replication),
          stability_pass = dplyr::if_else(
            stability_checked,
            replication >= stability_threshold,
            TRUE
          )
        ) |>
        dplyr::select(-replication.stb, -sd_replication.stb)
    }
  }

  evidence <- assigned |>
    dplyr::mutate(neg_log10_lfsr = -log10(pmax(lfsr, 10^-posterior_evidence_cap))) |>
    dplyr::group_by(program) |>
    dplyr::summarise(
      posterior_evidence = if (posterior_stat == "median") {
        stats::median(neg_log10_lfsr, na.rm = TRUE)
      } else {
        mean(neg_log10_lfsr, na.rm = TRUE)
      },
      .groups = "drop"
    )
  out <- out |>
    dplyr::left_join(evidence, by = "program")

  cor_mat <- NULL
  out$max_abs_factor_corr <- rep(NA_real_, nrow(out))
  out$redundant <- rep(FALSE, nrow(out))
  if (!is.null(fit) && fit$n_factors >= 2) {
    F_mat <- fit$F_pm
    if (!is.null(F_mat) && ncol(F_mat) >= 2) {
      cor_mat <- stats::cor(F_mat, use = "pairwise.complete.obs")
      diag(cor_mat) <- NA_real_
      max_abs <- apply(abs(cor_mat), 1, max, na.rm = TRUE)
      names(max_abs) <- as.character(seq_len(ncol(F_mat)))
      out$max_abs_factor_corr <- max_abs[as.character(out$program)]
      out$redundant <- !is.na(out$max_abs_factor_corr) &
        out$max_abs_factor_corr > redundancy_threshold
    }
  }

  pass_cols <- cbind(
    size = out$size_pass,
    internal = out$internal_pass,
    stability = out$stability_pass
  )
  fail_reason <- apply(!pass_cols, 1, function(f) {
    f[is.na(f)] <- TRUE
    failed <- names(f)[f]
    if (length(failed) == 0) {
      return("valid")
    }
    return(paste(failed, collapse = "; "))
  })
  out$fail_reason <- fail_reason
  out$status <- ifelse(fail_reason == "valid", "valid", fail_reason)

  out <- out |>
    dplyr::arrange(
      dplyr::desc(status == "valid"),
      dplyr::desc(factor_strength),
      dplyr::desc(n_snps)
    )

  return(list(
    programs = out,
    memberships = cal$memberships,
    assigned = assigned,
    factor_correlation = cor_mat,
    null_summary = cal$null_summary,
    settings = list(
      edge_threshold = edge_threshold,
      min_module_size = as.integer(min_module_size),
      min_mean_internal = min_mean_internal,
      min_connectedness = min_connectedness,
      n_null = n_null,
      alpha_membership = alpha_membership,
      n_candidate_tier = as.integer(n_candidate_tier),
      n_rep = n_rep,
      frac_traits = frac_traits,
      top_n = as.integer(top_n),
      stability_threshold = stability_threshold,
      posterior_stat = posterior_stat,
      posterior_evidence_cap = posterior_evidence_cap,
      redundancy_threshold = redundancy_threshold,
      seed = seed,
      cores = cores
    )
  ))
}


.program_internal_coherence <- function(s_matrix, assigned_snps, edge_threshold) {
  snps <- intersect(assigned_snps, colnames(s_matrix))
  n <- length(snps)
  if (n < 2) {
    return(data.frame(
      n_snps = n,
      mean_internal_similarity = NA_real_,
      median_internal_similarity = NA_real_,
      mean_external_similarity = NA_real_,
      separation = NA_real_,
      connectedness = NA_real_,
      n_internal_pairs = 0L,
      stringsAsFactors = FALSE
    ))
  }
  S <- s_matrix[snps, snps, drop = FALSE]
  diag(S) <- 1
  internal_vals <- S[upper.tri(S)]
  internal_vals <- internal_vals[is.finite(internal_vals)]
  n_pairs <- length(internal_vals)
  mean_int <- if (n_pairs > 0) mean(internal_vals) else NA_real_
  median_int <- if (n_pairs > 0) stats::median(internal_vals) else NA_real_
  connectedness <- if (n_pairs > 0) {
    mean(abs(internal_vals) >= edge_threshold)
  } else {
    NA_real_
  }
  rest <- setdiff(colnames(s_matrix), snps)
  mean_ext <- if (length(rest) > 0) {
    mean(s_matrix[snps, rest, drop = FALSE], na.rm = TRUE)
  } else {
    NA_real_
  }
  return(data.frame(
    n_snps = n,
    mean_internal_similarity = mean_int,
    median_internal_similarity = median_int,
    mean_external_similarity = mean_ext,
    separation = mean_int - mean_ext,
    connectedness = connectedness,
    n_internal_pairs = n_pairs,
    stringsAsFactors = FALSE
  ))
}
