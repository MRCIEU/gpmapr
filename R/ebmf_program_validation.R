#' @title Validate EBMF Programs Against Filtering and Evidence Metrics
#' @description One-stop summary of EBMF programs for filtering and reporting.
#' Folds the graph-reliability metrics (`mean_internal_similarity`,
#' `connectedness`) together with trait-subsampling **stability**, and the
#' descriptive **additional scores** (factor strength, aggregate posterior
#' evidence and factor distinctiveness) into a single per-program table with
#' pass/fail flags and an overall `status`.
#'
#' Filters applied (a program must pass all four to get `status == "valid"`):
#' \itemize{
#'   \item **Size** — at least `min_module_size` SNPs passing the
#'   lFSR/magnitude filter (`n_snps_filtered`). This is an interpretability
#'   floor, not a discriminating gate: a two-SNP program cannot be annotated.
#'   \item **Coherence** — the program's loading vector must be aligned with
#'   the SNP similarity graph more than a random reassignment of the same
#'   loadings would be (`coherence_q < coherence_q_threshold`, from
#'   `calibrate_program_coherence()`). This uses every SNP weighted by its
#'   squared posterior loading, so it does not depend on the lFSR/magnitude
#'   cutoff, and it is calibrated against the graph's own null, so it transfers
#'   between datasets whose similarity baselines differ. `size_pass` and
#'   `coherence_pass` are reported separately and `internal_pass` combines them.
#'   The hard-gated `mean_internal_similarity` and `connectedness` are still
#'   reported for comparison but no longer decide `status`.
#'   \item **Trait-subsampling stability** — the top-loading SNPs must be
#'   recovered with replication `>= stability_threshold` when `1 - frac_traits`
#'   of the traits are held out and flashier is refit. Because this refits
#'   flashier `n_rep` times, it is only run for programs that pass the
#'   internal-coherence checks; programs failing coherence are already rejected
#'   and are not re-checked.
#'   \item **Factor distinctiveness** — for each program, compare its assigned
#'   (lFSR/magnitude-filtered) SNP membership set to every other program via the
#'   reciprocal containment score `pair_redundancy(A,B) = min(shared / n_snps_A,
#'   shared / n_snps_B)`; `max_pair_redundancy` is the maximum such score across
#'   all other programs and `most_redundant_program` is the partner that
#'   achieves it. The program fails (`redundancy_pass == FALSE`) when
#'   `max_pair_redundancy >= snp_redundancy_threshold`.
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
#'   \item `max_pair_redundancy` / `most_redundant_program` — the reciprocal
#'   SNP-containment redundancy score and its partner program; reported for
#'   context alongside the gating `redundancy_pass` described above.
#' }
#' @param clustering_result Result of `run_univariate_clustering()`.
#' @param s_matrix SNP-by-SNP similarity matrix used for the coherence metrics.
#'   Defaults to `clustering_result$s_matrix`.
#' @param edge_threshold Absolute similarity used for connectedness. Defaults to
#'   `0.2`.
#' @param min_module_size Minimum assigned SNPs for the coherence metrics to be
#'   meaningful and for `size_pass`. Defaults to `5`.
#' @param min_mean_internal Minimum mean internal SNP similarity. Reported as
#'   `internal_similarity_pass` for comparison with earlier runs; it no longer
#'   gates `status`. Defaults to `0.3`.
#' @param min_connectedness Minimum pair-connectedness. Reported as
#'   `connectedness_pass` for comparison; it no longer gates `status`. Defaults
#'   to `0.5`.
#' @param coherence_n_perm Permutations per program for the loading-permutation
#'   coherence null. Defaults to `2000`.
#' @param coherence_q_threshold BH-adjusted q below which `coherence_pass` is
#'   TRUE. Defaults to `0.05`.
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
#' @param snp_redundancy_threshold Maximum reciprocal SNP-containment
#'   (`max_pair_redundancy`) before `redundancy_pass` fails. Defaults to `0.75`;
#'   at `0.9` the check never bound on real data, where the observed maximum
#'   containment was `0.5`.
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
                                    coherence_n_perm = 2000L,
                                    coherence_q_threshold = 0.05,
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
                                    snp_redundancy_threshold = 0.75,
                                    seed = 1,
                                    cores = 1,
                                    verbose = TRUE) {
  .assert_ebmf_result(clustering_result)
  params <- clustering_result$parameters
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
      internal_similarity_pass = is.finite(mean_internal_similarity) &
        mean_internal_similarity >= min_mean_internal,
      connectedness_pass = is.finite(connectedness) &
        connectedness >= min_connectedness
    )

  fit <- clustering_result$cluster_details$flash_fit

  # Loading-weighted coherence against a permutation null. This is the gate;
  # the hard-gated mean_internal_similarity / connectedness above are retained
  # for comparison but no longer decide status.
  coherence_cal <- calibrate_program_coherence(
    s_matrix,
    f_pm = if (is.null(fit)) NULL else fit$F_pm,
    n_perm = coherence_n_perm,
    edge_threshold = edge_threshold,
    seed = seed
  )
  coherence_columns <- c(
    "weighted_internal", "weighted_quorum", "n_eff",
    "coherence_emp_p", "coherence_q", "quorum_emp_p", "quorum_q"
  )
  if (nrow(coherence_cal) == 0) {
    for (column in coherence_columns) {
      out[[column]] <- rep(NA_real_, nrow(out))
    }
  } else {
    out <- out |>
      dplyr::left_join(coherence_cal, by = "program")
  }
  out$coherence_pass <- is.finite(out$coherence_q) &
    out$coherence_q < coherence_q_threshold
  out$quorum_pass <- is.finite(out$quorum_q) &
    out$quorum_q < coherence_q_threshold
  # The weighted mean and the weighted quorum measure different things --
  # average alignment versus how many pairs clear the edge threshold -- and once
  # null-calibrated they do not always agree. Only the mean gates; this flags
  # where the two readings diverge.
  out$shape_disagreement <- out$coherence_pass != out$quorum_pass
  out$internal_pass <- out$size_pass & out$coherence_pass
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

  redundancy_by_program <- .program_membership_redundancy(assigned)
  out <- out |>
    dplyr::left_join(redundancy_by_program, by = "program") |>
    dplyr::mutate(
      max_pair_redundancy = tidyr::replace_na(max_pair_redundancy, NA_real_),
      most_redundant_program = as.integer(
        tidyr::replace_na(most_redundant_program, NA_integer_)
      ),
      redundancy_pass = is.na(max_pair_redundancy) |
        max_pair_redundancy < snp_redundancy_threshold
    )

  pass_cols <- cbind(
    size = out$size_pass,
    coherence = out$coherence_pass,
    stability = out$stability_pass,
    redundancy = out$redundancy_pass
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
  out$status <- as.character(
    ifelse(fail_reason == "valid", "valid", fail_reason)
  )

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
      coherence_n_perm = as.integer(coherence_n_perm),
      coherence_q_threshold = coherence_q_threshold,
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
      snp_redundancy_threshold = snp_redundancy_threshold,
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


# Loading-weighted coherence of one program on the SNP similarity graph.
#
# Unlike .program_internal_coherence(), which first binarises the program with
# the lFSR/magnitude gate and then averages over the surviving SNPs, this uses
# *every* SNP, weighted by w (normally the squared posterior loading). The
# statistics are the weighted analogues of the hard-gated ones:
#
#   weighted_internal = sum_{i!=j} w_i w_j S_ij / sum_{i!=j} w_i w_j
#   weighted_quorum   = the same with S replaced by 1{|S| >= edge_threshold}
#
# With w normalised to sum 1 and the diagonal of S set to 1, both reduce to
# (w'Sw - sum w_i^2) / (1 - sum w_i^2).
#
# n_eff = 1 / sum(w_i^2) is the participation ratio: the weighted analogue of
# program size, and the number to read when asking how many SNPs a program is
# really resting on.
.program_loading_coherence <- function(s_matrix, w, edge_threshold) {
  empty <- data.frame(
    weighted_internal = NA_real_,
    weighted_quorum = NA_real_,
    n_eff = NA_real_,
    stringsAsFactors = FALSE
  )
  w <- .normalise_coherence_weights(w)
  if (is.null(w)) {
    return(empty)
  }
  sum_sq <- sum(w^2)
  denom <- 1 - sum_sq
  # n_eff stays reportable even when the coherence statistics are not: a
  # program whose weight sits on a single SNP has no off-diagonal pairs, and
  # n_eff = 1 alongside NA coherence says exactly that.
  if (!is.finite(denom) || denom <= 0) {
    empty$n_eff <- 1 / sum_sq
    return(empty)
  }
  s_diag1 <- s_matrix
  s_diag1[!is.finite(s_diag1)] <- 0
  diag(s_diag1) <- 1
  adjacency <- (abs(s_diag1) >= edge_threshold) * 1
  diag(adjacency) <- 1
  return(data.frame(
    weighted_internal = (sum(w * (s_diag1 %*% w)) - sum_sq) / denom,
    weighted_quorum = (sum(w * (adjacency %*% w)) - sum_sq) / denom,
    n_eff = 1 / sum_sq,
    stringsAsFactors = FALSE
  ))
}


# Non-negative weights summing to 1, or NULL when the program carries no
# usable loading mass.
.normalise_coherence_weights <- function(w) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  total <- sum(w)
  if (!is.finite(total) || total <= 0) {
    return(NULL)
  }
  return(w / total)
}


#' @title Permutation-Calibrate EBMF Program Coherence
#' @description Test whether each program's loading vector is aligned with the
#' SNP similarity graph more than a random reassignment of the same loadings
#' would be. For every program the weights `w = F_pm[, k]^2` are permuted across
#' the SNP index `n_perm` times and the loading-weighted coherence recomputed,
#' giving an empirical p-value and a BH-adjusted q-value.
#'
#' This replaces comparing `mean_internal_similarity` to an absolute constant.
#' Because permuting `w` preserves its distribution, the null is automatically
#' matched on program size and on the shape of the loading vector, so the same
#' threshold transfers between datasets whose similarity graphs have different
#' baselines — which absolute thresholds do not.
#' @param s_matrix SNP-by-SNP similarity matrix.
#' @param f_pm SNP-by-program posterior mean loading matrix (`flash_fit$F_pm`).
#' @param n_perm Number of permutations per program. Defaults to `2000`.
#' @param edge_threshold Absolute similarity defining an edge for the weighted
#'   quorum statistic. Defaults to `0.2`.
#' @param seed RNG seed.
#' @param chunk_size Permutations evaluated per matrix product, bounding peak
#'   memory. Defaults to `500`.
#' @return A dataframe with one row per program: `program`,
#'   `weighted_internal`, `weighted_quorum`, `n_eff`, `coherence_emp_p`,
#'   `coherence_q`, `quorum_emp_p` and `quorum_q`.
#' @export
calibrate_program_coherence <- function(s_matrix,
                                        f_pm,
                                        n_perm = 2000L,
                                        edge_threshold = 0.2,
                                        seed = 1,
                                        chunk_size = 500L) {
  empty <- data.frame(
    program = integer(0),
    weighted_internal = numeric(0),
    weighted_quorum = numeric(0),
    n_eff = numeric(0),
    coherence_emp_p = numeric(0),
    coherence_q = numeric(0),
    quorum_emp_p = numeric(0),
    quorum_q = numeric(0),
    stringsAsFactors = FALSE
  )
  if (is.null(s_matrix) || is.null(f_pm) || ncol(f_pm) == 0) {
    return(empty)
  }
  n_perm <- as.integer(n_perm)
  chunk_size <- max(1L, as.integer(chunk_size))

  snps <- intersect(rownames(f_pm), colnames(s_matrix))
  if (length(snps) < 3) {
    return(empty)
  }
  s_sub <- s_matrix[snps, snps, drop = FALSE]
  s_sub[!is.finite(s_sub)] <- 0
  diag(s_sub) <- 1
  adjacency <- (abs(s_sub) >= edge_threshold) * 1
  diag(adjacency) <- 1
  f_sub <- f_pm[snps, , drop = FALSE]

  n <- length(snps)
  set.seed(seed)
  rows <- lapply(seq_len(ncol(f_sub)), function(k) {
    w <- .normalise_coherence_weights(f_sub[, k]^2)
    if (is.null(w)) {
      return(data.frame(
        program = k, weighted_internal = NA_real_, weighted_quorum = NA_real_,
        n_eff = NA_real_, coherence_emp_p = NA_real_, quorum_emp_p = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    sum_sq <- sum(w^2)
    denom <- 1 - sum_sq
    if (!is.finite(denom) || denom <= 0) {
      return(data.frame(
        program = k, weighted_internal = NA_real_, weighted_quorum = NA_real_,
        n_eff = 1 / sum_sq, coherence_emp_p = NA_real_,
        quorum_emp_p = NA_real_, stringsAsFactors = FALSE
      ))
    }
    obs_int <- (sum(w * (s_sub %*% w)) - sum_sq) / denom
    obs_quo <- (sum(w * (adjacency %*% w)) - sum_sq) / denom

    # A permutation preserves sum(w^2), so the denominator is constant and the
    # null is evaluated in chunks of the permuted weight matrix.
    ge_int <- 0L
    ge_quo <- 0L
    remaining <- n_perm
    while (remaining > 0) {
      this_chunk <- min(chunk_size, remaining)
      perm <- matrix(
        w[unlist(lapply(seq_len(this_chunk), function(i) sample.int(n)))],
        nrow = n
      )
      null_int <- (colSums(perm * (s_sub %*% perm)) - sum_sq) / denom
      null_quo <- (colSums(perm * (adjacency %*% perm)) - sum_sq) / denom
      ge_int <- ge_int + sum(null_int >= obs_int, na.rm = TRUE)
      ge_quo <- ge_quo + sum(null_quo >= obs_quo, na.rm = TRUE)
      remaining <- remaining - this_chunk
    }
    return(data.frame(
      program = k,
      weighted_internal = obs_int,
      weighted_quorum = obs_quo,
      n_eff = 1 / sum_sq,
      coherence_emp_p = (1 + ge_int) / (1 + n_perm),
      quorum_emp_p = (1 + ge_quo) / (1 + n_perm),
      stringsAsFactors = FALSE
    ))
  })
  out <- dplyr::bind_rows(rows)
  out$coherence_q <- stats::p.adjust(out$coherence_emp_p, method = "BH")
  out$quorum_q <- stats::p.adjust(out$quorum_emp_p, method = "BH")
  out$program <- as.integer(out$program)
  return(out[, names(empty), drop = FALSE])
}


# Reciprocal SNP-containment redundancy between programs, computed on the
# assigned (lFSR/magnitude-gated) membership set. For every pair:
#   pair_redundancy(A, B) = min(|A n B| / |A|, |A n B| / |B|)
# and each program reports its maximum pair score plus the partner that
# achieves it. Programs with no assigned SNPs (or only one program) have no
# comparable partner and get NA / NA_integer_.
.program_membership_redundancy <- function(assigned) {
  empty <- data.frame(
    program = integer(0),
    max_pair_redundancy = numeric(0),
    most_redundant_program = integer(0),
    stringsAsFactors = FALSE
  )
  if (is.null(assigned) || nrow(assigned) == 0) {
    return(empty)
  }

  membership_wide <- as.matrix(table(assigned$snp_id, assigned$program))
  program_ids <- as.integer(colnames(membership_wide))
  if (length(program_ids) == 0) {
    return(empty)
  }
  if (length(program_ids) == 1) {
    return(data.frame(
      program = program_ids,
      max_pair_redundancy = NA_real_,
      most_redundant_program = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  M <- membership_wide
  sizes <- colSums(M)
  shared <- as.matrix(crossprod(M))

  n_prog <- ncol(M)
  max_pair_redundancy <- rep(NA_real_, n_prog)
  most_redundant_program <- rep(NA_integer_, n_prog)
  if (n_prog >= 2) {
    for (i in seq_len(n_prog)) {
      if (sizes[i] == 0) next
      row_vals <- rep(NA_real_, n_prog)
      for (j in seq_len(n_prog)) {
        if (i == j || sizes[j] == 0) next
        s <- shared[i, j]
        row_vals[j] <- min(s / sizes[i], s / sizes[j])
      }
      if (all(!is.finite(row_vals))) next
      best <- which.max(row_vals)
      max_pair_redundancy[i] <- row_vals[best]
      most_redundant_program[i] <- program_ids[best]
    }
  }

  return(data.frame(
    program = program_ids,
    max_pair_redundancy = max_pair_redundancy,
    most_redundant_program = most_redundant_program,
    stringsAsFactors = FALSE
  ))
}
