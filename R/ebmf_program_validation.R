#' @title Validate EBMF Programs Against Filtering and Evidence Metrics
#' @description One-stop summary of EBMF programs for filtering and reporting.
#' Every check is calibrated against the observed matrix itself, so the same
#' thresholds transfer between datasets whose size, density and similarity
#' baseline differ.
#'
#' Filters applied (a program must pass all of them to get `status == "valid"`):
#' \itemize{
#'   \item **Size** — at least `min_module_size` SNPs passing the
#'   lFSR/magnitude filter (`n_snps_filtered`). An interpretability floor, not
#'   a discriminating gate.
#'   \item **Similarity** — the program's member SNPs must be more similar to
#'   each other than random SNP sets of the same size drawn from the same
#'   matrix. The statistic is the mean pairwise similarity `S` among the
#'   assigned SNPs (`mean_internal_similarity`); its null comes from
#'   `similarity_n_perm` random same-size sets, giving `similarity_emp_p` and
#'   a BH-adjusted `similarity_q`. Because the null is size-matched, a tight
#'   five-SNP program is not penalised for being small.
#'   \item **Trait-subsampling stability** — the top-loading SNPs must be
#'   recovered with replication `>= stability_threshold` (mean best Jaccard
#'   overlap) when `1 - frac_traits` of the traits are held out and flashier is
#'   refit. Only run for programs passing size and similarity unless
#'   `stability_all_programs = TRUE`.
#'   \item **SNP redundancy** — off by default (`snp_redundancy_threshold =
#'   Inf`); see `max_pair_redundancy` below.
#' }
#' Reported, not gated:
#' \itemize{
#'   \item **Strength** — whether the factor explains more of the matrix than
#'   factors fitted to structureless data do. `factor_pve` is the proportion
#'   of variance explained by the factor (flashier's `pve`); its null is the
#'   pooled `pve` of every factor from the row-permuted null fits in
#'   `calibrate_ebmf_programs()`, giving `strength_emp_p`, a BH-adjusted
#'   `strength_q` and the `strength_pass` flag. If the null fits produce no
#'   factors at all, noise cannot match any program and `strength_pass` is
#'   TRUE. Not a gate: on the Suzuki T2D benchmark it removed clearly real
#'   small programs, because a factor's variance explained scales with how
#'   many SNPs and traits it spans.
#'   \item `mean_internal_similarity` and `connectedness` (share of member
#'   pairs with `|S| >= edge_threshold`), with `internal_similarity_pass` and
#'   `connectedness_pass` against the fixed `min_mean_internal` /
#'   `min_connectedness` for comparison with earlier runs.
#'   \item `parent_program` / `parent_containment` — when at least 60% of a
#'   program's SNPs lie within 500 kb of the SNPs of a larger program, that
#'   program is reported as its parent. Nested sub-programs are kept; the
#'   parent link says where they sit.
#'   \item `posterior_evidence` — aggregate (`median` or `mean`) of
#'   `-log10(lfsr)` across a program's assigned SNPs.
#'   \item `max_abs_factor_corr` / `redundant` — the maximum absolute
#'   correlation of the program's SNP loading vector (`F_pm`) with any other
#'   program; `redundant` flags programs whose maximum correlation exceeds
#'   `redundancy_threshold`.
#'   \item `max_pair_redundancy` / `most_redundant_program` — the reciprocal
#'   SNP-containment score `min(shared / n_snps_A, shared / n_snps_B)` and its
#'   partner program. It gates `redundancy_pass` only when
#'   `snp_redundancy_threshold` is finite.
#' }
#' @param clustering_result Result of `run_univariate_clustering()`.
#' @param s_matrix SNP-by-SNP similarity matrix used for the similarity metrics.
#'   Defaults to `clustering_result$s_matrix`, which excludes the target row.
#' @param edge_threshold Absolute similarity used for connectedness. Defaults to
#'   `0.2`.
#' @param min_module_size Minimum assigned SNPs for `size_pass`. Defaults to
#'   `5`.
#' @param min_mean_internal Fixed mean internal similarity reported as
#'   `internal_similarity_pass`; does not gate `status`. Defaults to `0.3`.
#' @param min_connectedness Fixed connectedness reported as
#'   `connectedness_pass`; does not gate `status`. Defaults to `0.5`.
#' @param similarity_n_perm Random same-size SNP sets drawn for the similarity
#'   null. Defaults to `1000`.
#' @param similarity_q_threshold BH-adjusted q below which `similarity_pass` is
#'   TRUE. Defaults to `0.05`.
#' @param strength_q_threshold BH-adjusted q below which the reported
#'   `strength_pass` flag is TRUE (not a gate). Defaults to `0.05`.
#' @param n_null Number of row-permuted null fits, used for membership
#'   calibration and the strength null.
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
#' @param stability_all_programs If `TRUE`, report `replication` (and gate
#'   `stability_pass`) for every program, not only those passing size and
#'   similarity. The refits are the same either way, so this costs nothing
#'   extra. Defaults to `FALSE`.
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
#'   (`max_pair_redundancy`) before `redundancy_pass` fails. Defaults to `Inf`,
#'   which turns the gate off.
#' @param seed RNG seed.
#' @param cores Number of cores for the parallel trait-subsample refits (passed
#'   to `stability_ebmf_programs()`). Defaults to `1` (serial).
#' @param verbose Print progress messages.
#' @return A list with:
#'   \itemize{
#'     \item programs: the per-program table (metrics, pass/fail flags,
#'       `fail_reason` and `status`)
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
                                    similarity_n_perm = 1000L,
                                    similarity_q_threshold = 0.05,
                                    strength_q_threshold = 0.05,
                                    n_null = 10,
                                    alpha_membership = 0.05,
                                    n_candidate_tier = 25L,
                                    n_rep = 10,
                                    frac_traits = 0.8,
                                    top_n = 25L,
                                    stability_threshold = 0.7,
                                    stability_all_programs = FALSE,
                                    posterior_stat = c("median", "mean"),
                                    posterior_evidence_cap = 50,
                                    redundancy_threshold = 0.9,
                                    snp_redundancy_threshold = Inf,
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

  # calibrate_ebmf_programs() reports its own L1 factor strength; the strength
  # check here uses the variance explained instead, so drop it to avoid two
  # competing columns.
  programs <- cal$programs |>
    dplyr::select(-dplyr::any_of(c("raw_factor_signal", "factor_strength")))
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

  similarity_null <- .program_similarity_null(
    s_matrix,
    assigned,
    programs$program,
    n_perm = similarity_n_perm,
    seed = seed
  )

  out <- programs |>
    dplyr::left_join(coherence, by = "program") |>
    dplyr::left_join(similarity_null, by = "program") |>
    dplyr::left_join(n_snps_by_program, by = "program") |>
    dplyr::left_join(n_snps_filtered_by_program, by = "program") |>
    dplyr::mutate(
      n_snps = as.integer(tidyr::replace_na(n_snps, 0L)),
      n_snps_filtered = as.integer(tidyr::replace_na(n_snps_filtered, 0L)),
      size_pass = n_snps_filtered >= min_module_size,
      similarity_pass = is.finite(similarity_q) &
        similarity_q < similarity_q_threshold,
      internal_similarity_pass = is.finite(mean_internal_similarity) &
        mean_internal_similarity >= min_mean_internal,
      connectedness_pass = is.finite(connectedness) &
        connectedness >= min_connectedness
    )
  out$internal_pass <- out$size_pass & out$similarity_pass

  fit <- clustering_result$cluster_details$flash_fit

  null_pve <- cal$null_summary$factor_pve
  out$factor_pve <- if (!is.null(fit$pve)) {
    as.numeric(fit$pve[out$program])
  } else {
    rep(NA_real_, nrow(out))
  }
  if (length(null_pve) > 0) {
    out$strength_emp_p <- vapply(out$factor_pve, function(x) {
      if (!is.finite(x)) {
        return(NA_real_)
      }
      return((1 + sum(null_pve >= x)) / (1 + length(null_pve)))
    }, numeric(1))
    out$strength_q <- stats::p.adjust(out$strength_emp_p, method = "BH")
    out$strength_pass <- is.finite(out$strength_q) &
      out$strength_q < strength_q_threshold
  } else {
    out$strength_emp_p <- rep(NA_real_, nrow(out))
    out$strength_q <- rep(NA_real_, nrow(out))
    out$strength_pass <- is.finite(out$factor_pve)
  }

  out$replication <- rep(NA_real_, nrow(out))
  out$sd_replication <- rep(NA_real_, nrow(out))
  out$stability_checked <- rep(FALSE, nrow(out))
  out$stability_pass <- rep(TRUE, nrow(out))
  stability_programs <- if (isTRUE(stability_all_programs)) {
    out$program
  } else {
    out$program[out$internal_pass %in% TRUE]
  }
  need_stability <- length(stability_programs) > 0 && n_rep > 0L
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
        dplyr::filter(program %in% stability_programs) |>
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

  out <- out |>
    dplyr::left_join(
      .program_parents(assigned, clustering_result$snp_info),
      by = "program"
    )

  # strength_pass is reported, not gated: on the Suzuki T2D benchmark it
  # removed clearly real small programs.
  pass_cols <- cbind(
    size = out$size_pass,
    similarity = out$similarity_pass,
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
      dplyr::desc(similarity_z),
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
      similarity_n_perm = as.integer(similarity_n_perm),
      similarity_q_threshold = similarity_q_threshold,
      strength_q_threshold = strength_q_threshold,
      n_null = n_null,
      alpha_membership = alpha_membership,
      n_candidate_tier = as.integer(n_candidate_tier),
      n_rep = n_rep,
      frac_traits = frac_traits,
      top_n = as.integer(top_n),
      stability_threshold = stability_threshold,
      stability_all_programs = isTRUE(stability_all_programs),
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


# Size-matched permutation null for each program's mean internal similarity.
#
# The observed statistic is the mean pairwise S among a program's assigned
# SNPs. Its null is the same statistic over `n_perm` random SNP sets of the
# same size, drawn from SNPs that have any non-zero similarity to another SNP
# (a SNP observed only for the dropped target row carries no profile). Nulls
# are shared between programs of the same size. Returns one row per program:
# similarity_lift (observed / null mean), similarity_z, similarity_emp_p and
# the BH-adjusted similarity_q; programs with fewer than two SNPs get NA.
.program_similarity_null <- function(s_matrix, assigned, programs, n_perm = 1000L, seed = 1) {
  empty_row <- function(pg) {
    return(data.frame(
      program = as.integer(pg), similarity_lift = NA_real_,
      similarity_z = NA_real_, similarity_emp_p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  if (length(programs) == 0) {
    return(data.frame(
      program = integer(0), similarity_lift = numeric(0),
      similarity_z = numeric(0), similarity_emp_p = numeric(0),
      similarity_q = numeric(0), stringsAsFactors = FALSE
    ))
  }
  S <- s_matrix
  S[!is.finite(S)] <- 0
  off_diag <- S
  diag(off_diag) <- 0
  pool <- colnames(S)[colSums(off_diag != 0) > 0]
  mean_upper <- function(ids) {
    sub <- S[ids, ids, drop = FALSE]
    return(mean(sub[upper.tri(sub)]))
  }

  n_perm <- as.integer(n_perm)
  set.seed(seed)
  null_by_size <- list()
  rows <- lapply(programs, function(pg) {
    snps <- intersect(unique(assigned$snp_id[assigned$program == pg]), colnames(S))
    n <- length(snps)
    if (n < 2 || length(pool) < n || n_perm < 1) {
      return(empty_row(pg))
    }
    key <- as.character(n)
    if (is.null(null_by_size[[key]])) {
      null_by_size[[key]] <<- vapply(
        seq_len(n_perm),
        function(i) mean_upper(sample(pool, n)),
        numeric(1)
      )
    }
    null <- null_by_size[[key]]
    obs <- mean_upper(snps)
    null_mean <- mean(null)
    null_sd <- stats::sd(null)
    return(data.frame(
      program = as.integer(pg),
      similarity_lift = if (null_mean > 0) obs / null_mean else NA_real_,
      similarity_z = if (is.finite(null_sd) && null_sd > 0) (obs - null_mean) / null_sd else NA_real_,
      similarity_emp_p = (1 + sum(null >= obs)) / (1 + n_perm),
      stringsAsFactors = FALSE
    ))
  })
  out <- dplyr::bind_rows(rows)
  out$similarity_q <- stats::p.adjust(out$similarity_emp_p, method = "BH")
  return(out)
}


# Parent of each program: the larger program that contains most of its loci.
#
# For program A and every program B with more assigned SNPs, containment is the
# share of A's SNPs lying within `window` bp of any of B's SNPs on the same
# chromosome (a shared SNP counts at distance 0). The B with the highest
# containment (ties to the smaller B) is reported as A's parent when the
# containment reaches `min_containment`. Needs snp_info with snp_id, chr and bp;
# otherwise every parent is NA.
.program_parents <- function(assigned, snp_info, window = 5e5, min_containment = 0.6) {
  program_ids <- sort(unique(as.integer(assigned$program)))
  out <- data.frame(
    program = program_ids,
    parent_program = rep(NA_integer_, length(program_ids)),
    parent_containment = rep(NA_real_, length(program_ids)),
    stringsAsFactors = FALSE
  )
  if (length(program_ids) < 2 || is.null(snp_info) ||
      !all(c("snp_id", "chr", "bp") %in% names(snp_info))) {
    return(out)
  }
  positions <- snp_info |>
    dplyr::transmute(snp_id = as.character(snp_id), chr = as.character(chr), bp = as.numeric(bp)) |>
    dplyr::distinct(snp_id, .keep_all = TRUE)
  members <- lapply(program_ids, function(pg) {
    ids <- unique(as.character(assigned$snp_id[assigned$program == pg]))
    return(positions[match(ids, positions$snp_id), , drop = FALSE] |>
             dplyr::filter(!is.na(snp_id)))
  })
  sizes <- vapply(members, nrow, integer(1))

  for (a in seq_along(program_ids)) {
    pa <- members[[a]]
    if (nrow(pa) == 0) {
      next
    }
    candidates <- which(sizes > sizes[a])
    if (length(candidates) == 0) {
      next
    }
    containment <- vapply(candidates, function(b) {
      pb <- members[[b]]
      near <- vapply(seq_len(nrow(pa)), function(i) {
        return(any(pb$chr == pa$chr[i] & abs(pb$bp - pa$bp[i]) <= window))
      }, logical(1))
      return(mean(near))
    }, numeric(1))
    best <- candidates[order(-containment, sizes[candidates])[1]]
    best_value <- max(containment)
    if (best_value >= min_containment) {
      out$parent_program[a] <- program_ids[best]
      out$parent_containment[a] <- best_value
    }
  }
  return(out)
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
