#' @title Summarise EBMF Programs with a Confidence Tier
#' @description One-stop summary of EBMF programs for reporting. Programs with
#' fewer than `min_module_size` member SNPs are dropped first (and listed in
#' `dropped_programs`). Each remaining program gets a `confidence_tier` of
#' `"high"`, `"medium"` or `"low"` from a size-matched similarity test and
#' connectedness, or `"single_trait"` when it passes the similarity test but
#' its trait loadings rest on essentially one trait.
#'
#' Member SNPs are those passing the lFSR/magnitude filter. Similarities are
#' plain signed `S` (not flipped by loading sign and not absolute): on the Suzuki
#' T2D benchmark, flipping by loading sign rewarded noise programs whose few
#' opposite-sign members happened to be anti-similar, and weak members' loading
#' signs are not reliable enough to flip on.
#'
#' Confidence tier, checked in this order:
#' \itemize{
#'   \item **low** — `similarity_q >= similarity_q_threshold` (or missing):
#'   the member SNPs are no more alike than random SNP groups of the same
#'   size. A program that fails this is low whatever its number of traits, so
#'   it can be set aside.
#'   \item **single_trait** — passes the similarity test, but
#'   `effective_n_traits < min_effective_traits`. A program is defined as
#'   multi-trait; a factor carried by one trait is just that trait's hits,
#'   which noise produces too, and a set of SNPs sharing one trait looks
#'   coherent to the similarity checks whether or not it is real. It needs a
#'   closer look rather than being discarded, so it sorts above low.
#'   \item **high** — `similarity_q < high_similarity_q` and
#'   `connectedness >= high_min_connectedness`.
#'   \item **medium** — `similarity_q < similarity_q_threshold`, but not high.
#' }
#' `effective_n_traits` is `1 / sum(p^2)`, where `p` are the program's squared
#' trait loadings (`L_pm`, target row excluded) as shares of their total: 1 for
#' a program on one trait, `n` for `n` equally loaded traits. Near-duplicate
#' traits each count, so remove them from the matrix beforehand.
#' The checks behind the tier:
#' \itemize{
#'   \item **Similarity** — `mean_internal_similarity`, the mean `S` over
#'   member pairs, is compared with `similarity_n_perm` random SNP sets of the
#'   same size from the same trait. `similarity_emp_p` is the share of random
#'   sets at least as similar and `similarity_q` its BH adjustment across the
#'   kept programs. This is a significance test, so it favours larger programs;
#'   the trait's background similarity (see `background`) sets how hard it is
#'   to pass.
#'   \item **Connectedness** — share of member pairs with
#'   `S >= edge_threshold`.
#' }
#' Reported, not used for the tier:
#' \itemize{
#'   \item `excess_similarity` / `excess_connectedness` — the member mean
#'   similarity and connectedness minus their means in the random sets: effect
#'   sizes that do not grow with program size. `similarity_lift` (observed /
#'   random mean; NA when the random mean is at most 0.01) and `similarity_z`
#'   come from the same random sets.
#'   \item `factor_pve` — the proportion of variance explained by the factor
#'   (flashier's `pve`).
#'   \item `replication` / `sd_replication` — mean best Jaccard overlap of the
#'   program's top-loading SNPs with a factor refit after holding out
#'   `1 - frac_traits` of the traits, over `n_rep` refits.
#'   \item `parent_program` / `parent_containment` — when at least 90% of a
#'   program's SNPs lie within 500 kb of the SNPs of a larger program, that
#'   program is reported as its parent. Nested sub-programs are kept; the
#'   parent link says where they sit, and `order_programs_by_confidence()`
#'   lists them under it.
#'   \item `posterior_evidence` — aggregate (`median` or `mean`) of
#'   `-log10(lfsr)` across a program's assigned SNPs.
#'   \item `max_abs_factor_corr` / `redundant` — the maximum absolute
#'   correlation of the program's SNP loading vector (`F_pm`) with any other
#'   program; `redundant` flags programs whose maximum correlation exceeds
#'   `redundancy_threshold`.
#'   \item `max_pair_redundancy` / `most_redundant_program` — the reciprocal
#'   SNP-containment score `min(shared / n_snps_A, shared / n_snps_B)` and its
#'   partner program.
#' }
#' @param clustering_result Result of `run_univariate_clustering()`.
#' @param s_matrix SNP-by-SNP similarity matrix used for the similarity metrics.
#'   Defaults to `clustering_result$s_matrix`, which excludes the target row.
#' @param edge_threshold Similarity at or above which a member pair counts as
#'   connected. Defaults to `0.2`.
#' @param min_module_size Programs with fewer member SNPs are dropped before any
#'   other check. Defaults to `5`.
#' @param similarity_n_perm Random same-size SNP sets drawn for the similarity
#'   test. Defaults to `10000`, so that the BH-adjusted q can fall below
#'   `high_similarity_q` even when only a few programs sit at the permutation
#'   floor.
#' @param similarity_q_threshold BH-adjusted q below which a program is at
#'   least medium confidence. Defaults to `0.05`.
#' @param high_similarity_q BH-adjusted q below which (with enough
#'   connectedness) a program is high confidence. Defaults to `0.01`.
#' @param high_min_connectedness Minimum connectedness for high confidence.
#'   Defaults to `0.4`.
#' @param min_effective_traits Programs that pass the similarity test but
#'   whose `effective_n_traits` is below this are `"single_trait"`. Defaults
#'   to `1.5`.
#' @param n_null Number of row-permuted null fits used for membership
#'   calibration.
#' @param alpha_membership Membership FDR level passed to
#'   `calibrate_ebmf_programs()` for `core` cell labelling.
#' @param n_candidate_tier Per program, how many non-core cells to label as the
#'   uncertified "candidate" tier.
#' @param n_rep Number of trait-subsample refits for `replication`. Set to `0`
#'   to skip it.
#' @param frac_traits Fraction of trait rows sampled per replication refit.
#'   Defaults to `0.8` (hold out 20% of traits).
#' @param top_n Size of each program's top-loading SNP set used for
#'   replication matching.
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
#'     \item programs: the per-program table (metrics, `effective_n_traits` and
#'       `confidence_tier`), ordered by `order_programs_by_confidence()`
#'     \item dropped_programs: programs below `min_module_size` (`program`,
#'       `n_snps_filtered`), removed before any check
#'     \item background: the trait's background similarity over all SNP pairs
#'       (`mean_similarity`, `frac_pairs_edge` = share with `S >=
#'       edge_threshold`, `frac_pairs_abs_edge` = share with `|S| >=
#'       edge_threshold`, `median_shared_traits` per pair)
#'     \item memberships: calibrated SNP x program posterior from
#'       `calibrate_ebmf_programs()`
#'     \item assigned: gated assigned memberships of the kept programs
#'     \item factor_correlation: K x K correlation matrix of SNP loading vectors
#'     \item null_summary: permutation-null summary
#'     \item settings: settings used
#'   }
#' @export
summarise_ebmf_programs <- function(clustering_result,
                                    s_matrix = NULL,
                                    edge_threshold = 0.2,
                                    min_module_size = 5L,
                                    similarity_n_perm = 10000L,
                                    similarity_q_threshold = 0.05,
                                    high_similarity_q = 0.01,
                                    high_min_connectedness = 0.4,
                                    min_effective_traits = 1.5,
                                    n_null = 10,
                                    alpha_membership = 0.05,
                                    n_candidate_tier = 25L,
                                    n_rep = 10,
                                    frac_traits = 0.8,
                                    top_n = 25L,
                                    posterior_stat = c("median", "mean"),
                                    posterior_evidence_cap = 50,
                                    redundancy_threshold = 0.9,
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

  # calibrate_ebmf_programs() reports its own L1 factor strength; the variance
  # explained is reported here instead, so drop it to avoid two competing
  # columns.
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

  # Programs below the size floor are dropped before any check, so they never
  # enter the random-set baseline, the BH correction, parent detection or the
  # reported tables.
  size_by_program <- n_snps_filtered_by_program$n_snps_filtered[
    match(programs$program, n_snps_filtered_by_program$program)
  ]
  size_by_program[is.na(size_by_program)] <- 0L
  keep <- size_by_program >= min_module_size
  dropped_programs <- data.frame(
    program = programs$program[!keep],
    n_snps_filtered = as.integer(size_by_program[!keep]),
    stringsAsFactors = FALSE
  )
  programs <- programs[keep, , drop = FALSE]
  assigned <- assigned[assigned$program %in% programs$program, , drop = FALSE]

  members_of <- function(pg) {
    return(unique(as.character(assigned$snp_id[assigned$program == pg])))
  }

  coherence_rows <- lapply(programs$program, function(pg) {
    return(.program_internal_coherence(s_matrix, members_of(pg), edge_threshold))
  })
  coherence <- dplyr::bind_rows(coherence_rows)
  coherence_columns <- c(
    "n_snps", "mean_internal_similarity", "median_internal_similarity",
    "connectedness", "n_internal_pairs"
  )
  for (column in coherence_columns) {
    if (!column %in% names(coherence)) {
      coherence[[column]] <- numeric(0)
    }
  }
  coherence$program <- programs$program
  coherence$n_snps <- NULL

  similarity_null <- .program_similarity_null(
    s_matrix,
    lapply(stats::setNames(programs$program, programs$program), members_of),
    edge_threshold = edge_threshold,
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
      n_snps_filtered = as.integer(tidyr::replace_na(n_snps_filtered, 0L))
    )
  fit <- clustering_result$cluster_details$flash_fit
  out$effective_n_traits <- .program_effective_traits(
    fit, out$program, target_trait_id = params$target_trait_id
  )
  out$confidence_tier <- .program_confidence_tier(
    similarity_q = out$similarity_q,
    connectedness = out$connectedness,
    similarity_q_threshold = similarity_q_threshold,
    high_similarity_q = high_similarity_q,
    high_min_connectedness = high_min_connectedness,
    effective_n_traits = out$effective_n_traits,
    min_effective_traits = min_effective_traits
  )

  out$factor_pve <- if (!is.null(fit$pve)) {
    as.numeric(fit$pve[out$program])
  } else {
    rep(NA_real_, nrow(out))
  }

  out$replication <- rep(NA_real_, nrow(out))
  out$sd_replication <- rep(NA_real_, nrow(out))
  if (n_rep > 0L && nrow(out) > 0) {
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
        dplyr::select(program, replication, sd_replication)
      out <- out |>
        dplyr::left_join(st, by = "program", suffix = c("", ".stb")) |>
        dplyr::mutate(
          replication = dplyr::coalesce(replication, replication.stb),
          sd_replication = dplyr::coalesce(sd_replication, sd_replication.stb)
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
      )
    )

  out <- out |>
    dplyr::left_join(
      .program_parents(assigned, clustering_result$snp_info),
      by = "program"
    )

  out <- order_programs_by_confidence(out)

  return(list(
    programs = out,
    dropped_programs = dropped_programs,
    background = .similarity_background(
      s_matrix, clustering_result$overlap_matrix, edge_threshold
    ),
    memberships = cal$memberships,
    assigned = assigned,
    factor_correlation = cor_mat,
    null_summary = cal$null_summary,
    settings = list(
      edge_threshold = edge_threshold,
      min_module_size = as.integer(min_module_size),
      similarity_n_perm = as.integer(similarity_n_perm),
      similarity_q_threshold = similarity_q_threshold,
      high_similarity_q = high_similarity_q,
      high_min_connectedness = high_min_connectedness,
      min_effective_traits = min_effective_traits,
      n_null = n_null,
      alpha_membership = alpha_membership,
      n_candidate_tier = as.integer(n_candidate_tier),
      n_rep = n_rep,
      frac_traits = frac_traits,
      top_n = as.integer(top_n),
      posterior_stat = posterior_stat,
      posterior_evidence_cap = posterior_evidence_cap,
      redundancy_threshold = redundancy_threshold,
      seed = seed,
      cores = cores
    )
  ))
}


# Confidence tier from the similarity test and connectedness:
#   low          = fails the similarity test (similarity_q >=
#                  similarity_q_threshold, or missing), whatever the traits
#   single_trait = passes it, but effective_n_traits < min_effective_traits
#                  (skipped when effective_n_traits is NULL or missing)
#   high         = similarity_q < high_similarity_q &
#                  connectedness >= high_min_connectedness
#   medium       = similarity_q < similarity_q_threshold (not high)
# Returned as an ordered factor, levels c("low", "single_trait", "medium",
# "high").
.program_confidence_tier <- function(similarity_q,
                                     connectedness,
                                     similarity_q_threshold = 0.05,
                                     high_similarity_q = 0.01,
                                     high_min_connectedness = 0.4,
                                     effective_n_traits = NULL,
                                     min_effective_traits = 1.5) {
  medium <- is.finite(similarity_q) & similarity_q < similarity_q_threshold
  high <- is.finite(similarity_q) & similarity_q < high_similarity_q &
    is.finite(connectedness) & connectedness >= high_min_connectedness
  tier <- ifelse(high, "high", ifelse(medium, "medium", "low"))
  if (!is.null(effective_n_traits)) {
    single <- is.finite(effective_n_traits) &
      effective_n_traits < min_effective_traits
    tier[single & tier != "low"] <- "single_trait"
  }
  return(factor(
    tier,
    levels = c("low", "single_trait", "medium", "high"),
    ordered = TRUE
  ))
}


# Effective number of traits each program loads on: 1 / sum(p^2), where p are
# the program's squared trait loadings (L_pm, target row excluded) as shares of
# their total. One trait gives 1; n equally loaded traits give n. NA when the
# fit has no trait loadings for the program.
.program_effective_traits <- function(fit, programs, target_trait_id = NULL) {
  L <- if (!is.null(fit)) fit$L_pm else NULL
  if (is.null(L) || length(programs) == 0) {
    return(rep(NA_real_, length(programs)))
  }
  keep <- rep(TRUE, nrow(L))
  if (!is.null(target_trait_id) && !is.null(rownames(L))) {
    keep <- !rownames(L) %in% as.character(target_trait_id)
  }
  return(vapply(programs, function(k) {
    if (is.na(k) || k > ncol(L)) {
      return(NA_real_)
    }
    l2 <- L[keep, k]^2
    total <- sum(l2, na.rm = TRUE)
    if (!is.finite(total) || total <= 0) {
      return(NA_real_)
    }
    return(1 / sum((l2 / total)^2, na.rm = TRUE))
  }, numeric(1)))
}


# Similarity among a program's member SNPs: mean and median pairwise S, and
# connectedness (share of pairs with S >= edge_threshold).
.program_internal_coherence <- function(s_matrix, assigned_snps, edge_threshold) {
  snps <- intersect(assigned_snps, colnames(s_matrix))
  n <- length(snps)
  if (n < 2) {
    return(data.frame(
      n_snps = n,
      mean_internal_similarity = NA_real_,
      median_internal_similarity = NA_real_,
      connectedness = NA_real_,
      n_internal_pairs = 0L,
      stringsAsFactors = FALSE
    ))
  }
  S <- s_matrix[snps, snps, drop = FALSE]
  internal_vals <- S[upper.tri(S)]
  internal_vals <- internal_vals[is.finite(internal_vals)]
  n_pairs <- length(internal_vals)
  return(data.frame(
    n_snps = n,
    mean_internal_similarity = if (n_pairs > 0) mean(internal_vals) else NA_real_,
    median_internal_similarity = if (n_pairs > 0) stats::median(internal_vals) else NA_real_,
    connectedness = if (n_pairs > 0) mean(internal_vals >= edge_threshold) else NA_real_,
    n_internal_pairs = n_pairs,
    stringsAsFactors = FALSE
  ))
}


# Size-matched random-set test for each program's similarity and connectedness.
#
# `members_by_program` is a named list (names = program ids) of member SNP ids.
# For a program of n SNPs, each of `n_perm` draws takes n random SNPs (from
# those with any non-zero similarity to another SNP), so the baseline is
# matched on size and drawn from the trait's own similarity graph. Programs of
# the same size share one baseline. Returns one row per program:
# excess_similarity and excess_connectedness (observed minus random mean),
# similarity_lift (observed / random mean, NA when that mean is <= 0.01),
# similarity_z, similarity_emp_p and the BH-adjusted similarity_q.
.program_similarity_null <- function(s_matrix,
                                     members_by_program,
                                     edge_threshold = 0.2,
                                     n_perm = 10000L,
                                     seed = 1) {
  columns <- c(
    "excess_similarity", "excess_connectedness", "similarity_lift",
    "similarity_z", "similarity_emp_p"
  )
  empty_row <- function(pg) {
    row <- data.frame(program = as.integer(pg), stringsAsFactors = FALSE)
    row[columns] <- NA_real_
    return(row)
  }
  if (length(members_by_program) == 0) {
    out <- data.frame(program = integer(0), stringsAsFactors = FALSE)
    out[c(columns, "similarity_q")] <- list(numeric(0))
    return(out)
  }
  S <- s_matrix
  S[!is.finite(S)] <- 0
  off_diag <- S
  diag(off_diag) <- 0
  pool <- colnames(S)[colSums(off_diag != 0) > 0]
  pair_stats <- function(ids) {
    sub <- S[ids, ids, drop = FALSE]
    vals <- sub[upper.tri(sub)]
    return(c(mean(vals), mean(vals >= edge_threshold)))
  }

  n_perm <- as.integer(n_perm)
  set.seed(seed)
  null_by_size <- list()
  rows <- lapply(names(members_by_program), function(pg) {
    snps <- intersect(members_by_program[[pg]], colnames(S))
    n <- length(snps)
    if (n < 2 || length(pool) < n || n_perm < 1) {
      return(empty_row(pg))
    }
    key <- as.character(n)
    if (is.null(null_by_size[[key]])) {
      null_by_size[[key]] <<- vapply(
        seq_len(n_perm),
        function(i) pair_stats(sample(pool, n)),
        numeric(2)
      )
    }
    null <- null_by_size[[key]]
    obs <- pair_stats(snps)
    null_mean <- mean(null[1, ])
    null_sd <- stats::sd(null[1, ])
    return(data.frame(
      program = as.integer(pg),
      excess_similarity = obs[1] - null_mean,
      excess_connectedness = obs[2] - mean(null[2, ]),
      similarity_lift = if (null_mean > 0.01) obs[1] / null_mean else NA_real_,
      similarity_z = if (is.finite(null_sd) && null_sd > 0) (obs[1] - null_mean) / null_sd else NA_real_,
      similarity_emp_p = (1 + sum(null[1, ] >= obs[1])) / (1 + n_perm),
      stringsAsFactors = FALSE
    ))
  })
  out <- dplyr::bind_rows(rows)
  out$similarity_q <- stats::p.adjust(out$similarity_emp_p, method = "BH")
  return(out)
}


# Background similarity of a trait's SNPs: mean S over all pairs of SNPs with
# any profile, the share of pairs with S >= edge_threshold and with |S| >=
# edge_threshold, and the median number of traits two SNPs are both observed on
# (from the overlap matrix, when supplied).
.similarity_background <- function(s_matrix, overlap_matrix, edge_threshold) {
  S <- s_matrix
  S[!is.finite(S)] <- 0
  off_diag <- S
  diag(off_diag) <- 0
  pool <- colnames(S)[colSums(off_diag != 0) > 0]
  if (length(pool) < 2) {
    return(data.frame(
      n_snps = length(pool), mean_similarity = NA_real_,
      frac_pairs_edge = NA_real_, frac_pairs_abs_edge = NA_real_,
      median_shared_traits = NA_real_
    ))
  }
  sub <- S[pool, pool, drop = FALSE]
  vals <- sub[upper.tri(sub)]
  shared <- NA_real_
  if (!is.null(overlap_matrix) && all(pool %in% colnames(overlap_matrix))) {
    ov <- overlap_matrix[pool, pool, drop = FALSE]
    shared <- stats::median(ov[upper.tri(ov)])
  }
  return(data.frame(
    n_snps = length(pool),
    mean_similarity = mean(vals),
    frac_pairs_edge = mean(vals >= edge_threshold),
    frac_pairs_abs_edge = mean(abs(vals) >= edge_threshold),
    median_shared_traits = shared
  ))
}


# Parent of each program: the larger program that contains nearly all of its
# loci.
#
# For program A and every program B with more assigned SNPs, containment is the
# share of A's SNPs lying within `window` bp of any of B's SNPs on the same
# chromosome (a shared SNP counts at distance 0). The B with the highest
# containment (ties to the smaller B) is reported as A's parent when the
# containment reaches `min_containment`. Needs snp_info with snp_id, chr and bp;
# otherwise every parent is NA.
.program_parents <- function(assigned, snp_info, window = 5e5, min_containment = 0.9) {
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
