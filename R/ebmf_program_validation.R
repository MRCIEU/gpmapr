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
#'   \item **high** — `similarity_q < high_similarity_q`, and either
#'   `connectedness >= high_min_connectedness` or `similarity_z >=
#'   high_min_similarity_z`. On real data large programs are looser (on the
#'   Suzuki T2D benchmark connectedness fell as programs grew), so
#'   connectedness alone demotes them; a large `similarity_z` lets a program
#'   that is far above its size-matched random sets reach high without it.
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
#'   \item **similarity_z** — `(mean_internal_similarity - random mean) /
#'   random SD`, from the same random sets as `similarity_q`. Unlike
#'   `similarity_emp_p`, it has no floor at `1 / similarity_n_perm`, so it keeps
#'   ranking programs that beat every random set.
#' }
#' Reported, not used for the tier:
#' \itemize{
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
#'   \item `posterior_evidence` — median of `-log10(lfsr)` across a program's
#'   assigned SNPs, with lFSR floored at `1e-50` because flashier's lFSR can
#'   underflow towards 0.
#' }
#' @param clustering_result Result of `run_univariate_clustering()`.
#' @param s_matrix SNP-by-SNP similarity matrix used for the similarity metrics.
#'   Defaults to `clustering_result$s_matrix`, which excludes the target row.
#' @param edge_threshold Similarity at or above which a member pair counts as
#'   connected. Defaults to `0.2`.
#' @param min_module_size Programs with fewer member SNPs are dropped before any
#'   other check. Defaults to `3`.
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
#' @param high_min_similarity_z `similarity_z` at or above which a program
#'   that passes `high_similarity_q` is high confidence even below
#'   `high_min_connectedness`. Defaults to `6`, just above every multi-trait
#'   null program that passed `high_similarity_q` in 50 null simulations
#'   (maximum 5.3); `Inf` restores the connectedness-only rule.
#' @param min_effective_traits Programs that pass the similarity test but
#'   whose `effective_n_traits` is below this are `"single_trait"`. Defaults
#'   to `1.5`.
#' @param n_rep Number of trait-subsample refits for `replication`. Set to `0`
#'   to skip it.
#' @param frac_traits Fraction of trait rows sampled per replication refit.
#'   Defaults to `0.8` (hold out 20% of traits).
#' @param top_n Size of each program's top-loading SNP set used for
#'   replication matching. Defaults to `20`.
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
#'     \item assigned: gated assigned memberships of the kept programs
#'     \item settings: settings used
#'   }
#' @export
summarise_ebmf_programs <- function(clustering_result,
                                    s_matrix = NULL,
                                    edge_threshold = 0.2,
                                    min_module_size = 3L,
                                    similarity_n_perm = 10000L,
                                    similarity_q_threshold = 0.05,
                                    high_similarity_q = 0.01,
                                    high_min_connectedness = 0.4,
                                    high_min_similarity_z = 6,
                                    min_effective_traits = 1.5,
                                    n_rep = 10,
                                    frac_traits = 0.8,
                                    top_n = 20L,
                                    seed = 1,
                                    cores = 1,
                                    verbose = TRUE) {
  .assert_ebmf_result(clustering_result)
  params <- clustering_result$parameters

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

  posterior <- ebmf_posterior_table(clustering_result)
  ebmf_memberships <- posterior |>
    dplyr::filter(is.finite(abs_loading), abs_loading > 0)
  assigned <- ebmf_memberships |>
    dplyr::filter(
      !is.na(lfsr),
      lfsr < lfsr_threshold,
      is.na(mag_threshold) | abs_loading > mag_threshold
    )

  programs <- data.frame(
    program = sort(unique(as.integer(posterior$program))),
    stringsAsFactors = FALSE
  )

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
    "n_snps", "mean_internal_similarity", "connectedness", "n_internal_pairs"
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
    similarity_z = out$similarity_z,
    high_min_similarity_z = high_min_similarity_z,
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

  # flashier's lFSR can underflow towards 0, so -log10(lfsr) is capped at 50.
  evidence <- assigned |>
    dplyr::mutate(neg_log10_lfsr = -log10(pmax(lfsr, 1e-50))) |>
    dplyr::group_by(program) |>
    dplyr::summarise(
      posterior_evidence = stats::median(neg_log10_lfsr, na.rm = TRUE),
      .groups = "drop"
    )
  out <- out |>
    dplyr::left_join(evidence, by = "program")

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
    assigned = assigned,
    settings = list(
      edge_threshold = edge_threshold,
      min_module_size = as.integer(min_module_size),
      similarity_n_perm = as.integer(similarity_n_perm),
      similarity_q_threshold = similarity_q_threshold,
      high_similarity_q = high_similarity_q,
      high_min_connectedness = high_min_connectedness,
      high_min_similarity_z = high_min_similarity_z,
      min_effective_traits = min_effective_traits,
      n_rep = n_rep,
      frac_traits = frac_traits,
      top_n = as.integer(top_n),
      seed = seed,
      cores = cores
    )
  ))
}


#' @title Distinctness of a Trait's EBMF Programs
#' @description How alike the programs of one EBMF fit are, on both axes of the
#' fit. For every pair of programs:
#' \itemize{
#'   \item `phi` — Tucker's congruence coefficient between the two programs'
#'   trait loadings (`L_pm`), target row left out.
#'   \item `alignment` — the squared uncentred correlation of their SNP
#'   loadings (`F_pm`).
#' }
#' Loadings are weighted by their confidence, \eqn{w(\mathrm{lfsr})}, as in
#' [compare_program_pairs_loadings()], so the values sit on the same scale as
#' the cross-trait link tiers (Tucker's conventions: 0.85 fair similarity, 0.95
#' equal).
#'
#' `effective_rank` summarises the whole set: the exponential of the entropy
#' of the normalised eigenvalues of the programs' Gram matrix, whose entries
#' are \eqn{\phi \cdot \cos_F}, the cosine between two programs' rank-one fits
#' (Roy & Vetterli 2007). It equals the number of programs when they are
#' mutually orthogonal and falls as programs overlap: two exact copies of one
#' program count as a single program of double weight, which the entropy
#' scores as somewhat less than two. A pair whose cosine is undefined (a
#' program with no confident loading) counts as orthogonal.
#' @param clustering_result Result of `run_univariate_clustering()`.
#' @param programs Optional program (factor) ids to compare, typically
#'   `summarise_ebmf_programs()$programs$program`. Defaults to every factor.
#' @return A list with:
#'   \itemize{
#'     \item pairs: one row per program pair with `program_a`, `program_b`,
#'       `phi` and `alignment`
#'     \item summary: one row with `n_programs`, `max_abs_phi`,
#'       `n_pairs_phi_0.7` (pairs with \eqn{|\phi| \ge 0.7}), `max_alignment`
#'       and `effective_rank`
#'   }
#' @export
program_distinctness <- function(clustering_result, programs = NULL) {
  .assert_ebmf_result(clustering_result)
  fit <- clustering_result$cluster_details$flash_fit
  n_factors <- if (is.null(fit)) 0L else fit$n_factors
  if (is.null(programs)) {
    programs <- seq_len(n_factors)
  }
  programs <- sort(unique(as.integer(programs)))
  if (any(is.na(programs) | programs < 1L | programs > n_factors)) {
    stop("programs must be factor ids between 1 and ", n_factors)
  }

  L <- fit$L_pm[, programs, drop = FALSE]
  F <- fit$F_pm[, programs, drop = FALSE]
  if (!is.null(fit$L_lfsr)) {
    L <- L * .confidence_factor(fit$L_lfsr[, programs, drop = FALSE])
  }
  if (!is.null(fit$F_lfsr)) {
    F <- F * .confidence_factor(fit$F_lfsr[, programs, drop = FALSE])
  }
  target_id <- clustering_result$parameters$target_trait_id
  if (!is.null(target_id) && !is.null(rownames(L))) {
    L <- L[rownames(L) != as.character(target_id), , drop = FALSE]
  }

  cosine <- function(A) {
    norms <- sqrt(colSums(A^2))
    out <- crossprod(A) / outer(norms, norms)
    out[!is.finite(out)] <- NA_real_
    return(out)
  }
  phi <- cosine(L)
  cos_f <- cosine(F)

  n <- length(programs)
  pairs <- if (n >= 2) {
    idx <- t(utils::combn(n, 2))
    data.frame(
      program_a = programs[idx[, 1]],
      program_b = programs[idx[, 2]],
      phi = phi[idx],
      alignment = cos_f[idx]^2,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      program_a = integer(0), program_b = integer(0),
      phi = numeric(0), alignment = numeric(0)
    )
  }

  effective_rank <- NA_real_
  if (n > 0) {
    gram <- phi * cos_f
    gram[!is.finite(gram)] <- 0
    diag(gram) <- 1
    eigenvalues <- pmax(eigen(gram, symmetric = TRUE, only.values = TRUE)$values, 0)
    p <- eigenvalues[eigenvalues > 0] / sum(eigenvalues)
    effective_rank <- exp(-sum(p * log(p)))
  }
  safe_max <- function(x) {
    return(if (any(is.finite(x))) max(x, na.rm = TRUE) else NA_real_)
  }

  summary <- data.frame(
    n_programs = n,
    max_abs_phi = safe_max(abs(pairs$phi)),
    n_pairs_phi_0.7 = sum(abs(pairs$phi) >= 0.7, na.rm = TRUE),
    max_alignment = safe_max(pairs$alignment),
    effective_rank = effective_rank,
    stringsAsFactors = FALSE
  )
  return(list(pairs = pairs, summary = summary))
}


# Confidence tier from the similarity test and connectedness:
#   low          = fails the similarity test (similarity_q >=
#                  similarity_q_threshold, or missing), whatever the traits
#   single_trait = passes it, but effective_n_traits < min_effective_traits
#                  (skipped when effective_n_traits is NULL or missing)
#   high         = similarity_q < high_similarity_q &
#                  (connectedness >= high_min_connectedness |
#                   similarity_z >= high_min_similarity_z)
#   medium       = similarity_q < similarity_q_threshold (not high)
# Returned as an ordered factor, levels c("low", "single_trait", "medium",
# "high").
.program_confidence_tier <- function(similarity_q,
                                     connectedness,
                                     similarity_q_threshold = 0.05,
                                     high_similarity_q = 0.01,
                                     high_min_connectedness = 0.4,
                                     effective_n_traits = NULL,
                                     min_effective_traits = 1.5,
                                     similarity_z = NULL,
                                     high_min_similarity_z = Inf) {
  medium <- is.finite(similarity_q) & similarity_q < similarity_q_threshold
  connected <- is.finite(connectedness) & connectedness >= high_min_connectedness
  far_above_random <- if (is.null(similarity_z)) {
    rep(FALSE, length(similarity_q))
  } else {
    is.finite(similarity_z) & similarity_z >= high_min_similarity_z
  }
  high <- is.finite(similarity_q) & similarity_q < high_similarity_q &
    (connected | far_above_random)
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


# Similarity among a program's member SNPs: mean pairwise S, and connectedness
# (share of pairs with S >= edge_threshold).
.program_internal_coherence <- function(s_matrix, assigned_snps, edge_threshold) {
  snps <- intersect(assigned_snps, colnames(s_matrix))
  n <- length(snps)
  if (n < 2) {
    return(data.frame(
      n_snps = n,
      mean_internal_similarity = NA_real_,
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
    connectedness = if (n_pairs > 0) mean(internal_vals >= edge_threshold) else NA_real_,
    n_internal_pairs = n_pairs,
    stringsAsFactors = FALSE
  ))
}


# Size-matched random-set test for each program's mean internal similarity.
#
# `members_by_program` is a named list (names = program ids) of member SNP ids.
# For a program of n SNPs, each of `n_perm` draws takes n random SNPs (from
# those with any non-zero similarity to another SNP), so the baseline is
# matched on size and drawn from the trait's own similarity graph. Programs of
# the same size share one baseline. Returns one row per program:
# similarity_z, similarity_emp_p and the BH-adjusted similarity_q.
.program_similarity_null <- function(s_matrix,
                                     members_by_program,
                                     n_perm = 10000L,
                                     seed = 1) {
  columns <- c("similarity_z", "similarity_emp_p")
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
  mean_similarity <- function(ids) {
    sub <- S[ids, ids, drop = FALSE]
    return(mean(sub[upper.tri(sub)]))
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
        function(i) mean_similarity(sample(pool, n)),
        numeric(1)
      )
    }
    null <- null_by_size[[key]]
    obs <- mean_similarity(snps)
    null_mean <- mean(null)
    null_sd <- stats::sd(null)
    return(data.frame(
      program = as.integer(pg),
      similarity_z = if (is.finite(null_sd) && null_sd > 0) (obs - null_mean) / null_sd else NA_real_,
      similarity_emp_p = (1 + sum(null >= obs)) / (1 + n_perm),
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
