#' @title Stability of EBMF Programs Across Trait Subsamples
#' @description Refit the EBMF pipeline on trait subsamples of the observed
#' matrix and score how consistently each program's filtered member set
#' reappears. For each replicate, a random fraction of trait rows is held out,
#' the fit is repeated with the stored pipeline parameters, and each factor's
#' SNPs are filtered using the same lFSR and magnitude thresholds. A reference
#' program's replication score is the mean across replicates of its best
#' Jaccard overlap with any replicate factor set, both sets built by the same
#' top-loading rule. Jaccard is symmetric, so a large replicate factor that
#' merely contains a small program's few top SNPs does not count as a
#' replication.
#' Descriptive: scores are not p-values. Use them to prioritise programmes whose
#' membership survives resampling.
#'
#' The replicate refits are independent, so they can be run in parallel with
#' `cores > 1` (fork-based `parallel::mclapply()`; serial on platforms without
#' fork support). Each replicate sets its own seed, so results are reproducible
#' regardless of `cores`. flashier itself is single-threaded, so the speedup
#' scales with the number of worker processes.
#' @param clustering_result Result of `run_univariate_clustering()`.
#' @param n_rep Number of subsample replicates.
#' @param frac_traits Fraction of trait rows sampled per replicate.
#' @param top_n Size of each program's member set used for matching.
#' @param seed RNG seed.
#' @param cores Number of cores for parallel replicate refits. Defaults to `1`
#'   (serial).
#' @param verbose Print progress messages.
#' @return A dataframe with columns `program`, `n_ref` (reference member-set
#'   size), `replication` (mean best Jaccard overlap), `sd_replication`.
#' @keywords internal
stability_ebmf_programs <- function(clustering_result,
                                    n_rep = 10,
                                    frac_traits = 0.8,
                                    top_n = 20,
                                    seed = 1,
                                    cores = 1,
                                    verbose = TRUE) {
  .assert_ebmf_result(clustering_result)
  params <- clustering_result$parameters
  posterior <- ebmf_posterior_table(clustering_result)
  if (nrow(posterior) == 0) {
    return(data.frame(program = integer(0), n_ref = integer(0),
                      replication = numeric(0), sd_replication = numeric(0)))
  }

  lfsr_threshold <- params$ebmf_lfsr_threshold
  magnitude_threshold <- params$ebmf_magnitude_threshold
  keep_membership <- function(lfsr, abs_loading) {
    lfsr_ok <- is.na(lfsr_threshold) |
      (!is.na(lfsr) & lfsr < lfsr_threshold)
    support_ok <- is.finite(abs_loading) & abs_loading > 0
    magnitude_ok <- is.na(magnitude_threshold) |
      (support_ok & abs_loading > magnitude_threshold)
    lfsr_ok & support_ok & magnitude_ok
  }

  ref_sets <- lapply(sort(unique(posterior$program)), function(pg) {
    rows <- posterior$program == pg &
      keep_membership(posterior$lfsr, posterior$abs_loading)
    s <- posterior$snp_id[rows]
    loading <- posterior$abs_loading[rows]
    s[order(-loading)][seq_len(.stability_ref_size(top_n, length(s)))]
  })
  names(ref_sets) <- sort(unique(posterior$program))

  ebmf_input <- clustering_result$x_star

  replicate_sets <- function(i) {
    set.seed(seed + i)
    keep_rows <- sample(nrow(ebmf_input),
                        max(2L, floor(frac_traits * nrow(ebmf_input))))
    input <- ebmf_input[keep_rows, , drop = FALSE]
    fit <- .cluster_snp_profiles_ebmf(
      input,
      greedy_Kmax = params$ebmf_greedy_Kmax,
      lfsr_threshold = params$ebmf_lfsr_threshold,
      magnitude_threshold = params$ebmf_magnitude_threshold,
      # Must match the main fit. With drop_global = FALSE the refits keep a
      # global factor the reference fit had removed, so programs resembling it
      # lose their top-loading SNPs to that factor and their replication is
      # depressed for a structural reason rather than a stability one.
      drop_global = isTRUE(params$ebmf_drop_global),
      prior = params$ebmf_prior,
      backfit = params$ebmf_backfit,
      missing_se = .ebmf_missing_se(params)
    )$details$flash_fit
    if (is.null(fit) || fit$n_factors == 0) {
      return(list())
    }
    lapply(seq_len(fit$n_factors), function(k) {
      abs_loading <- abs(fit$F_pm[, k])
      lfsr <- if (!is.null(fit$F_lfsr)) {
        fit$F_lfsr[, k]
      } else {
        rep(NA_real_, length(abs_loading))
      }
      keep <- keep_membership(lfsr, abs_loading)
      ids <- rownames(fit$F_pm)
      if (is.null(ids)) {
        ids <- as.character(seq_len(nrow(fit$F_pm)))
      }
      ids <- ids[keep]
      loading <- abs_loading[keep]
      ids[order(-loading)][seq_len(.stability_ref_size(top_n, length(ids)))]
    })
  }

  run_reps <- function(i) {
    if (verbose) {
      message("Stability replicate ", i, "/", n_rep)
    }
    tryCatch(replicate_sets(i), error = function(e) list())
  }

  n_cores <- max(1L, as.integer(cores))
  if (n_cores > n_rep) {
    n_cores <- n_rep
  }
  if (n_cores > 1) {
    replicate_sets_all <- parallel::mclapply(
      seq_len(n_rep),
      run_reps,
      mc.cores = n_cores
    )
  } else {
    replicate_sets_all <- lapply(seq_len(n_rep), run_reps)
  }

  scores <- do.call(rbind, lapply(names(ref_sets), function(pg) {
    ref <- ref_sets[[pg]]
    per_rep <- vapply(replicate_sets_all, function(sets) {
      return(.stability_best_overlap(ref, sets))
    }, numeric(1))
    data.frame(
      program = as.integer(pg),
      n_ref = length(ref),
      replication = mean(per_rep),
      sd_replication = stats::sd(per_rep)
    )
  }))

  return(scores[order(-scores$replication), , drop = FALSE])
}


# Best Jaccard overlap between a reference SNP set and any replicate factor set
# (0 when there is nothing to compare). Jaccard, not |ref n s| / |ref|: with the
# one-sided score a small program was "replicated" by any large replicate factor
# that happened to contain its few top SNPs.
.stability_best_overlap <- function(ref, sets) {
  if (length(sets) == 0 || length(ref) == 0) {
    return(0)
  }
  return(max(vapply(sets, function(s) {
    return(length(intersect(ref, s)) / length(union(ref, s)))
  }, numeric(1))))
}


# Size of the top-loading reference set used for stability matching. A flat
# top_n measures programs on different fractions of themselves -- a 7-SNP
# program is matched on all 7 of itself while a 265-SNP program is matched on
# its top 20 -- which makes replication track program size rather than
# stability. Capping at half the program's SNPs keeps the fraction comparable.
.stability_ref_size <- function(top_n, n_available) {
  if (n_available <= 0) {
    return(0L)
  }
  return(as.integer(max(1L, min(as.integer(top_n), ceiling(n_available / 2)))))
}


# Refits must give unobserved cells the main fit's SE. Results saved before
# ebmf_missing_se existed were fit with SE 1.
.ebmf_missing_se <- function(params) {
  return(if (is.null(params$ebmf_missing_se)) 1 else params$ebmf_missing_se)
}
