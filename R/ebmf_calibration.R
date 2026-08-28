#' @title Extract Ungated EBMF Posterior Table
#' @description Return the full SNP x program posterior from an EBMF clustering
#' result without applying any membership gates. Every cell of the factor
#' loading matrix becomes one row, so downstream reporting can apply its own
#' calibrated thresholds (see [calibrate_ebmf_programs()]) instead of the
#' pipeline's fixed lFSR / magnitude gates.
#' @param clustering_result Result of `run_univariate_clustering()` with
#'   `cluster_type = "ebmf"`.
#' @return A dataframe with one row per SNP x program:
#'   \itemize{
#'     \item snp_id, program
#'     \item loading, abs_loading: posterior mean factor loading (F_pm)
#'     \item lfsr: local false sign rate where available
#'   }
#' @export
ebmf_posterior_table <- function(clustering_result) {
  if (!is.null(clustering_result$parameters) &&
      !identical(clustering_result$parameters$cluster_type, "ebmf")) {
    stop("clustering_result must come from cluster_type = 'ebmf'")
  }
  fit <- clustering_result$cluster_details$flash_fit
  if (is.null(fit) || fit$n_factors == 0) {
    return(data.frame(
      snp_id = character(0), program = integer(0),
      loading = numeric(0), abs_loading = numeric(0), lfsr = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  F_pm <- fit$F_pm
  F_lfsr <- fit$F_lfsr
  snp_ids <- rownames(F_pm)
  out <- do.call(rbind, lapply(seq_len(fit$n_factors), function(k) {
    data.frame(
      snp_id = snp_ids,
      program = k,
      loading = F_pm[, k],
      abs_loading = abs(F_pm[, k]),
      lfsr = if (!is.null(F_lfsr)) F_lfsr[, k] else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  return(out[order(out$program, -out$abs_loading), , drop = FALSE])
}


.ebmf_factor_signal <- function(fit) {
  if (is.null(fit) || fit$n_factors == 0) {
    return(numeric(0))
  }
  L_mat <- fit$L_pm
  F_mat <- fit$F_pm
  if (is.null(L_mat) || is.null(F_mat)) {
    return(numeric(0))
  }
  raw <- colSums(abs(L_mat), na.rm = TRUE) * colSums(abs(F_mat), na.rm = TRUE)
  return(unname(as.numeric(raw)))
}


#' @title Calibrate EBMF Programs Against Permutation Nulls
#' @description Post-hoc, permutation-calibrated reporting layer for EBMF
#' results. Generates `n_null` structureless versions of the *observed* matrix
#' by independently shuffling each trait row across SNPs — preserving every
#' trait's exact sparsity, effect distribution, and noise regime — refits the
#' EBMF pipeline on each permuted matrix, and reports, for every fitted program:
#'
#' \describe{
#'   \item{Factor strength}{The proportion of the observed matrix signal
#'   attributable to the factor, computed from its fitted reconstruction:
#'   \eqn{\sum_i |L_{ik}| \sum_j |F_{jk}|} (raw factor signal) over the
#'   observed matrix signal \eqn{\sum |X|}.}
#'   \item{SNP membership (emp_p / q)}{Each SNP x program loading receives an
#'   empirical p-value against the pooled null loading distribution, BH-corrected
#'   into q-values; cells with `q <= alpha_membership` are labelled `core`. This
#'   is a membership gate only — it no longer drives a program-level
#'   detect/certify decision (see `summarise_ebmf_programs()` for the
#'   program-level filters).}
#' }
#'
#' No thresholds are applied inside the pipeline; this function is the
#' reporting layer. The Louvain path is unaffected.
#' @param clustering_result Result of `run_univariate_clustering()` with
#'   `cluster_type = "ebmf"`. Its recorded parameters are reused verbatim for
#'   every null replicate.
#' @param n_null Number of permutation replicates. More replicates sharpen the
#'   empirical quantiles; 20 is a reasonable minimum.
#' @param alpha_membership Membership FDR level for `core` cells.
#' @param n_candidate_tier Per program, how many non-core cells to label as the
#'   uncertified "candidate" tier (top-loaded first).
#' @param seed RNG seed for the permutations.
#' @param verbose Print progress messages.
#' @return A list with:
#'   \itemize{
#'     \item programs: dataframe (program, loading_mass, n_candidate,
#'       raw_factor_signal, factor_strength)
#'     \item memberships: dataframe (snp_id, program, loading, abs_loading,
#'       lfsr, emp_p, q, core, tier)
#'     \item null_summary: list with per-replicate max masses and pooled loading
#'       quantiles
#'     \item settings: calibration settings used
#'   }
#' @export
calibrate_ebmf_programs <- function(clustering_result,
                                    n_null = 20,
                                    alpha_membership = 0.05,
                                    n_candidate_tier = 25L,
                                    seed = 1,
                                    verbose = TRUE) {
  params <- clustering_result$parameters
  if (is.null(params) || !identical(params$cluster_type, "ebmf")) {
    stop("clustering_result must come from cluster_type = 'ebmf'")
  }

  fit_null <- function(i) {
    set.seed(seed + i)
    X <- clustering_result$x_matrix
    perm <- lapply(seq_len(nrow(X)), function(j) sample(ncol(X)))
    if (params$ebmf_se_mode == "matrix" &&
        !is.null(clustering_result$beta_matrix)) {
      beta <- clustering_result$beta_matrix
      se <- clustering_result$se_matrix
      beta_perm <- beta
      se_perm <- se
      for (j in seq_len(nrow(beta))) {
        idx <- perm[[j]]
        beta_perm[j, ] <- beta[j, idx]
        se_perm[j, ] <- se[j, idx]
      }
      signs <- orient_pleiotropy_matrix(
        X, target_trait_id = params$target_trait_id
      )$target_signs
      input <- sweep(beta_perm, 2, signs, `*`)
      se_input <- se_perm
      if (identical(params$ebmf_beta_scale, "trait")) {
        rs <- sqrt(rowMeans(input^2, na.rm = TRUE))
        rs[!is.finite(rs) | rs <= 0] <- 1
        input <- sweep(input, 1, rs, `/`)
        se_input <- sweep(se_input, 1, rs, `/`)
      }
      .cluster_snp_profiles_ebmf(
        input,
        greedy_Kmax = params$ebmf_greedy_Kmax,
        lfsr_threshold = params$ebmf_lfsr_threshold,
        magnitude_threshold = params$ebmf_magnitude_threshold,
        drop_global = FALSE,
        prior = params$ebmf_prior,
        backfit = params$ebmf_backfit,
        observed_se_matrix = se_input
      )
    } else {
      X_perm <- X
      for (j in seq_len(nrow(X))) {
        X_perm[j, ] <- X[j, perm[[j]]]
      }
      oriented <- orient_pleiotropy_matrix(
        X_perm, target_trait_id = params$target_trait_id
      )
      X_star <- compress_effect_matrix(
        oriented$x_matrix,
        method = params$compress_method,
        asinh_scale = params$compress_scale
      )
      .cluster_snp_profiles_ebmf(
        X_star,
        greedy_Kmax = params$ebmf_greedy_Kmax,
        lfsr_threshold = params$ebmf_lfsr_threshold,
        magnitude_threshold = params$ebmf_magnitude_threshold,
        drop_global = FALSE,
        prior = params$ebmf_prior,
        backfit = params$ebmf_backfit
      )
    }
  }

  factor_masses <- function(fit) {
    if (is.null(fit) || fit$n_factors == 0) {
      return(numeric(0))
    }
    colSums(abs(fit$F_pm), na.rm = TRUE)
  }

  null_max_mass <- numeric(n_null)
  null_cells <- vector("list", n_null)
  for (i in seq_len(n_null)) {
    if (verbose) {
      message("Null replicate ", i, "/", n_null)
    }
    fit <- tryCatch(fit_null(i)$details$flash_fit, error = function(e) NULL)
    masses <- factor_masses(fit)
    null_max_mass[i] <- if (length(masses)) max(masses) else 0
    if (!is.null(fit) && fit$n_factors > 0) {
      null_cells[[i]] <- as.numeric(abs(fit$F_pm))
    }
  }
  null_cells <- unlist(null_cells)
  null_cells <- null_cells[is.finite(null_cells)]

  obs_fit <- clustering_result$cluster_details$flash_fit
  obs_signal <- .ebmf_factor_signal(obs_fit)
  obs_input <- clustering_result$ebmf_input
  if (is.null(obs_input)) {
    obs_input <- clustering_result$x_star
  }
  total_signal <- sum(abs(obs_input), na.rm = TRUE)

  posterior <- ebmf_posterior_table(clustering_result)
  min_core <- max(
    as.integer(if (is.null(params$min_module_size)) 3L else params$min_module_size),
    2L
  )
  if (nrow(posterior) == 0) {
    return(list(
      programs = data.frame(program = integer(0), loading_mass = numeric(0),
                            n_candidate = integer(0),
                            raw_factor_signal = numeric(0),
                            factor_strength = numeric(0),
                            stringsAsFactors = FALSE),
      memberships = cbind(posterior, emp_p = numeric(0), q = numeric(0),
                          core = logical(0), tier = character(0)),
      null_summary = list(
        max_masses = null_max_mass,
        loading_quantiles = stats::quantile(null_cells,
                                            c(.5, .9, .95, .99))
      ),
      settings = list(n_null = n_null, alpha_membership = alpha_membership,
                      seed = seed,
                      min_core_members = min_core,
                      n_candidate_tier = n_candidate_tier)
    ))
  }

  programs <- aggregate(abs_loading ~ program, data = posterior, FUN = sum)
  names(programs)[names(programs) == "abs_loading"] <- "loading_mass"
  programs$program <- as.integer(as.character(programs$program))
  programs$raw_factor_signal <- obs_signal[programs$program]
  programs$factor_strength <- if (is.finite(total_signal) && total_signal > 0) {
    programs$raw_factor_signal / total_signal
  } else {
    NA_real_
  }

  if (length(null_cells) == 0) {
    warning(
      "Permutation nulls produced no fitted factors; membership cannot be ",
      "calibrated (emp_p/q set to NA). Increase n_null.",
      call. = FALSE
    )
    posterior$emp_p <- NA_real_
    posterior$q <- NA_real_
    posterior$core <- FALSE
  } else {
    posterior$emp_p <- vapply(posterior$abs_loading, function(x) {
      (1 + sum(null_cells >= x)) / (1 + length(null_cells))
    }, numeric(1))
    posterior$q <- stats::p.adjust(posterior$emp_p, method = "BH")
    posterior$core <- posterior$q <= alpha_membership
  }

  # Two-tier reporting: calibrated core plus an explicitly uncertified
  # candidate tier. Candidate cells are retained for every fitted program so
  # a real-data report can show plausible biology even when no program passes
  # the stricter membership gate.
  posterior$tier <- NA_character_
  posterior$tier[posterior$core] <- "core"
  n_candidate <- integer(nrow(programs))
  for (i in seq_len(nrow(programs))) {
    pg <- programs$program[i]
    idx <- which(posterior$program == pg & is.na(posterior$tier))
    top <- order(-posterior$abs_loading[idx])[seq_len(min(n_candidate_tier,
                                                          length(idx)))]
    posterior$tier[idx[top]] <- "candidate"
    n_candidate[i] <- length(top)
  }
  programs$n_candidate <- n_candidate

  return(list(
    programs = programs[order(-programs$factor_strength, -programs$n_candidate,
                              na.last = TRUE), ],
    memberships = posterior[, c("snp_id", "program", "loading", "abs_loading",
                                "lfsr", "emp_p", "q", "core", "tier")],
    null_summary = list(
      max_masses = null_max_mass,
      loading_quantiles = stats::quantile(null_cells, c(.5, .9, .95, .99)),
      n_pooled_cells = length(null_cells)
    ),
    settings = list(n_null = n_null, alpha_membership = alpha_membership,
                    seed = seed, min_core_members = min_core,
                    n_candidate_tier = n_candidate_tier)
  ))
}


#' @title Stability of EBMF Programs Across Trait Subsamples
#' @description Refit the EBMF pipeline on trait subsamples of the observed
#' matrix and score how consistently each program's filtered member set
#' reappears. For each replicate, a random fraction of trait rows is held out,
#' the fit is repeated with the stored pipeline parameters, and each factor's
#' SNPs are filtered using the same lFSR and magnitude thresholds. A reference
#' program's replication score is the mean across replicates of its best overlap
#' with any replicate factor set.
#' Descriptive: scores are not p-values. Use them to prioritise programmes whose
#' membership survives resampling.
#'
#' The replicate refits are independent, so they can be run in parallel with
#' `cores > 1` (fork-based `parallel::mclapply()`; serial on platforms without
#' fork support). Each replicate sets its own seed, so results are reproducible
#' regardless of `cores`. flashier itself is single-threaded, so the speedup
#' scales with the number of worker processes.
#' @param clustering_result Result of `run_univariate_clustering()` with
#'   `cluster_type = "ebmf"`.
#' @param n_rep Number of subsample replicates.
#' @param frac_traits Fraction of trait rows sampled per replicate.
#' @param top_n Size of each program's member set used for matching.
#' @param seed RNG seed.
#' @param cores Number of cores for parallel replicate refits. Defaults to `1`
#'   (serial).
#' @param verbose Print progress messages.
#' @return A dataframe with columns `program`, `n_ref` (reference member-set
#'   size), `replication` (mean best-overlap fraction), `sd_replication`.
#' @export
stability_ebmf_programs <- function(clustering_result,
                                    n_rep = 10,
                                    frac_traits = 0.8,
                                    top_n = 20,
                                    seed = 1,
                                    cores = 1,
                                    verbose = TRUE) {
  params <- clustering_result$parameters
  if (is.null(params) || !identical(params$cluster_type, "ebmf")) {
    stop("clustering_result must come from cluster_type = 'ebmf'")
  }
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
    s[order(-loading)][seq_len(min(top_n, length(s)))]
  })
  names(ref_sets) <- sort(unique(posterior$program))

  ebmf_input <- clustering_result$x_star
  se_input <- NULL
  if (identical(params$ebmf_se_mode, "matrix") &&
      !is.null(clustering_result$beta_matrix)) {
    signs <- orient_pleiotropy_matrix(
      clustering_result$x_matrix,
      target_trait_id = params$target_trait_id
    )$target_signs
    beta_o <- sweep(clustering_result$beta_matrix[
      rownames(clustering_result$x_star), , drop = FALSE
    ], 2, signs[colnames(clustering_result$x_star)], `*`)
    ebmf_input <- beta_o
    se_input <- clustering_result$se_matrix[
      rownames(beta_o), , drop = FALSE
    ]
    if (identical(params$ebmf_beta_scale, "trait")) {
      rs <- sqrt(rowMeans(ebmf_input^2, na.rm = TRUE))
      rs[!is.finite(rs) | rs <= 0] <- 1
      ebmf_input <- sweep(ebmf_input, 1, rs, `/`)
      se_input <- sweep(se_input, 1, rs, `/`)
    }
  }

  replicate_sets <- function(i) {
    set.seed(seed + i)
    keep_rows <- sample(nrow(ebmf_input),
                        max(2L, floor(frac_traits * nrow(ebmf_input))))
    input <- ebmf_input[keep_rows, , drop = FALSE]
    se_sub <- if (!is.null(se_input)) se_input[keep_rows, , drop = FALSE] else NULL
    fit <- .cluster_snp_profiles_ebmf(
      input,
      greedy_Kmax = params$ebmf_greedy_Kmax,
      lfsr_threshold = params$ebmf_lfsr_threshold,
      magnitude_threshold = params$ebmf_magnitude_threshold,
      drop_global = FALSE,
      prior = params$ebmf_prior,
      backfit = params$ebmf_backfit,
      observed_se_matrix = se_sub
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
      ids[order(-loading)][seq_len(min(top_n, length(ids)))]
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
      if (length(sets) == 0 || length(ref) == 0) {
        return(0)
      }
      max(vapply(sets, function(s) {
        length(intersect(ref, s)) / length(ref)
      }, numeric(1)))
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
