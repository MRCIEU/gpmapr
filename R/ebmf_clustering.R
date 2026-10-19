# Fit EBMF to a features x SNPs matrix and turn the factors into programs:
# drop at most one global mega-factor, then gate SNP membership on lFSR and
# loading magnitude. Returns the fit details run_univariate_clustering()
# stores as cluster_details (on a failed fit, the error and no programs).
.cluster_snp_profiles_ebmf <- function(x_matrix,
                                       greedy_Kmax = 50L,
                                       lfsr_threshold = 0.05,
                                       magnitude_threshold = 0.25,
                                       drop_global = TRUE,
                                       prior = "point_normal",
                                       backfit = TRUE,
                                       missing_se = 1) {
  fit_error <- NULL
  flash_fit <- tryCatch(
    run_ebmf(
      beta_matrix = x_matrix,
      missing_se = missing_se,
      greedy_Kmax = greedy_Kmax,
      backfit = backfit,
      ebnm_fn = .resolve_ebnm_fn(prior),
      verbose = 0L
    ),
    error = function(e) {
      fit_error <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(flash_fit)) {
    warning("EBMF fit failed: ", fit_error, call. = FALSE)
    return(list(
      method = "ebmf",
      details = list(
        flash_fit = NULL,
        membership = NULL,
        n_programs = 0L,
        n_multi_program = 0L,
        dropped_global_factors = integer(0),
        error = fit_error
      )
    ))
  }

  # Recorded before any global factor is removed, so the cap check below asks
  # what the greedy search itself produced.
  n_greedy_factors <- flash_fit$n_factors

  dropped_global <- integer(0)
  if (drop_global && flash_fit$n_factors > 0) {
    preliminary <- extract_ebmf_clusters(
      flash_fit,
      lfsr_threshold = lfsr_threshold,
      magnitude_threshold = magnitude_threshold
    )
    dropped_global <- identify_ebmf_global_factors(
      flash_fit,
      membership = preliminary$membership
    )
    if (length(dropped_global) > 0) {
      flash_fit <- remove_ebmf_factors(flash_fit, kset = dropped_global)
    }
  }

  extracted <- extract_ebmf_clusters(
    flash_fit,
    lfsr_threshold = lfsr_threshold,
    magnitude_threshold = magnitude_threshold
  )

  membership <- extracted$membership

  return(list(
    method = "ebmf",
    details = list(
      flash_fit = flash_fit,
      membership = membership,
      n_programs = extracted$n_programs,
      n_multi_program = extracted$n_multi_program,
      dropped_global_factors = dropped_global,
      greedy_Kmax = as.integer(greedy_Kmax),
      # TRUE when the greedy search stopped because it hit the cap rather than
      # because it ran out of factors worth adding -- i.e. K was set by
      # greedy_Kmax, not by the data. Raise the cap and refit if so.
      greedy_Kmax_reached = isTRUE(
        (n_greedy_factors + length(dropped_global)) >= as.integer(greedy_Kmax)
      )
    )
  ))
}


#' @title Run EBMF Factorization
#' @description Fit an Empirical Bayes Matrix Factorization model using
#' \pkg{flashier}. The model decomposes the features x SNPs z-score matrix
#' \eqn{Y = LF' + E} into latent factors (\eqn{F}, SNP programs) and
#' loadings (\eqn{L}, feature contributions), with sparse priors estimated
#' from the data.
#'
#' Observed cells are z-scores, so each enters with standard error 1. An
#' unobserved (`NA`) cell, meaning no colocalising association, is fit as an
#' observed 0 with standard error \code{missing_se}, so the pattern of observed
#' cells informs the fit. Because \code{S} carries the noise, no residual
#' variance is estimated.
#'
#' @param beta_matrix Features x SNPs matrix of signed z-scores, `NA` where
#'   unobserved, typically \code{run_univariate_clustering()$ebmf_input}. Rows
#'   and columns that are entirely missing or zero are removed before fitting.
#' @param missing_se Standard error of an unobserved cell, which is fit as an
#'   observed 0. \code{1} (default) treats it as confidently zero; a larger
#'   value (e.g. \code{2.5}) says only that the effect is small, with unknown
#'   sign.
#' @param greedy_Kmax Maximum number of factors to add greedily. Factors are
#'   only added while they improve the variational lower bound. Defaults to 50.
#' @param backfit Logical; if \code{TRUE} (default), all factors are cyclically
#'   updated after the greedy phase until convergence.
#' @param ebnm_fn EBNM function for priors on loadings and factors. Defaults to
#'   \code{ebnm::ebnm_point_normal} (point-normal prior). Pass a list of two
#'   functions to use different priors for loadings (features) and factors (SNPs).
#' @param verbose Verbosity: 0 = silent, 1 = summary, 2 = ELBO updates,
#'   3 = per-iteration.
#' @return A \code{flash} object from \pkg{flashier}. Key elements:
#'   \itemize{
#'     \item n_factors: number of discovered factors (programs)
#'     \item L_pm, L_lfsr: posterior means and lFSR for feature loadings (n x K)
#'     \item F_pm, F_lfsr: posterior means and lFSR for SNP factors (p x K)
#'     \item pve: proportion of variance explained per factor
#'     \item elbo: variational lower bound
#'   }
#' @keywords internal
run_ebmf <- function(beta_matrix,
                     missing_se = 1,
                     greedy_Kmax = 50L,
                     backfit = TRUE,
                     ebnm_fn = NULL,
                     verbose = 1L) {
  if (!requireNamespace("flashier", quietly = TRUE)) {
    stop("Package 'flashier' is required for EBMF clustering", call. = FALSE)
  }
  if (!requireNamespace("ebnm", quietly = TRUE)) {
    stop("Package 'ebnm' is required for EBMF clustering", call. = FALSE)
  }
  if (!is.numeric(missing_se) || length(missing_se) != 1 ||
        !is.finite(missing_se) || missing_se <= 0) {
    stop("missing_se must be a single positive number")
  }
  if (!is.matrix(beta_matrix)) {
    stop("beta_matrix must be a matrix")
  }
  if (is.null(ebnm_fn)) {
    ebnm_fn <- ebnm::ebnm_point_normal
  }

  prepared <- .prepare_ebmf_flash_inputs(beta_matrix, missing_se = missing_se)

  fit <- flashier::flash(
    data = prepared$data,
    S = prepared$S,
    ebnm_fn = ebnm_fn,
    var_type = NULL,
    greedy_Kmax = greedy_Kmax,
    backfit = backfit,
    nullcheck = TRUE,
    verbose = verbose
  )

  return(fit)
}


#' @title Identify Global EBMF Factors To Drop
#' @description Suggest factor indices that look like a shared/global axis,
#' but only when the fit has non-trivial total PVE. Factor 1 is not assumed
#' global. Returns an empty integer vector for near-null fits. By default at
#' most one factor is suggested (highest PVE share among candidates).
#' @param flash_fit A \code{flash} object from \code{run_ebmf()}.
#' @param membership Optional SNP x program logical membership matrix from
#'   \code{extract_ebmf_clusters()$membership}.
#' @param min_total_pve Require this total PVE before dropping anything.
#'   Defaults to \code{0.01}.
#' @param min_pve_share Candidate if the factor explains at least this share of
#'   total PVE, regardless of how many SNPs it soft-assigns. Defaults to
#'   \code{0.7}.
#' @param min_soft_frac Candidate if the factor soft-assigns at least this
#'   fraction of SNPs. Defaults to \code{0.25}.
#' @param max_drop Maximum number of factors to return. Defaults to 1.
#' @return Integer vector of 1-based factor indices to pass to
#'   \code{remove_ebmf_factors()}, possibly empty.
#' @keywords internal
identify_ebmf_global_factors <- function(flash_fit,
                                         membership = NULL,
                                         min_total_pve = 0.01,
                                         min_pve_share = 0.7,
                                         min_soft_frac = 0.25,
                                         max_drop = 1L) {
  if (!inherits(flash_fit, "flash")) {
    stop("flash_fit must be a flash object from run_ebmf()")
  }

  K <- flash_fit$n_factors
  if (K == 0 || max_drop < 1) {
    return(integer(0))
  }

  pve <- flash_fit$pve
  total_pve <- sum(pve, na.rm = TRUE)
  if (!is.finite(total_pve) || total_pve < min_total_pve) {
    return(integer(0))
  }

  share <- pve / total_pve
  if (is.null(membership)) {
    soft_frac <- rep(0, K)
  } else {
    if (!is.matrix(membership) || ncol(membership) != K) {
      stop("membership must be a logical matrix with one column per factor")
    }
    soft_frac <- as.numeric(colSums(membership) / nrow(membership))
  }

  # A mega-factor shows up either as very broad membership or as a factor that
  # carries a dominant share of the fitted signal. The PVE-share branch is
  # deliberately independent of soft_frac: the real BMI mega-factor held 40% of
  # the factor strength on only 27% of SNPs, so requiring broad membership as
  # well as high PVE let it through.
  candidates <- which(soft_frac >= min_soft_frac | share >= min_pve_share)
  if (length(candidates) == 0) {
    return(integer(0))
  }

  ord <- candidates[order(
    share[candidates],
    soft_frac[candidates],
    decreasing = TRUE
  )]
  return(as.integer(utils::head(ord, max_drop)))
}


#' @title Extract EBMF SNP Clusters
#' @description Extract overlapping SNP clusters from a fitted \code{flash}
#' object. For each discovered factor (program), SNPs are selected if their
#' local false sign rate (lFSR) is below \code{lfsr_threshold} \emph{and}
#' their absolute posterior mean factor loading exceeds
#' \code{magnitude_threshold}. SNPs may belong to multiple programs
#' (overlapping clusters).
#'
#' @param flash_fit A \code{flash} object from \code{run_ebmf()}.
#' @param lfsr_threshold lFSR threshold for including a SNP in a program.
#'   Defaults to 0.05.
#' @param magnitude_threshold Minimum absolute posterior mean factor loading.
#'   Defaults to 0.25.
#' @return A list with:
#'   \itemize{
#'     \item clusters: dataframe with columns \code{program}, \code{snp_id},
#'       \code{loading}, \code{lfsr}
#'     \item membership: logical matrix (SNPs x programs)
#'     \item n_programs: number of discovered programs
#'     \item n_assigned: number of SNPs assigned to at least one program
#'     \item n_multi_program: number of SNPs assigned to more than one program
#'     \item lfsr_threshold: threshold used
#'     \item magnitude_threshold: threshold used
#'   }
#' @keywords internal
extract_ebmf_clusters <- function(flash_fit,
                                  lfsr_threshold = 0.05,
                                  magnitude_threshold = 0.25) {
  if (!inherits(flash_fit, "flash")) {
    stop("flash_fit must be a flash object from run_ebmf()")
  }

  K <- flash_fit$n_factors
  if (K == 0) {
    return(list(
      clusters = data.frame(
        program = integer(0), snp_id = character(0),
        loading = numeric(0), lfsr = numeric(0),
        stringsAsFactors = FALSE
      ),
      membership = matrix(FALSE, nrow = 0, ncol = 0),
      n_programs = 0L,
      n_assigned = 0L,
      n_multi_program = 0L,
      lfsr_threshold = lfsr_threshold,
      magnitude_threshold = magnitude_threshold
    ))
  }

  F_pm <- flash_fit$F_pm
  F_lfsr <- flash_fit$F_lfsr

  snp_ids <- rownames(F_pm)
  if (is.null(snp_ids)) {
    snp_ids <- as.character(seq_len(nrow(F_pm)))
  }

  membership <- matrix(
    FALSE,
    nrow = nrow(F_pm),
    ncol = K,
    dimnames = list(snp_ids, paste0("program_", seq_len(K)))
  )

  cluster_rows <- vector("list", K)
  for (k in seq_len(K)) {
    if (is.na(lfsr_threshold)) {
      pass_lfsr <- TRUE
    } else {
      pass_lfsr <- !is.na(F_lfsr[, k]) & F_lfsr[, k] < lfsr_threshold
    }
    if (is.na(magnitude_threshold)) {
      pass_mag <- TRUE
    } else {
      pass_mag <- !is.na(F_pm[, k]) & abs(F_pm[, k]) > magnitude_threshold
    }
    selected <- pass_lfsr & pass_mag
    membership[, k] <- selected

    idx <- which(selected)
    if (length(idx) > 0) {
      cluster_rows[[k]] <- data.frame(
        program = k,
        snp_id = snp_ids[idx],
        loading = F_pm[idx, k],
        lfsr = F_lfsr[idx, k],
        stringsAsFactors = FALSE
      )
    }
  }

  clusters <- dplyr::bind_rows(cluster_rows)
  row_sums <- rowSums(membership)

  return(list(
    clusters = clusters,
    membership = membership,
    n_programs = K,
    n_assigned = sum(row_sums > 0),
    n_multi_program = sum(row_sums > 1),
    lfsr_threshold = lfsr_threshold,
    magnitude_threshold = magnitude_threshold
  ))
}


#' @title Remove EBMF Factors
#' @description Drop one or more factors from a fitted \code{flash} object.
#' Remaining factors are re-indexed. Prefer
#' \code{identify_ebmf_global_factors()} over assuming factor 1 is global.
#' An empty \code{kset} returns \code{flash_fit} unchanged.
#' @param flash_fit A \code{flash} object from \code{run_ebmf()}.
#' @param kset Integer vector of factor indices to remove.
#' @return A \code{flash} object with \code{kset} removed.
#' @keywords internal
remove_ebmf_factors <- function(flash_fit, kset = integer(0)) {
  if (!requireNamespace("flashier", quietly = TRUE)) {
    stop("Package 'flashier' is required to remove EBMF factors", call. = FALSE)
  }
  if (!inherits(flash_fit, "flash")) {
    stop("flash_fit must be a flash object from run_ebmf()")
  }
  if (length(kset) == 0) {
    return(flash_fit)
  }
  if (any(is.na(kset))) {
    stop("kset must not contain NA")
  }
  kset <- as.integer(kset)
  if (flash_fit$n_factors == 0) {
    return(flash_fit)
  }
  if (any(kset < 1L) || any(kset > flash_fit$n_factors)) {
    stop("kset must be within 1:n_factors")
  }
  return(flashier::flash_factors_remove(flash_fit, kset = kset))
}


#' @title Extract Ungated EBMF Posterior Table
#' @description Return the full SNP x program posterior from an EBMF clustering
#' result without applying any membership gates. Every cell of the factor
#' loading matrix becomes one row, so callers can apply the lFSR / magnitude
#' gate themselves or use the loadings directly.
#' @param clustering_result Result of `run_univariate_clustering()`.
#' @return A dataframe with one row per SNP x program:
#'   \itemize{
#'     \item snp_id, program
#'     \item loading, abs_loading: posterior mean factor loading (F_pm)
#'     \item lfsr: local false sign rate where available
#'   }
#' @keywords internal
ebmf_posterior_table <- function(clustering_result) {
  .assert_ebmf_result(clustering_result)
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


.ebmf_slice_has_signal <- function(x) {
  return(any(!is.na(x) & x != 0))
}


# flashier inputs from a features x SNPs z-score matrix: rows and columns with
# no observed non-zero value are removed (they carry nothing to factorise), and
# every unobserved cell enters as a 0 with standard error missing_se, observed
# cells with standard error 1. An unobserved cell means no colocalising
# association, so it is evidence of a small effect; leaving it out (NA) would
# discard the pattern of observed cells that, on a matrix this sparse, is most
# of what EBMF has to separate programs by.
.prepare_ebmf_flash_inputs <- function(beta_matrix, missing_se = 1) {
  if (!is.matrix(beta_matrix)) {
    stop("beta_matrix must be a matrix")
  }
  row_keep <- apply(beta_matrix, 1, .ebmf_slice_has_signal)
  col_keep <- apply(beta_matrix, 2, .ebmf_slice_has_signal)
  removed_rows <- sum(!row_keep)
  removed_cols <- sum(!col_keep)
  if (removed_rows > 0 || removed_cols > 0) {
    warning(
      "Removed ", removed_rows, " all-zero/all-missing row(s) and ",
      removed_cols, " all-zero/all-missing column(s) for EBMF",
      call. = FALSE
    )
  }
  data <- beta_matrix[row_keep, col_keep, drop = FALSE]
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("No rows or columns remain after removing all-zero/all-missing slices")
  }

  missing <- is.na(data)
  data[missing] <- 0
  S <- 1
  if (any(missing) && missing_se != 1) {
    S <- matrix(1, nrow = nrow(data), ncol = ncol(data), dimnames = dimnames(data))
    S[missing] <- missing_se
  }
  return(list(data = data, S = S))
}


.resolve_ebnm_fn <- function(ebnm_fn) {
  if (is.function(ebnm_fn)) {
    return(ebnm_fn)
  }
  if (!is.character(ebnm_fn) || length(ebnm_fn) != 1L) {
    stop("ebnm_fn must be a function or character string")
  }

  resolved <- switch(
    ebnm_fn,
    point_normal = ebnm::ebnm_point_normal,
    point_laplace = ebnm::ebnm_point_laplace,
    stop("Unknown ebnm_fn: ", ebnm_fn, call. = FALSE)
  )

  return(resolved)
}
