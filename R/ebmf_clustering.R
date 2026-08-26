#' @title Build EBMF Feature Matrix
#' @description Construct a features x SNPs matrix for Empirical Bayes Matrix
#' Factorization. Cell values are signed z-scores (\eqn{\beta / \mathrm{SE}}).
#' Phenotypic traits and molecular QTL studies are retained as individual
#' study-level rows. Sparse features observed in fewer than
#' \code{min_snp_signals} SNPs are dropped. SNP columns are oriented by
#' \eqn{\mathrm{sign}(z_{\mathrm{target}})}. Matching standard errors for
#' flashier are 1 for observed cells and missing for unobserved cells (filled
#' with a large value inside \code{run_ebmf()}).
#' @param trait_id Numeric ID of the target trait whose associated SNPs define
#'   the columns.
#' @param p_threshold P-value threshold for including a target-trait SNP.
#'   Defaults to 5e-8.
#' @param snp_key Column used to name SNP columns: \code{"variant_id"},
#'   \code{"display_snp"}, or \code{"coloc_group_id"}. Defaults to
#'   \code{"variant_id"}.
#' @param min_snp_signals Minimum non-missing SNP count required to keep a
#'   phenotypic trait or gene row. Defaults to 5.
#' @param compress_method Soft compression applied after orientation via
#'   \code{compress_effect_matrix()}. Defaults to \code{"asinh"} so extreme /
#'   corrupt z-scores do not dominate flashier sum-of-squares without hard
#'   clipping. Use \code{"none"} to disable.
#' @return A list with:
#'   \itemize{
#'     \item beta_matrix: features x SNPs matrix of (optionally compressed)
#'       signed z-scores
#'     \item se_matrix: features x SNPs matrix of SEs for flashier (\code{S})
#'     \item trait_info: feature metadata (\code{trait_id}, \code{trait_name},
#'       \code{feature_type}, \code{trait_category}, \code{trait_label})
#'     \item verification_trait_matrix: SNPs x phenotypic traits (oriented,
#'       compressed), for post-clustering trait verification helpers
#'     \item verification_trait_info: metadata for verification matrix columns
#'     \item snp_info, z_target, coloc_groups, target_trait_id, prep_summary
#'   }
#' @export
build_ebmf_matrix <- function(trait_id,
                              p_threshold = NULL,
                              snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                              min_snp_signals = 5L,
                              compress_method = c("asinh", "none")) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }
  snp_key <- match.arg(snp_key)
  compress_method <- match.arg(compress_method)
  if (!is.numeric(min_snp_signals) || min_snp_signals < 1) {
    stop("min_snp_signals must be a positive number")
  }

  pert <- build_perturbation_matrices(
    trait_id = trait_id,
    p_threshold = p_threshold,
    snp_key = snp_key
  )
  filtered <- filter_perturbation_features(
    trait_matrix = pert$trait_matrix,
    gene_matrix = pert$gene_matrix,
    min_snps = min_snp_signals
  )
  oriented <- orient_perturbation_matrices(
    trait_matrix = filtered$trait_matrix,
    gene_matrix = filtered$gene_matrix,
    z_target = pert$z_target
  )

  compressed <- compress_perturbation_matrices(
    trait_matrix = oriented$trait_matrix,
    gene_matrix = oriented$gene_matrix,
    method = compress_method
  )
  trait_mat <- compressed$trait_matrix
  gene_mat <- compressed$gene_matrix

  trait_info <- .ebmf_feature_info_from_perturbation(
    feature_info = pert$trait_info,
    kept_ids = colnames(trait_mat),
    feature_type = "trait",
    coloc_groups = pert$coloc_groups
  )
  gene_info <- .ebmf_feature_info_from_perturbation(
    feature_info = pert$gene_info,
    kept_ids = colnames(gene_mat),
    feature_type = "gene",
    coloc_groups = pert$coloc_groups
  )

  y_traits <- t(trait_mat)
  y_genes <- t(gene_mat)
  if (nrow(y_genes) > 0) {
    rownames(y_genes) <- paste0("gene:", rownames(y_genes))
  }
  if (nrow(gene_info) > 0) {
    gene_info$trait_id <- paste0("gene:", gene_info$trait_id)
  }

  beta_matrix <- rbind(y_traits, y_genes)
  se_matrix <- matrix(
    1,
    nrow = nrow(beta_matrix),
    ncol = ncol(beta_matrix),
    dimnames = dimnames(beta_matrix)
  )
  se_matrix[is.na(beta_matrix)] <- NA_real_

  cleaned <- .filter_ebmf_matrices(
    beta_matrix = beta_matrix,
    se_matrix = se_matrix
  )

  trait_info <- dplyr::bind_rows(trait_info, gene_info) |>
    dplyr::filter(trait_id %in% rownames(cleaned$beta_matrix))

  snp_info <- pert$snp_info |>
    dplyr::filter(as.character(snp_id) %in% colnames(cleaned$beta_matrix))

  z_target <- pert$z_target[colnames(cleaned$beta_matrix)]

  prep_summary <- data.frame(
    n_features = nrow(cleaned$beta_matrix),
    n_traits = sum(trait_info$feature_type == "trait"),
    n_genes = sum(trait_info$feature_type == "gene"),
    n_snps = ncol(cleaned$beta_matrix),
    n_traits_dropped = filtered$n_traits_dropped,
    n_genes_dropped = filtered$n_genes_dropped,
    min_snp_signals = as.integer(min_snp_signals),
    compress_method = compress_method,
    max_abs_raw = max(
      c(abs(oriented$trait_matrix), abs(oriented$gene_matrix)),
      na.rm = TRUE
    ),
    max_abs_used = max(c(abs(trait_mat), abs(gene_mat)), na.rm = TRUE),
    fraction_observed = mean(!is.na(beta_matrix)),
    stringsAsFactors = FALSE
  )

  return(list(
    beta_matrix = cleaned$beta_matrix,
    se_matrix = cleaned$se_matrix,
    trait_info = trait_info,
    verification_trait_matrix = trait_mat[
      colnames(cleaned$beta_matrix),
      ,
      drop = FALSE
    ],
    verification_trait_info = trait_info |>
      dplyr::filter(feature_type == "trait"),
    snp_info = snp_info,
    z_target = z_target,
    coloc_groups = pert$coloc_groups,
    target_trait_id = trait_id,
    prep_summary = prep_summary
  ))
}


#' @title Run EBMF Factorization
#' @description Fit an Empirical Bayes Matrix Factorization model using
#' \pkg{flashier}. The model decomposes the features x SNPs z-score matrix
#' \eqn{Y = LF' + E} into latent factors (\eqn{F}, SNP programs) and
#' loadings (\eqn{L}, feature contributions), with sparse priors estimated
#' from the data.
#'
#' Known per-observation standard errors are passed to flashier via \code{S}.
#' For z-score matrices, \code{se_mode = "unit"} (default) is recommended:
#' missing cells stay \code{NA} and flashier uses \code{S = 1}. The legacy
#' \code{se_mode = "matrix"} path fills missing values with 0 and large SEs.
#'
#' @param beta_matrix Features x SNPs matrix of signed z-scores, typically
#'   \code{build_ebmf_matrix()$beta_matrix}. Rows and columns that are entirely
#'   missing or zero are removed before fitting.
#' @param se_matrix Features x SNPs matrix of standard errors aligned with
#'   \code{beta_matrix}. Used to mark missing cells when \code{se_mode = "unit"}
#'   and passed as flashier \code{S} when \code{se_mode = "matrix"}.
#' @param se_mode \code{"unit"} (default): data keeps \code{NA} missingness and
#'   flashier uses \code{S = 1}. \code{"matrix"}: pass \code{se_matrix} with
#'   missing cells as \code{Y = 0}, \code{S = 1e6}.
#' @param greedy_Kmax Maximum number of factors to add greedily. Factors are
#'   only added while they improve the variational lower bound. Defaults to 50.
#' @param backfit Logical; if \code{TRUE} (default), all factors are cyclically
#'   updated after the greedy phase until convergence.
#' @param ebnm_fn EBNM function for priors on loadings and factors. Defaults to
#'   \code{ebnm::ebnm_point_normal} (point-normal prior). Pass a list of two
#'   functions to use different priors for loadings (features) and factors (SNPs).
#' @param var_type Residual variance structure. Defaults to \code{NULL} when
#'   \code{se_mode = "unit"} (S accounts for residual variance) and \code{1L}
#'   (per-row) when \code{se_mode = "matrix"}. Use \code{0} for constant,
#'   \code{1} for per-row, \code{2} for per-column, or \code{c(1, 2)} for
#'   Kronecker.
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
#' @export
run_ebmf <- function(beta_matrix,
                     se_matrix,
                     se_mode = c("unit", "matrix"),
                     greedy_Kmax = 50L,
                     backfit = TRUE,
                     ebnm_fn = NULL,
                     var_type = NULL,
                     verbose = 1L) {
  if (!requireNamespace("flashier", quietly = TRUE)) {
    stop("Package 'flashier' is required for EBMF clustering", call. = FALSE)
  }
  if (!requireNamespace("ebnm", quietly = TRUE)) {
    stop("Package 'ebnm' is required for EBMF clustering", call. = FALSE)
  }

  se_mode <- match.arg(se_mode)
  if (!is.matrix(beta_matrix)) {
    stop("beta_matrix must be a matrix")
  }
  if (!is.matrix(se_matrix)) {
    stop("se_matrix must be a matrix")
  }
  if (!identical(dim(beta_matrix), dim(se_matrix))) {
    stop("beta_matrix and se_matrix must have the same dimensions")
  }

  if (is.null(ebnm_fn)) {
    ebnm_fn <- ebnm::ebnm_point_normal
  }
  if (missing(var_type)) {
    var_type <- NULL
  }

  prepared <- .prepare_ebmf_flash_inputs(
    beta_matrix = beta_matrix,
    se_matrix = se_matrix,
    se_mode = se_mode
  )

  fit <- flashier::flash(
    data = prepared$data,
    S = prepared$S,
    ebnm_fn = ebnm_fn,
    var_type = var_type,
    greedy_Kmax = greedy_Kmax,
    backfit = backfit,
    nullcheck = TRUE,
    verbose = verbose
  )

  return(fit)
}


#' @title Run EBMF Comparison Grid
#' @description Fit a grid of EBMF models varying EBNM prior and residual
#' variance structure. The feature matrix is built once via
#' \code{build_ebmf_matrix()} and reused across runs. Defaults favour
#' \code{se_mode = "unit"} with \code{var_type} in \code{NULL} / \code{1L}
#' rather than Kronecker variance, which previously produced soft mega-clusters
#' on near-null fits.
#' @inheritParams build_ebmf_matrix
#' @inheritParams run_ebmf
#' @param ebnm_fns Character vector of EBNM priors: \code{"point_normal"} and/or
#'   \code{"point_laplace"}. Defaults to both.
#' @param se_modes Character vector of \code{se_mode} values. Defaults to
#'   \code{"unit"}.
#' @param var_types List of \code{var_type} values passed to \code{run_ebmf()}.
#'   Defaults to \code{list(NULL, 1L)}.
#' @param lfsr_threshold Passed to \code{extract_ebmf_clusters()}.
#' @param magnitude_threshold Passed to \code{extract_ebmf_clusters()}.
#'   Defaults to \code{0.25}.
#' @param save_path Optional path to save results as \code{.rds}.
#' @return A list with \code{summary} (one-row-per-run dataframe) and
#'   \code{results} (named list of per-run outputs).
#' @export
run_ebmf_comparison <- function(trait_id,
                                p_threshold = NULL,
                                snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                min_snp_signals = 5L,
                                compress_method = c("asinh", "none"),
                                ebnm_fns = c("point_normal", "point_laplace"),
                                se_modes = c("unit"),
                                var_types = list(NULL, 1L),
                                greedy_Kmax = 20L,
                                backfit = TRUE,
                                lfsr_threshold = 0.05,
                                magnitude_threshold = 0.25,
                                save_path = NULL,
                                verbose = 0L) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }

  snp_key <- match.arg(snp_key)
  compress_method <- match.arg(compress_method)
  se_modes <- match.arg(se_modes, c("unit", "matrix"), several.ok = TRUE)

  configs <- expand.grid(
    ebnm_fn = ebnm_fns,
    se_mode = se_modes,
    var_type_idx = seq_along(var_types),
    stringsAsFactors = FALSE
  )
  configs$var_type <- vapply(
    configs$var_type_idx,
    function(i) return(.ebmf_var_type_label(var_types[[i]])),
    character(1)
  )
  configs$run_id <- paste0(
    "ebnm_", configs$ebnm_fn,
    "__se_", configs$se_mode,
    "__var_", configs$var_type
  )

  message("Building EBMF matrix once for comparison grid")
  ebmf_data <- build_ebmf_matrix(
    trait_id = trait_id,
    p_threshold = p_threshold,
    snp_key = snp_key,
    min_snp_signals = min_snp_signals,
    compress_method = compress_method
  )

  n_runs <- nrow(configs)
  results <- vector("list", n_runs)
  names(results) <- configs$run_id
  summary_rows <- vector("list", n_runs)

  for (i in seq_len(n_runs)) {
    cfg <- configs[i, , drop = FALSE]
    run_id <- cfg$run_id
    message("EBMF run ", i, "/", n_runs, ": ", run_id)

    run_result <- tryCatch({
      flash_fit <- run_ebmf(
        beta_matrix = ebmf_data$beta_matrix,
        se_matrix = ebmf_data$se_matrix,
        se_mode = cfg$se_mode,
        greedy_Kmax = greedy_Kmax,
        backfit = backfit,
        ebnm_fn = .resolve_ebnm_fn(cfg$ebnm_fn),
        var_type = var_types[[cfg$var_type_idx]],
        verbose = verbose
      )

      clusters <- extract_ebmf_clusters(
        flash_fit,
        lfsr_threshold = lfsr_threshold,
        magnitude_threshold = magnitude_threshold
      )

      drivers <- summarise_ebmf_program_drivers(
        flash_fit = flash_fit,
        trait_info = ebmf_data$trait_info
      )

      list(
        ebnm_fn = cfg$ebnm_fn,
        se_mode = cfg$se_mode,
        var_type = var_types[[cfg$var_type_idx]],
        ebmf_data = ebmf_data,
        flash_fit = flash_fit,
        ebmf_clusters = clusters,
        program_drivers = drivers
      )
    }, error = function(e) {
      warning("Run '", run_id, "' failed: ", conditionMessage(e), call. = FALSE)
      return(list(
        ebnm_fn = cfg$ebnm_fn,
        se_mode = cfg$se_mode,
        var_type = var_types[[cfg$var_type_idx]],
        error = conditionMessage(e)
      ))
    })

    results[[run_id]] <- run_result
    summary_rows[[i]] <- .ebmf_run_summary_row(run_id, run_result)
  }

  summary_df <- dplyr::bind_rows(summary_rows)
  output <- list(summary = summary_df, results = results, ebmf_data = ebmf_data)

  if (!is.null(save_path)) {
    saveRDS(output, save_path)
    message("Saved to ", save_path)
  }

  return(output)
}


#' @title Summarise EBMF Comparison Results
#' @description Rebuild or return the summary table from a comparison result.
#' @param comparison Output list from \code{run_ebmf_comparison()}.
#' @return A dataframe with one row per run.
#' @export
summarise_ebmf_comparison <- function(comparison) {
  if (!is.null(comparison$summary)) {
    return(comparison$summary)
  }
  summary_rows <- lapply(names(comparison$results), function(run_id) {
    return(.ebmf_run_summary_row(run_id, comparison$results[[run_id]]))
  })
  return(dplyr::bind_rows(summary_rows))
}


#' @title Select Best EBMF Comparison Run
#' @description Prefer successful runs with total PVE at least
#' \code{min_total_pve}, then highest total PVE, then highest ELBO. Tiny
#' positive PVE values (e.g. \code{1e-8}) are treated as null fits.
#' @param comparison Output list from \code{run_ebmf_comparison()}.
#' @param min_total_pve Minimum acceptable total PVE. Defaults to \code{0.01}.
#' @return Character \code{run_id}, or \code{NULL} if none meet the PVE floor.
#' @export
select_ebmf_comparison_run <- function(comparison, min_total_pve = 0.01) {
  summary_df <- summarise_ebmf_comparison(comparison)
  ok <- summary_df |>
    dplyr::filter(status == "ok", !is.na(n_factors), n_factors > 0)
  if (nrow(ok) == 0) {
    return(NULL)
  }

  if (!"total_pve" %in% names(ok)) {
    ok$total_pve <- NA_real_
  }
  if (!"elbo" %in% names(ok)) {
    ok$elbo <- NA_real_
  }

  candidates <- ok |>
    dplyr::filter(!is.na(total_pve), total_pve >= min_total_pve)
  if (nrow(candidates) == 0) {
    best_null <- max(ok$total_pve, na.rm = TRUE)
    warning(
      "No EBMF run reached min_total_pve = ", min_total_pve,
      " (best total_pve = ", signif(best_null, 3), "). ",
      "Soft membership from near-null fits should not be trusted.",
      call. = FALSE
    )
    return(NULL)
  }

  candidates <- candidates |>
    dplyr::arrange(
      dplyr::desc(total_pve),
      dplyr::desc(elbo),
      dplyr::desc(n_factors)
    )

  return(candidates$run_id[1])
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
#'   total PVE and soft-assigns at least 50\% of SNPs. Defaults to \code{0.4}.
#' @param min_soft_frac Candidate if the factor soft-assigns at least this
#'   fraction of SNPs. Defaults to \code{0.85}.
#' @param max_drop Maximum number of factors to return. Defaults to 1.
#' @return Integer vector of 1-based factor indices to pass to
#'   \code{remove_ebmf_factors()}, possibly empty.
#' @export
identify_ebmf_global_factors <- function(flash_fit,
                                         membership = NULL,
                                         min_total_pve = 0.01,
                                         min_pve_share = 0.4,
                                         min_soft_frac = 0.85,
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

  candidates <- which(
    soft_frac >= min_soft_frac |
      (share >= min_pve_share & soft_frac >= 0.5)
  )
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


#' @title Summarise EBMF Program Feature Drivers
#' @description Rank feature loadings for each discovered EBMF program.
#' @param flash_fit A \code{flash} object from \code{run_ebmf()}.
#' @param trait_info Optional feature metadata with \code{trait_id} and
#'   \code{trait_name}.
#' @param lfsr_threshold lFSR threshold for feature loadings. Defaults to
#'   \code{0.05}.
#' @param magnitude_threshold Minimum absolute loading. Defaults to \code{0.05}.
#' @return A dataframe of significant feature loadings per program.
#' @export
summarise_ebmf_program_drivers <- function(flash_fit,
                                           trait_info = NULL,
                                           lfsr_threshold = 0.05,
                                           magnitude_threshold = 0.05) {
  if (!inherits(flash_fit, "flash")) {
    stop("flash_fit must be a flash object from run_ebmf()")
  }

  K <- flash_fit$n_factors
  if (K == 0) {
    return(data.frame(
      program = integer(0),
      trait_id = character(0),
      loading = numeric(0),
      lfsr = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  L_pm <- flash_fit$L_pm
  L_lfsr <- flash_fit$L_lfsr
  trait_ids <- rownames(L_pm)
  if (is.null(trait_ids)) {
    trait_ids <- as.character(seq_len(nrow(L_pm)))
  }

  summaries <- lapply(seq_len(K), function(k) {
    sig <- !is.na(L_lfsr[, k]) &
      L_lfsr[, k] < lfsr_threshold &
      abs(L_pm[, k]) > magnitude_threshold
    if (sum(sig) == 0) {
      return(NULL)
    }
    out <- data.frame(
      program = k,
      trait_id = trait_ids[sig],
      loading = L_pm[sig, k],
      lfsr = L_lfsr[sig, k],
      stringsAsFactors = FALSE
    ) |>
      dplyr::arrange(dplyr::desc(abs(loading)))
    return(out)
  })

  result <- dplyr::bind_rows(summaries)

  if (!is.null(trait_info) && nrow(result) > 0) {
    trait_info <- trait_info |>
      dplyr::mutate(trait_id = as.character(trait_id))
    join_cols <- c("trait_id", "trait_name")
    optional_cols <- c("trait_label", "trait_category", "feature_type", "is_qtl")
    join_cols <- unique(c(join_cols, intersect(optional_cols, names(trait_info))))
    result <- result |>
      dplyr::left_join(trait_info[, join_cols, drop = FALSE], by = "trait_id")
  }

  return(result)
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
#' @export
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
#' @export
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


.ebmf_slice_has_signal <- function(x) {
  return(any(!is.na(x) & x != 0))
}


.ebmf_run_summary_row <- function(run_id, run_result) {
  se_mode <- if (is.null(run_result$se_mode)) {
    NA_character_
  } else {
    as.character(run_result$se_mode)
  }
  var_type <- .ebmf_var_type_label(run_result$var_type)

  if (!is.null(run_result$error)) {
    ebnm_fn <- if (is.null(run_result$ebnm_fn)) {
      NA_character_
    } else {
      as.character(run_result$ebnm_fn)
    }
    return(data.frame(
      run_id = run_id,
      ebnm_fn = ebnm_fn,
      se_mode = se_mode,
      var_type = var_type,
      n_factors = NA_integer_,
      total_pve = NA_real_,
      pve = NA_character_,
      elbo = NA_real_,
      n_programs_with_snps = NA_integer_,
      n_assigned = NA_integer_,
      n_multi_program = NA_integer_,
      program_sizes = NA_character_,
      status = "failed",
      error = run_result$error,
      stringsAsFactors = FALSE
    ))
  }

  clusters <- run_result$ebmf_clusters
  program_sizes <- integer(0)
  if (!is.null(clusters$membership) && ncol(clusters$membership) > 0) {
    program_sizes <- sort(colSums(clusters$membership), decreasing = TRUE)
  }

  pve <- run_result$flash_fit$pve
  total_pve <- if (is.null(pve) || length(pve) == 0) {
    NA_real_
  } else {
    sum(pve, na.rm = TRUE)
  }
  pve_str <- if (is.null(pve) || length(pve) == 0) {
    NA_character_
  } else {
    paste(round(pve, 4), collapse = ", ")
  }
  elbo <- if (is.null(run_result$flash_fit$elbo)) {
    NA_real_
  } else {
    as.numeric(run_result$flash_fit$elbo)
  }

  return(data.frame(
    run_id = run_id,
    ebnm_fn = as.character(run_result$ebnm_fn),
    se_mode = se_mode,
    var_type = var_type,
    n_factors = run_result$flash_fit$n_factors,
    total_pve = total_pve,
    pve = pve_str,
    elbo = elbo,
    n_programs_with_snps = sum(program_sizes > 0),
    n_assigned = clusters$n_assigned,
    n_multi_program = clusters$n_multi_program,
    program_sizes = paste(program_sizes, collapse = ", "),
    status = "ok",
    error = NA_character_,
    stringsAsFactors = FALSE
  ))
}


.filter_ebmf_matrices <- function(beta_matrix, se_matrix) {
  if (!is.matrix(beta_matrix)) {
    stop("beta_matrix must be a matrix")
  }
  if (!is.matrix(se_matrix)) {
    stop("se_matrix must be a matrix")
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

  beta_matrix <- beta_matrix[row_keep, col_keep, drop = FALSE]
  se_matrix <- se_matrix[row_keep, col_keep, drop = FALSE]

  if (nrow(beta_matrix) == 0 || ncol(beta_matrix) == 0) {
    stop("No rows or columns remain after removing all-zero/all-missing slices")
  }

  if (anyNA(beta_matrix)) {
    beta_matrix[is.na(beta_matrix)] <- 0
  }
  if (anyNA(se_matrix)) {
    se_matrix[is.na(se_matrix)] <- 1e6
  }

  return(list(
    beta_matrix = beta_matrix,
    se_matrix = se_matrix
  ))
}


.prepare_ebmf_flash_inputs <- function(beta_matrix, se_matrix, se_mode) {
  filtered <- .filter_ebmf_matrices(
    beta_matrix = beta_matrix,
    se_matrix = se_matrix
  )

  # Restore NA missingness; flashier skips NA cells natively. In "matrix"
  # mode the observed per-cell SEs are passed as flashier's S (with NAs at
  # unobserved cells), replacing the legacy Y = 0 / S = 1e6 fill which
  # degraded fits.
  data <- filtered$beta_matrix
  S <- if (se_mode == "matrix") filtered$se_matrix else 1
  missing <- is.na(data) |
    (!is.matrix(S) && FALSE) | (is.matrix(S) & is.na(S)) |
    (is.matrix(S) & S >= 1e6)
  data[missing] <- NA_real_
  if (is.matrix(S)) {
    S[is.na(S)] <- NA_real_
    S[S >= 1e6] <- NA_real_
  }

  return(list(data = data, S = S))
}


.ebmf_var_type_label <- function(var_type) {
  if (is.null(var_type)) {
    return("null")
  }
  return(paste(var_type, collapse = ","))
}


.ebmf_feature_info_from_perturbation <- function(feature_info,
                                                 kept_ids,
                                                 feature_type,
                                                 coloc_groups) {
  if (is.null(feature_info) || nrow(feature_info) == 0 || length(kept_ids) == 0) {
    return(data.frame(
      trait_id = character(0),
      trait_name = character(0),
      feature_type = character(0),
      trait_category = character(0),
      trait_label = character(0),
      is_qtl = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  info <- feature_info |>
    dplyr::mutate(
      trait_id = as.character(feature_id),
      trait_name = as.character(feature_name)
    ) |>
    dplyr::filter(trait_id %in% as.character(kept_ids))

  if (feature_type == "trait") {
    cats <- coloc_groups |>
      dplyr::mutate(trait_id = as.character(trait_id)) |>
      dplyr::filter(
        trait_id %in% info$trait_id,
        !is.na(trait_category),
        trait_category != ""
      ) |>
      dplyr::group_by(trait_id) |>
      dplyr::summarise(
        trait_category = dplyr::first(trait_category),
        .groups = "drop"
      )
    info <- info |>
      dplyr::left_join(cats, by = "trait_id") |>
      dplyr::mutate(
        feature_type = "trait",
        trait_category = dplyr::if_else(
          is.na(trait_category),
          "Unknown",
          trait_category
        ),
        trait_label = trait_category,
        is_qtl = FALSE
      )
  } else {
    info <- info |>
      dplyr::mutate(
        feature_type = "gene",
        trait_category = "Gene",
        trait_label = dplyr::if_else(
          is.na(trait_name) | trait_name == "",
          trait_id,
          trait_name
        ),
        is_qtl = TRUE
      )
  }

  return(info |>
    dplyr::select(
      trait_id, trait_name, feature_type, trait_category, trait_label, is_qtl
    ) |>
    dplyr::arrange(trait_id))
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
