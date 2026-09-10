#' @title Run Univariate SNP Program Discovery
#' @description Run the full univariate analysis pipeline (pleiotropy matrix ->
#' sparse/ubiquitous trait filter -> orientation -> compression -> SNP cosine
#' similarity -> Empirical Bayes Matrix Factorization programs) from a single
#' entry point. The function takes a trait object as returned by
#' `trait(trait_id, include_associations = TRUE)` — or a simulated object with
#' the same shape — so no API calls are made when `coloc_groups` (and, for
#' `associations = "full"`, `full_associations`) are supplied.
#'
#' The simulated `coloc_groups` dataframe must include at least: `coloc_group_id`,
#' `variant_id`, `display_snp`, `chr`, `bp`, `trait_id`, `trait_name`, `min_p`,
#' `beta`, and `se` (plus `snp_key` if not `"variant_id"`). Rows where
#' `trait_id == target_trait_id` define the target SNPs (columns); all other
#' rows become background trait profiles. Optional annotation columns
#' (`gene_id`, `gene`, `tissue`, `trait_category`) enrich downstream
#' interpretation functions.
#' @param trait_object A trait result list (from `trait()` or a simulation) with
#'   a `coloc_groups` dataframe. For `associations = "full"` it must also contain
#'   `full_associations`.
#' @param target_trait_id Target trait id. Defaults to `trait_object$trait$id`.
#' @param associations Which association set backs the matrix: `"coloc"`
#'   (default; colocalisation rows only) or `"full"` (all full associations).
#' @param p_threshold Optional p-value threshold for including target SNPs.
#' @param snp_key Column used to name SNP columns.
#' @param include_trans If `FALSE` (default), rows flagged `cis_trans == "trans"`
#'   are removed before matrix construction (requires a `cis_trans` column;
#'   ignored otherwise).
#' @param trait_subset Which trait rows back the trait x SNP matrix:
#'   `"all"` (default) uses every available background trait, while
#'   `"phenotypic"` keeps only phenotype traits (rows whose `feature_type` is
#'   `"phenotypic"` or missing) plus the target trait, dropping molecular /
#'   gene-annotated rows. Molecular traits are much sparser, so clustering on
#'   phenotype traits only can recover structure that sparse molecular rows
#'   dilute. Requires `trait_info` with a `feature_type` column (present for
#'   `trait()` and simulated objects).
#' @param min_snp_signals Minimum non-`NA` SNPs per background trait row.
#' @param max_snp_fraction Maximum non-`NA` fraction per background trait row.
#' @param compress_method Effect compression passed to `compress_effect_matrix()`.
#' @param compress_scale Asinh scale passed to `compress_effect_matrix()`.
#' @param ebmf_greedy_Kmax Maximum number of EBMF factors.
#' @param ebmf_lfsr_threshold lFSR threshold for EBMF program membership.
#' @param ebmf_magnitude_threshold Minimum absolute posterior mean loading for
#'   EBMF program membership.
#' @param ebmf_drop_global If `TRUE` (default), drop at most one global
#'   mega-factor (via `identify_ebmf_global_factors()`) before extracting EBMF
#'   programs. This guards against the near-null "one giant cluster" failure
#'   mode.
#' @param ebmf_prior Prior family for EBMF loadings/factors:
#'   `"point_normal"` (default) or `"point_laplace"` (heavier tails, often
#'   better when true effects are large relative to noise).
#' @param ebmf_backfit If `TRUE` (default), cyclically refit all factors after
#'   the greedy phase. `FALSE` is faster but leaves greedy-order artifacts.
#' @param ebmf_se_mode Noise model for EBMF. `"unit"` (default) treats every
#'   observed z-score as having standard error 1. `"matrix"` passes the
#'   observed per-cell standard errors to flashier, so imprecise estimates are
#'   down-weighted and winner's-curse inflation of tiny-but-precise effects is
#'   avoided. In `"matrix"` mode the EBMF input is the oriented *beta* matrix
#'   (compression is skipped, as it would violate the noise model), so
#'   `ebmf_magnitude_threshold` reads raw-beta-scale loadings — consider `NA`
#'   to gate membership on lFSR alone. Caveat: betas from different traits are
#'   on study-native scales, so factors can be dominated by large-unit traits;
#'   z-scores (`"unit"`) harmonise units at the cost of discarding precision.
#' @param ebmf_beta_scale Scale normalisation for `ebmf_se_mode = "matrix"`:
#'   `"none"` uses study-native betas; `"trait"` (recommended) divides each
#'   trait row's betas and SEs by that row's root-mean-square beta, which
#'   harmonises units across traits while preserving within-trait precision
#'   information (z-scores remain unchanged by this since beta/se is
#'   scale-invariant).
#' @param min_module_size,min_mean_internal,min_connectedness Recorded pipeline
#'   settings for the downstream program-validation layer and used (for
#'   `min_module_size`) by `calibrate_ebmf_programs()`. They do not affect the
#'   EBMF fit itself; pass them explicitly to `summarise_ebmf_programs()` to
#'   gate programs.
#' @return A list with:
#'   \itemize{
#'     \item x_matrix: raw traits x SNPs z-score matrix after trait filtering
#'     \item x_star: oriented, compressed traits x SNPs matrix used for similarity
#'     \item trait_matrix: SNP x trait version of `x_star`
#'     \item s_matrix: SNP-by-SNP cosine similarity matrix
#'     \item overlap_matrix, eligible_matrix: joint-observation counts / eligibility
#'     \item cluster_membership: SNP x program logical membership matrix from EBMF
#'     \item cluster_details: EBMF fit details (including the `flash` object)
#'     \item beta_matrix, se_matrix: trait x SNP beta / SE matrices when available
#'     \item ebmf_input: the features x SNPs matrix actually passed to flashier
#'     \item trait_info, snp_info: row / column metadata
#'     \item dropped_trait_ids: background traits removed by the feature-type
#'       (`trait_subset`) and/or sparse/ubiquitous filter
#'     \item coloc_groups: coloc groups actually used (after trans filtering)
#'     \item parameters: list of settings used
#'   }
#' @export
run_univariate_clustering <- function(trait_object,
                                      target_trait_id = NULL,
                                      associations = c("coloc", "full"),
                                      p_threshold = NULL,
                                      snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                      include_trans = FALSE,
                                      trait_subset = c("all", "phenotypic"),
                                      min_snp_signals = 5L,
                                      max_snp_fraction = 0.8,
                                      compress_method = c("none", "asinh"),
                                      compress_scale = 2,
                                      ebmf_greedy_Kmax = 50L,
                                      ebmf_lfsr_threshold = 0.05,
                                      ebmf_magnitude_threshold = 0.25,
                                      ebmf_drop_global = TRUE,
                                      ebmf_prior = c("point_normal", "point_laplace"),
                                      ebmf_backfit = TRUE,
                                      ebmf_se_mode = c("unit", "matrix"),
                                      ebmf_beta_scale = c("none", "trait"),
                                      min_module_size = 3L,
                                      min_mean_internal = 0.3,
                                      min_connectedness = 0.5) {
  associations <- match.arg(associations)
  snp_key <- match.arg(snp_key)
  trait_subset <- match.arg(trait_subset)
  compress_method <- match.arg(compress_method)
  ebmf_prior <- match.arg(ebmf_prior)
  ebmf_se_mode <- match.arg(ebmf_se_mode)
  ebmf_beta_scale <- match.arg(ebmf_beta_scale)

  if (is.null(trait_object)) {
    stop("trait_object is required")
  }
  coloc_groups <- trait_object$coloc_groups
  if (!is.data.frame(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("trait_object$coloc_groups must be a non-empty dataframe")
  }

  if (is.null(target_trait_id)) {
    if (!is.null(trait_object$trait) && !is.null(trait_object$trait$id)) {
      target_trait_id <- trait_object$trait$id
    } else {
      stop("target_trait_id is required (could not infer from trait_object$trait$id)")
    }
  }
  target_id <- as.character(target_trait_id)

  if (!include_trans && "cis_trans" %in% names(coloc_groups)) {
    coloc_groups <- coloc_groups |>
      dplyr::filter(is.na(cis_trans) | tolower(cis_trans) != "trans")
    if (nrow(coloc_groups) == 0) {
      stop("No coloc_groups rows left after removing trans markers")
    }
  }

  if (associations == "full") {
    full_associations <- trait_object$full_associations
    if (!is.data.frame(full_associations) || nrow(full_associations) == 0) {
      stop("associations = 'full' requires trait_object$full_associations")
    }
    locus_data <- .prepare_dense_locus_data(
      trait_id = target_trait_id,
      coloc_groups = coloc_groups,
      full_associations = full_associations,
      p_threshold = p_threshold,
      snp_key = snp_key
    )
  } else {
    locus_data <- .prepare_pleiotropy_locus_data(
      trait_id = target_trait_id,
      coloc_groups = coloc_groups,
      p_threshold = p_threshold,
      snp_key = snp_key
    )
  }
  pleiotropy <- .finalize_pleiotropy_from_locus(locus_data, target_id)

  X <- pleiotropy$x_matrix
  if (!target_id %in% rownames(X)) {
    stop("target trait row missing from the built pleiotropy matrix")
  }

  subset_dropped <- character(0)
  if (trait_subset == "phenotypic") {
    info <- pleiotropy$trait_info
    if (is.null(info) || !"feature_type" %in% names(info)) {
      stop("trait_subset = 'phenotypic' requires trait_info with a feature_type column")
    }
    keep_ids <- as.character(info$trait_id[
      as.character(info$trait_id) == target_id |
        is.na(info$feature_type) |
        info$feature_type == "phenotypic"
    ])
    subset_dropped <- setdiff(rownames(X), keep_ids)
    if (length(subset_dropped) > 0) {
      X <- X[!rownames(X) %in% subset_dropped, , drop = FALSE]
      for (mat in c("beta_matrix", "se_matrix")) {
        if (!is.null(pleiotropy[[mat]])) {
          pleiotropy[[mat]] <- pleiotropy[[mat]][
            !rownames(pleiotropy[[mat]]) %in% subset_dropped, , drop = FALSE
          ]
        }
      }
      pleiotropy$trait_info <- info |>
        dplyr::filter(!as.character(trait_id) %in% subset_dropped)
    }
  }

  trait_snp_counts <- rowSums(!is.na(X))
  trait_snp_frac <- trait_snp_counts / ncol(X)
  sparse_trait_ids <- names(trait_snp_counts)[trait_snp_counts < min_snp_signals]
  ubiquitous_trait_ids <- names(trait_snp_frac)[trait_snp_frac > max_snp_fraction]
  drop_trait_ids <- setdiff(
    union(sparse_trait_ids, ubiquitous_trait_ids),
    target_id
  )
  if (length(drop_trait_ids) > 0) {
    X <- X[!rownames(X) %in% drop_trait_ids, , drop = FALSE]
    for (mat in c("beta_matrix", "se_matrix")) {
      if (!is.null(pleiotropy[[mat]])) {
        pleiotropy[[mat]] <- pleiotropy[[mat]][
          !rownames(pleiotropy[[mat]]) %in% drop_trait_ids, , drop = FALSE
        ]
      }
    }
    pleiotropy$trait_info <- pleiotropy$trait_info |>
      dplyr::filter(!as.character(trait_id) %in% drop_trait_ids)
  }

  oriented <- orient_pleiotropy_matrix(X, target_trait_id = target_trait_id)
  X_star <- compress_effect_matrix(
    oriented$x_matrix,
    method = compress_method,
    asinh_scale = compress_scale
  )
  trait_matrix <- t(X_star)

  similarity <- snp_similarity_matrix(X_star)

  ebmf_x_input <- X_star
  ebmf_se_input <- NULL
  if (ebmf_se_mode == "matrix") {
    if (is.null(pleiotropy$beta_matrix) || is.null(pleiotropy$se_matrix)) {
      stop("ebmf_se_mode = 'matrix' requires beta and se columns in coloc_groups")
    }
    shared_rows <- intersect(
      rownames(oriented$x_matrix),
      rownames(pleiotropy$beta_matrix)
    )
    beta_oriented <- sweep(
      pleiotropy$beta_matrix[shared_rows, colnames(X), drop = FALSE],
      2, oriented$target_signs[colnames(X)], `*`
    )
    ebmf_x_input <- beta_oriented[rownames(X), , drop = FALSE]
    ebmf_se_input <- pleiotropy$se_matrix[
      rownames(ebmf_x_input), colnames(X), drop = FALSE
    ]
    if (ebmf_beta_scale == "trait") {
      row_scale <- sqrt(rowMeans(ebmf_x_input^2, na.rm = TRUE))
      row_scale[!is.finite(row_scale) | row_scale <= 0] <- 1
      ebmf_x_input <- sweep(ebmf_x_input, 1, row_scale, `/`)
      ebmf_se_input <- sweep(ebmf_se_input, 1, row_scale, `/`)
    }
  }

  ebmf_fit <- .cluster_snp_profiles_ebmf(
    ebmf_x_input,
    greedy_Kmax = ebmf_greedy_Kmax,
    lfsr_threshold = ebmf_lfsr_threshold,
    magnitude_threshold = ebmf_magnitude_threshold,
    drop_global = ebmf_drop_global,
    prior = ebmf_prior,
    backfit = ebmf_backfit,
    observed_se_matrix = ebmf_se_input
  )

  return(list(
    x_matrix = X,
    x_star = X_star,
    trait_matrix = trait_matrix,
    s_matrix = similarity$s_matrix,
    overlap_matrix = similarity$overlap_matrix,
    eligible_matrix = similarity$eligible_matrix,
    cluster_membership = ebmf_fit$details$membership,
    cluster_details = ebmf_fit$details,
    beta_matrix = if (!is.null(pleiotropy$beta_matrix)) {
      pleiotropy$beta_matrix[rownames(X), , drop = FALSE]
    } else {
      NULL
    },
    se_matrix = if (!is.null(pleiotropy$se_matrix)) {
      pleiotropy$se_matrix[rownames(X), , drop = FALSE]
    } else {
      NULL
    },
    ebmf_input = ebmf_x_input,
    trait_info = pleiotropy$trait_info,
    snp_info = pleiotropy$snp_info,
    dropped_trait_ids = union(drop_trait_ids, subset_dropped),
    coloc_groups = coloc_groups,
    parameters = list(
      method = "ebmf",
      target_trait_id = target_id,
      associations = associations,
      p_threshold = p_threshold,
      snp_key = snp_key,
      include_trans = include_trans,
      trait_subset = trait_subset,
      min_snp_signals = min_snp_signals,
      max_snp_fraction = max_snp_fraction,
      compress_method = compress_method,
      compress_scale = compress_scale,
      ebmf_greedy_Kmax = ebmf_greedy_Kmax,
      ebmf_lfsr_threshold = ebmf_lfsr_threshold,
      ebmf_magnitude_threshold = ebmf_magnitude_threshold,
      ebmf_drop_global = ebmf_drop_global,
      ebmf_prior = ebmf_prior,
      ebmf_backfit = ebmf_backfit,
      ebmf_se_mode = ebmf_se_mode,
      ebmf_beta_scale = ebmf_beta_scale,
      min_module_size = min_module_size,
      min_mean_internal = min_mean_internal,
      min_connectedness = min_connectedness
    )
  ))
}


.assert_ebmf_result <- function(clustering_result) {
  params <- clustering_result$parameters
  if (is.null(params) || !identical(params$method, "ebmf")) {
    stop("clustering_result must come from run_univariate_clustering()")
  }
  return(invisible(TRUE))
}


.cluster_snp_profiles_ebmf <- function(x_matrix,
                                       greedy_Kmax = 50L,
                                       lfsr_threshold = 0.05,
                                       magnitude_threshold = 0.25,
                                       drop_global = TRUE,
                                       prior = "point_normal",
                                       backfit = TRUE,
                                       observed_se_matrix = NULL) {
  beta_matrix <- x_matrix
  se_matrix <- matrix(
    1,
    nrow = nrow(beta_matrix),
    ncol = ncol(beta_matrix),
    dimnames = dimnames(beta_matrix)
  )
  se_mode <- "unit"
  if (!is.null(observed_se_matrix)) {
    stopifnot(identical(dim(beta_matrix), dim(observed_se_matrix)))
    se_matrix[!is.na(beta_matrix)] <- observed_se_matrix[!is.na(beta_matrix)]
    se_matrix[is.na(beta_matrix)] <- NA_real_
    se_mode <- "matrix"
  } else {
    se_matrix[is.na(beta_matrix)] <- NA_real_
  }

  fit_error <- NULL
  flash_fit <- tryCatch(
    run_ebmf(
      beta_matrix = beta_matrix,
      se_matrix = se_matrix,
      se_mode = se_mode,
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
      dropped_global_factors = dropped_global
    )
  ))
}
