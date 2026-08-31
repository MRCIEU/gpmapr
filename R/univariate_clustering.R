#' @title Run Univariate SNP Clustering Pipeline
#' @description Run the full univariate analysis pipeline (pleiotropy matrix ->
#' sparse/ubiquitous trait filter -> orientation -> compression -> SNP cosine
#' similarity -> signed Louvain clustering -> module reliability) from a single
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
#' @param similarity_threshold Edge threshold applied before Louvain clustering.
#' @param cluster_type Clustering method: `"louvain"` (default, signed Louvain
#'   on the SNP similarity graph), `"spectral"` (normalised Laplacian embedding
#'   plus k-means), or `"ebmf"` (empirical Bayes matrix factorization via
#'   flashier; soft factors are converted to hard clusters by assigning each
#'   SNP to the factor with the largest absolute posterior loading).
#' @param spectral_k Number of clusters for `cluster_type = "spectral"`.
#' @param ebmf_greedy_Kmax Maximum number of EBMF factors for
#'   `cluster_type = "ebmf"`.
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
#' @param ebmf_hard_assignment If `TRUE`, additionally create a legacy hard
#'   SNP-to-program assignment by selecting the strongest surviving factor.
#'   Defaults to `FALSE`: EBMF results remain soft and overlapping.
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
#' @param louvain_gamma Signed Louvain resolution parameter.
#' @param seed RNG seed for Louvain node-order randomisation.
#' @param min_module_size Minimum SNPs for a module to be reliable.
#' @param min_mean_internal Minimum mean internal similarity for reliability.
#' @param min_connectedness Minimum pair-connectedness for reliability.
#' @param compute_specific_traits If `TRUE`, also run
#'   `summarise_module_specific_traits()` on the reliable modules using the
#'   compressed trait matrix (no API calls). Defaults to `FALSE`.
#' @param min_specificity Soft specificity threshold for module-specific traits.
#' @param top_n_traits Number of top traits returned per module.
#' @return A list with:
#'   \itemize{
#'     \item x_matrix: raw traits x SNPs z-score matrix after trait filtering
#'     \item x_star: oriented, compressed traits x SNPs matrix used for similarity
#'     \item trait_matrix: SNP x trait version of `x_star`
#'     \item s_matrix: SNP-by-SNP cosine similarity matrix
#'     \item overlap_matrix, eligible_matrix: joint-observation counts / eligibility
#'     \item clusters: named Louvain cluster assignment per SNP
#'     \item clusters_reliable: subset of `clusters` restricted to reliable modules
#'     \item module_quality: reliability metrics dataframe
#'     \item ebmf_input: the features x SNPs matrix actually passed to flashier
#'       when `cluster_type = "ebmf"` (`NULL` otherwise)
#'     \item trait_info, snp_info: row / column metadata
#'     \item dropped_trait_ids: background traits removed by the feature-type
#'       (`trait_subset`) and/or sparse/ubiquitous filter
#'     \item group_traits: module-specific traits (when `compute_specific_traits = TRUE`)
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
                                      similarity_threshold = 0,
                                      cluster_type = c("louvain", "spectral", "ebmf"),
                                      spectral_k = 3L,
                                      ebmf_greedy_Kmax = 50L,
                                      ebmf_lfsr_threshold = 0.05,
                                      ebmf_magnitude_threshold = 0.25,
                                      ebmf_drop_global = TRUE,
                                      ebmf_prior = c("point_normal", "point_laplace"),
                                      ebmf_backfit = TRUE,
                                      ebmf_se_mode = c("unit", "matrix"),
                                      ebmf_beta_scale = c("none", "trait"),
                                      ebmf_hard_assignment = FALSE,
                                      louvain_gamma = 2,
                                      seed = 1L,
                                      min_module_size = 3L,
                                      min_mean_internal = 0.3,
                                      min_connectedness = 0.5,
                                      compute_specific_traits = FALSE,
                                      min_specificity = 1.25,
                                      top_n_traits = 10L) {
  associations <- match.arg(associations)
  snp_key <- match.arg(snp_key)
  trait_subset <- match.arg(trait_subset)
  compress_method <- match.arg(compress_method)
  cluster_type <- match.arg(cluster_type)
  ebmf_prior <- match.arg(ebmf_prior)
  ebmf_se_mode <- match.arg(ebmf_se_mode)
  ebmf_beta_scale <- match.arg(ebmf_beta_scale)
  if (!is.logical(ebmf_hard_assignment) || length(ebmf_hard_assignment) != 1L ||
      is.na(ebmf_hard_assignment)) {
    stop("ebmf_hard_assignment must be TRUE or FALSE")
  }

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
  if (cluster_type == "ebmf" && ebmf_se_mode == "matrix") {
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

  clusters <- switch(
    cluster_type,
    louvain = cluster_snp_profiles_louvain(
      similarity$s_matrix,
      similarity_threshold = similarity_threshold,
      gamma = louvain_gamma,
      seed = seed
    ),
    spectral = cluster_snp_profiles_spectral(
      similarity$s_matrix,
      k = spectral_k,
      similarity_threshold = similarity_threshold
    ),
    ebmf = .cluster_snp_profiles_ebmf(
      ebmf_x_input,
      greedy_Kmax = ebmf_greedy_Kmax,
      lfsr_threshold = ebmf_lfsr_threshold,
      magnitude_threshold = ebmf_magnitude_threshold,
      drop_global = ebmf_drop_global,
      prior = ebmf_prior,
      backfit = ebmf_backfit,
      hard_assignment = ebmf_hard_assignment,
      observed_se_matrix = ebmf_se_input
    )
  )

  if (length(clusters$cluster) >= 2) {
    module_quality <- summarise_snp_module_quality(
      similarity$s_matrix,
      clusters$cluster,
      edge_threshold = similarity_threshold,
      min_module_size = min_module_size,
      min_mean_internal = min_mean_internal,
      min_connectedness = min_connectedness
    )
  } else {
    module_quality <- .empty_module_quality()
  }
  reliable_ids <- module_quality$cluster[module_quality$reliable]
  clusters_reliable <- clusters$cluster[clusters$cluster %in% reliable_ids]

  group_traits <- NULL
  if (compute_specific_traits && length(clusters_reliable) > 0) {
    group_traits <- summarise_module_specific_traits(
      trait_matrix = trait_matrix,
      groups = clusters_reliable,
      trait_info = pleiotropy$trait_info,
      coloc_groups = coloc_groups,
      exclude_trait_ids = target_trait_id,
      min_group_size = min_module_size,
      top_n = top_n_traits,
      n_categories = 3L,
      min_specificity = min_specificity
    )
  }

  return(list(
    x_matrix = X,
    x_star = X_star,
    trait_matrix = trait_matrix,
    s_matrix = similarity$s_matrix,
    overlap_matrix = similarity$overlap_matrix,
    eligible_matrix = similarity$eligible_matrix,
    clusters = clusters$cluster,
    clusters_reliable = clusters_reliable,
    cluster_membership = if (!is.null(clusters$details$membership)) {
      clusters$details$membership
    } else {
      NULL
    },
    cluster_details = clusters$details,
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
    ebmf_input = if (cluster_type == "ebmf") ebmf_x_input else NULL,
    module_quality = module_quality,
    trait_info = pleiotropy$trait_info,
    snp_info = pleiotropy$snp_info,
    dropped_trait_ids = union(drop_trait_ids, subset_dropped),
    group_traits = group_traits,
    coloc_groups = coloc_groups,
    parameters = list(
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
      similarity_threshold = similarity_threshold,
      cluster_type = cluster_type,
      spectral_k = spectral_k,
      ebmf_greedy_Kmax = ebmf_greedy_Kmax,
      ebmf_lfsr_threshold = ebmf_lfsr_threshold,
      ebmf_magnitude_threshold = ebmf_magnitude_threshold,
      ebmf_drop_global = ebmf_drop_global,
      ebmf_prior = ebmf_prior,
      ebmf_backfit = ebmf_backfit,
      ebmf_se_mode = ebmf_se_mode,
      ebmf_beta_scale = ebmf_beta_scale,
      ebmf_hard_assignment = ebmf_hard_assignment,
      louvain_gamma = louvain_gamma,
      seed = seed,
      min_module_size = min_module_size,
      min_mean_internal = min_mean_internal,
      min_connectedness = min_connectedness
    )
  ))
}


.empty_module_quality <- function() {
  return(data.frame(
    cluster = integer(0),
    n_snps = integer(0),
    mean_internal_similarity = numeric(0),
    mean_external_similarity = numeric(0),
    separation = numeric(0),
    connectedness = numeric(0),
    n_internal_edges = integer(0),
    n_internal_pairs = integer(0),
    n_components = integer(0),
    largest_component_frac = numeric(0),
    mean_silhouette = numeric(0),
    reliable = logical(0),
    stringsAsFactors = FALSE
  ))
}


.cluster_snp_profiles_ebmf <- function(x_matrix,
                                       greedy_Kmax = 50L,
                                       lfsr_threshold = 0.05,
                                       magnitude_threshold = 0.25,
                                       drop_global = TRUE,
                                       prior = "point_normal",
                                       backfit = TRUE,
                                       hard_assignment = TRUE,
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
      cluster = stats::setNames(integer(0), character(0)),
      method = "ebmf",
      n_clusters = 0L,
      similarity_threshold = NULL,
      gamma = NA_real_,
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

  snp_ids <- colnames(x_matrix)
  cluster <- setNames(rep(NA_integer_, length(snp_ids)), snp_ids)
  membership <- extracted$membership
  if (hard_assignment && !is.null(membership) && ncol(membership) > 0 &&
      nrow(membership) > 0) {
    loadings <- abs(flash_fit$F_pm[rownames(membership), , drop = FALSE])
    loadings[!membership] <- -Inf
    best <- max.col(loadings, ties.method = "first")
    has_program <- rowSums(membership) > 0
    assigned <- rownames(membership)[has_program]
    cluster[assigned] <- best[has_program]
  }

  cluster <- cluster[!is.na(cluster)]
  cluster_ids <- sort(unique(cluster))
  cluster_program_map <- data.frame(
    cluster = integer(0),
    program = integer(0),
    stringsAsFactors = FALSE
  )
  if (length(cluster_ids) > 0) {
    reindex <- setNames(seq_along(cluster_ids), as.character(cluster_ids))
    cluster_program_map <- data.frame(
      cluster = seq_along(cluster_ids),
      program = as.integer(cluster_ids),
      stringsAsFactors = FALSE
    )
    cluster <- stats::setNames(
      as.integer(reindex[as.character(cluster)]),
      names(cluster)
    )
  }

  return(list(
    cluster = cluster,
    method = "ebmf",
    n_clusters = length(unique(cluster)),
    similarity_threshold = NULL,
    gamma = NA_real_,
    details = list(
      flash_fit = flash_fit,
      membership = membership,
      n_programs = extracted$n_programs,
      n_multi_program = extracted$n_multi_program,
      hard_assignment = hard_assignment,
      cluster_program_map = cluster_program_map,
      dropped_global_factors = dropped_global
    )
  ))
}
