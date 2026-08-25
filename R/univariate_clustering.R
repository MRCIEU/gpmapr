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
#' @param min_snp_signals Minimum non-`NA` SNPs per background trait row.
#' @param max_snp_fraction Maximum non-`NA` fraction per background trait row.
#' @param compress_method Effect compression passed to `compress_effect_matrix()`.
#' @param compress_scale Asinh scale passed to `compress_effect_matrix()`.
#' @param similarity_threshold Edge threshold applied before Louvain clustering.
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
#'     \item trait_info, snp_info: row / column metadata
#'     \item dropped_trait_ids: background traits removed by the sparse/ubiquitous filter
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
                                      min_snp_signals = 5L,
                                      max_snp_fraction = 0.8,
                                      compress_method = c("none", "asinh"),
                                      compress_scale = 2,
                                      similarity_threshold = 0,
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
  compress_method <- match.arg(compress_method)

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

  clusters <- cluster_snp_profiles_louvain(
    similarity$s_matrix,
    similarity_threshold = similarity_threshold,
    gamma = louvain_gamma,
    seed = seed
  )

  module_quality <- summarise_snp_module_quality(
    similarity$s_matrix,
    clusters$cluster,
    edge_threshold = similarity_threshold,
    min_module_size = min_module_size,
    min_mean_internal = min_mean_internal,
    min_connectedness = min_connectedness
  )
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
    module_quality = module_quality,
    trait_info = pleiotropy$trait_info,
    snp_info = pleiotropy$snp_info,
    dropped_trait_ids = drop_trait_ids,
    group_traits = group_traits,
    coloc_groups = coloc_groups,
    parameters = list(
      target_trait_id = target_id,
      associations = associations,
      p_threshold = p_threshold,
      snp_key = snp_key,
      include_trans = include_trans,
      min_snp_signals = min_snp_signals,
      max_snp_fraction = max_snp_fraction,
      compress_method = compress_method,
      compress_scale = compress_scale,
      similarity_threshold = similarity_threshold,
      louvain_gamma = louvain_gamma,
      seed = seed,
      min_module_size = min_module_size,
      min_mean_internal = min_mean_internal,
      min_connectedness = min_connectedness
    )
  ))
}
