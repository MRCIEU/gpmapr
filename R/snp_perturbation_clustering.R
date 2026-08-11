#' @title Build SNP Perturbation Feature Matrices
#' @description Construct the two feature matrices from the SNP perturbation
#' clustering framework: a SNP x complex-trait matrix of signed z-scores
#' (\eqn{\beta / \mathrm{SE}}) and a SNP x gene matrix of signed gene-perturbation
#' scores. Gene scores collapse all molecular QTL rows for each SNP-gene pair
#' with a weighted Stouffer combination
#' \eqn{Z = \sum_i w_i Z_i / \sqrt{\sum_i w_i^2}} where \eqn{w_i = 1/\mathrm{SE}_i}.
#' Colocalisation probability is not currently returned by GPMap association
#' tables, so inverse-SE weights are used instead.
#' @param trait_id Numeric ID of the target trait whose GWAS SNPs define the rows.
#' @param p_threshold P-value threshold for including a target-trait SNP.
#'   Defaults to 5e-8.
#' @param snp_key Column used to name SNP rows: `"variant_id"`, `"display_snp"`,
#'   or `"coloc_group_id"`. Defaults to `"variant_id"`.
#' @return A list with:
#'   \itemize{
#'     \item trait_matrix: SNPs x complex traits (signed z)
#'     \item gene_matrix: SNPs x genes (weighted Stouffer z)
#'     \item trait_info, gene_info, snp_info
#'     \item z_target: named target-trait z-scores for orientation
#'     \item coloc_groups, target_trait_id
#'   }
#' @export
build_perturbation_matrices <- function(trait_id,
                                        p_threshold = NULL,
                                        snp_key = c(
                                          "variant_id",
                                          "display_snp",
                                          "coloc_group_id"
                                        )) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }
  snp_key <- match.arg(snp_key)
  target_id <- trait_id

  trait_data <- trait(target_id, include_full_associations = TRUE)
  coloc_groups <- trait_data$coloc_groups
  full_associations <- trait_data$full_associations

  locus_data <- .prepare_dense_locus_data(
    trait_id = target_id,
    coloc_groups = coloc_groups,
    full_associations = full_associations,
    p_threshold = p_threshold,
    snp_key = snp_key
  )

  snp_ids <- as.character(unique(locus_data$target_snps$snp_id))
  cg <- locus_data$cg
  molecular <- .is_molecular_gene_row(cg)

  trait_matrix <- .snp_feature_matrix_from_long(
    cg = cg[!molecular, , drop = FALSE],
    snp_ids = snp_ids,
    feature_id_col = "trait_id",
    feature_name_col = "trait_name",
    collapse = "best_p"
  )
  gene_matrix <- .snp_feature_matrix_from_long(
    cg = cg[molecular, , drop = FALSE],
    snp_ids = snp_ids,
    feature_id_col = "gene_id",
    feature_name_col = "gene",
    collapse = "stouffer"
  )

  snp_info <- locus_data$target_snps |>
    dplyr::filter(snp_id %in% snp_ids) |>
    dplyr::distinct()

  z_target <- stats::setNames(
    rep(NA_real_, length(snp_ids)),
    snp_ids
  )
  if (nrow(locus_data$z_target) > 0) {
    z_idx <- match(snp_ids, locus_data$z_target$snp_id)
    matched <- !is.na(z_idx)
    z_target[matched] <- locus_data$z_target$z[z_idx[matched]]
  }

  return(list(
    trait_matrix = trait_matrix$x_matrix,
    gene_matrix = gene_matrix$x_matrix,
    trait_info = trait_matrix$feature_info,
    gene_info = gene_matrix$feature_info,
    snp_info = snp_info,
    z_target = z_target,
    coloc_groups = coloc_groups,
    target_trait_id = target_id
  ))
}


#' @title Filter Sparse Perturbation Features
#' @description Drop trait/gene columns observed in fewer than `min_snps` SNPs.
#' Sparse gene columns dominate missingness patterns and collapse Leiden into
#' singletons when left unfiltered.
#' @param trait_matrix SNPs x traits matrix.
#' @param gene_matrix SNPs x genes matrix.
#' @param min_snps Minimum non-missing SNP count per feature column. Defaults to 5.
#' @return A list with filtered `trait_matrix`, `gene_matrix`, and counts of
#'   columns retained / dropped.
#' @export
filter_perturbation_features <- function(trait_matrix,
                                         gene_matrix,
                                         min_snps = 5L) {
  if (!is.matrix(trait_matrix) || !is.matrix(gene_matrix)) {
    stop("trait_matrix and gene_matrix must be matrices")
  }
  if (!identical(rownames(trait_matrix), rownames(gene_matrix))) {
    stop("trait_matrix and gene_matrix must have identical rownames (SNP ids)")
  }
  if (!is.numeric(min_snps) || min_snps < 1) {
    stop("min_snps must be a positive number")
  }

  trait_keep <- colSums(!is.na(trait_matrix)) >= min_snps
  gene_keep <- colSums(!is.na(gene_matrix)) >= min_snps

  return(list(
    trait_matrix = trait_matrix[, trait_keep, drop = FALSE],
    gene_matrix = gene_matrix[, gene_keep, drop = FALSE],
    n_traits_kept = sum(trait_keep),
    n_genes_kept = sum(gene_keep),
    n_traits_dropped = sum(!trait_keep),
    n_genes_dropped = sum(!gene_keep),
    min_snps = as.integer(min_snps)
  ))
}


#' @title Orient Perturbation Matrices To Target Trait
#' @description Multiply each SNP row by \eqn{\mathrm{sign}(z_{\mathrm{target}})}
#' so profiles reflect increasing target-trait liability (same role as
#' `orient_pleiotropy_matrix()` in the trait-matrix workflows).
#' @param trait_matrix SNPs x traits matrix.
#' @param gene_matrix SNPs x genes matrix.
#' @param z_target Named numeric vector of target-trait z-scores (names = SNP ids).
#' @return A list with oriented `trait_matrix`, `gene_matrix`, and `target_signs`.
#' @export
orient_perturbation_matrices <- function(trait_matrix, gene_matrix, z_target) {
  if (!is.matrix(trait_matrix) || !is.matrix(gene_matrix)) {
    stop("trait_matrix and gene_matrix must be matrices")
  }
  if (!identical(rownames(trait_matrix), rownames(gene_matrix))) {
    stop("trait_matrix and gene_matrix must have identical rownames (SNP ids)")
  }
  snp_ids <- rownames(trait_matrix)
  if (is.null(names(z_target))) {
    stop("z_target must be a named vector of SNP ids")
  }
  z <- as.numeric(z_target[snp_ids])
  signs <- sign(z)
  signs[!is.finite(signs) | signs == 0] <- 1
  names(signs) <- snp_ids

  trait_oriented <- sweep(trait_matrix, 1, signs, "*")
  gene_oriented <- sweep(gene_matrix, 1, signs, "*")
  if (!is.matrix(trait_oriented)) {
    trait_oriented <- matrix(
      trait_oriented,
      nrow = nrow(trait_matrix),
      dimnames = dimnames(trait_matrix)
    )
  }
  if (!is.matrix(gene_oriented)) {
    gene_oriented <- matrix(
      gene_oriented,
      nrow = nrow(gene_matrix),
      dimnames = dimnames(gene_matrix)
    )
  }

  return(list(
    trait_matrix = trait_oriented,
    gene_matrix = gene_oriented,
    target_signs = signs
  ))
}


#' @title Compress Heavy-Tailed Effect Scores
#' @description Soft-compress signed effect / z-score matrices so extreme values
#' (e.g. corrupt \eqn{|z| \approx 1000} or ultra-large-N GWAS hits) do not
#' dominate cosine similarity or mean-absolute rankings. Unlike hard clipping,
#' \code{asinh} preserves sign and rank order while smoothly shrinking tails:
#' \eqn{\mathrm{asinh}(z) = \log(z + \sqrt{z^2 + 1})}.
#' @param x_matrix Numeric matrix of signed effects / z-scores.
#' @param method Compression method. \code{"asinh"} (default) or \code{"none"}.
#' @return Matrix of the same dimensions / dimnames with finite entries
#'   transformed (NAs retained).
#' @export
compress_effect_matrix <- function(x_matrix, method = c("asinh", "none")) {
  method <- match.arg(method)
  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }
  if (method == "none") {
    return(x_matrix)
  }

  out <- x_matrix
  obs <- is.finite(out)
  out[obs] <- asinh(out[obs])
  return(out)
}


#' @title Compress Perturbation Feature Matrices
#' @description Apply \code{compress_effect_matrix()} to trait and gene matrices.
#' @inheritParams compress_effect_matrix
#' @param trait_matrix SNPs x traits matrix.
#' @param gene_matrix SNPs x genes matrix.
#' @return A list with compressed \code{trait_matrix}, \code{gene_matrix}, and
#'   \code{method}.
#' @export
compress_perturbation_matrices <- function(trait_matrix,
                                           gene_matrix,
                                           method = c("asinh", "none")) {
  method <- match.arg(method)
  if (!is.matrix(trait_matrix) || !is.matrix(gene_matrix)) {
    stop("trait_matrix and gene_matrix must be matrices")
  }
  if (!identical(rownames(trait_matrix), rownames(gene_matrix))) {
    stop("trait_matrix and gene_matrix must have identical rownames (SNP ids)")
  }
  return(list(
    trait_matrix = compress_effect_matrix(trait_matrix, method = method),
    gene_matrix = compress_effect_matrix(gene_matrix, method = method),
    method = method
  ))
}


#' @title Column-Wise Feature Scaling
#' @description Scale each feature column across SNPs. By default this **does not
#' centre** columns: mean-centring after filling missing associations with zero
#' turns shared missingness into strong artificial similarity and destroys
#' Leiden community structure on sparse GPMap matrices. Missing values are
#' treated as zero after computing the scale from observed entries.
#' @param x_matrix Numeric matrix (SNPs x features).
#' @param center Logical; if `TRUE`, subtract the column mean (including zeros
#'   for missing cells). Defaults to `FALSE`.
#' @return The column-scaled matrix with the same dimensions and dimnames.
#' @export
scale_feature_columns <- function(x_matrix, center = FALSE) {
  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }
  if (ncol(x_matrix) == 0 || nrow(x_matrix) == 0) {
    out <- x_matrix
    out[is.na(out)] <- 0
    return(out)
  }

  col_sds <- apply(x_matrix, 2, function(col) {
    vals <- col[is.finite(col)]
    if (length(vals) < 2) {
      return(1)
    }
    s <- stats::sd(vals)
    if (!is.finite(s) || s == 0) {
      return(1)
    }
    return(s)
  })

  x <- x_matrix
  x[is.na(x)] <- 0
  if (isTRUE(center)) {
    x <- sweep(x, 2, colMeans(x), "-")
  }
  out <- sweep(x, 2, col_sds, "/")
  if (!is.matrix(out)) {
    out <- matrix(out, nrow = nrow(x_matrix), dimnames = dimnames(x_matrix))
  }
  return(out)
}


#' @title Column-Wise Z-Score Normalisation
#' @description Alias for `scale_feature_columns(x_matrix, center = TRUE)`.
#' Prefer `scale_feature_columns(center = FALSE)` for sparse perturbation
#' matrices.
#' @inheritParams scale_feature_columns
#' @return The column-standardised matrix.
#' @export
zscore_feature_columns <- function(x_matrix) {
  return(scale_feature_columns(x_matrix, center = TRUE))
}


#' @title Weighted Trait And Gene SNP Similarity
#' @description Compute cosine similarity among SNP rows of the trait and gene
#' matrices separately, then combine:
#' \deqn{S = w_{\mathrm{gene}} S_{\mathrm{gene}} + w_{\mathrm{trait}} S_{\mathrm{trait}}}
#' @param trait_matrix SNPs x traits matrix (typically after
#'   `scale_feature_columns()`).
#' @param gene_matrix SNPs x genes matrix (typically after
#'   `scale_feature_columns()`).
#' @param w_gene Weight for gene similarity. Defaults to 0.3.
#' @param w_trait Weight for trait similarity. Defaults to 0.7.
#' @return A list with `s_matrix`, `s_trait`, `s_gene`, `w_gene`, and `w_trait`.
#' @export
weighted_snp_similarity <- function(trait_matrix,
                                    gene_matrix,
                                    w_gene = 0.3,
                                    w_trait = 0.7) {
  if (!is.matrix(trait_matrix) || !is.matrix(gene_matrix)) {
    stop("trait_matrix and gene_matrix must be matrices")
  }
  if (!identical(rownames(trait_matrix), rownames(gene_matrix))) {
    stop("trait_matrix and gene_matrix must have identical rownames (SNP ids)")
  }
  if (!is.numeric(w_gene) || !is.numeric(w_trait) || w_gene < 0 || w_trait < 0) {
    stop("w_gene and w_trait must be non-negative numbers")
  }
  if (w_gene + w_trait <= 0) {
    stop("w_gene + w_trait must be positive")
  }

  s_trait <- .row_cosine_similarity(trait_matrix)
  s_gene <- .row_cosine_similarity(gene_matrix)
  s_matrix <- w_gene * s_gene + w_trait * s_trait
  dimnames(s_matrix) <- list(rownames(trait_matrix), rownames(trait_matrix))

  return(list(
    s_matrix = s_matrix,
    s_trait = s_trait,
    s_gene = s_gene,
    w_gene = w_gene,
    w_trait = w_trait
  ))
}


#' @title Cluster SNPs From A Perturbation Similarity Matrix
#' @description Cluster SNPs using the combined SNP x SNP similarity matrix from
#' `weighted_snp_similarity()`. Supports hierarchical clustering and Leiden
#' community detection via \pkg{igraph} (`igraph::cluster_leiden()`).
#' @param s_matrix Symmetric SNP x SNP similarity matrix.
#' @param method `"hierarchical"` or `"leiden"`.
#' @param k Number of clusters for hierarchical clustering. Ignored for Leiden.
#' @param linkage Linkage method for `hclust()`. Defaults to `"average"`.
#' @param similarity_threshold Minimum similarity for Leiden edges. Weaker
#'   edges are dropped. Defaults to 0.2.
#' @param resolution Leiden resolution parameter. Defaults to 0.5.
#' @return A list with `cluster` (named SNP assignments), `method`, `n_clusters`,
#'   and method-specific `details`.
#' @export
cluster_perturbation_snps <- function(s_matrix,
                                      method = c("hierarchical", "leiden"),
                                      k = 5L,
                                      linkage = "average",
                                      similarity_threshold = 0.2,
                                      resolution = 0.5) {
  method <- match.arg(method)
  if (!is.matrix(s_matrix)) {
    stop("s_matrix must be a matrix")
  }
  if (nrow(s_matrix) != ncol(s_matrix)) {
    stop("s_matrix must be square")
  }

  snp_ids <- rownames(s_matrix)
  if (is.null(snp_ids)) {
    snp_ids <- as.character(seq_len(nrow(s_matrix)))
  }

  if (method == "hierarchical") {
    if (is.null(k) || k < 1) {
      stop("k must be a positive integer for hierarchical clustering")
    }
    dist_matrix <- stats::as.dist(1 - s_matrix)
    hc <- stats::hclust(dist_matrix, method = linkage)
    cluster <- stats::cutree(hc, k = k)
    names(cluster) <- snp_ids
    return(list(
      cluster = cluster,
      method = "hierarchical",
      n_clusters = length(unique(cluster)),
      k = as.integer(k),
      details = list(hclust = hc, linkage = linkage)
    ))
  }

  return(.cluster_perturbation_leiden(
    s_matrix = s_matrix,
    snp_ids = snp_ids,
    similarity_threshold = similarity_threshold,
    resolution = resolution
  ))
}


.snp_feature_matrix_from_long <- function(cg,
                                          snp_ids,
                                          feature_id_col,
                                          feature_name_col,
                                          collapse = c("best_p", "stouffer")) {
  collapse <- match.arg(collapse)
  snp_ids <- as.character(snp_ids)
  empty <- list(
    x_matrix = matrix(
      numeric(0),
      nrow = length(snp_ids),
      ncol = 0,
      dimnames = list(snp_ids, NULL)
    ),
    feature_info = data.frame(
      feature_id = character(0),
      feature_name = character(0),
      stringsAsFactors = FALSE
    )
  )

  if (is.null(cg) || nrow(cg) == 0) {
    return(empty)
  }
  if (!all(c(feature_id_col, feature_name_col, "snp_id", "z") %in% names(cg))) {
    stop("cg is missing required columns for feature matrix construction")
  }

  feature_ids <- cg[[feature_id_col]]
  feature_names <- cg[[feature_name_col]]
  long <- cg |>
    dplyr::mutate(
      snp_id = as.character(snp_id),
      feature_id = as.character(feature_ids),
      feature_name = as.character(feature_names),
      z = as.numeric(z),
      se = if ("se" %in% names(cg)) as.numeric(se) else NA_real_,
      min_p = if ("min_p" %in% names(cg)) as.numeric(min_p) else NA_real_
    ) |>
    dplyr::filter(
      snp_id %in% snp_ids,
      !is.na(feature_id),
      feature_id != "",
      !is.na(z)
    )

  if (nrow(long) == 0) {
    return(empty)
  }

  if (collapse == "best_p") {
    scored <- long |>
      dplyr::group_by(snp_id, feature_id, feature_name) |>
      dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::select(snp_id, feature_id, feature_name, z)
  } else {
    scored <- long |>
      dplyr::filter(!is.na(se), se > 0) |>
      dplyr::group_by(snp_id, feature_id, feature_name) |>
      dplyr::summarise(
        z = {
          w <- 1 / se
          sum(w * z) / sqrt(sum(w * w))
        },
        .groups = "drop"
      )
  }

  if (nrow(scored) == 0) {
    return(empty)
  }

  feature_info <- scored |>
    dplyr::group_by(feature_id) |>
    dplyr::summarise(
      feature_name = dplyr::first(feature_name[!is.na(feature_name)]),
      .groups = "drop"
    ) |>
    dplyr::arrange(feature_id)

  wide <- scored |>
    dplyr::select(snp_id, feature_id, z) |>
    tidyr::pivot_wider(
      names_from = feature_id,
      values_from = z,
      values_fn = mean
    )

  x_matrix <- matrix(
    NA_real_,
    nrow = length(snp_ids),
    ncol = nrow(feature_info),
    dimnames = list(snp_ids, feature_info$feature_id)
  )
  row_idx <- match(wide$snp_id, snp_ids)
  feature_cols <- setdiff(names(wide), "snp_id")
  col_idx <- match(feature_cols, feature_info$feature_id)
  keep <- !is.na(row_idx)
  if (any(keep)) {
    x_matrix[row_idx[keep], col_idx] <- as.matrix(
      wide[keep, feature_cols, drop = FALSE]
    )
  }

  return(list(x_matrix = x_matrix, feature_info = feature_info))
}


.row_cosine_similarity <- function(x_matrix) {
  x <- x_matrix
  x[is.na(x)] <- 0
  if (ncol(x) == 0) {
    n <- nrow(x)
    out <- matrix(0, nrow = n, ncol = n, dimnames = list(rownames(x), rownames(x)))
    diag(out) <- 1
    return(out)
  }

  row_norms <- sqrt(rowSums(x * x))
  zero_rows <- row_norms == 0
  row_norms_div <- row_norms
  row_norms_div[zero_rows] <- 1
  x_norm <- sweep(x, 1, row_norms_div, "/")
  if (!is.matrix(x_norm)) {
    x_norm <- matrix(x_norm, nrow = nrow(x), dimnames = dimnames(x))
  }
  if (any(zero_rows)) {
    x_norm[zero_rows, ] <- 0
  }
  s_matrix <- tcrossprod(x_norm)
  if (any(zero_rows)) {
    s_matrix[zero_rows, ] <- 0
    s_matrix[, zero_rows] <- 0
    diag(s_matrix)[zero_rows] <- 1
  }
  return(s_matrix)
}


.cluster_perturbation_leiden <- function(s_matrix,
                                         snp_ids,
                                         similarity_threshold,
                                         resolution) {
  if (!is.null(similarity_threshold) && similarity_threshold > 0) {
    w <- s_matrix
    w[abs(w) < similarity_threshold] <- 0
  } else {
    w <- s_matrix
  }
  w[w < 0] <- 0
  diag(w) <- 0

  graph <- igraph::graph_from_adjacency_matrix(
    w,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  if (igraph::ecount(graph) == 0) {
    cluster <- stats::setNames(seq_along(snp_ids), snp_ids)
    return(list(
      cluster = cluster,
      method = "leiden",
      n_clusters = length(cluster),
      details = list(
        note = "No edges after thresholding; each SNP is its own cluster"
      )
    ))
  }

  communities <- igraph::cluster_leiden(
    graph,
    resolution = resolution
  )
  membership <- igraph::membership(communities)
  member_names <- names(membership)
  if (is.null(member_names)) {
    member_names <- snp_ids
  }
  cluster <- stats::setNames(as.integer(membership), member_names)

  return(list(
    cluster = cluster,
    method = "leiden",
    n_clusters = length(unique(cluster)),
    details = list(
      communities = communities,
      resolution = resolution,
      similarity_threshold = similarity_threshold
    )
  ))
}
