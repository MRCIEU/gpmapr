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


