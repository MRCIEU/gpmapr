#' @title Orient Pleiotropy Matrix Relative to Target Trait
#' @description Re-scale each SNP column of a pleiotropy matrix so that profiles
#' reflect what happens to the phenome when the variant increases liability for
#' the target trait. Allele-swapping noise is removed by multiplying each column
#' by the sign of the target trait's z-score at that SNP:
#' \deqn{X^*_{tj} = \mathrm{sign}(z_j) \cdot X_{tj}}
#' @param x_matrix A numeric matrix (traits x SNPs), such as
#'   `run_univariate_clustering()$x_matrix`.
#' @param target_trait_id Numeric ID of the target trait used to define column signs.
#' @param z_target Optional named numeric vector of target-trait z-scores, one per
#'   SNP column. If NULL, extracted from the row of `x_matrix` matching
#'   `target_trait_id`.
#' @return A list with:
#'   \itemize{
#'     \item x_matrix: the oriented matrix \eqn{X^*}
#'     \item z_target: target-trait z-scores per SNP column
#'     \item target_signs: directional signs (\code{+1}, \code{-1}, or \code{0})
#'       applied to each column
#'   }
#' @keywords internal
orient_pleiotropy_matrix <- function(x_matrix, target_trait_id, z_target = NULL) {
  if (missing(target_trait_id) || is.null(target_trait_id)) {
    stop("target_trait_id is required")
  }

  if (is.null(z_target)) {
    target_row <- as.character(target_trait_id)
    if (!target_row %in% rownames(x_matrix)) {
      stop("target_trait_id not found in x_matrix rownames")
    }
    z_target <- x_matrix[target_row, , drop = TRUE]
  }

  if (length(z_target) != ncol(x_matrix)) {
    stop("length of z_target must match ncol(x_matrix)")
  }

  if (!is.null(colnames(x_matrix))) {
    if (is.null(names(z_target))) {
      names(z_target) <- colnames(x_matrix)
    } else if (!identical(names(z_target), colnames(x_matrix))) {
      stop("names of z_target must match colnames(x_matrix)")
    }
  }

  target_signs <- sign(z_target)
  x_oriented <- sweep(x_matrix, 2, target_signs, `*`)

  return(list(
    x_matrix = x_oriented,
    z_target = z_target,
    target_signs = target_signs
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
#' @param asinh_scale Positive scale for the \code{asinh} transformation. The
#'   transformation is \code{asinh_scale * asinh(x / asinh_scale)}; larger
#'   values apply gentler compression. Ignored when \code{method = "none"}.
#' @return Matrix of the same dimensions / dimnames with finite entries
#'   transformed (NAs retained).
#' @export
compress_effect_matrix <- function(x_matrix,
                                   method = c("asinh", "none"),
                                   asinh_scale = 1) {
  method <- match.arg(method)
  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }
  if (!is.numeric(asinh_scale) || length(asinh_scale) != 1L ||
      !is.finite(asinh_scale) || asinh_scale <= 0) {
    stop("asinh_scale must be one positive finite number")
  }
  if (method == "none") {
    return(x_matrix)
  }

  out <- x_matrix
  obs <- is.finite(out)
  out[obs] <- asinh_scale * asinh(out[obs] / asinh_scale)
  return(out)
}


.finalize_pleiotropy_from_locus <- function(locus_data, target_id) {
  cg <- locus_data$cg
  if (!"tissue" %in% names(cg)) {
    cg$tissue <- NA_character_
  }
  if (!"gene" %in% names(cg)) {
    cg$gene <- NA_character_
  }

  # Sorted (C-locale) so the matrix does not depend on API row order: EBMF's
  # greedy initialisation is order-dependent.
  snp_ids <- sort(as.character(unique(locus_data$target_snps$snp_id)), method = "radix")
  built <- .build_study_level_pleiotropy_from_locus(cg, snp_ids)

  snp_info <- locus_data$target_snps |>
    dplyr::filter(snp_id %in% colnames(built$x_matrix)) |>
    dplyr::distinct() |>
    dplyr::arrange(snp_id)

  return(list(
    x_matrix = built$x_matrix,
    beta_matrix = if (!is.null(built$beta_matrix)) built$beta_matrix else NULL,
    se_matrix = if (!is.null(built$se_matrix)) built$se_matrix else NULL,
    trait_info = built$trait_info,
    snp_info = snp_info,
    target_trait_id = target_id
  ))
}


.build_study_level_pleiotropy_from_locus <- function(cg, snp_ids) {
  snp_ids <- as.character(snp_ids)
  has_beta_se <- all(c("beta", "se") %in% names(cg))
  z_long <- cg |>
    dplyr::mutate(snp_id = as.character(snp_id)) |>
    dplyr::group_by(trait_id, snp_id) |>
    dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select("trait_id", "trait_name", "snp_id", "z")

  z_wide <- z_long |>
    tidyr::pivot_wider(
      names_from = "snp_id",
      values_from = "z"
    ) |>
    dplyr::arrange(trait_id)

  if (!"gene_id" %in% names(cg)) {
    cg$gene_id <- NA_integer_
  }

  meta <- cg |>
    dplyr::mutate(
      trait_id = as.character(trait_id),
      is_molecular = !is.na(gene_id)
    ) |>
    dplyr::group_by(trait_id) |>
    dplyr::summarise(
      feature_type = if (any(is_molecular)) "molecular" else "phenotypic",
      gene = {
        vals <- as.character(gene)
        vals <- vals[!is.na(vals) & vals != ""]
        if (length(vals) == 0) NA_character_ else vals[[1]]
      },
      tissue = {
        vals <- as.character(tissue)
        vals <- vals[!is.na(vals) & vals != ""]
        if (length(vals) == 0) NA_character_ else vals[[1]]
      },
      .groups = "drop"
    )

  trait_info <- z_wide |>
    dplyr::transmute(
      trait_id = as.character(trait_id),
      trait_name = trait_name
    ) |>
    dplyr::distinct(trait_id, .keep_all = TRUE) |>
    dplyr::left_join(meta, by = "trait_id") |>
    dplyr::mutate(
      feature_type = dplyr::coalesce(feature_type, "phenotypic")
    ) |>
    dplyr::select("trait_id", "trait_name", "feature_type", "gene", "tissue")
  snp_cols <- intersect(snp_ids, setdiff(names(z_wide), c("trait_id", "trait_name")))
  x_matrix <- matrix(numeric(0), nrow = 0, ncol = length(snp_ids))
  dimnames(x_matrix) <- list(character(0), snp_ids)

  if (length(snp_cols) > 0) {
    x_matrix <- as.matrix(z_wide[, snp_cols, drop = FALSE])
    rownames(x_matrix) <- as.character(z_wide$trait_id)
    x_matrix <- .ensure_matrix_columns(x_matrix, snp_ids)
  }

  out <- list(x_matrix = x_matrix, trait_info = trait_info)

  if (has_beta_se && nrow(z_wide) > 0 && length(snp_cols) > 0) {
    value_long <- cg |>
      dplyr::mutate(snp_id = as.character(snp_id)) |>
      dplyr::group_by(trait_id, snp_id) |>
      dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::select("trait_id", "snp_id", "beta", "se") |>
      dplyr::distinct()
    beta_wide <- value_long |>
      dplyr::select("trait_id", "snp_id", "beta") |>
      tidyr::pivot_wider(names_from = "snp_id", values_from = "beta") |>
      dplyr::arrange(trait_id)
    se_wide <- value_long |>
      dplyr::select("trait_id", "snp_id", "se") |>
      tidyr::pivot_wider(names_from = "snp_id", values_from = "se") |>
      dplyr::arrange(trait_id)
    beta_cols <- intersect(snp_ids, setdiff(names(beta_wide), "trait_id"))
    beta_matrix <- as.matrix(beta_wide[, beta_cols, drop = FALSE])
    rownames(beta_matrix) <- as.character(beta_wide$trait_id)
    se_matrix <- as.matrix(se_wide[, beta_cols, drop = FALSE])
    rownames(se_matrix) <- as.character(se_wide$trait_id)
    out$beta_matrix <- .ensure_matrix_columns(beta_matrix, snp_ids)
    out$se_matrix <- .ensure_matrix_columns(se_matrix, snp_ids)
  }

  return(out)
}


.ensure_matrix_columns <- function(x_matrix, snp_ids) {
  snp_ids <- as.character(snp_ids)
  missing_cols <- setdiff(snp_ids, colnames(x_matrix))
  if (length(missing_cols) > 0) {
    na_cols <- matrix(
      NA_real_,
      nrow = nrow(x_matrix),
      ncol = length(missing_cols),
      dimnames = list(NULL, missing_cols)
    )
    x_matrix <- cbind(x_matrix, na_cols)
  }
  return(x_matrix[, snp_ids, drop = FALSE])
}


# Target SNPs (columns) are every target-trait coloc row with a usable beta/SE,
# named by variant_id; background rows are the other traits colocalising at
# those loci.
.prepare_pleiotropy_locus_data <- function(trait_id, coloc_groups) {
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("No coloc_groups data available")
  }

  if (!all(c("beta", "se") %in% names(coloc_groups))) {
    stop("coloc_groups must include beta and se columns (use include_associations = TRUE)")
  }

  if (!"variant_id" %in% names(coloc_groups)) {
    stop("coloc_groups must include column: variant_id")
  }

  # Bind locally so dplyr::filter does not treat both sides as the column
  target_id <- as.character(trait_id)
  coloc_groups <- coloc_groups |>
    dplyr::mutate(trait_id = as.character(trait_id))

  target_snps <- coloc_groups |>
    dplyr::filter(
      trait_id == target_id,
      !is.na(beta),
      !is.na(se),
      se > 0
    ) |>
    dplyr::mutate(snp_id = as.character(variant_id)) |>
    dplyr::distinct(coloc_group_id, snp_id, variant_id, display_snp, chr, bp)

  if (nrow(target_snps) == 0) {
    stop("No target-trait SNPs after filtering")
  }

  cg <- coloc_groups |>
    dplyr::filter(coloc_group_id %in% target_snps$coloc_group_id) |>
    dplyr::filter(!is.na(beta), !is.na(se), se > 0) |>
    dplyr::mutate(z = beta / se) |>
    dplyr::inner_join(
      target_snps |> dplyr::select("coloc_group_id", "snp_id"),
      by = "coloc_group_id"
    )

  z_target <- cg |>
    dplyr::filter(trait_id == target_id) |>
    dplyr::group_by(snp_id) |>
    dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select("snp_id", "z")

  return(list(
    target_snps = target_snps,
    cg = cg,
    z_target = z_target
  ))
}


.normalize_pleiotropy_columns <- function(x_matrix) {
  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }

  x <- x_matrix
  x[is.na(x)] <- 0

  col_norms <- sqrt(colSums(x^2))
  zero_cols <- col_norms == 0
  if (any(zero_cols)) {
    warning(
      "One or more SNP columns have zero norm; those columns are treated as zero vectors"
    )
  }

  col_norms_div <- col_norms
  col_norms_div[zero_cols] <- 1
  x_norm <- sweep(x, 2, col_norms_div, "/")
  if (!is.matrix(x_norm)) {
    x_norm <- matrix(x_norm, nrow = nrow(x), dimnames = dimnames(x))
  }
  if (any(zero_cols)) {
    col_mask <- rep(1, ncol(x_norm))
    col_mask[zero_cols] <- 0
    x_norm <- sweep(x_norm, 2, col_mask, "*")
  }

  return(list(
    x_norm = x_norm,
    col_norms = setNames(col_norms, colnames(x_matrix))
  ))
}


#' @title SNP–SNP Cosine Similarity Matrix
#' @description Compute pairwise cosine similarity between SNP pleiotropy
#' profiles. Missing trait-SNP values are treated as zero effects: each column
#' of \eqn{X^*} is normalised to unit length and similarity is obtained by
#' cross-product.
#' @param x_matrix A numeric matrix (traits x SNPs), typically the oriented matrix
#'   from `orient_pleiotropy_matrix()$x_matrix`.
#' @return A list with:
#'   \itemize{
#'     \item s_matrix: symmetric SNP-by-SNP cosine similarity matrix
#'     \item col_norms: Euclidean norm of each SNP column before normalisation
#'     \item overlap_matrix: number of jointly observed traits for each SNP pair
#'     \item eligible_matrix: whether both SNPs of a pair have a non-zero profile
#'   }
#' @export
snp_similarity_matrix <- function(x_matrix) {
  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }

  observed <- !is.na(x_matrix)
  overlap_matrix <- crossprod(observed * 1L)
  snp_ids <- colnames(x_matrix)
  dimnames(overlap_matrix) <- list(snp_ids, snp_ids)

  norm <- .normalize_pleiotropy_columns(x_matrix)
  s_matrix <- crossprod(norm$x_norm)
  col_norms <- norm$col_norms
  nonzero_cols <- col_norms > 0
  eligible_matrix <- outer(nonzero_cols, nonzero_cols, `&`)
  dimnames(eligible_matrix) <- list(snp_ids, snp_ids)

  if (!is.null(snp_ids)) {
    dimnames(s_matrix) <- list(snp_ids, snp_ids)
  }

  return(list(
    s_matrix = s_matrix,
    col_norms = col_norms,
    overlap_matrix = overlap_matrix,
    eligible_matrix = eligible_matrix
  ))
}
