#' @title Build Pleiotropy Matrix
#' @description Construct a traits x SNPs matrix of signed z-scores (beta / se). Each
#' column is a target-trait SNP; each row is a background trait with an effect at that
#' SNP's locus.
#'
#' @param trait_id Numeric ID of the target trait whose associated SNPs define the columns.
#' @param coloc_groups A dataframe of coloc_groups with `beta` and `se` columns
#'   (from `trait(..., include_associations = TRUE)` or `traits(...)`). If NULL,
#'   will fetch via `trait(trait_id, include_associations = TRUE)`.
#' @param p_threshold P-value threshold for including a target-trait SNP. Defaults to 5e-8.
#' @param snp_key Column used to name SNP columns: `"variant_id"`, `"display_snp"`, or
#'   `"coloc_group_id"`. Defaults to `"variant_id"`.
#' @param associations Which association set to use for the matrix entries:
#'   `"coloc"` (default) uses colocalisation rows only (each trait-SNP pair keeps its
#'   smallest `min_p`), while `"full"` uses all full associations
#'   (`/associations-full`), keeping every observed effect per study-SNP pair.
#' @return A list with:
#'   \itemize{
#'     \item x_matrix: numeric matrix (traits x SNPs) of z-scores; `NA` where a trait
#'       has no effect at that SNP's locus
#'     \item trait_info: dataframe mapping row indices to `trait_id` and `trait_name`
#'       (plus `feature_type`, and `gene`/`tissue` when available)
#'     \item snp_info: dataframe mapping column names to `coloc_group_id`, `variant_id`,
#'       `display_snp`, `chr`, and `bp`
#'     \item target_trait_id: the target trait ID used
#'   }
#' @export
build_pleiotropy_matrix <- function(trait_id,
                                    coloc_groups = NULL,
                                    p_threshold = NULL,
                                    snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                    associations = c("coloc", "full")) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }

  snp_key <- match.arg(snp_key)
  associations <- match.arg(associations)
  target_id <- trait_id

  if (is.null(coloc_groups)) {
    coloc_groups <- trait(target_id, include_associations = TRUE)$coloc_groups
  }

  if (associations == "full") {
    full_associations <- trait(target_id, include_full_associations = TRUE)$full_associations
    locus_data <- .prepare_dense_locus_data(
      trait_id = target_id,
      coloc_groups = coloc_groups,
      full_associations = full_associations,
      p_threshold = p_threshold,
      snp_key = snp_key
    )
  } else {
    locus_data <- .prepare_pleiotropy_locus_data(
      trait_id = target_id,
      coloc_groups = coloc_groups,
      p_threshold = p_threshold,
      snp_key = snp_key
    )
  }

  return(.finalize_pleiotropy_from_locus(locus_data, target_id))
}


#' @title Orient Pleiotropy Matrix Relative to Target Trait
#' @description Re-scale each SNP column of a pleiotropy matrix so that profiles
#' reflect what happens to the phenome when the variant increases liability for
#' the target trait. Allele-swapping noise is removed by multiplying each column
#' by the sign of the target trait's z-score at that SNP:
#' \deqn{X^*_{tj} = \mathrm{sign}(z_j) \cdot X_{tj}}
#' @param x_matrix A numeric matrix (traits x SNPs) as returned by
#'   `build_pleiotropy_matrix()$x_matrix`.
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
#' @export
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


.finalize_pleiotropy_from_locus <- function(locus_data, target_id) {
  cg <- locus_data$cg
  if (!"tissue" %in% names(cg)) {
    cg$tissue <- NA_character_
  }
  if (!"gene" %in% names(cg)) {
    cg$gene <- NA_character_
  }

  snp_ids <- as.character(unique(locus_data$target_snps$snp_id))
  built <- .build_study_level_pleiotropy_from_locus(cg, snp_ids)

  snp_info <- locus_data$target_snps |>
    dplyr::filter(snp_id %in% colnames(built$x_matrix)) |>
    dplyr::distinct()

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
    )

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
      tidyr::pivot_wider(names_from = "snp_id", values_from = "beta")
    se_wide <- value_long |>
      dplyr::select("trait_id", "snp_id", "se") |>
      tidyr::pivot_wider(names_from = "snp_id", values_from = "se")
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


.is_molecular_gene_row <- function(cg) {
  if (!"gene_id" %in% names(cg)) {
    return(rep(FALSE, nrow(cg)))
  }
  return(!is.na(cg$gene_id))
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


.prepare_pleiotropy_locus_data <- function(trait_id, coloc_groups, p_threshold, snp_key) {
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("No coloc_groups data available")
  }

  if (!all(c("beta", "se") %in% names(coloc_groups))) {
    stop("coloc_groups must include beta and se columns (use include_associations = TRUE)")
  }

  if (!snp_key %in% names(coloc_groups)) {
    stop("coloc_groups must include column: ", snp_key)
  }

  # Bind locally so dplyr::filter does not treat both sides as the column
  target_id <- as.character(trait_id)
  coloc_groups <- coloc_groups |>
    dplyr::mutate(trait_id = as.character(trait_id))

  target_snps <- coloc_groups |>
    dplyr::filter(
      trait_id == target_id,
      if (!is.null(p_threshold)) min_p <= p_threshold else TRUE,
      !is.na(beta),
      !is.na(se),
      se > 0
    ) |>
    dplyr::mutate(snp_id = as.character(.data[[snp_key]])) |>
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


.prepare_dense_locus_data <- function(trait_id, coloc_groups, full_associations, p_threshold, snp_key) {
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("No coloc_groups data available")
  }
  if (is.null(full_associations) || nrow(full_associations) == 0) {
    stop("full_associations is required (use include_full_associations = TRUE)")
  }
  if (!snp_key %in% names(coloc_groups)) {
    stop("coloc_groups must include column: ", snp_key)
  }

  target_id <- as.character(trait_id)
  coloc_groups <- coloc_groups |>
    dplyr::mutate(trait_id = as.character(trait_id))

  target_snps <- coloc_groups |>
    dplyr::filter(
      trait_id == target_id,
      if (!is.null(p_threshold)) min_p <= p_threshold else TRUE
    ) |>
    dplyr::mutate(snp_id = as.character(.data[[snp_key]])) |>
    dplyr::distinct(coloc_group_id, snp_id, variant_id, display_snp, chr, bp)

  if (nrow(target_snps) == 0) {
    stop("No target-trait SNPs after filtering")
  }

  if (!"tissue" %in% names(coloc_groups)) {
    coloc_groups$tissue <- NA_character_
  }

  study_map <- coloc_groups |>
    dplyr::filter(!is.na(study_id)) |>
    dplyr::distinct(study_id, trait_id, trait_name, gene_id, gene, tissue)

  if ("existing_study_id" %in% names(coloc_groups)) {
    existing_map <- coloc_groups |>
      dplyr::filter(!is.na(existing_study_id)) |>
      dplyr::transmute(
        study_id = existing_study_id,
        trait_id = trait_id,
        trait_name = trait_name,
        gene_id = gene_id,
        gene = gene,
        tissue = tissue
      ) |>
      dplyr::distinct()
    study_map <- dplyr::bind_rows(study_map, existing_map) |>
      dplyr::distinct(study_id, trait_id, trait_name, gene_id, gene, tissue)
  }

  cg <- full_associations |>
    dplyr::filter(!is.na(beta), !is.na(se), se > 0) |>
    dplyr::inner_join(study_map, by = "study_id") |>
    dplyr::inner_join(
      target_snps |> dplyr::select("variant_id", "snp_id", "coloc_group_id"),
      by = "variant_id"
    ) |>
    dplyr::mutate(z = beta / se, min_p = p)

  z_target <- cg |>
    dplyr::filter(trait_id == target_id) |>
    dplyr::mutate(snp_id = as.character(snp_id)) |>
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


.normalize_pleiotropy_columns <- function(x_matrix, na_as_zero = TRUE) {
  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }

  x <- x_matrix
  if (na_as_zero) {
    x[is.na(x)] <- 0
  } else if (anyNA(x)) {
    stop("x_matrix contains NA; set na_as_zero = TRUE or impute missing values")
  }

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


#' @title Summarise Phenotypic vs Molecular Matrix Contribution
#' @description Report how much of a traits-x-SNPs matrix is carried by phenotypic
#' versus molecular rows. Useful after sparse filtering to see whether clustering
#' profiles are driven by complex traits or QTL features.
#' @param x_matrix Numeric traits x SNPs matrix.
#' @param trait_info Dataframe with `trait_id` matching `rownames(x_matrix)` and
#'   `feature_type` values such as `"phenotypic"` / `"molecular"`.
#' @return A dataframe with one row per `feature_type`:
#'   \itemize{
#'     \item n_features: number of rows of that type
#'     \item n_nonzero_cells: non-`NA` cells
#'     \item mean_abs_z: mean absolute z among non-`NA` cells
#'     \item squared_mass: sum of squared z (NA treated as 0)
#'     \item frac_squared_mass: share of total squared mass
#'     \item mean_snp_frac_mass: mean over SNPs of that type's share of the
#'       column's squared mass (SNPs with zero total mass skipped)
#'   }
#' @export
summarise_feature_type_contribution <- function(x_matrix, trait_info) {
  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }
  if (is.null(rownames(x_matrix))) {
    stop("x_matrix must have rownames matching trait_info$trait_id")
  }
  if (is.null(trait_info) || nrow(trait_info) == 0) {
    stop("trait_info is required")
  }
  if (!all(c("trait_id", "feature_type") %in% names(trait_info))) {
    stop("trait_info must include trait_id and feature_type")
  }

  info <- trait_info |>
    dplyr::mutate(
      trait_id = as.character(trait_id),
      feature_type = as.character(feature_type)
    ) |>
    dplyr::distinct(trait_id, .keep_all = TRUE)

  row_ids <- rownames(x_matrix)
  types <- info$feature_type[match(row_ids, info$trait_id)]
  types[is.na(types)] <- "unknown"

  x_sq <- x_matrix
  x_sq[is.na(x_sq)] <- 0
  x_sq <- x_sq^2
  total_mass <- sum(x_sq)
  col_mass <- colSums(x_sq)

  type_levels <- unique(types)
  rows <- lapply(type_levels, function(ftype) {
    idx <- which(types == ftype)
    sub <- x_matrix[idx, , drop = FALSE]
    sub_sq <- x_sq[idx, , drop = FALSE]
    type_mass <- sum(sub_sq)
    type_col_mass <- colSums(sub_sq)
    keep_cols <- col_mass > 0
    mean_snp_frac <- if (any(keep_cols)) {
      mean(type_col_mass[keep_cols] / col_mass[keep_cols])
    } else {
      NA_real_
    }
    return(data.frame(
      feature_type = ftype,
      n_features = length(idx),
      n_nonzero_cells = sum(!is.na(sub)),
      mean_abs_z = mean(abs(sub), na.rm = TRUE),
      squared_mass = type_mass,
      frac_squared_mass = if (total_mass > 0) type_mass / total_mass else NA_real_,
      mean_snp_frac_mass = mean_snp_frac,
      stringsAsFactors = FALSE
    ))
  })

  out <- dplyr::bind_rows(rows) |>
    dplyr::arrange(dplyr::desc(frac_squared_mass))
  return(out)
}


#' @title SNP–SNP Cosine Similarity Matrix
#' @description Compute pairwise cosine similarity between SNP pleiotropy profiles.
#' With zero filling, each column of \eqn{X^*} is normalised to unit length and
#' similarity is obtained by cross-product. Pairwise-complete similarity instead
#' normalises each SNP pair over its jointly observed trait rows.
#' @param x_matrix A numeric matrix (traits x SNPs), typically the oriented matrix
#'   from `orient_pleiotropy_matrix()$x_matrix`.
#' @param na_as_zero Logical; treat `NA` entries as zero before normalisation.
#'   Used only when `missing_method = "zero_fill"`. Defaults to `TRUE`.
#' @param missing_method How to handle missing trait-SNP values. `"zero_fill"`
#'   (default) preserves the original behaviour and treats missing values as zero
#'   effects. `"pairwise_complete"` calculates each SNP-pair cosine using only
#'   traits observed for both SNPs.
#' @param min_overlap Minimum number of jointly observed traits required when
#'   `missing_method = "pairwise_complete"`. Pairs below this threshold receive
#'   similarity zero (an absent graph edge). Defaults to `1`.
#' @param overlap_shrinkage Non-negative shrinkage strength for pairwise-complete
#'   similarity. Each cosine is multiplied by
#'   \eqn{n_{\mathrm{overlap}} / (n_{\mathrm{overlap}} + \lambda)}. Defaults to
#'   `0` (no shrinkage).
#' @return A list with:
#'   \itemize{
#'     \item s_matrix: symmetric SNP-by-SNP cosine similarity matrix
#'     \item col_norms: Euclidean norm of each SNP column before normalisation
#'     \item overlap_matrix: number of jointly observed traits for each SNP pair
#'     \item eligible_matrix: whether each pair met the overlap and norm criteria
#'     \item missing_method: missing-value method used
#'     \item min_overlap, overlap_shrinkage: pairwise-complete settings used
#'   }
#' @export
snp_similarity_matrix <- function(x_matrix,
                                  na_as_zero = TRUE,
                                  missing_method = c("zero_fill", "pairwise_complete"),
                                  min_overlap = 1L,
                                  overlap_shrinkage = 0) {
  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }

  missing_method <- match.arg(missing_method)
  if (
    length(min_overlap) != 1L ||
      !is.numeric(min_overlap) ||
      !is.finite(min_overlap) ||
      min_overlap < 1 ||
      min_overlap != floor(min_overlap) ||
      min_overlap > .Machine$integer.max
  ) {
    stop("min_overlap must be a positive integer")
  }
  if (
    length(overlap_shrinkage) != 1L ||
      !is.numeric(overlap_shrinkage) ||
      !is.finite(overlap_shrinkage) ||
      overlap_shrinkage < 0
  ) {
    stop("overlap_shrinkage must be a non-negative number")
  }

  observed <- !is.na(x_matrix)
  overlap_matrix <- crossprod(observed * 1L)
  snp_ids <- colnames(x_matrix)
  dimnames(overlap_matrix) <- list(snp_ids, snp_ids)

  if (missing_method == "zero_fill") {
    norm <- .normalize_pleiotropy_columns(
      x_matrix,
      na_as_zero = na_as_zero
    )
    s_matrix <- crossprod(norm$x_norm)
    col_norms <- norm$col_norms
    nonzero_cols <- col_norms > 0
    eligible_matrix <- outer(nonzero_cols, nonzero_cols, `&`)
    dimnames(eligible_matrix) <- list(snp_ids, snp_ids)
  } else {
    x_zero <- x_matrix
    x_zero[is.na(x_zero)] <- 0
    observed_numeric <- observed * 1
    numerator <- crossprod(x_zero)
    shared_sum_squares <- crossprod(x_zero^2, observed_numeric)
    denominator <- sqrt(shared_sum_squares * t(shared_sum_squares))
    eligible_matrix <- overlap_matrix >= as.integer(min_overlap) &
      denominator > 0
    s_matrix <- matrix(
      0,
      nrow = ncol(x_matrix),
      ncol = ncol(x_matrix),
      dimnames = list(snp_ids, snp_ids)
    )
    s_matrix[eligible_matrix] <- numerator[eligible_matrix] /
      denominator[eligible_matrix]

    if (overlap_shrinkage > 0) {
      shrinkage <- overlap_matrix[eligible_matrix] /
        (overlap_matrix[eligible_matrix] + overlap_shrinkage)
      s_matrix[eligible_matrix] <- s_matrix[eligible_matrix] * shrinkage
    }

    s_matrix[] <- pmax(-1, pmin(1, s_matrix))
    col_norms <- sqrt(colSums(x_zero^2))
    names(col_norms) <- snp_ids
    nonzero_cols <- col_norms > 0
    diag(s_matrix)[nonzero_cols] <- 1
    diag(eligible_matrix) <- nonzero_cols
  }

  if (!is.null(snp_ids)) {
    dimnames(s_matrix) <- list(snp_ids, snp_ids)
  }

  return(list(
    s_matrix = s_matrix,
    col_norms = col_norms,
    overlap_matrix = overlap_matrix,
    eligible_matrix = eligible_matrix,
    missing_method = missing_method,
    min_overlap = as.integer(min_overlap),
    overlap_shrinkage = overlap_shrinkage
  ))
}
