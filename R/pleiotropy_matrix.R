#' @title Build Pleiotropy Matrix
#' @description Construct a traits x SNPs matrix of signed z-scores (beta / se) from
#' colocalisation data. Each column is a target-trait SNP; each row is a background
#' trait with a colocalisation signal at that SNP's locus.
#'
#' @param trait_id Numeric ID of the target trait whose associated SNPs define the columns.
#' @param coloc_groups A dataframe of coloc_groups with `beta` and `se` columns
#'   (from `trait(..., include_associations = TRUE)` or `traits(...)`). If NULL,
#'   will fetch via `trait(trait_id, include_associations = TRUE)`.
#' @param p_threshold P-value threshold for including a target-trait SNP. Defaults to 5e-8.
#' @param snp_key Column used to name SNP columns: `"variant_id"`, `"display_snp"`, or
#'   `"coloc_group_id"`. Defaults to `"variant_id"`.
#' @param collapse_gene_tissue If `TRUE`, molecular QTL rows are collapsed by
#'   gene x tissue (Stouffer over cis/trans and QTL types). If `FALSE` (default),
#'   each QTL study remains its own trait row.
#' @return A list with:
#'   \itemize{
#'     \item x_matrix: numeric matrix (traits x SNPs) of z-scores; `NA` where a trait
#'       has no colocalisation signal at that SNP's locus
#'     \item trait_info: dataframe mapping row indices to `trait_id` and `trait_name`
#'       (plus `feature_type`, and `gene`/`tissue` when available)
#'     \item snp_info: dataframe mapping column names to `coloc_group_id`, `variant_id`,
#'       `display_snp`, `chr`, and `bp`
#'     \item target_trait_id: the target trait ID used
#'     \item collapse_gene_tissue: whether gene x tissue collapse was applied
#'   }
#' @export
build_pleiotropy_matrix <- function(trait_id,
                                    coloc_groups = NULL,
                                    p_threshold = NULL,
                                    snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                    collapse_gene_tissue = FALSE) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }

  snp_key <- match.arg(snp_key)
  if (!is.logical(collapse_gene_tissue) || length(collapse_gene_tissue) != 1L) {
    stop("collapse_gene_tissue must be a single logical")
  }
  target_id <- trait_id

  if (is.null(coloc_groups)) {
    coloc_groups <- trait(target_id, include_associations = TRUE)$coloc_groups
  }

  locus_data <- .prepare_pleiotropy_locus_data(
    trait_id = target_id,
    coloc_groups = coloc_groups,
    p_threshold = p_threshold,
    snp_key = snp_key
  )

  return(.finalize_pleiotropy_from_locus(
    locus_data,
    target_id,
    collapse_gene_tissue = collapse_gene_tissue
  ))
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


#' @title Build Bivariate Pleiotropy Matrices
#' @description Construct two pleiotropy matrices for a trait pair, with columns
#' defined by each trait's independent SNPs and rows aligned to the union of
#' background traits pleiotropically associated with either target.
#' @param trait_id_1 Numeric ID of the first target trait (columns of \eqn{X_1}).
#' @param trait_id_2 Numeric ID of the second target trait (columns of \eqn{X_2}).
#' @param coloc_groups A dataframe of coloc_groups with `beta` and `se` columns
#'   covering both traits. If NULL, fetched via
#'   `traits(c(trait_id_1, trait_id_2), include_associations = TRUE)` (coloc mode)
#'   or per-trait `trait(..., include_full_associations = TRUE)` (full mode).
#' @param p_threshold P-value threshold for including target-trait SNPs. If NULL,
#'   all SNPs with valid z-scores are retained.
#' @param snp_key Column used to name SNP columns. Defaults to `"variant_id"`.
#' @param association_source `"coloc"` (default): sparse matrix from coloc-merged
#'   associations only. `"full"`: denser matrix from `/associations-full` at the
#'   same target SNP columns (still defined by coloc loci).
#' @param full_associations_1,full_associations_2 Optional full-association
#'   tables for each target when `association_source = "full"`. Fetched if NULL.
#' @param collapse_gene_tissue If `TRUE`, molecular QTL rows are collapsed by
#'   gene x tissue (Stouffer over cis/trans and QTL types). If `FALSE` (default),
#'   each QTL study remains its own trait row.
#' @return A list with:
#'   \itemize{
#'     \item x1_matrix, x2_matrix: aligned numeric matrices (traits x SNPs)
#'     \item trait_info: dataframe of shared background trait metadata
#'       (`trait_id`, `trait_name`, `feature_type`, optional `gene`/`tissue`)
#'     \item snp_info_1, snp_info_2: SNP metadata for each matrix
#'     \item trait_id_1, trait_id_2: target trait IDs
#'     \item association_source: the source used
#'     \item collapse_gene_tissue: whether gene x tissue collapse was applied
#'   }
#' @export
build_bivariate_pleiotropy_matrices <- function(trait_id_1,
                                                trait_id_2,
                                                coloc_groups = NULL,
                                                p_threshold = NULL,
                                                snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                                association_source = c("coloc", "full"),
                                                full_associations_1 = NULL,
                                                full_associations_2 = NULL,
                                                collapse_gene_tissue = FALSE) {
  if (missing(trait_id_1) || is.null(trait_id_1) ||
      missing(trait_id_2) || is.null(trait_id_2)) {
    stop("trait_id_1 and trait_id_2 are required")
  }

  snp_key <- match.arg(snp_key)
  association_source <- match.arg(association_source)
  if (!is.logical(collapse_gene_tissue) || length(collapse_gene_tissue) != 1L) {
    stop("collapse_gene_tissue must be a single logical")
  }

  if (identical(association_source, "full")) {
    if (is.null(coloc_groups) || is.null(full_associations_1) ||
        is.null(full_associations_2)) {
      trait_data_1 <- trait(
        trait_id_1,
        include_associations = TRUE,
        include_full_associations = TRUE
      )
      trait_data_2 <- trait(
        trait_id_2,
        include_associations = TRUE,
        include_full_associations = TRUE
      )
      if (is.null(coloc_groups)) {
        coloc_groups <- dplyr::bind_rows(
          trait_data_1$coloc_groups,
          trait_data_2$coloc_groups
        )
      }
      if (is.null(full_associations_1)) {
        full_associations_1 <- trait_data_1$full_associations
      }
      if (is.null(full_associations_2)) {
        full_associations_2 <- trait_data_2$full_associations
      }
    }
    pleiotropy_1 <- .build_pleiotropy_matrix_dense(
      trait_id = trait_id_1,
      coloc_groups = coloc_groups,
      full_associations = full_associations_1,
      p_threshold = p_threshold,
      snp_key = snp_key,
      collapse_gene_tissue = collapse_gene_tissue
    )
    pleiotropy_2 <- .build_pleiotropy_matrix_dense(
      trait_id = trait_id_2,
      coloc_groups = coloc_groups,
      full_associations = full_associations_2,
      p_threshold = p_threshold,
      snp_key = snp_key,
      collapse_gene_tissue = collapse_gene_tissue
    )
  } else {
    if (is.null(coloc_groups)) {
      coloc_groups <- traits(
        c(trait_id_1, trait_id_2),
        include_associations = TRUE
      )$coloc_groups
    }
    pleiotropy_1 <- build_pleiotropy_matrix(
      trait_id = trait_id_1,
      coloc_groups = coloc_groups,
      p_threshold = p_threshold,
      snp_key = snp_key,
      collapse_gene_tissue = collapse_gene_tissue
    )
    pleiotropy_2 <- build_pleiotropy_matrix(
      trait_id = trait_id_2,
      coloc_groups = coloc_groups,
      p_threshold = p_threshold,
      snp_key = snp_key,
      collapse_gene_tissue = collapse_gene_tissue
    )
  }

  shared_trait_ids <- union(
    rownames(pleiotropy_1$x_matrix),
    rownames(pleiotropy_2$x_matrix)
  )
  shared_trait_ids <- union(
    shared_trait_ids,
    as.character(c(trait_id_1, trait_id_2))
  )

  trait_info <- dplyr::bind_rows(
    pleiotropy_1$trait_info,
    pleiotropy_2$trait_info
  ) |>
    dplyr::mutate(trait_id = as.character(trait_id)) |>
    dplyr::distinct(trait_id, .keep_all = TRUE) |>
    dplyr::filter(trait_id %in% shared_trait_ids)

  missing_trait_ids <- setdiff(shared_trait_ids, trait_info$trait_id)
  if (length(missing_trait_ids) > 0) {
    trait_info <- dplyr::bind_rows(
      trait_info,
      data.frame(
        trait_id = as.character(missing_trait_ids),
        trait_name = NA_character_,
        feature_type = "phenotypic",
        gene = NA_character_,
        tissue = NA_character_,
        stringsAsFactors = FALSE
      )
    )
  }

  return(list(
    x1_matrix = .align_pleiotropy_rows(pleiotropy_1$x_matrix, shared_trait_ids),
    x2_matrix = .align_pleiotropy_rows(pleiotropy_2$x_matrix, shared_trait_ids),
    trait_info = trait_info,
    snp_info_1 = pleiotropy_1$snp_info,
    snp_info_2 = pleiotropy_2$snp_info,
    trait_id_1 = trait_id_1,
    trait_id_2 = trait_id_2,
    association_source = association_source,
    collapse_gene_tissue = collapse_gene_tissue,
    coloc_groups = coloc_groups
  ))
}


.build_pleiotropy_matrix_dense <- function(trait_id,
                                           coloc_groups,
                                           full_associations,
                                           p_threshold,
                                           snp_key,
                                           collapse_gene_tissue = FALSE) {
  locus_data <- .prepare_dense_locus_data(
    trait_id = trait_id,
    coloc_groups = coloc_groups,
    full_associations = full_associations,
    p_threshold = p_threshold,
    snp_key = snp_key
  )
  return(.finalize_pleiotropy_from_locus(
    locus_data,
    trait_id,
    collapse_gene_tissue = collapse_gene_tissue
  ))
}


.finalize_pleiotropy_from_locus <- function(locus_data,
                                            target_id,
                                            collapse_gene_tissue = FALSE) {
  cg <- locus_data$cg
  if (!"tissue" %in% names(cg)) {
    cg$tissue <- NA_character_
  }
  if (!"gene" %in% names(cg)) {
    cg$gene <- NA_character_
  }

  snp_ids <- as.character(unique(locus_data$target_snps$snp_id))
  if (isTRUE(collapse_gene_tissue)) {
    built <- .build_gene_tissue_pleiotropy_from_locus(cg, snp_ids)
    x_matrix <- built$x_matrix
    trait_info <- built$trait_info
  } else {
    built <- .build_study_level_pleiotropy_from_locus(cg, snp_ids)
    x_matrix <- built$x_matrix
    trait_info <- built$trait_info
  }

  snp_info <- locus_data$target_snps |>
    dplyr::filter(snp_id %in% colnames(x_matrix)) |>
    dplyr::distinct()

  return(list(
    x_matrix = x_matrix,
    trait_info = trait_info,
    snp_info = snp_info,
    target_trait_id = target_id,
    collapse_gene_tissue = collapse_gene_tissue
  ))
}


.build_study_level_pleiotropy_from_locus <- function(cg, snp_ids) {
  snp_ids <- as.character(snp_ids)
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

  return(list(x_matrix = x_matrix, trait_info = trait_info))
}


.build_gene_tissue_pleiotropy_from_locus <- function(cg, snp_ids) {
  snp_ids <- as.character(snp_ids)
  molecular <- .is_molecular_gene_row(cg)
  trait_cg <- cg[!molecular, , drop = FALSE]
  mol_cg <- cg[molecular, , drop = FALSE]

  trait_part <- .phenotypic_trait_z_to_matrix(trait_cg, snp_ids)
  gene_tissue_part <- .gene_tissue_z_to_matrix(
    .collapse_gene_tissue_z_scores(mol_cg),
    snp_ids
  )

  parts <- list()
  if (nrow(trait_part$x_matrix) > 0) {
    parts[[length(parts) + 1L]] <- trait_part
  }
  if (nrow(gene_tissue_part$x_matrix) > 0) {
    parts[[length(parts) + 1L]] <- gene_tissue_part
  }
  if (length(parts) == 0) {
    stop("No phenotypic or gene x tissue rows could be constructed")
  }

  x_matrix <- do.call(rbind, lapply(parts, function(p) {
    return(.ensure_matrix_columns(p$x_matrix, snp_ids))
  }))
  trait_info <- dplyr::bind_rows(lapply(parts, `[[`, "trait_info"))

  return(list(x_matrix = x_matrix, trait_info = trait_info))
}


.phenotypic_trait_z_to_matrix <- function(cg, snp_ids) {
  snp_ids <- as.character(snp_ids)
  empty <- list(
    x_matrix = matrix(
      numeric(0),
      nrow = 0,
      ncol = length(snp_ids),
      dimnames = list(character(0), snp_ids)
    ),
    trait_info = data.frame(
      trait_id = character(0),
      trait_name = character(0),
      feature_type = character(0),
      gene = character(0),
      tissue = character(0),
      stringsAsFactors = FALSE
    )
  )
  if (is.null(cg) || nrow(cg) == 0) {
    return(empty)
  }

  z_long <- cg |>
    dplyr::mutate(snp_id = as.character(snp_id)) |>
    dplyr::group_by(trait_id, trait_name, snp_id) |>
    dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select("trait_id", "trait_name", "snp_id", "z")

  z_wide <- z_long |>
    tidyr::pivot_wider(
      names_from = "snp_id",
      values_from = "z",
      values_fn = mean
    )

  trait_info <- z_wide |>
    dplyr::transmute(
      trait_id = as.character(trait_id),
      trait_name = trait_name,
      feature_type = "phenotypic",
      gene = NA_character_,
      tissue = NA_character_
    ) |>
    dplyr::distinct(trait_id, .keep_all = TRUE)

  snp_cols <- intersect(snp_ids, setdiff(names(z_wide), c("trait_id", "trait_name")))
  x_matrix <- empty$x_matrix
  if (length(snp_cols) > 0) {
    x_matrix <- as.matrix(z_wide[, snp_cols, drop = FALSE])
    rownames(x_matrix) <- as.character(z_wide$trait_id)
    x_matrix <- .ensure_matrix_columns(x_matrix, snp_ids)
  }

  return(list(x_matrix = x_matrix, trait_info = trait_info))
}


.collapse_gene_tissue_z_scores <- function(cg) {
  if (is.null(cg) || nrow(cg) == 0) {
    return(data.frame(
      gene_id = integer(0),
      gene = character(0),
      tissue = character(0),
      snp_id = character(0),
      z = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  if (!"tissue" %in% names(cg)) {
    cg$tissue <- NA_character_
  }
  if (!"se" %in% names(cg)) {
    cg$se <- NA_real_
  }
  if (!"min_p" %in% names(cg)) {
    cg$min_p <- NA_real_
  }

  long <- cg |>
    dplyr::filter(!is.na(gene_id), !is.na(z)) |>
    dplyr::mutate(
      snp_id = as.character(snp_id),
      gene = ifelse(
        is.na(gene) | gene == "",
        as.character(gene_id),
        as.character(gene)
      ),
      tissue = ifelse(
        is.na(tissue) | tissue == "",
        "unknown",
        as.character(tissue)
      ),
      z = as.numeric(z),
      se = as.numeric(se),
      min_p = as.numeric(min_p)
    )

  if (nrow(long) == 0) {
    return(long[, c("gene_id", "gene", "tissue", "snp_id", "z"), drop = FALSE])
  }

  with_se <- long |> dplyr::filter(!is.na(se), se > 0)
  without_se <- long |> dplyr::filter(is.na(se) | se <= 0)

  stouffer_part <- data.frame(
    gene_id = integer(0),
    gene = character(0),
    tissue = character(0),
    snp_id = character(0),
    z = numeric(0),
    stringsAsFactors = FALSE
  )
  if (nrow(with_se) > 0) {
    stouffer_part <- with_se |>
      dplyr::group_by(gene_id, gene, tissue, snp_id) |>
      dplyr::summarise(
        z = {
          w <- 1 / se
          sum(w * z) / sqrt(sum(w * w))
        },
        .groups = "drop"
      )
  }

  best_p_part <- data.frame(
    gene_id = integer(0),
    gene = character(0),
    tissue = character(0),
    snp_id = character(0),
    z = numeric(0),
    stringsAsFactors = FALSE
  )
  if (nrow(without_se) > 0) {
    best_p_part <- without_se |>
      dplyr::group_by(gene_id, gene, tissue, snp_id) |>
      dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::select("gene_id", "gene", "tissue", "snp_id", "z")
  }

  # Prefer Stouffer when SE exists; fill remaining gene x tissue x snp from best-p
  covered <- stouffer_part |>
    dplyr::transmute(key = paste(gene_id, tissue, snp_id, sep = "\r"))
  best_p_part <- best_p_part |>
    dplyr::mutate(key = paste(gene_id, tissue, snp_id, sep = "\r")) |>
    dplyr::filter(!key %in% covered$key) |>
    dplyr::select(-"key")

  return(dplyr::bind_rows(stouffer_part, best_p_part))
}


.gene_tissue_z_to_matrix <- function(gene_tissue_z, snp_ids) {
  snp_ids <- as.character(snp_ids)
  empty <- list(
    x_matrix = matrix(
      numeric(0),
      nrow = 0,
      ncol = length(snp_ids),
      dimnames = list(character(0), snp_ids)
    ),
    trait_info = data.frame(
      trait_id = character(0),
      trait_name = character(0),
      feature_type = character(0),
      gene = character(0),
      tissue = character(0),
      stringsAsFactors = FALSE
    )
  )
  if (is.null(gene_tissue_z) || nrow(gene_tissue_z) == 0) {
    return(empty)
  }

  long <- gene_tissue_z |>
    dplyr::mutate(
      snp_id = as.character(snp_id),
      feature_id = paste0("gene:", gene_id, "|", tissue),
      feature_name = paste0(gene, " | ", tissue),
      feature_type = "molecular",
      gene = as.character(gene),
      tissue = as.character(tissue)
    )

  feature_wide <- long |>
    tidyr::pivot_wider(
      id_cols = c("feature_id", "feature_name", "feature_type", "gene", "tissue"),
      names_from = "snp_id",
      values_from = "z",
      values_fn = mean
    )

  trait_info <- feature_wide |>
    dplyr::transmute(
      trait_id = feature_id,
      trait_name = feature_name,
      feature_type = feature_type,
      gene = gene,
      tissue = tissue
    )

  snp_cols <- intersect(
    snp_ids,
    setdiff(
      names(feature_wide),
      c("feature_id", "feature_name", "feature_type", "gene", "tissue")
    )
  )
  x_matrix <- empty$x_matrix
  if (length(snp_cols) > 0) {
    x_matrix <- as.matrix(feature_wide[, snp_cols, drop = FALSE])
    rownames(x_matrix) <- feature_wide$feature_id
    x_matrix <- .ensure_matrix_columns(x_matrix, snp_ids)
  }

  return(list(x_matrix = x_matrix, trait_info = trait_info))
}


.is_molecular_gene_row <- function(cg) {
  if (!"gene_id" %in% names(cg)) {
    return(rep(FALSE, nrow(cg)))
  }
  !is.na(cg$gene_id)
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
  target_id <- trait_id

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

  target_id <- trait_id

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


.align_pleiotropy_rows <- function(x_matrix, trait_ids) {
  aligned <- matrix(
    NA_real_,
    nrow = length(trait_ids),
    ncol = ncol(x_matrix),
    dimnames = list(trait_ids, colnames(x_matrix))
  )
  shared_rows <- intersect(rownames(x_matrix), trait_ids)
  if (length(shared_rows) > 0) {
    aligned[shared_rows, ] <- x_matrix[shared_rows, , drop = FALSE]
  }
  return(aligned)
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


#' @title Cross-Trait SNP Cosine Similarity Matrix
#' @description Compute rectangular cosine similarity between SNP profiles from two
#' oriented pleiotropy matrices sharing the same background trait rows.
#' \deqn{S_{jk} = \frac{{x^*_{1,j}}^{\top} x^*_{2,k}}{\|x^*_{1,j}\| \|x^*_{2,k}\|}}
#' Columns of both matrices are L2-normalised, then
#' \eqn{S = X_1^{*\top} X_2^*}.
#' @param x1_matrix Oriented pleiotropy matrix for trait 1 (traits x SNPs).
#' @param x2_matrix Oriented pleiotropy matrix for trait 2 (traits x SNPs).
#' @param na_as_zero Logical; treat `NA` entries as zero before normalisation.
#' @return A list with:
#'   \itemize{
#'     \item s_matrix: rectangular Trait-1 SNPs x Trait-2 SNPs similarity matrix
#'     \item col_norms_1, col_norms_2: column norms before normalisation
#'   }
#' @export
cross_snp_similarity_matrix <- function(x1_matrix, x2_matrix, na_as_zero = TRUE) {
  if (!identical(rownames(x1_matrix), rownames(x2_matrix))) {
    stop("x1_matrix and x2_matrix must have identical rownames")
  }
  if (nrow(x1_matrix) != nrow(x2_matrix)) {
    stop("x1_matrix and x2_matrix must have the same number of rows")
  }

  norm_1 <- .normalize_pleiotropy_columns(x1_matrix, na_as_zero = na_as_zero)
  norm_2 <- .normalize_pleiotropy_columns(x2_matrix, na_as_zero = na_as_zero)

  s_matrix <- crossprod(norm_1$x_norm, norm_2$x_norm)
  if (!is.null(colnames(x1_matrix)) && !is.null(colnames(x2_matrix))) {
    dimnames(s_matrix) <- list(colnames(x1_matrix), colnames(x2_matrix))
  }

  return(list(
    s_matrix = s_matrix,
    col_norms_1 = norm_1$col_norms,
    col_norms_2 = norm_2$col_norms
  ))
}


#' @title Cluster Cross-Trait SNPs From a Similarity Matrix
#' @description Cluster Trait-1 and Trait-2 SNPs independently from a rectangular
#' cross-similarity matrix (`S^*`). Hierarchical clustering uses Euclidean distance
#' across cross-trait similarity profiles. Signed Louvain clusters cosine similarity
#' between those profiles (rows of `S^*` for Trait-1 SNPs; columns for Trait-2 SNPs).
#' @param s_matrix Rectangular cross-trait similarity matrix from
#'   `cross_snp_similarity_matrix()$s_matrix`.
#' @param method Clustering method: `"hierarchical"` or `"louvain"`.
#' @param linkage Linkage method passed to `stats::hclust()` for hierarchical
#'   clustering and for dendrogram ordering when `method = "louvain"`.
#'   Defaults to `"average"`.
#' @param k1 Optional number of row (Trait-1 SNP) clusters when
#'   `method = "hierarchical"`.
#' @param k2 Optional number of column (Trait-2 SNP) clusters when
#'   `method = "hierarchical"`.
#' @param similarity_threshold Minimum absolute profile similarity for off-diagonal
#'   SNP pairs when `method = "louvain"`. Weaker edges are zeroed before signed
#'   Louvain. Set to `NULL` or `0` to use the full matrix. Defaults to `NULL`.
#' @param gamma Resolution parameter for signed Louvain when `method = "louvain"`.
#'   Values greater than 1 tend to yield more, smaller modules. Defaults to `1`.
#' @return A list with:
#'   \itemize{
#'     \item method: clustering method used
#'     \item hclust_rows: hierarchical clustering of Trait-1 SNPs (ordering)
#'     \item hclust_cols: hierarchical clustering of Trait-2 SNPs (ordering)
#'     \item clusters_rows: row cluster assignments (`k1` for hierarchical;
#'       always returned for Louvain)
#'     \item clusters_cols: column cluster assignments (`k2` for hierarchical;
#'       always returned for Louvain)
#'     \item profile_similarity_rows, profile_similarity_cols: square profile
#'       similarity matrices when `method = "louvain"`
#'     \item louvain_rows, louvain_cols: full output from `cluster_snp_profiles()`
#'       when `method = "louvain"`
#'   }
#' @export
cluster_cross_trait_snps <- function(s_matrix,
                                     method = c("hierarchical", "louvain"),
                                     linkage = "average",
                                     k1 = NULL,
                                     k2 = NULL,
                                     similarity_threshold = NULL,
                                     gamma = 1) {
  method <- match.arg(method)

  if (!is.matrix(s_matrix)) {
    stop("s_matrix must be a matrix")
  }
  if (!is.null(similarity_threshold) && similarity_threshold < 0) {
    stop("similarity_threshold must be NULL, 0, or a non-negative value")
  }
  if (gamma <= 0) {
    stop("gamma must be positive")
  }

  if (method == "hierarchical") {
    return(.cluster_cross_trait_hierarchical(
      s_matrix = s_matrix,
      linkage = linkage,
      k1 = k1,
      k2 = k2
    ))
  }

  profile_rows <- .profile_similarity_matrix(s_matrix, margin = "rows")
  profile_cols <- .profile_similarity_matrix(s_matrix, margin = "cols")

  louvain_rows <- cluster_snp_profiles(
    profile_rows,
    method = "community",
    community_algorithm = "louvain",
    similarity_threshold = similarity_threshold,
    gamma = gamma
  )
  louvain_cols <- cluster_snp_profiles(
    profile_cols,
    method = "community",
    community_algorithm = "louvain",
    similarity_threshold = similarity_threshold,
    gamma = gamma
  )

  dist_rows <- stats::dist(1 - profile_rows, method = "euclidean")
  dist_cols <- stats::dist(1 - profile_cols, method = "euclidean")
  hclust_rows <- stats::hclust(dist_rows, method = linkage)
  hclust_cols <- stats::hclust(dist_cols, method = linkage)

  return(list(
    method = "louvain",
    hclust_rows = hclust_rows,
    hclust_cols = hclust_cols,
    clusters_rows = louvain_rows$cluster,
    clusters_cols = louvain_cols$cluster,
    profile_similarity_rows = profile_rows,
    profile_similarity_cols = profile_cols,
    louvain_rows = louvain_rows,
    louvain_cols = louvain_cols,
    similarity_threshold = similarity_threshold,
    gamma = gamma
  ))
}


.cluster_cross_trait_hierarchical <- function(s_matrix, linkage, k1, k2) {
  s <- s_matrix
  s[!is.finite(s)] <- 0

  dist_rows <- stats::dist(s, method = "euclidean")
  hclust_rows <- stats::hclust(dist_rows, method = linkage)

  dist_cols <- stats::dist(t(s), method = "euclidean")
  hclust_cols <- stats::hclust(dist_cols, method = linkage)

  out <- list(
    method = "hierarchical",
    hclust_rows = hclust_rows,
    hclust_cols = hclust_cols
  )

  if (!is.null(k1)) {
    out$clusters_rows <- stats::cutree(hclust_rows, k = k1)
  }
  if (!is.null(k2)) {
    out$clusters_cols <- stats::cutree(hclust_cols, k = k2)
  }

  return(out)
}


.normalize_pleiotropy_rows <- function(x_matrix, na_as_zero = TRUE) {
  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }

  x <- x_matrix
  if (na_as_zero) {
    x[is.na(x)] <- 0
  } else if (anyNA(x)) {
    stop("x_matrix contains NA; set na_as_zero = TRUE or impute missing values")
  }

  row_norms <- sqrt(rowSums(x^2))
  zero_rows <- row_norms == 0
  if (any(zero_rows)) {
    warning(
      "One or more rows have zero norm; those rows are treated as zero vectors"
    )
  }

  row_norms_div <- row_norms
  row_norms_div[zero_rows] <- 1
  x_norm <- sweep(x, 1, row_norms_div, "/")
  if (!is.matrix(x_norm)) {
    x_norm <- matrix(x_norm, nrow = nrow(x), dimnames = dimnames(x))
  }
  if (any(zero_rows)) {
    row_mask <- rep(1, nrow(x_norm))
    row_mask[zero_rows] <- 0
    x_norm <- sweep(x_norm, 1, row_mask, "*")
  }

  return(list(
    x_norm = x_norm,
    row_norms = stats::setNames(row_norms, rownames(x_matrix))
  ))
}


.profile_similarity_matrix <- function(x_matrix,
                                       margin = c("rows", "cols"),
                                       na_as_zero = TRUE) {
  margin <- match.arg(margin)

  if (margin == "rows") {
    norm <- .normalize_pleiotropy_rows(x_matrix, na_as_zero = na_as_zero)
    s_matrix <- tcrossprod(norm$x_norm)
    snp_ids <- rownames(x_matrix)
  } else {
    norm <- .normalize_pleiotropy_columns(x_matrix, na_as_zero = na_as_zero)
    s_matrix <- crossprod(norm$x_norm)
    snp_ids <- colnames(x_matrix)
  }

  if (!is.null(snp_ids)) {
    dimnames(s_matrix) <- list(snp_ids, snp_ids)
  }

  return(s_matrix)
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


#' @title Summarise Cross-Trait SNP Modules
#' @description Summarise rectangular blocks from a clustered cross-trait similarity
#' matrix and classify each block as concordant, discordant, or trait-specific.
#' @param s_matrix Rectangular cross-trait similarity matrix (`S^*`).
#' @param clusters_rows Named integer vector of Trait-1 SNP cluster assignments.
#' @param clusters_cols Named integer vector of Trait-2 SNP cluster assignments.
#' @param concordant_threshold Mean similarity above which a block is classified as
#'   shared concordant (Block A). Defaults to `0.5`.
#' @param discordant_threshold Mean similarity below which a block is classified as
#'   shared discordant (Block B). Defaults to `-0.5`.
#' @param specific_threshold Absolute mean similarity below which a block is classified
#'   as trait-specific (Block C). Defaults to `0.2`.
#' @return A dataframe with one row per row-cluster/column-cluster block, including
#'   `mean_similarity`, SNP counts, and `block_type`.
#' @export
summarise_cross_trait_modules <- function(s_matrix,
                                          clusters_rows,
                                          clusters_cols,
                                          concordant_threshold = 0.5,
                                          discordant_threshold = -0.5,
                                          specific_threshold = 0.2) {
  if (length(clusters_rows) != nrow(s_matrix)) {
    stop("length(clusters_rows) must match nrow(s_matrix)")
  }
  if (length(clusters_cols) != ncol(s_matrix)) {
    stop("length(clusters_cols) must match ncol(s_matrix)")
  }

  if (is.null(names(clusters_rows)) && !is.null(rownames(s_matrix))) {
    names(clusters_rows) <- rownames(s_matrix)
  }
  if (is.null(names(clusters_cols)) && !is.null(colnames(s_matrix))) {
    names(clusters_cols) <- colnames(s_matrix)
  }

  row_cluster_ids <- sort(unique(clusters_rows))
  col_cluster_ids <- sort(unique(clusters_cols))

  blocks <- lapply(row_cluster_ids, function(row_cluster) {
    row_snps <- names(clusters_rows)[clusters_rows == row_cluster]
    lapply(col_cluster_ids, function(col_cluster) {
      col_snps <- names(clusters_cols)[clusters_cols == col_cluster]
      vals <- s_matrix[row_snps, col_snps, drop = FALSE]
      mean_sim <- mean(vals, na.rm = TRUE)

      if (mean_sim >= concordant_threshold) {
        block_type <- "concordant"
      } else if (mean_sim <= discordant_threshold) {
        block_type <- "discordant"
      } else if (abs(mean_sim) < specific_threshold) {
        block_type <- "trait_specific"
      } else {
        block_type <- "mixed"
      }

      data.frame(
        row_cluster = row_cluster,
        col_cluster = col_cluster,
        n_snps_trait_1 = length(row_snps),
        n_snps_trait_2 = length(col_snps),
        mean_similarity = mean_sim,
        min_similarity = min(vals, na.rm = TRUE),
        max_similarity = max(vals, na.rm = TRUE),
        block_type = block_type,
        stringsAsFactors = FALSE
      )
    })
  })

  result <- dplyr::bind_rows(unlist(blocks, recursive = FALSE)) |>
    dplyr::arrange(dplyr::desc(abs(mean_similarity)))

  return(result)
}


#' @title SNP–SNP Cosine Similarity Matrix
#' @description Compute pairwise cosine similarity between SNP pleiotropy profiles.
#' Each column of the oriented matrix \eqn{X^*} is normalised to unit length, then
#' \deqn{S_{jk} = \frac{{x^*_j}^{\top} x^*_k}{\|x^*_j\| \|x^*_k\|}}
#' is obtained as \eqn{{X^*}^\top X^*} on the unit-normalised columns.
#' @param x_matrix A numeric matrix (traits x SNPs), typically the oriented matrix
#'   from `orient_pleiotropy_matrix()$x_matrix`.
#' @param na_as_zero Logical; treat `NA` entries as zero before normalisation.
#'   Defaults to `TRUE` (no colocalisation signal at a trait is treated as zero effect).
#' @return A list with:
#'   \itemize{
#'     \item s_matrix: symmetric SNP-by-SNP cosine similarity matrix
#'     \item col_norms: Euclidean norm of each SNP column before normalisation
#'   }
#' @export
snp_similarity_matrix <- function(x_matrix, na_as_zero = TRUE) {
  norm <- .normalize_pleiotropy_columns(x_matrix, na_as_zero = na_as_zero)
  s_matrix <- crossprod(norm$x_norm)
  if (!is.null(colnames(x_matrix))) {
    dimnames(s_matrix) <- list(colnames(x_matrix), colnames(x_matrix))
  }

  return(list(
    s_matrix = s_matrix,
    col_norms = norm$col_norms
  ))
}


#' @title Cluster SNPs by Pleiotropy Profile Similarity
#' @description Group SNPs into functional classes based on geometric parallelism of
#' their multi-trait profiles. Supports hierarchical clustering on \eqn{1 - S},
#' spectral clustering on the similarity graph, community detection, and Gaussian
#' mixture models in a latent space derived from the oriented pleiotropy matrix.
#' @param s_matrix Symmetric SNP-by-SNP cosine similarity matrix from
#'   `snp_similarity_matrix()$s_matrix`.
#' @param k Number of clusters (biological modules). For signed Louvain community
#'   detection, `k` is ignored and the number of modules is data-driven.
#' @param method Clustering method: `"hierarchical"`, `"spectral"`, `"community"`,
#'   or `"gmm"`.
#' @param x_matrix Oriented pleiotropy matrix (traits x SNPs). Required for
#'   `method = "gmm"`; typically `orient_pleiotropy_matrix()$x_matrix`.
#' @param linkage Linkage method passed to `stats::hclust()` when
#'   `method = "hierarchical"`. Defaults to `"average"`.
#' @param community_algorithm Community detection algorithm when
#'   `method = "community"`: `"leading_eigen"` (recursive splits via `steps = k - 1`
#'   on the positive similarity subgraph; requires **igraph**) or `"louvain"`
#'   (signed Louvain maximising Gomez signed modularity on the full cosine
#'   similarity matrix, including negative anti-parallel edges).
#' @param similarity_threshold Minimum absolute similarity for off-diagonal SNP pairs.
#'   Weaker edges are zeroed before clustering to reduce background connectivity.
#'   Set to `NULL` or `0` to use the full matrix. Defaults to `NULL`.
#' @param gamma Resolution parameter for signed Louvain (`method = "community"`,
#'   `community_algorithm = "louvain"`). Values greater than 1 tend to yield more,
#'   smaller modules. Defaults to `1`.
#' @param na_as_zero Treat `NA` entries as zero when building latent coordinates
#'   for GMM clustering.
#' @param n_latent Number of latent dimensions for GMM clustering. Defaults to
#'   `min(k + 1, n_traits - 1, n_snps - 1)`.
#' @return A list with:
#'   \itemize{
#'     \item cluster: named integer vector of cluster assignments per SNP
#'     \item method: clustering method used
#'     \item k: requested number of clusters (`NA` for Louvain)
#'     \item n_clusters: number of clusters returned
#'     \item similarity_threshold: edge threshold applied before clustering
#'     \item gamma: Louvain resolution used (NA for other methods)
#'     \item details: method-specific objects (e.g. `hclust`, `kmeans`, `mclust`)
#'   }
#' @export
cluster_snp_profiles <- function(s_matrix,
                                 k = 3,
                                 method = c("hierarchical", "spectral", "community", "gmm"),
                                 x_matrix = NULL,
                                 linkage = "average",
                                 community_algorithm = c("leading_eigen", "louvain"),
                                 similarity_threshold = NULL,
                                 gamma = 1,
                                 na_as_zero = TRUE,
                                 n_latent = NULL) {
  method <- match.arg(method)
  community_algorithm <- match.arg(community_algorithm)

  if (k < 2 && method != "community") {
    stop("k must be at least 2")
  }
  if (method == "gmm" && is.null(x_matrix)) {
    stop("x_matrix is required when method = 'gmm'")
  }
  if (!is.null(similarity_threshold) && similarity_threshold < 0) {
    stop("similarity_threshold must be NULL, 0, or a non-negative value")
  }
  if (gamma <= 0) {
    stop("gamma must be positive")
  }

  snp_ids <- colnames(s_matrix)
  if (is.null(snp_ids)) {
    snp_ids <- as.character(seq_len(ncol(s_matrix)))
  }

  S <- s_matrix
  dimnames(S) <- list(snp_ids, snp_ids)
  diag(S) <- 1
  S[!is.finite(S)] <- 0

  if (!is.null(similarity_threshold) && similarity_threshold > 0) {
    S <- .sparsify_similarity_matrix(S, threshold = similarity_threshold)
  }

  if (method == "hierarchical") {
    result <- .cluster_snps_hierarchical(S, k = k, linkage = linkage)
  } else if (method == "spectral") {
    result <- .cluster_snps_spectral(S, k = k)
  } else if (method == "community") {
    result <- .cluster_snps_community(
      S,
      k = k,
      algorithm = community_algorithm,
      gamma = gamma
    )
  } else {
    result <- .cluster_snps_gmm(
      x_matrix = x_matrix,
      k = k,
      na_as_zero = na_as_zero,
      n_latent = n_latent
    )
  }

  clusters <- result$cluster
  names(clusters) <- snp_ids

  return(list(
    cluster = clusters,
    method = method,
    k = if (method == "community" && community_algorithm == "louvain") NA_integer_ else k,
    n_clusters = length(unique(clusters)),
    similarity_threshold = similarity_threshold,
    gamma = if (method == "community" && community_algorithm == "louvain") gamma else NA_real_,
    details = result$details
  ))
}


.sparsify_similarity_matrix <- function(s_matrix, threshold) {
  n <- nrow(s_matrix)
  out <- s_matrix
  off_diag <- matrix(TRUE, n, n)
  diag(off_diag) <- FALSE
  out[off_diag & abs(out) < threshold] <- 0
  diag(out) <- 1
  return(out)
}


.cluster_snps_hierarchical <- function(s_matrix, k, linkage) {
  dist_matrix <- 1 - s_matrix
  dist_obj <- stats::as.dist(dist_matrix)
  hc <- stats::hclust(dist_obj, method = linkage)
  cluster <- stats::cutree(hc, k = k)

  return(list(
    cluster = cluster,
    details = list(
      hclust = hc,
      dist = dist_obj
    )
  ))
}


.cluster_snps_spectral <- function(s_matrix, k) {
  affinity <- (s_matrix + 1) / 2
  affinity[affinity < 0] <- 0
  diag(affinity) <- 0

  degree <- rowSums(affinity)
  degree[degree == 0] <- 1
  d_inv_sqrt <- diag(1 / sqrt(degree))
  laplacian <- diag(nrow(affinity)) - d_inv_sqrt %*% affinity %*% d_inv_sqrt

  eig <- eigen(laplacian, symmetric = TRUE)
  embedding <- eig$vectors[, seq_len(k), drop = FALSE]
  row_norms <- sqrt(rowSums(embedding^2))
  row_norms[row_norms == 0] <- 1
  embedding <- embedding / row_norms

  km <- stats::kmeans(embedding, centers = k, nstart = 25)

  return(list(
    cluster = km$cluster,
    details = list(
      embedding = embedding,
      kmeans = km,
      eigenvalues = eig$values[seq_len(k)]
    )
  ))
}


.igraph_from_signed_adjacency <- function(affinity) {
  adj <- affinity
  diag(adj) <- 0
  adj[!is.finite(adj)] <- 0

  presence <- abs(adj) > 0
  graph <- igraph::graph_from_adjacency_matrix(
    presence * 1,
    mode = "undirected",
    weighted = FALSE,
    diag = FALSE
  )

  if (igraph::ecount(graph) > 0) {
    edge_ends <- igraph::as_edgelist(graph)
    edge_weights <- adj[cbind(edge_ends[, 1], edge_ends[, 2])]
    igraph::E(graph)$weight <- edge_weights
  }

  graph
}


.cluster_louvain_signed <- function(s_matrix, gamma = 1, qtype = "sta", seed = NULL) {
  W <- s_matrix
  diag(W) <- 0
  W[!is.finite(W)] <- 0

  n <- nrow(W)
  W0 <- W * (W > 0)
  W1 <- -W * (W < 0)
  s0 <- sum(W0)
  s1 <- sum(W1)

  qtype <- match.arg(qtype, c("sta", "pos", "smp", "gja", "neg"))

  if (qtype == "smp") {
    d0 <- if (s0 > 0) 1 / s0 else 0
    d1 <- if (s1 > 0) 1 / s1 else 0
  } else if (qtype == "gja") {
    denom <- s0 + s1
    d0 <- if (denom > 0) 1 / denom else 0
    d1 <- d0
  } else if (qtype == "sta") {
    d0 <- if (s0 > 0) 1 / s0 else 0
    d1 <- if ((s0 + s1) > 0) 1 / (s0 + s1) else 0
  } else if (qtype == "pos") {
    d0 <- if (s0 > 0) 1 / s0 else 0
    d1 <- 0
  } else if (qtype == "neg") {
    d0 <- 0
    d1 <- if (s1 > 0) 1 / s1 else 0
  }

  if (s0 == 0) {
    s0 <- 1
    d0 <- 0
  }
  if (s1 == 0) {
    s1 <- 1
    d1 <- 0
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  h <- 1L
  nh <- n
  ci <- list(NULL, seq_len(n))
  q <- c(-1, 0)

  while (q[h + 1L] - q[h] > 1e-10) {
    if (h > 300L) {
      stop("Signed Louvain exceeded maximum hierarchy depth", call. = FALSE)
    }

    kn0 <- colSums(W0)
    kn1 <- colSums(W1)
    km0 <- kn0
    km1 <- kn1
    knm0 <- W0
    knm1 <- W1

    m <- seq_len(nh)
    flag <- TRUE
    it <- 0L

    while (flag) {
      it <- it + 1L
      if (it > 1000L) {
        stop("Signed Louvain iteration limit exceeded", call. = FALSE)
      }
      flag <- FALSE

      for (u in sample(nh)) {
        ma <- m[u]
        dQ0 <- (knm0[u, ] + W0[u, u] - knm0[u, ma]) -
          gamma * kn0[u] * (km0 + kn0[u] - km0[ma]) / s0
        dQ1 <- (knm1[u, ] + W1[u, u] - knm1[u, ma]) -
          gamma * kn1[u] * (km1 + kn1[u] - km1[ma]) / s1

        dQ <- d0 * dQ0 - d1 * dQ1
        dQ[ma] <- 0

        if (max(dQ) > 1e-10) {
          flag <- TRUE
          mb <- which.max(dQ)

          knm0[, mb] <- knm0[, mb] + W0[, u]
          knm0[, ma] <- knm0[, ma] - W0[, u]
          knm1[, mb] <- knm1[, mb] + W1[, u]
          knm1[, ma] <- knm1[, ma] - W1[, u]
          km0[mb] <- km0[mb] + kn0[u]
          km0[ma] <- km0[ma] - kn0[u]
          km1[mb] <- km1[mb] + kn1[u]
          km1[ma] <- km1[ma] - kn1[u]
          m[u] <- mb
        }
      }
    }

    h <- h + 1L
    ci[[h + 1L]] <- numeric(n)
    m_factor <- as.integer(as.factor(m))

    for (u in seq_len(nh)) {
      ci[[h + 1L]][ci[[h]] == u] <- m_factor[u]
    }

    nh <- max(m_factor)
    wn0 <- matrix(0, nh, nh)
    wn1 <- matrix(0, nh, nh)
    for (u in seq_len(nh)) {
      for (v in u:nh) {
        idx_u <- m_factor == u
        idx_v <- m_factor == v
        val0 <- sum(W0[idx_u, idx_v, drop = FALSE])
        val1 <- sum(W1[idx_u, idx_v, drop = FALSE])
        wn0[u, v] <- val0
        wn0[v, u] <- val0
        wn1[u, v] <- val1
        wn1[v, u] <- val1
      }
    }
    W0 <- wn0
    W1 <- wn1

    q0 <- sum(diag(W0)) - sum(W0 %*% W0) / s0
    q1 <- sum(diag(W1)) - sum(W1 %*% W1) / s1
    q <- c(q, d0 * q0 - d1 * q1)
  }

  membership <- as.integer(as.factor(ci[[h + 1L]]))

  list(
    cluster = membership,
    modularity = q[length(q)],
    qtype = qtype
  )
}


.cluster_snps_community <- function(s_matrix, k, algorithm, gamma = 1) {
  affinity <- s_matrix
  diag(affinity) <- 0
  affinity[!is.finite(affinity)] <- 0

  if (algorithm == "louvain") {
    signed <- .cluster_louvain_signed(affinity, gamma = gamma, seed = 1L)
    graph <- .igraph_from_signed_adjacency(affinity)

    return(list(
      cluster = as.integer(signed$cluster),
      details = list(
        igraph = graph,
        modularity = signed$modularity,
        qtype = signed$qtype,
        gamma = gamma,
        algorithm = "signed_louvain"
      )
    ))
  }

  affinity_pos <- affinity
  affinity_pos[affinity_pos < 0] <- 0

  graph <- igraph::graph_from_adjacency_matrix(
    affinity_pos,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  communities <- igraph::cluster_leading_eigen(
    graph,
    steps = k - 1L,
    weights = igraph::E(graph)$weight
  )
  membership <- igraph::membership(communities)

  return(list(
    cluster = as.integer(membership),
    details = list(
      igraph = graph,
      communities = communities
    )
  ))
}


.cluster_snps_gmm <- function(x_matrix, k, na_as_zero, n_latent) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required for GMM clustering", call. = FALSE)
  }

  x <- x_matrix
  if (na_as_zero) {
    x[is.na(x)] <- 0
  } else if (anyNA(x)) {
    stop("x_matrix contains NA; set na_as_zero = TRUE or impute missing values")
  }

  col_norms <- sqrt(colSums(x^2))
  col_norms[col_norms == 0] <- 1
  x_norm <- sweep(x, 2, col_norms, "/")

  snp_coords <- t(x_norm)
  if (!is.null(colnames(x_matrix))) {
    rownames(snp_coords) <- colnames(x_matrix)
  }

  max_latent <- min(nrow(snp_coords) - 1L, ncol(snp_coords) - 1L, k + 1L)
  if (max_latent < 1L) {
    stop("Not enough dimensions for latent-space GMM clustering")
  }

  if (is.null(n_latent)) {
    n_latent <- max_latent
  } else {
    n_latent <- min(as.integer(n_latent), max_latent)
  }

  pca <- stats::prcomp(snp_coords, center = TRUE, scale. = FALSE)
  latent <- pca$x[, seq_len(n_latent), drop = FALSE]

  if (!"package:mclust" %in% search()) {
    attachNamespace("mclust")
  }
  fit <- mclust::Mclust(latent, G = k, verbose = FALSE)
  if (is.null(fit)) {
    warning("mclust::Mclust() failed; falling back to k-means in latent space")
    km <- stats::kmeans(latent, centers = k, nstart = 25)
    return(list(
      cluster = km$cluster,
      details = list(
        kmeans = km,
        latent = latent,
        pca = pca,
        mclust = NULL
      )
    ))
  }

  return(list(
    cluster = fit$classification,
    details = list(
      mclust = fit,
      latent = latent,
      pca = pca
    )
  ))
}


#' @title Summarise SNP Module Phenotype Drivers
#' @description For each SNP cluster, compute the mean oriented pleiotropy profile
#' across member SNPs and rank background traits by absolute effect magnitude to
#' identify functional drivers of each biological module.
#' @param x_matrix Oriented pleiotropy matrix (traits x SNPs), typically
#'   `orient_pleiotropy_matrix()$x_matrix`.
#' @param cluster Named integer vector of cluster assignments per SNP column
#'   (as returned by `cluster_snp_profiles()$cluster`).
#' @param trait_info Optional dataframe with `trait_id` and `trait_name` columns
#'   (as returned by `build_pleiotropy_matrix()$trait_info`).
#' @param exclude_trait_id Optional trait ID to exclude from summaries (e.g. the
#'   target trait when interpreting background drivers).
#' @param na_as_zero Treat `NA` entries as zero before averaging. Defaults to `TRUE`.
#' @param top_n If not `NULL`, return only the top `n` traits per cluster ranked by
#'   `abs_mean_z`.
#' @return A dataframe with columns: `cluster`, `trait_id`, `trait_name` (if
#'   available), `mean_z`, `abs_mean_z`, `n_snps`.
#' @export
summarise_snp_modules <- function(x_matrix,
                                  cluster,
                                  trait_info = NULL,
                                  exclude_trait_id = NULL,
                                  na_as_zero = TRUE,
                                  top_n = NULL) {
  if (length(cluster) != ncol(x_matrix)) {
    stop("length(cluster) must match ncol(x_matrix)")
  }

  if (is.null(names(cluster)) && !is.null(colnames(x_matrix))) {
    names(cluster) <- colnames(x_matrix)
  }

  x <- x_matrix
  if (na_as_zero) {
    x[is.na(x)] <- 0
  }

  cluster_ids <- sort(unique(cluster))
  summaries <- lapply(cluster_ids, function(cl) {
    snp_cols <- names(cluster)[cluster == cl]
    module_matrix <- x[, snp_cols, drop = FALSE]
    mean_z <- rowMeans(module_matrix, na.rm = TRUE)

    out <- data.frame(
      cluster = cl,
      trait_id = rownames(x_matrix),
      mean_z = mean_z,
      abs_mean_z = abs(mean_z),
      n_snps = length(snp_cols),
      stringsAsFactors = FALSE
    )

    if (!is.null(exclude_trait_id)) {
      out <- out[out$trait_id != as.character(exclude_trait_id), , drop = FALSE]
    }

    out <- out[order(-out$abs_mean_z), , drop = FALSE]

    if (!is.null(top_n)) {
      out <- utils::head(out, top_n)
    }

    return(out)
  })

  result <- dplyr::bind_rows(summaries)

  if (!is.null(trait_info)) {
    trait_info <- trait_info |>
      dplyr::mutate(trait_id = as.character(trait_id))
    result <- result |>
      dplyr::left_join(trait_info, by = "trait_id")
  }

  return(result)
}
