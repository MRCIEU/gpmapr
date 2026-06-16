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
#' @return A list with:
#'   \itemize{
#'     \item x_matrix: numeric matrix (traits x SNPs) of z-scores; `NA` where a trait
#'       has no colocalisation signal at that SNP's locus
#'     \item trait_info: dataframe mapping row indices to `trait_id` and `trait_name`
#'     \item snp_info: dataframe mapping column names to `coloc_group_id`, `variant_id`,
#'       `display_snp`, `chr`, and `bp`
#'     \item target_trait_id: the target trait ID used
#'   }
#' @export
build_pleiotropy_matrix <- function(trait_id,
                                    coloc_groups = NULL,
                                    p_threshold = NULL,
                                    snp_key = c("variant_id", "display_snp", "coloc_group_id")) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }

  snp_key <- match.arg(snp_key)
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

  z_long <- locus_data$cg |>
    dplyr::group_by(trait_id, snp_id) |>
    dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select("trait_id", "trait_name", "snp_id", "z")

  z_wide <- z_long |>
    tidyr::pivot_wider(
      names_from = "snp_id",
      values_from = "z"
    )

  trait_info <- z_wide |>
    dplyr::select("trait_id", "trait_name") |>
    dplyr::distinct()

  snp_ids <- setdiff(names(z_wide), c("trait_id", "trait_name"))
  x_matrix <- as.matrix(z_wide[, snp_ids, drop = FALSE])
  rownames(x_matrix) <- as.character(z_wide$trait_id)

  snp_info <- locus_data$target_snps |>
    dplyr::filter(snp_id %in% snp_ids) |>
    dplyr::distinct()

  return(list(
    x_matrix = x_matrix,
    trait_info = trait_info,
    snp_info = snp_info,
    target_trait_id = target_id
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
#'   `traits(c(trait_id_1, trait_id_2), include_associations = TRUE)`.
#' @param p_threshold P-value threshold for including target-trait SNPs. If NULL,
#'   all SNPs with valid z-scores are retained.
#' @param snp_key Column used to name SNP columns. Defaults to `"variant_id"`.
#' @return A list with:
#'   \itemize{
#'     \item x1_matrix, x2_matrix: aligned numeric matrices (traits x SNPs)
#'     \item trait_info: dataframe of shared background trait metadata
#'     \item snp_info_1, snp_info_2: SNP metadata for each matrix
#'     \item trait_id_1, trait_id_2: target trait IDs
#'   }
#' @export
build_bivariate_pleiotropy_matrices <- function(trait_id_1,
                                                trait_id_2,
                                                coloc_groups = NULL,
                                                p_threshold = NULL,
                                                snp_key = c("variant_id", "display_snp", "coloc_group_id")) {
  if (missing(trait_id_1) || is.null(trait_id_1) ||
      missing(trait_id_2) || is.null(trait_id_2)) {
    stop("trait_id_1 and trait_id_2 are required")
  }

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
    snp_key = snp_key
  )
  pleiotropy_2 <- build_pleiotropy_matrix(
    trait_id = trait_id_2,
    coloc_groups = coloc_groups,
    p_threshold = p_threshold,
    snp_key = snp_key
  )

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
    dplyr::distinct(trait_id, .keep_all = TRUE) |>
    dplyr::filter(as.character(trait_id) %in% shared_trait_ids)

  missing_trait_ids <- setdiff(shared_trait_ids, as.character(trait_info$trait_id))
  if (length(missing_trait_ids) > 0) {
    trait_info <- dplyr::bind_rows(
      trait_info,
      data.frame(
        trait_id = as.integer(missing_trait_ids),
        trait_name = NA_character_,
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
    trait_id_2 = trait_id_2
  ))
}


#' @title Build Collapsed-Trait Pleiotropy Matrix
#' @description Construct a pleiotropy matrix whose rows are `collapsed_trait`
#' groups rather than individual background traits. Gene-linked QTL signals are
#' first merged by `gene_id` (across QTL types and coloc groups), then mapped to
#' enriched pathways. Genes outside enriched pathways remain as gene-level
#' collapsed traits.
#' @inheritParams build_pleiotropy_matrix
#' @param collapsed_trait_map Optional gene-to-collapsed-trait mapping from
#'   `build_pathway_collapsed_trait_map()`. If NULL, it is built from genes at
#'   the target SNP loci.
#' @param pathway_source Optional pathway source passed to `pathway_enrichment()`.
#' @param pathway_p_value_threshold FDR threshold for pathway enrichment.
#' @param minimum_count_in_network Minimum gene overlap per pathway term.
#' @return A list with:
#'   \itemize{
#'     \item x_matrix: collapsed traits x SNPs matrix
#'     \item collapsed_trait_info: metadata for each collapsed_trait row
#'     \item z_target: target-trait z-scores per SNP column (for orientation)
#'     \item snp_info, target_trait_id, pathway_enrichment, collapsed_trait_map
#'   }
#' @export
build_collapsed_pleiotropy_matrix <- function(trait_id,
                                               coloc_groups = NULL,
                                               p_threshold = NULL,
                                               snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                               collapsed_trait_map = NULL,
                                               pathway_source = NULL,
                                               pathway_p_value_threshold = 0.05,
                                               minimum_count_in_network = NULL) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }

  snp_key <- match.arg(snp_key)
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

  if (!"gene_id" %in% names(locus_data$cg)) {
    stop("coloc_groups must include gene_id for collapsed_trait matrices")
  }

  gene_z <- .collapse_gene_z_scores(locus_data$cg)
  if (nrow(gene_z) == 0) {
    stop("No gene-linked colocalisation signals available for collapsed_trait matrix")
  }

  if (is.null(collapsed_trait_map)) {
    map_result <- build_pathway_collapsed_trait_map(
      genes = unique(gene_z$gene_id),
      source = pathway_source,
      p_value_threshold = pathway_p_value_threshold,
      minimum_count_in_network = minimum_count_in_network
    )
    collapsed_trait_map <- map_result$collapsed_trait_map
    pathway_enrichment <- map_result$pathway_enrichment
  } else {
    pathway_enrichment <- NULL
  }

  matrix_result <- .collapsed_z_to_matrix(
    gene_z = gene_z,
    collapsed_trait_map = collapsed_trait_map,
    snp_ids = as.character(unique(locus_data$target_snps$snp_id))
  )

  z_target <- stats::setNames(
    locus_data$z_target$z,
    as.character(locus_data$z_target$snp_id)
  )
  z_target <- z_target[colnames(matrix_result$x_matrix)]

  snp_info <- locus_data$target_snps |>
    dplyr::filter(snp_id %in% colnames(matrix_result$x_matrix)) |>
    dplyr::distinct()

  return(list(
    x_matrix = matrix_result$x_matrix,
    collapsed_trait_info = matrix_result$collapsed_trait_info,
    z_target = z_target,
    snp_info = snp_info,
    target_trait_id = target_id,
    pathway_enrichment = pathway_enrichment,
    collapsed_trait_map = collapsed_trait_map
  ))
}


#' @title Build Pathway Collapsed-Trait Map
#' @description Map genes to enriched pathway terms (`collapsed_trait`) and retain
#' gene-level collapsed traits for genes not assigned to any enriched pathway.
#' @param genes Numeric gene IDs or gene names accepted by `pathway_enrichment()`.
#' @param source Optional pathway source: `"Reactome"`, `"KEGG"`, or `"HP"`.
#' @param p_value_threshold FDR threshold for enriched pathways.
#' @param minimum_count_in_network Minimum overlap per pathway term.
#' @return A list with `collapsed_trait_map` (gene_id, collapsed_trait_id,
#'   collapsed_trait_name, source) and `pathway_enrichment`.
#' @export
build_pathway_collapsed_trait_map <- function(genes,
                                              source = NULL,
                                              p_value_threshold = 0.05,
                                              minimum_count_in_network = NULL) {
  if (is.null(genes) || length(genes) == 0) {
    stop("genes is required")
  }

  pathway_enrichment <- pathway_enrichment(
    genes = genes,
    source = source,
    p_value_threshold = p_value_threshold,
    minimum_count_in_network = minimum_count_in_network
  )

  collapsed_trait_map <- data.frame(
    gene_id = integer(0),
    collapsed_trait_id = character(0),
    collapsed_trait_name = character(0),
    source = character(0),
    stringsAsFactors = FALSE
  )

  if (is.data.frame(pathway_enrichment$results) && nrow(pathway_enrichment$results) > 0) {
    pathway_rows <- lapply(seq_len(nrow(pathway_enrichment$results)), function(i) {
      row <- pathway_enrichment$results[i, , drop = FALSE]
      overlap_genes <- row$gene_ids[[1]]
      if (length(overlap_genes) == 0) {
        return(NULL)
      }
      data.frame(
        gene_id = overlap_genes,
        collapsed_trait_id = paste0(row$source, ":", row$term_id),
        collapsed_trait_name = row$description,
        source = row$source,
        stringsAsFactors = FALSE
      )
    })
    collapsed_trait_map <- dplyr::bind_rows(pathway_rows)
  }

  return(list(
    collapsed_trait_map = collapsed_trait_map,
    pathway_enrichment = pathway_enrichment
  ))
}


#' @title Build Bivariate Collapsed-Trait Pleiotropy Matrices
#' @description Like `build_bivariate_pleiotropy_matrices()`, but rows are shared
#' `collapsed_trait` groups derived from pathway enrichment on the union of
#' gene-linked signals at both target SNP sets.
#' @inheritParams build_bivariate_pleiotropy_matrices
#' @inheritParams build_collapsed_pleiotropy_matrix
#' @return A list with aligned `x1_matrix`, `x2_matrix`, `collapsed_trait_info`,
#'   `z_target_1`, `z_target_2`, SNP metadata, and the shared `collapsed_trait_map`.
#' @export
build_bivariate_collapsed_pleiotropy_matrices <- function(trait_id_1,
                                                          trait_id_2,
                                                          coloc_groups = NULL,
                                                          p_threshold = NULL,
                                                          snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                                          pathway_source = NULL,
                                                          pathway_p_value_threshold = 0.05,
                                                          minimum_count_in_network = NULL) {
  if (missing(trait_id_1) || is.null(trait_id_1) ||
      missing(trait_id_2) || is.null(trait_id_2)) {
    stop("trait_id_1 and trait_id_2 are required")
  }

  snp_key <- match.arg(snp_key)

  if (is.null(coloc_groups)) {
    coloc_groups <- traits(
      c(trait_id_1, trait_id_2),
      include_associations = TRUE
    )$coloc_groups
  }

  locus_data_1 <- .prepare_pleiotropy_locus_data(
    trait_id = trait_id_1,
    coloc_groups = coloc_groups,
    p_threshold = p_threshold,
    snp_key = snp_key
  )
  locus_data_2 <- .prepare_pleiotropy_locus_data(
    trait_id = trait_id_2,
    coloc_groups = coloc_groups,
    p_threshold = p_threshold,
    snp_key = snp_key
  )

  gene_z_1 <- .collapse_gene_z_scores(locus_data_1$cg)
  gene_z_2 <- .collapse_gene_z_scores(locus_data_2$cg)

  if (!"gene_id" %in% names(locus_data_1$cg)) {
    stop("coloc_groups must include gene_id for collapsed_trait matrices")
  }

  shared_genes <- union(gene_z_1$gene_id, gene_z_2$gene_id)

  if (length(shared_genes) == 0) {
    stop("No gene-linked colocalisation signals available for collapsed_trait matrices")
  }

  map_result <- build_pathway_collapsed_trait_map(
    genes = shared_genes,
    source = pathway_source,
    p_value_threshold = pathway_p_value_threshold,
    minimum_count_in_network = minimum_count_in_network
  )
  collapsed_trait_map <- map_result$collapsed_trait_map

  matrix_1 <- .collapsed_z_to_matrix(
    gene_z = gene_z_1,
    collapsed_trait_map = collapsed_trait_map,
    snp_ids = as.character(unique(locus_data_1$target_snps$snp_id))
  )
  matrix_2 <- .collapsed_z_to_matrix(
    gene_z = gene_z_2,
    collapsed_trait_map = collapsed_trait_map,
    snp_ids = as.character(unique(locus_data_2$target_snps$snp_id))
  )

  shared_collapsed_trait_ids <- union(
    rownames(matrix_1$x_matrix),
    rownames(matrix_2$x_matrix)
  )

  collapsed_trait_info <- dplyr::bind_rows(
    matrix_1$collapsed_trait_info,
    matrix_2$collapsed_trait_info
  ) |>
    dplyr::distinct(collapsed_trait_id, .keep_all = TRUE) |>
    dplyr::filter(collapsed_trait_id %in% shared_collapsed_trait_ids)

  return(list(
    x1_matrix = .align_pleiotropy_rows(matrix_1$x_matrix, shared_collapsed_trait_ids),
    x2_matrix = .align_pleiotropy_rows(matrix_2$x_matrix, shared_collapsed_trait_ids),
    collapsed_trait_info = collapsed_trait_info,
    z_target_1 = stats::setNames(locus_data_1$z_target$z, as.character(locus_data_1$z_target$snp_id))[
      colnames(matrix_1$x_matrix)
    ],
    z_target_2 = stats::setNames(locus_data_2$z_target$z, as.character(locus_data_2$z_target$snp_id))[
      colnames(matrix_2$x_matrix)
    ],
    snp_info_1 = locus_data_1$target_snps |> dplyr::distinct(),
    snp_info_2 = locus_data_2$target_snps |> dplyr::distinct(),
    trait_id_1 = trait_id_1,
    trait_id_2 = trait_id_2,
    collapsed_trait_map = collapsed_trait_map,
    pathway_enrichment = map_result$pathway_enrichment
  ))
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

  target_snps <- coloc_groups |>
    dplyr::filter(
      trait_id == trait_id,
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
    dplyr::filter(trait_id == trait_id) |>
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


.collapse_gene_z_scores <- function(cg) {
  cg |>
    dplyr::filter(!is.na(gene_id)) |>
    dplyr::group_by(gene_id, gene, snp_id) |>
    dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select("gene_id", "gene", "snp_id", "z")
}


.collapsed_z_to_matrix <- function(gene_z, collapsed_trait_map, snp_ids) {
  snp_ids <- as.character(snp_ids)
  gene_z <- gene_z |> dplyr::mutate(snp_id = as.character(snp_id))
  if (nrow(collapsed_trait_map) > 0) {
    pathway_z <- gene_z |>
      dplyr::inner_join(collapsed_trait_map, by = "gene_id") |>
      dplyr::group_by(collapsed_trait_id, collapsed_trait_name, source, snp_id) |>
      dplyr::summarise(z = mean(z), .groups = "drop")
  } else {
    pathway_z <- data.frame(
      collapsed_trait_id = character(0),
      collapsed_trait_name = character(0),
      source = character(0),
      snp_id = character(0),
      z = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  mapped_genes <- if (nrow(collapsed_trait_map) > 0) {
    unique(collapsed_trait_map$gene_id)
  } else {
    integer(0)
  }

  gene_only_z <- gene_z |>
    dplyr::filter(!gene_id %in% mapped_genes) |>
    dplyr::mutate(
      collapsed_trait_id = paste0("gene:", gene_id),
      collapsed_trait_name = paste0("gene:", gene),
      source = "gene"
    ) |>
    dplyr::select("collapsed_trait_id", "collapsed_trait_name", "source", "snp_id", "z")

  collapsed_long <- dplyr::bind_rows(pathway_z, gene_only_z)

  if (nrow(collapsed_long) == 0) {
    stop("No collapsed_trait rows could be constructed from gene-linked signals")
  }

  collapsed_wide <- collapsed_long |>
    tidyr::pivot_wider(
      id_cols = c("collapsed_trait_id", "collapsed_trait_name", "source"),
      names_from = "snp_id",
      values_from = "z",
      values_fn = mean
    )

  collapsed_trait_info <- collapsed_wide |>
    dplyr::select("collapsed_trait_id", "collapsed_trait_name", "source")

  snp_cols <- intersect(snp_ids, setdiff(names(collapsed_wide), c(
    "collapsed_trait_id", "collapsed_trait_name", "source"
  )))

  if (length(snp_cols) == 0) {
    stop("No collapsed_trait rows could be constructed from gene-linked signals")
  }

  x_matrix <- as.matrix(collapsed_wide[, snp_cols, drop = FALSE])
  rownames(x_matrix) <- collapsed_wide$collapsed_trait_id

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
  x_matrix <- x_matrix[, snp_ids, drop = FALSE]

  return(list(
    x_matrix = x_matrix,
    collapsed_trait_info = collapsed_trait_info
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
  if (any(col_norms == 0)) {
    warning("One or more SNP columns have zero norm; similarities involving them are set to NA")
    col_norms[col_norms == 0] <- NA_real_
  }

  x_norm <- sweep(x, 2, col_norms, "/")

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


#' @title Dual Hierarchical Clustering of Cross-Trait SNP Similarity
#' @description Cluster Trait-1 and Trait-2 SNPs independently from a rectangular
#' cross-similarity matrix. Row SNPs are clustered on Euclidean distance across
#' similarity profiles; column SNPs are clustered on the transposed matrix.
#' @param s_matrix Rectangular cross-trait similarity matrix from
#'   `cross_snp_similarity_matrix()$s_matrix`.
#' @param linkage Linkage method passed to `stats::hclust()`. Defaults to `"average"`.
#' @param k1 Optional number of row (Trait-1 SNP) clusters.
#' @param k2 Optional number of column (Trait-2 SNP) clusters.
#' @return A list with:
#'   \itemize{
#'     \item hclust_rows: hierarchical clustering of Trait-1 SNPs
#'     \item hclust_cols: hierarchical clustering of Trait-2 SNPs
#'     \item clusters_rows: row cluster assignments if `k1` is set
#'     \item clusters_cols: column cluster assignments if `k2` is set
#'   }
#' @export
cluster_cross_trait_snps <- function(s_matrix,
                                     linkage = "average",
                                     k1 = NULL,
                                     k2 = NULL) {
  if (!is.matrix(s_matrix)) {
    stop("s_matrix must be a matrix")
  }

  s <- s_matrix
  s[is.na(s)] <- 0

  dist_rows <- stats::dist(s, method = "euclidean")
  hclust_rows <- stats::hclust(dist_rows, method = linkage)

  dist_cols <- stats::dist(t(s), method = "euclidean")
  hclust_cols <- stats::hclust(dist_cols, method = linkage)

  out <- list(
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
#' @param k Number of clusters (biological modules). For Louvain community detection,
#'   `k` is ignored and the number of modules is data-driven.
#' @param method Clustering method: `"hierarchical"`, `"spectral"`, `"community"`,
#'   or `"gmm"`.
#' @param x_matrix Oriented pleiotropy matrix (traits x SNPs). Required for
#'   `method = "gmm"`; typically `orient_pleiotropy_matrix()$x_matrix`.
#' @param linkage Linkage method passed to `stats::hclust()` when
#'   `method = "hierarchical"`. Defaults to `"average"`.
#' @param community_algorithm Community detection algorithm when
#'   `method = "community"`: `"leading_eigen"` (recursive splits via `steps = k - 1`)
#'   or `"louvain"` (data-driven module count; requires the **igraph** package).
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
#'     \item details: method-specific objects (e.g. `hclust`, `kmeans`, `mclust`)
#'   }
#' @export
cluster_snp_profiles <- function(s_matrix,
                                 k = 3,
                                 method = c("hierarchical", "spectral", "community", "gmm"),
                                 x_matrix = NULL,
                                 linkage = "average",
                                 community_algorithm = c("leading_eigen", "louvain"),
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

  snp_ids <- colnames(s_matrix)
  if (is.null(snp_ids)) {
    snp_ids <- as.character(seq_len(ncol(s_matrix)))
  }

  S <- s_matrix
  dimnames(S) <- list(snp_ids, snp_ids)
  diag(S) <- 1
  S[is.na(S)] <- 0

  if (method == "hierarchical") {
    result <- .cluster_snps_hierarchical(S, k = k, linkage = linkage)
  } else if (method == "spectral") {
    result <- .cluster_snps_spectral(S, k = k)
  } else if (method == "community") {
    result <- .cluster_snps_community(S, k = k, algorithm = community_algorithm)
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
    details = result$details
  ))
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


.cluster_snps_community <- function(s_matrix, k, algorithm) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for community detection", call. = FALSE)
  }

  affinity <- s_matrix
  diag(affinity) <- 0
  affinity[affinity < 0] <- 0

  graph <- igraph::graph_from_adjacency_matrix(
    affinity,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  if (algorithm == "louvain") {
    communities <- igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight)
    membership <- igraph::membership(communities)
  } else {
    communities <- igraph::cluster_leading_eigen(
      graph,
      steps = k - 1L,
      weights = igraph::E(graph)$weight
    )
    membership <- igraph::membership(communities)
  }

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
