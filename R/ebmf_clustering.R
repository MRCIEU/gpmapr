#' @title Build EBMF Matrix with Annotations
#' @description Construct a flat traits x SNPs z-score matrix and annotation vectors
#' for Empirical Bayes Matrix Factorization (EBMF). Rows represent traits (uncollapsed),
#' columns represent variants. Cell values are signed z-scores (`beta / se`).
#' Rows and columns that are entirely missing or zero are removed before the
#' matrix is returned.
#'
#' Annotations encode biological metadata for each row and column:
#' \itemize{
#'   \item Trait (row): \code{Tissue__Gene__EnrichedPathway__TraitCategory}
#'   \item SNP (column): \code{ColocGroup\{ID\}__VariantConsequence}
#' }
#' @inheritParams build_pleiotropy_matrix
#' @param pathway_source Optional pathway source for enrichment annotation.
#' @param pathway_p_value_threshold FDR threshold for pathway enrichment annotation.
#'   Defaults to 0.05.
#' @param minimum_count_in_network Minimum gene overlap per pathway term.
#' @return A list with:
#'   \itemize{
#'     \item x_matrix: traits x SNPs z-score matrix
#'     \item trait_annotations: named character vector of row annotations
#'     \item snp_annotations: named character vector of column annotations
#'     \item trait_info: dataframe of trait metadata
#'     \item snp_info: dataframe of SNP metadata
#'     \item target_trait_id: target trait ID
#'     \item pathway_enrichment: pathway enrichment results (if computed)
#'   }
#' @export
build_ebmf_matrix <- function(trait_id,
                              coloc_groups = NULL,
                              p_threshold = NULL,
                              snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                              pathway_source = NULL,
                              pathway_p_value_threshold = 0.05,
                              minimum_count_in_network = NULL) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }

  snp_key <- match.arg(snp_key)

  if (is.null(coloc_groups)) {
    coloc_groups <- trait(trait_id, include_associations = TRUE)$coloc_groups
  }

  pleiotropy <- build_pleiotropy_matrix(
    trait_id = trait_id,
    coloc_groups = coloc_groups,
    p_threshold = p_threshold,
    snp_key = snp_key
  )

  filtered <- .filter_ebmf_zero_slices(
    x_matrix = pleiotropy$x_matrix
  )

  trait_annot <- .build_ebmf_trait_annotations(
    trait_ids = rownames(filtered$x_matrix),
    coloc_groups = coloc_groups,
    pathway_source = pathway_source,
    pathway_p_value_threshold = pathway_p_value_threshold,
    minimum_count_in_network = minimum_count_in_network
  )

  snp_annot <- .build_ebmf_snp_annotations(
    snp_ids = colnames(filtered$x_matrix),
    snp_info = pleiotropy$snp_info,
    coloc_groups = coloc_groups,
    target_trait_id = trait_id
  )

  trait_info <- pleiotropy$trait_info |>
    dplyr::filter(as.character(trait_id) %in% rownames(filtered$x_matrix))

  snp_info <- pleiotropy$snp_info |>
    dplyr::filter(as.character(snp_id) %in% colnames(filtered$x_matrix))

  return(list(
    x_matrix = filtered$x_matrix,
    trait_annotations = trait_annot$annotations,
    snp_annotations = snp_annot,
    trait_info = trait_info,
    snp_info = snp_info,
    target_trait_id = trait_id,
    pathway_enrichment = trait_annot$pathway_enrichment
  ))
}


#' @title Run EBMF Factorization
#' @description Fit an Empirical Bayes Matrix Factorization model using
#' \pkg{flashier}. The model decomposes the traits x SNPs z-score matrix
#' \eqn{Y = LF' + E} into latent factors (\eqn{F}, SNP programs) and
#' loadings (\eqn{L}, trait contributions), with sparse priors estimated
#' from the data.
#'
#' By default a Kronecker (two-way) residual variance structure is used,
#' estimating per-row and per-column variance parameters. This naturally
#' down-weights noisy rows (singleton/unannotated traits) and columns
#' (intergenic SNPs), acting as an automatic noise shield.
#'
#' @param x_matrix Traits x SNPs z-score matrix, typically from
#'   \code{build_ebmf_matrix()$x_matrix}. Rows and columns that are entirely
#'   missing or zero are removed before fitting. \code{NA} entries are treated
#'   as missing data.
#' @param greedy_Kmax Maximum number of factors to add greedily. Factors are
#'   only added while they improve the variational lower bound. Defaults to 50.
#' @param backfit Logical; if \code{TRUE} (default), all factors are cyclically
#'   updated after the greedy phase until convergence.
#' @param ebnm_fn EBNM function for priors on loadings and factors. Defaults to
#'   \code{ebnm::ebnm_point_normal} (point-normal prior). Pass a list of two
#'   functions to use different priors for loadings (traits) and factors (SNPs).
#' @param var_type Residual variance structure. Defaults to \code{c(1, 2)} for
#'   Kronecker (two-way) estimation. Use \code{0} for constant variance,
#'   \code{1} for per-row, or \code{2} for per-column.
#' @param verbose Verbosity: 0 = silent, 1 = summary, 2 = ELBO updates,
#'   3 = per-iteration.
#' @return A \code{flash} object from \pkg{flashier}. Key elements:
#'   \itemize{
#'     \item n_factors: number of discovered factors (programs)
#'     \item L_pm, L_lfsr: posterior means and lFSR for trait loadings (n x K)
#'     \item F_pm, F_lfsr: posterior means and lFSR for SNP factors (p x K)
#'     \item pve: proportion of variance explained per factor
#'     \item elbo: variational lower bound
#'   }
#' @export
run_ebmf <- function(x_matrix,
                     greedy_Kmax = 50L,
                     backfit = TRUE,
                     ebnm_fn = NULL,
                     var_type = c(1L, 2L),
                     verbose = 1L) {
  if (!requireNamespace("flashier", quietly = TRUE)) {
    stop("Package 'flashier' is required for EBMF clustering", call. = FALSE)
  }
  if (!requireNamespace("ebnm", quietly = TRUE)) {
    stop("Package 'ebnm' is required for EBMF clustering", call. = FALSE)
  }

  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }

  if (is.null(ebnm_fn)) {
    ebnm_fn <- ebnm::ebnm_point_normal
  }

  filtered <- .filter_ebmf_zero_slices(x_matrix = x_matrix)

  fit <- flashier::flash(
    data = filtered$x_matrix,
    S = NULL,
    ebnm_fn = ebnm_fn,
    var_type = var_type,
    greedy_Kmax = greedy_Kmax,
    backfit = backfit,
    nullcheck = TRUE,
    verbose = verbose
  )

  return(fit)
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
#'   Defaults to 0.10.
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
                                  magnitude_threshold = 0.10) {
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


.ebmf_slice_has_signal <- function(x) {
  return(any(!is.na(x) & x != 0))
}


.filter_ebmf_zero_slices <- function(x_matrix, z_target = NULL) {
  if (!is.matrix(x_matrix)) {
    stop("x_matrix must be a matrix")
  }

  row_keep <- apply(x_matrix, 1, .ebmf_slice_has_signal)
  col_keep <- apply(x_matrix, 2, .ebmf_slice_has_signal)

  removed_rows <- sum(!row_keep)
  removed_cols <- sum(!col_keep)
  if (removed_rows > 0 || removed_cols > 0) {
    warning(
      "Removed ", removed_rows, " all-zero/all-missing row(s) and ",
      removed_cols, " all-zero/all-missing column(s) for EBMF",
      call. = FALSE
    )
  }

  x_matrix <- x_matrix[row_keep, col_keep, drop = FALSE]

  if (nrow(x_matrix) == 0 || ncol(x_matrix) == 0) {
    stop("No rows or columns remain after removing all-zero/all-missing slices")
  }

  if (!is.null(z_target)) {
    z_target <- z_target[colnames(x_matrix)]
  }

  return(list(
    x_matrix = x_matrix,
    z_target = z_target
  ))
}


.build_ebmf_trait_annotations <- function(trait_ids,
                                          coloc_groups,
                                          pathway_source = NULL,
                                          pathway_p_value_threshold = 0.05,
                                          minimum_count_in_network = NULL) {
  trait_meta <- coloc_groups |>
    dplyr::filter(as.character(trait_id) %in% trait_ids) |>
    dplyr::select(trait_id, tissue, gene, gene_id, trait_category) |>
    dplyr::group_by(trait_id) |>
    dplyr::summarise(
      tissue = dplyr::first(tissue[!is.na(tissue) & tissue != ""]),
      gene = dplyr::first(gene[!is.na(gene) & gene != ""]),
      gene_id = dplyr::first(gene_id[!is.na(gene_id)]),
      trait_category = dplyr::first(trait_category[!is.na(trait_category) & trait_category != ""]),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      trait_id = as.character(trait_id),
      tissue = dplyr::if_else(is.na(tissue), "None", tissue),
      gene = dplyr::if_else(is.na(gene), "None", gene),
      trait_category = dplyr::if_else(is.na(trait_category), "Unknown", trait_category)
    )

  # Pathway enrichment for gene annotation
  gene_ids <- unique(stats::na.omit(trait_meta$gene_id))
  gene_pathway <- data.frame(
    gene_id = integer(0),
    enriched_pathway = character(0),
    stringsAsFactors = FALSE
  )
  pathway_enrichment <- NULL

  if (length(gene_ids) > 0) {
    tryCatch({
      map_result <- build_pathway_feature_map(
        genes = gene_ids,
        source = pathway_source,
        p_value_threshold = pathway_p_value_threshold,
        minimum_count_in_network = minimum_count_in_network
      )
      pathway_enrichment <- map_result$pathway_enrichment

      if (nrow(map_result$feature_map) > 0) {
        gene_pathway <- map_result$feature_map |>
          dplyr::group_by(gene_id) |>
          dplyr::summarise(
            enriched_pathway = dplyr::first(feature_name),
            .groups = "drop"
          )
      }
    }, error = function(e) {
      warning("Pathway enrichment failed: ", conditionMessage(e), call. = FALSE)
      return(NULL)
    })
  }

  if (nrow(gene_pathway) > 0) {
    trait_meta <- trait_meta |>
      dplyr::left_join(gene_pathway, by = "gene_id")
  } else {
    trait_meta$enriched_pathway <- NA_character_
  }
  trait_meta <- trait_meta |>
    dplyr::mutate(
      enriched_pathway = dplyr::if_else(
        is.na(enriched_pathway), "Unannotated", enriched_pathway
      )
    )

  # Build annotation strings preserving row order
  meta_lookup <- stats::setNames(
    paste(
      trait_meta$tissue, trait_meta$gene,
      trait_meta$enriched_pathway, trait_meta$trait_category,
      sep = "__"
    ),
    trait_meta$trait_id
  )
  annotations <- dplyr::if_else(
    trait_ids %in% names(meta_lookup),
    meta_lookup[trait_ids],
    "None__None__Unannotated__Unknown"
  )
  names(annotations) <- trait_ids

  return(list(
    annotations = annotations,
    pathway_enrichment = pathway_enrichment
  ))
}


.build_ebmf_snp_annotations <- function(snp_ids,
                                        snp_info,
                                        coloc_groups,
                                        target_trait_id) {
  # Get coloc_group and gene info for target-trait SNPs
  target_gene_info <- coloc_groups |>
    dplyr::filter(trait_id == target_trait_id) |>
    dplyr::select(coloc_group_id, gene, gene_id) |>
    dplyr::distinct(coloc_group_id, .keep_all = TRUE)

  snp_meta <- snp_info |>
    dplyr::left_join(target_gene_info, by = "coloc_group_id") |>
    dplyr::mutate(
      consequence = dplyr::if_else(
        !is.na(gene) & gene != "", "Genic", "Intergenic"
      )
    )

  meta_lookup <- stats::setNames(
    paste0("ColocGroup", snp_meta$coloc_group_id, "__", snp_meta$consequence),
    as.character(snp_meta$snp_id)
  )

  annotations <- dplyr::if_else(
    snp_ids %in% names(meta_lookup),
    meta_lookup[snp_ids],
    paste0("ColocGroupNA__Unknown")
  )
  names(annotations) <- snp_ids

  return(annotations)
}
