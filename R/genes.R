#' @title All genes
#' @description Get all genes from the API
#' @return A dataframe containing all genes with the following columns:
#'   \itemize{
#'     \item id: the id of the gene
#'     \item gene: the name of the gene
#'     \item description: the description of the gene
#'     \item gene_biotype: the gene biotype
#'     \item chr: the chromosome of the gene
#'     \item start: the start position of the gene
#'     \item stop: the end position of the gene
#'     \item strand: the strand of the gene
#'     \item source: the source of the gene
#'     \item distinct_trait_categories: the number of trait categories that the gene is associated with via coloc groups
#'     \item distinct_protein_coding_genes: the number of genes that the gene is associated with via coloc groups
#'     \item num_study_extractions: the number of study extractions for this gene
#'     \item num_coloc_groups: the number of coloc groups for this gene
#'     \item num_coloc_studies: the number of studies that have coloc results for this gene
#'     \item num_rare_groups: the number of rare groups for this gene
#'   }
#' @export
all_genes <- function() {
  genes <- genes_api()
  return(genes$genes)
}

#' @title Gene
#' @description A collection of studies that are associated with a particular gene.
#' @param gene_id A numeric value specifying the gene id
#' @param include_associations A logical value specifying whether to include associations
#' (BETA, SE, P), defaults to FALSE
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs, defaults to FALSE
#' @param include_trans A logical value specifying whether to include trans genetic effects, defaults to TRUE
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs, defaults to 0.8
#' @return A list which contains the following elements:
#' \itemize{
#'   \item gene: A list containing metadata about the gene, including region, and neighboring genes.
#'   \item coloc_groups: a dataframe containing information about which studies have coloc results for this gene.
#' See below for details.
#'   \item study_extractions: a list of dataframes containing the study extractions for this trait.
#' See below for details.
#'   \item rare_results: (optional) a list of dataframes containing the rare results for this trait
#'   \item coloc_pairs: (optional) a dataframe containing all pairwise coloc results for this trait.
#'   \item variants: a dataframe containing the variants for each associated coloc group or rare group.
#' }
#' See below for details.
#' @details
#' The dataframes returned by this function are as follows:
#' @inheritSection coloc_groups_doc coloc_groups_dataframe
#' @inheritSection study_extractions_doc study_extractions_dataframe
#' @inheritSection rare_results_doc rare_results_dataframe
#' @inheritSection coloc_pairs_doc coloc_pairs_dataframe
#' @inheritSection variants_doc variants_dataframe
#' @export
gene <- function(gene_id,
  include_associations = FALSE,
  include_coloc_pairs = FALSE,
  include_trans = TRUE,
  h4_threshold = 0.8
) {
  if (is.null(gene_id)) {
    stop("gene_id is required")
  }
  gene_info <- gene_api(gene_id, include_associations, include_coloc_pairs, include_trans, h4_threshold)
  gene_info$tissues <- NULL

  gene_info <- cleanup_api_object(gene_info)
  if (include_associations) {
    new_groups <- merge_associations(gene_info$coloc_groups, gene_info$rare_results, gene_info$associations)
    gene_info$coloc_groups <- new_groups$coloc_groups
    gene_info$rare_results <- new_groups$rare_results
  }

  return(gene_info)
}

#' @title Genes
#' @description Get specific genes from the API. The API returns collapsed/combined data for all requested genes.
#' @param gene_ids A vector of gene ids (1 or more)
#' @param include_associations A logical value specifying whether to include associations
#' (BETA, SE, P), defaults to FALSE
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs, defaults to FALSE
#' @param include_trans A logical value specifying whether to include trans genetic effects, defaults to TRUE
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs, defaults to 0.8
#' @return A list which contains the following elements:
#' \itemize{
#'   \item genes: gene metadata for the requested genes
#'   \item coloc_groups: a dataframe containing information about which studies have coloc results for all genes
#'   \item study_extractions: a dataframe containing the study extractions for all genes
#'   \item rare_results: a dataframe containing the rare results for all genes
#' }
#' @details
#' The dataframes returned by this function are as follows:
#' @inheritSection coloc_groups_doc coloc_groups_dataframe
#' @inheritSection study_extractions_doc study_extractions_dataframe
#' @inheritSection rare_results_doc rare_results_dataframe
#' @inheritSection coloc_pairs_doc coloc_pairs_dataframe
#' @export
genes <- function(gene_ids,
  include_associations = FALSE,
  include_coloc_pairs = FALSE,
  include_trans = TRUE,
  h4_threshold = 0.8
) {
  if (is.null(gene_ids) || length(gene_ids) == 0) stop("gene_ids is required")
  if (any(is.na(gene_ids))) stop("gene_ids must not contain NA values")
  if (length(gene_ids) > 10) stop("gene_ids must contain at most 10 gene ids")

  gene_ids <- unique(gene_ids)

  result <- specific_genes_api(
    gene_ids,
    include_associations = include_associations,
    include_coloc_pairs = include_coloc_pairs,
    include_trans = include_trans,
    h4_threshold = h4_threshold
  )

  if (include_associations && !is.null(result$associations)) {
    new_groups <- merge_associations(result$coloc_groups, result$rare_results, result$associations)
    result$coloc_groups <- new_groups$coloc_groups
    result$rare_results <- new_groups$rare_results
  }

  return(result)
}
