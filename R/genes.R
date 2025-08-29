#' @title Gene
#' @description A collection of studies that are associated with a particular gene.
#' @param gene_id A numeric value specifying the gene id
#' @param include_associations A logical value specifying whether to include associations
#' (BETA, SE, P), defaults to FALSE
#' @param coloc_group_threshold A character value specifying the group threshold for coloc groups, defaults to 'strong'
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs, defaults to FALSE
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
                 coloc_group_threshold = "strong",
                 include_coloc_pairs = FALSE,
                 h4_threshold = 0.8) {
  if (is.null(gene_id)) {
    stop("gene_id is required")
  }
  if (!coloc_group_threshold %in% c("strong", "moderate")) {
    stop("coloc_group_threshold must be either 'strong' or 'moderate'")
  }
  gene_info <- gene_api(gene_id, include_associations, include_coloc_pairs, h4_threshold)
  gene_info$coloc_groups <- dplyr::filter(gene_info$coloc_groups, group_threshold == coloc_group_threshold)
  gene_info$tissues <- NULL

  gene_info <- cleanup_api_object(gene_info)
  if (include_associations) gene_info <- merge_associations(gene_info)

  return(gene_info)
}