#' @title Variant
#' @description A collection of studies that are associated with a particular variant.
#' @param snp_id A character string specifying the SNP ID
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs
#' @param h4_threshold A numeric value specifying the cutoff for included coloc pairs, defaults to 0.8.
#'      Only used if include_coloc_pairs is TRUE.
#' @param include_summary_stats A logical value specifying whether to include summary stats
#' @return A list which contains the following elements:
#'  * variant: named list containing the variant information
#'  * coloc_groups: a dataframe containing information about which studies have coloc results for this variant
#'  * rare_results: a list of dataframes containing the rare variants
#'  * study_extractions: a list of dataframes containing the study extractions
#'  * summary_stats (optional): a list of dataframes containing the summary stats for each study,
#'      where the name of the dataframe is the study_extraction_id
#'  * coloc_pairs (optional): a dataframe containing information about which studies have coloc pairs for this variant
#'      where the study_extraction_a_id and study_extraction_b_id are the study_extraction_ids of the two studies.
#'      h4_threshold is the cutoff for included coloc pairs, defaults to 0.8
#' @details
#' The dataframes returned by this function are as follows:
#' @inheritSection coloc_groups_doc coloc_groups_dataframe
#' @inheritSection rare_results_doc rare_results_dataframe
#' @inheritSection study_extractions_doc study_extractions_dataframe
#' @inheritSection summary_statistics_doc summary_statistics_dataframe
#' @inheritSection coloc_pairs_doc coloc_pairs_dataframe
#' @export
variant <- function(snp_id,
                    include_coloc_pairs = FALSE,
                    h4_threshold = 0.8,
                    include_summary_stats = FALSE) {
  if (is.null(snp_id)) {
    stop("snp_id is required")
  }

  variant_info <- variant_api(snp_id, include_coloc_pairs = include_coloc_pairs, h4_threshold = h4_threshold)

  if (include_summary_stats) {
    summary_stats <- variant_summary_stats_api(snp_id)
    variant_info$summary_stats <- summary_stats
  }

  variant_info <- cleanup_api_object(variant_info)
  new_groups <- merge_associations(variant_info$coloc_groups, variant_info$rare_results, variant_info$associations)
  variant_info$coloc_groups <- new_groups$coloc_groups
  variant_info$rare_results <- new_groups$rare_results

  return(variant_info)
}

#' @title Variants
#' @description Get specific variants from the API. The API returns collapsed/combined data for all requested variants.
#' @param snp_ids A vector of SNP ids
#' @param include_associations A logical value specifying whether to include associations (BETA, SE, P)
#' @param p_value_threshold A numeric value specifying the p-value threshold for associations
#' @param expand Logical or character vector. FALSE (default) returns minimal data. TRUE for expanded data.
#' @return A list which contains the following elements:
#' \itemize{
#'   \item variants: a dataframe containing the variants for all requested SNPs
#'   \item coloc_groups: (if expanded) a dataframe containing the coloc groups for all variants
#'   \item study_extractions: (if expanded) a dataframe containing the study extractions for all variants
#'   \item rare_results: (if expanded) a dataframe containing the rare results for all variants
#' }
#' @details
#' The dataframes returned by this function are as follows:
#' @inheritSection coloc_groups_doc coloc_groups_dataframe
#' @inheritSection study_extractions_doc study_extractions_dataframe
#' @inheritSection rare_results_doc rare_results_dataframe
#' @inheritSection summary_statistics_doc summary_statistics_dataframe
#' @inheritSection coloc_pairs_doc coloc_pairs_dataframe
#' @export
variants <- function(snp_ids,
  include_associations = FALSE,
  p_value_threshold = NULL,
  expand = FALSE
) {
  if (is.null(snp_ids) || length(snp_ids) == 0) stop("snp_ids is required")
  if (any(is.na(snp_ids))) stop("snp_ids must not contain NA values")
  if (length(snp_ids) > 50) stop("snp_ids must contain at most 50 SNP ids")

  snp_ids <- unique(snp_ids)

  result <- variants_by_snp_id_api(
    snp_ids,
    include_associations = include_associations,
    p_value_threshold = p_value_threshold,
    expand = expand
  )
  return(result)
}