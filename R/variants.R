#' @title Variant
#' @description A collection of studies that are associated with a particular variant.
#' @param snp_id A character string specifying the SNP ID
#' @param coloc_group_threshold A character value specifying the group threshold for coloc groups, defaults to 'strong'
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
                    coloc_group_threshold = "strong",
                    include_coloc_pairs = FALSE,
                    h4_threshold = 0.8,
                    include_summary_stats = FALSE) {
  if (is.null(snp_id)) {
    stop("snp_id is required")
  }
  if (!coloc_group_threshold %in% c("strong", "moderate")) {
    stop("coloc_group_threshold must be either 'strong' or 'moderate'")
  }

  variant_info <- variant_api(snp_id, include_coloc_pairs = include_coloc_pairs, h4_threshold = h4_threshold)
  variant_info$coloc_groups <- dplyr::filter(variant_info$coloc_groups, group_threshold == coloc_group_threshold)

  if (include_summary_stats) {
    summary_stats <- variant_summary_stats_api(snp_id)
    variant_info$summary_stats <- summary_stats
  }

  variant_info <- cleanup_api_object(variant_info)
  variant_info <- merge_associations(variant_info)

  return(variant_info)
}