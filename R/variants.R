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
#'      where the name of each element is the study_id. Column names are uppercase (e.g. SNP, BP, BETA, SE, LBF_1).
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
#' @description Get specific variants from the API. The API accepts variant identifiers (snp_ids, rsids, or strings)
#' and returns collapsed/combined data. The API distinguishes between identifier types automatically.
#' Max 10 variants when expand=TRUE.
#' @param variants A vector of variant identifiers (snp_ids, rsids, or strings)
#' @param expand Logical. FALSE (default) returns minimal data. TRUE returns full VariantResponse (max 10)
#' @param include_associations Logical. Whether to include associations (BETA, SE, P). Only when expand=TRUE
#' @param include_coloc_pairs Logical. Whether to include coloc pairs. Only when expand=TRUE
#' @param h4_threshold Numeric. H4 threshold for coloc pairs, defaults to 0.8
#' @return A list which contains the following elements:
#' \itemize{
#'   \item variants: a dataframe containing the variants for all requested variants
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
variants <- function(
  variants,
  expand = FALSE,
  include_associations = FALSE,
  include_coloc_pairs = FALSE,
  h4_threshold = 0.8
) {
  if (is.null(variants) || length(variants) == 0) stop("variants is required")
  if (any(is.na(variants))) stop("variants must not contain NA values")
  if (expand && length(variants) > 10) stop("variants must contain at most 10 when expand=TRUE")
  if (!expand && length(variants) > 50) stop("variants must contain at most 50 when expand=FALSE")

  variants <- unique(variants)

  result <- variants_api(
    variants,
    expand = expand,
    include_associations = include_associations,
    include_coloc_pairs = include_coloc_pairs,
    h4_threshold = h4_threshold
  )

  if (expand && include_associations && is.data.frame(result$associations) && nrow(result$associations) > 0) {
    new_groups <- merge_associations(result$coloc_groups, result$rare_results, result$associations)
    result$coloc_groups <- new_groups$coloc_groups
    result$rare_results <- new_groups$rare_results
  }

  return(cleanup_api_object(result))
}