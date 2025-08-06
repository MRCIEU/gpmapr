#' Trait: a collection of studies that are associated with a particular phenotype.
#' A trait will include a common study and occasionally a rare study.
#' @param trait_id A numeric value specifying the trait id
#' @param include_associations A logical value specifying whether to include associations
#' (BETA, SE, P), defaults to FALSE
#' @param coloc_group_threshold A character value specifying the group threshold for coloc groups, defaults to 'strong'
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs, defaults to FALSE
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs, defaults to 0.8
#' @return A list which contains the following elements:
#' \itemize{
#'   \item trait: A list containing metadata about the trait,
#' including common and rare studies associated with the trait
#'   \item coloc_groups: a dataframe containing information about which studies have coloc results for this trait.
#' See below for details.
#'   \item study_extractions: a list of dataframes containing the study extractions for this trait.
#' See below for details.
#'   \item rare_results: (optional) a list of dataframes containing the rare results for this trait
#'   \item coloc_pairs: (optional) a dataframe containing all pairwise coloc results for this trait.
#' See below for details.
#' @details
#' The dataframes returned by this function are as follows:
#' @inheritSection coloc_groups_doc coloc_groups dataframe
#' @inheritSection study_extractions_doc study_extractions dataframe
#' @inheritSection rare_results_doc rare_results dataframe
#' @inheritSection coloc_pairs_doc coloc_pairs dataframe
#' @export
trait <- function(trait_id,
                  include_associations = FALSE,
                  coloc_group_threshold = "strong",
                  include_coloc_pairs = FALSE,
                  h4_threshold = 0.8) {

  if (!coloc_group_threshold %in% c("strong", "moderate")) {
    stop("coloc_group_threshold must be either 'strong' or 'moderate'")
  }
  trait_info <- trait_api(trait_id, include_associations, include_coloc_pairs, h4_threshold)
  trait_info$upload_study_extractions <- NULL
  trait_info$coloc_groups <- dplyr::filter(trait_info$coloc_groups, group_threshold == coloc_group_threshold)
  return(trait_info)
}