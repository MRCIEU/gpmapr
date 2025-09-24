#' @title Trait
#' @description A collection of studies that are associated with a particular phenotype.
#' A trait will include a common study and occasionally a rare study.
#' @param trait_id A numeric value specifying the trait id
#' @param include_associations A logical value specifying whether to include associations
#' (BETA, SE, P), defaults to FALSE
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
trait <- function(trait_id,
                  include_associations = FALSE,
                  include_coloc_pairs = FALSE,
                  h4_threshold = 0.8) {
  if (is.null(trait_id)) {
    stop("trait_id is required")
  }

  trait_info <- trait_api(trait_id, include_associations, include_coloc_pairs, h4_threshold)

  trait_info <- cleanup_api_object(trait_info)
  if (include_associations) trait_info <- merge_associations(trait_info)

  return(trait_info)
}