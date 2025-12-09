#' @title Traits
#' @description Get all traits from the API
#' @return A dataframe containing all traits with the following columns:
#'   \itemize{
#'     \item id: the id of the trait
#'     \item data_type: the data type of the trait
#'     \item trait: the internal string id of the trait
#'     \item trait_name: the name of the trait
#'     \item trait_category: the trait category of the trait
#'     \item variant_type: the type of variant
#'     \item sample_size: the sample size of the trait
#'     \item category: the category of the trait (continuous, categorical)
#'     \item ancestry: the ancestry of the trait
#'     \item heritability: the LDSC heritability score of the trait
#'     \item heritability_se: the standard error of the LDSC heritability score of the trait
#'     \item num_study_extractions: the number of study extractions for this trait
#'     \item num_coloc_groups: the number of coloc groups for this trait
#'     \item num_coloc_studies: the number of studies that have coloc results for this trait
#'     \item num_rare_results: the number of rare results for this trait
#'   }
#' @export
traits <- function() {
  traits <- traits_api()
  return(traits$traits)
}

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

  trait_info <- trait_api(trait_id, include_associations)

  if (include_coloc_pairs) {
    response <- trait_coloc_pairs_api(trait_id, h4_threshold)
    trait_info$coloc_pairs <- response
  }

  trait_info$trait$source_url <- create_source_url(trait_info$trait$trait)

  trait_info <- cleanup_api_object(trait_info)
  if (include_associations) trait_info <- merge_associations(trait_info)

  return(trait_info)
}

