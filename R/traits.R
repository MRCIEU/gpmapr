#' @title All traits
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
all_traits <- function() {
  traits <- traits_api()
  return(traits$traits)
}

#' @title Trait
#' @description A collection of studies that are associated with a particular phenotype.
#' A trait will include a common study and occasionally a rare study.
#' When trait_id is a GUID (from GWAS upload), fetches the upload result instead.
#' @param trait_id A numeric value or GUID (from GWAS upload) specifying the trait id
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

  if (is_guid(trait_id)) {
    trait_info <- get_gwas_api(trait_id, include_summary_stats = FALSE)
    if (include_associations && !is.null(trait_info$associations)) {
      new_groups <- merge_associations(trait_info$coloc_groups, trait_info$rare_results, trait_info$associations)
      trait_info$coloc_groups <- new_groups$coloc_groups
      trait_info$rare_results <- new_groups$rare_results
    }
    coloc_trait_id <- if (!is.null(trait_info$trait$id)) trait_info$trait$id else trait_id
    if (include_coloc_pairs && !is.null(coloc_trait_id) && !is_guid(coloc_trait_id)) {
      response <- trait_coloc_pairs_api(coloc_trait_id, h4_threshold)
      trait_info$coloc_pairs <- response
    }
  } else {
    trait_info <- trait_api(trait_id, include_associations)

    if (include_coloc_pairs) {
      response <- trait_coloc_pairs_api(trait_id, h4_threshold)
      trait_info$coloc_pairs <- response
    }

    if (include_associations) {
      new_groups <- merge_associations(trait_info$coloc_groups, trait_info$rare_results, trait_info$associations)
      trait_info$coloc_groups <- new_groups$coloc_groups
      trait_info$rare_results <- new_groups$rare_results
    }
  }

  if (!is.null(trait_info$trait) && !is.null(trait_info$trait$trait)) {
    trait_info$trait$source_url <- create_source_url(trait_info$trait$trait)
  }

  trait_info <- cleanup_api_object(trait_info)
  return(trait_info)
}


#' @title Traits
#' @description Get specific traits from the API. The API returns collapsed/combined data for all requested traits.
#' When a trait ID is a GUID (from GWAS upload), fetches the upload result instead.
#' @param trait_ids A vector of trait ids (numeric) or GUIDs (from GWAS upload)
#' @param include_associations A logical value specifying whether to include associations
#' (BETA, SE, P), defaults to FALSE
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs, defaults to FALSE.
#' Coloc pairs are fetched from a separate endpoint per trait.
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs, defaults to 0.8
#' @return A list which contains the following elements:
#' \itemize{
#'   \item traits: trait metadata for the requested traits
#'   \item coloc_groups: a dataframe containing information about which studies have coloc results for all traits.
#' See below for details.
#'   \item study_extractions: a dataframe containing the study extractions for all traits.
#' See below for details.
#'   \item rare_results: a dataframe containing the rare results for all traits
#'   \item coloc_pairs: (optional) a dataframe containing all pairwise coloc results for all traits.
#' }
#' @details
#' The dataframes returned by this function are as follows:
#' @inheritSection coloc_groups_doc coloc_groups_dataframe
#' @inheritSection study_extractions_doc study_extractions_dataframe
#' @inheritSection rare_results_doc rare_results_dataframe
#' @inheritSection coloc_pairs_doc coloc_pairs_dataframe
#' @export
traits <- function(trait_ids,
  include_associations = FALSE,
  include_coloc_pairs = FALSE,
  h4_threshold = 0.8
) {
  if (is.null(trait_ids) || length(trait_ids) == 0) stop("trait_ids is required")
  if (any(is.na(trait_ids))) stop("trait_ids must not contain NA values")
  if (length(trait_ids) > 10) stop("trait_ids must contain at most 10 trait ids")

  trait_ids <- unique(trait_ids)
  numeric_ids <- trait_ids[!is_guid(trait_ids)]
  guid_ids <- trait_ids[is_guid(trait_ids)]

  if (length(numeric_ids) > 0) {
    result <- specific_traits_api(
      numeric_ids,
      include_associations = include_associations,
      h4_threshold = h4_threshold
    )
    if (include_associations && is.data.frame(result$associations) && nrow(result$associations) > 0) {
      new_groups <- merge_associations(result$coloc_groups, result$rare_results, result$associations)
      result$coloc_groups <- new_groups$coloc_groups
      result$rare_results <- new_groups$rare_results
    }
    if (include_coloc_pairs && length(numeric_ids) > 0) {
      result$coloc_pairs <- lapply(numeric_ids, trait_coloc_pairs_api, h4_threshold = h4_threshold) |>
        dplyr::bind_rows() |>
        dplyr::distinct()
    }
  } else {
    result <- list(traits = NULL, coloc_groups = NULL, study_extractions = NULL, rare_results = NULL)
  }

  if (length(guid_ids) > 0) {
    upload_results <- lapply(guid_ids, get_gwas_api)
    for (upload in upload_results) {
      result <- merge_trait_result(result, upload)
    }
  }

  return(result)
}

#' Merge a trait result into a traits result
#'
#' @param result A list of trait results
#' @param upload A list of upload results
#' @return A list of trait results with the upload results merged in
#' @keywords internal
#' @noRd
merge_trait_result <- function(result, upload) {
  if (is.null(upload)) return(result)
  if (!is.null(upload$trait)) {
    trait_df <- as.data.frame(t(upload$trait))
    result$traits <- dplyr::bind_rows(result$traits, trait_df)
  }
  if (is.data.frame(upload$coloc_groups) && nrow(upload$coloc_groups) > 0) {
    result$coloc_groups <- dplyr::bind_rows(result$coloc_groups, upload$coloc_groups)
  }
  if (is.data.frame(upload$study_extractions) && nrow(upload$study_extractions) > 0) {
    result$study_extractions <- dplyr::bind_rows(result$study_extractions, upload$study_extractions)
  }
  if (is.data.frame(upload$rare_results) && nrow(upload$rare_results) > 0) {
    result$rare_results <- dplyr::bind_rows(result$rare_results, upload$rare_results)
  }
  return(result)
}
