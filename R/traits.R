#' @title All traits
#' @description Get all traits from the API
#' @return A dataframe containing all traits with the following columns:
#'   \itemize{
#'     \item id: the id of the trait
#'     \item data_type: the data type of the trait
#'     \item pipeline_id: the string id the GPMap pipeline uses for the trait, e.g. its OpenGWAS dataset id
#'     \item trait_name: the name of the trait
#'     \item trait_category: the trait category of the trait
#'     \item variant_type: the type of variant
#'     \item sample_size: the sample size of the trait
#'     \item variable_type: the type of variable the trait is (continuous, categorical)
#'     \item ancestry: the ancestry of the trait
#'     \item heritability: the LDSC heritability score of the trait
#'     \item heritability_se: the standard error of the LDSC heritability score of the trait
#'     \item num_study_extractions: the number of study extractions for this trait
#'     \item num_coloc_groups: the number of coloc groups for this trait
#'     \item num_coloc_studies: the number of studies that have coloc results for this trait
#'     \item num_rare_results: the number of rare results for this trait
#'   }
#' @examplesIf gpmap_available()
#' traits <- all_traits()
#' head(traits[, c("id", "trait_name")])
#' @export
all_traits <- function() {
  traits <- traits_api()
  return(traits$traits)
}

#' @title Trait duplicates
#' @description Get all traits that have been marked as a duplicate of another trait
#' @return A dataframe in one of two shapes, depending on the API version:
#'   \itemize{
#'     \item one row per duplicated trait, with `trait_id`, `trait_name`,
#'       `parent_trait_id` and `parent_trait_name` (the trait it duplicates)
#'     \item one row per trait in a duplicate group, with `trait_id`,
#'       `trait_name`, `duplicate_group` and `is_representative` (the trait
#'       the group's other members duplicate), plus `match_rules` and the
#'       proxy-group columns
#'   }
#'   [remove_duplicate_traits()] accepts either.
#' @examplesIf gpmap_available()
#' head(trait_duplicates())
#' @export
trait_duplicates <- function() {
  duplicates <- trait_duplicates_api()
  return(duplicates)
}

#' @title Remove duplicate traits
#' @description Drop rows for traits marked as a duplicate of another trait, keeping the parent trait.
#' A trait in `keep_trait_ids` (e.g. the target trait) is never removed; if it is itself a duplicate,
#' its parent is removed instead, so only one copy of the measurement remains.
#' @param data A dataframe with a `trait_id` column, e.g. `trait()$coloc_groups`
#' @param duplicates The output of `trait_duplicates()`, in either of its shapes: duplicate
#'   rows with a `parent_trait_id`, or `duplicate_group` rows where each group's
#'   `is_representative` trait is the parent of the others. Fetched from the API if not supplied.
#' @param keep_trait_ids Trait ids that must not be removed
#' @return `data` without the duplicate traits' rows
#' @export
remove_duplicate_traits <- function(data, duplicates = trait_duplicates(), keep_trait_ids = NULL) {
  if (!is.data.frame(data) || !"trait_id" %in% names(data)) {
    stop("data must be a dataframe with a trait_id column")
  }
  keep_trait_ids <- as.character(keep_trait_ids)
  pairs <- .duplicate_parent_pairs(duplicates)
  duplicate_ids <- pairs$trait_id
  parent_ids <- pairs$parent_trait_id

  drop_ids <- union(duplicate_ids, parent_ids[duplicate_ids %in% keep_trait_ids])
  drop_ids <- setdiff(drop_ids, keep_trait_ids)
  return(data[!as.character(data$trait_id) %in% drop_ids, , drop = FALSE])
}

# Duplicate -> parent pairs from either shape of trait_duplicates(): rows that
# name a parent_trait_id, or duplicate_group rows whose is_representative trait
# is the parent of every other member of the group.
.duplicate_parent_pairs <- function(duplicates) {
  if (all(c("trait_id", "parent_trait_id") %in% names(duplicates))) {
    return(data.frame(
      trait_id = as.character(duplicates$trait_id),
      parent_trait_id = as.character(duplicates$parent_trait_id),
      stringsAsFactors = FALSE
    ))
  }
  if (!all(c("trait_id", "duplicate_group", "is_representative") %in% names(duplicates))) {
    stop(
      "duplicates must have trait_id and parent_trait_id, or trait_id, ",
      "duplicate_group and is_representative"
    )
  }
  representative <- duplicates$is_representative %in% TRUE
  parent_of_group <- stats::setNames(
    as.character(duplicates$trait_id[representative]),
    as.character(duplicates$duplicate_group[representative])
  )
  members <- duplicates[!representative & !is.na(duplicates$duplicate_group), , drop = FALSE]
  parent <- unname(parent_of_group[as.character(members$duplicate_group)])
  return(data.frame(
    trait_id = as.character(members$trait_id[!is.na(parent)]),
    parent_trait_id = parent[!is.na(parent)],
    stringsAsFactors = FALSE
  ))
}

#' @title Trait
#' @description A collection of studies that are associated with a particular phenotype.
#' A trait will include a common study and occasionally a rare study.
#' When trait_id is a GUID (from GWAS upload), fetches the upload result instead, as [get_gwas()] does: see
#' its documentation for how upload results mark the upload's own rows.
#' @param trait_id A numeric value or GUID (from GWAS upload) specifying the trait id
#' @param include_associations A logical value specifying whether to add each study's association at its
#' variant (beta, se, p, eaf, imputed) to the coloc_groups, rare_results and study_extractions, defaults to FALSE
#' @param include_full_associations A logical value specifying whether to include the stored associations
#' of the trait's studies and the studies that colocalise with it, at the trait's coloc group, rare result and
#' study extraction variants (see "Which associations are stored"), defaults to FALSE
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs, defaults to FALSE
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs, defaults to 0.8
#' @return A list which contains the following elements:
#' \itemize{
#'   \item trait: A list containing metadata about the trait, including its `pipeline_id`, `source_url`,
#' and the common and rare studies associated with it
#'   \item coloc_groups: a dataframe containing information about which studies have coloc results for this trait.
#' See below for details.
#'   \item study_extractions: a list of dataframes containing the study extractions for this trait.
#' See below for details.
#'   \item rare_results: (optional) a list of dataframes containing the rare results for this trait
#'   \item coloc_pairs: (optional) a dataframe containing all pairwise coloc results for this trait.
#'   \item full_associations: (optional) a dataframe of associations (`variant_id`, `study_id`, `beta`, `se`,
#' `p`, ...): see "Which associations are stored" below for which variants and studies it covers.
#' }
#' See below for details.
#' @details
#' The dataframes returned by this function are as follows:
#' @inheritSection coloc_groups_doc coloc_groups_dataframe
#' @inheritSection study_extractions_doc study_extractions_dataframe
#' @inheritSection rare_results_doc rare_results_dataframe
#' @inheritSection coloc_pairs_doc coloc_pairs_dataframe
#' @inheritSection variants_doc variants_dataframe
#' @inheritSection associations_doc Which associations are stored
#' @examplesIf gpmap_available()
#' result <- trait(5020)
#' head(result$coloc_groups)
#' @export
trait <- function(trait_id,
                  include_associations = FALSE,
                  include_full_associations = FALSE,
                  include_coloc_pairs = FALSE,
                  h4_threshold = 0.8) {
  if (is.null(trait_id)) {
    stop("trait_id is required")
  }

  if (is_guid(trait_id)) {
    trait_info <- get_gwas_api(trait_id, include_summary_stats = FALSE, include_associations = include_associations)
    trait_info <- .decorate_gwas_upload(trait_info, guid = trait_id)
  } else {
    trait_info <- trait_api(trait_id, include_associations)

    if (include_full_associations) {
      trait_info$full_associations <- mask_padded_associations(trait_associations_full_api(trait_id))
    }

    trait_info <- fold_associations(trait_info)
  }

  trait_info <- cleanup_api_object(trait_info)
  if (!is.null(trait_info$trait$pipeline_id)) {
    trait_info$trait$source_url <- create_source_url(trait_info$trait$pipeline_id)
  }

  # Added after the clean up: the trait coloc pairs endpoint doesn't return the study extractions
  # its pairs refer to, so their details can't be joined on (see add_coloc_pair_details()).
  if (!is_guid(trait_id) && include_coloc_pairs) {
    coloc_pairs <- trait_coloc_pairs_api(trait_id, h4_threshold)
    if (length(coloc_pairs) > 0) trait_info$coloc_pairs <- coloc_pairs
  }
  return(trait_info)
}


#' @title Traits
#' @description Get specific traits from the API. The API returns collapsed/combined data for all requested traits.
#' When a trait ID is a GUID (from GWAS upload), fetches the upload result instead.
#' @param trait_ids A vector of trait ids (numeric) or GUIDs (from GWAS upload)
#' @param include_associations A logical value specifying whether to add each study's association at its
#' variant (beta, se, p, eaf, imputed) to the coloc_groups, rare_results and study_extractions, defaults to FALSE
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
#' @examplesIf gpmap_available()
#' result <- traits(c(4405, 4872))
#' result$traits
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
    result <- fold_associations(result)
  } else {
    result <- list(traits = NULL, coloc_groups = NULL, study_extractions = NULL, rare_results = NULL)
  }

  if (length(guid_ids) > 0) {
    upload_results <- lapply(guid_ids, function(guid) {
      return(.decorate_gwas_upload(get_gwas_api(guid, include_associations = include_associations), guid = guid))
    })
    for (upload in upload_results) {
      result <- merge_trait_result(result, upload)
    }
  }

  result <- cleanup_api_object(result)

  # Added after the clean up, as in trait(): the trait coloc pairs endpoint doesn't return the study
  # extractions its pairs refer to, so their details can't be joined on.
  if (include_coloc_pairs && length(numeric_ids) > 0) {
    coloc_pairs <- lapply(numeric_ids, trait_coloc_pairs_api, h4_threshold = h4_threshold) |>
      dplyr::bind_rows() |>
      dplyr::distinct()
    if (length(coloc_pairs) > 0) result$coloc_pairs <- coloc_pairs
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
  if (is.list(upload$trait)) {
    # One row per upload, from its single-value fields, named like the other traits.
    is_single_value <- vapply(upload$trait, function(x) is.atomic(x) && length(x) == 1, logical(1))
    trait_df <- as.data.frame(upload$trait[is_single_value], stringsAsFactors = FALSE)
    if ("name" %in% names(trait_df)) trait_df$trait_name <- trait_df$name
    result$traits <- dplyr::bind_rows(result$traits, trait_df)
  }
  for (table in c("coloc_groups", "study_extractions", "rare_results")) {
    if (!is.data.frame(upload[[table]]) || nrow(upload[[table]]) == 0) next
    if (is.data.frame(result[[table]]) && nrow(result[[table]]) > 0 && !"source" %in% names(result[[table]])) {
      result[[table]]$source <- "existing"
    }
    # An upload's own coloc_groups rows carry its GUID as trait_id (.decorate_gwas_upload()), so trait_id is
    # character there: make it character on both sides before binding.
    if (is.character(upload[[table]]$trait_id) && is.data.frame(result[[table]]) &&
          "trait_id" %in% names(result[[table]])) {
      result[[table]]$trait_id <- as.character(result[[table]]$trait_id)
    }
    result[[table]] <- dplyr::bind_rows(result[[table]], upload[[table]])
  }
  return(result)
}
