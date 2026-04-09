api_to_package_version <- list(
  "1.0.0" = c("0.0.0.9000"),
  "1.0.1" = c("0.0.1.0")
)

#' Check if a string is a GUID
#' @param id A string to check
#' @return A logical value indicating whether the string is a GUID (8-4-4-4-12 hex format)
#' @keywords internal
#' @noRd
is_guid <- function(id) {
  id <- as.character(id)
  result <- grepl("^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$", id)
  return(result)
}

#' Merge associations into coloc info
#'
#' @param coloc_groups A dataframe of coloc groups
#' @param rare_results A dataframe of rare results
#' @param associations A dataframe of associations
#' @return An updated list of coloc info with associations merged in
#' @keywords internal
#' @noRd
merge_associations <- function(coloc_groups, rare_results, associations) {
  if (is.null(associations)) return(coloc_groups = coloc_groups, rare_results = rare_results)

  join_keys <- c("study_id", "variant_id")
  if (!all(join_keys %in% colnames(associations))) {
    return(list(coloc_groups = coloc_groups, rare_results = rare_results))
  }
  associations <- dplyr::distinct(associations, dplyr::across(dplyr::all_of(join_keys)), .keep_all = TRUE)

  if (is.data.frame(coloc_groups) && nrow(coloc_groups) > 0) {
    coloc_groups <- dplyr::left_join(coloc_groups, associations, by = join_keys)
  }
  if (is.data.frame(rare_results) && nrow(rare_results) > 0) {
    rare_results <- dplyr::left_join(rare_results, associations, by = join_keys)
  }

  return(list(coloc_groups = coloc_groups, rare_results = rare_results))
}

#' Merge associations into coloc info for GWAS upload results
#'
#' Joins associations using study_id when study_id is not NA, or existing_study_id
#' when existing_study_id is not NA. GWAS uploads do not have rare results.
#'
#' @param coloc_groups A dataframe of coloc groups
#' @param associations A dataframe of associations
#' @return Coloc groups with associations merged in
#' @keywords internal
#' @noRd
merge_gwas_upload_associations <- function(coloc_groups, associations) {
  if (is.null(associations)) return(coloc_groups)
  if (!is.data.frame(coloc_groups) || nrow(coloc_groups) == 0) return(coloc_groups)

  study_rows <- coloc_groups[!is.na(coloc_groups[["study_id"]]), , drop = FALSE]
  existing_rows <- coloc_groups[
    is.na(coloc_groups[["study_id"]]) & !is.na(coloc_groups[["existing_study_id"]]),
    , drop = FALSE
  ]
  other_rows <- coloc_groups[
    is.na(coloc_groups[["study_id"]]) & is.na(coloc_groups[["existing_study_id"]]),
    , drop = FALSE
  ]
  parts <- list()

  assoc_cols <- c("beta", "se", "p", "eaf", "imputed")
  assoc_cols <- intersect(assoc_cols, colnames(associations))

  if (nrow(study_rows) > 0 && "study_id" %in% colnames(associations)) {
    assoc <- dplyr::distinct(associations, dplyr::across(dplyr::all_of(c("study_id", "variant_id"))), .keep_all = TRUE)
    assoc <- dplyr::select(assoc, dplyr::all_of(c("study_id", "variant_id", assoc_cols)))
    parts[[length(parts) + 1L]] <- dplyr::left_join(study_rows, assoc, by = c("study_id", "variant_id"))
  }
  if (nrow(existing_rows) > 0 && "existing_study_id" %in% colnames(associations)) {
    assoc <- dplyr::distinct(
      associations,
      dplyr::across(dplyr::all_of(c("existing_study_id", "variant_id"))),
      .keep_all = TRUE
    )
    assoc <- dplyr::select(assoc, dplyr::all_of(c("existing_study_id", "variant_id", assoc_cols)))
    parts[[length(parts) + 1L]] <- dplyr::left_join(existing_rows, assoc, by = c("existing_study_id", "variant_id"))
  }
  if (nrow(other_rows) > 0) {
    parts[[length(parts) + 1L]] <- other_rows
  }

  if (length(parts) > 0) {
    return(dplyr::bind_rows(parts))
  }
  return(coloc_groups)
}

#' Create a source url for a study
#'
#' @param study_name character string specifying the study name
#' @return A character string specifying the source url
#' @keywords internal
#' @noRd
create_source_url <- function(study_name) {
  replace_except_first_two <- function(x) {
    dash_positions <- gregexpr("-", x)[[1]]
    if (length(dash_positions) <= 2) {
      return(x)
    }
    result <- x
    for (i in 3:length(dash_positions)) {
      pos <- dash_positions[i]
      substr(result, pos, pos) <- "_"
    }
    return(result)
  }
  study_name <- replace_except_first_two(study_name)

  return(paste0("https://opengwas.io/datasets/", study_name))
}

#' Clean up an api object
#'
#' @param api_object A list of api info, as returned by the api
#' @return An updated list of api info with null, empty, and NA elements removed
#' @keywords internal
#' @noRd
cleanup_api_object <- function(api_object) {
  api_object <- api_object[!sapply(api_object, is.null)]
  api_object <- api_object[!sapply(api_object, length) == 0]
  return(api_object)
}