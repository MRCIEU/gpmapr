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


#' Attach GWAS-upload associations and a stable trait_id
#'
#' Upload coloc rows are keyed by `gwas_upload_id` / own `study_id`, not the GUID
#' used to fetch them. Stamp that GUID onto the upload's own rows so downstream
#' filters of the form `trait_id == <user id>` work for numeric traits and uploads.
#'
#' @param upload A GWAS upload result from `get_gwas_api()`.
#' @param guid The GUID used to fetch the upload.
#' @param include_associations Whether associations were requested.
#' @return The upload with associations merged and `trait_id` stamped.
#' @keywords internal
#' @noRd
.decorate_gwas_upload <- function(upload, guid, include_associations = FALSE) {
  if (is.null(upload)) {
    return(upload)
  }
  if (isTRUE(include_associations) && !is.null(upload$associations)) {
    upload$coloc_groups <- merge_gwas_upload_associations(
      upload$coloc_groups,
      upload$associations
    )
  }

  upload_id <- NULL
  trait_name <- NULL
  if (!is.null(upload$trait)) {
    upload_id <- upload$trait$id
    trait_name <- upload$trait$name
    upload$trait$guid <- guid
    if (is.null(upload$trait$trait_name) && !is.null(trait_name)) {
      upload$trait$trait_name <- trait_name
    }
  }

  upload$coloc_groups <- .stamp_gwas_upload_trait_id(
    coloc_groups = upload$coloc_groups,
    lookup_id = guid,
    upload_id = upload_id,
    trait_name = trait_name
  )
  return(upload)
}


#' Stamp a lookup id onto a GWAS upload's own coloc_groups rows
#'
#' @param coloc_groups A coloc_groups dataframe.
#' @param lookup_id The id callers use (usually the upload GUID).
#' @param upload_id Optional numeric upload id (`trait$id` / `gwas_upload_id`).
#' @param trait_name Optional display name for own rows missing `trait_name`.
#' @return `coloc_groups` with `trait_id` set on the upload's own rows.
#' @keywords internal
#' @noRd
.stamp_gwas_upload_trait_id <- function(coloc_groups,
                                        lookup_id,
                                        upload_id = NULL,
                                        trait_name = NULL) {
  if (!is.data.frame(coloc_groups) || nrow(coloc_groups) == 0) {
    return(coloc_groups)
  }
  if (!"trait_id" %in% names(coloc_groups)) {
    coloc_groups$trait_id <- NA_character_
  }
  coloc_groups$trait_id <- as.character(coloc_groups$trait_id)

  own_rows <- rep(FALSE, nrow(coloc_groups))
  if (!is.null(upload_id)) {
    upload_key <- as.character(upload_id)
    own_rows <- own_rows |
      (!is.na(coloc_groups$trait_id) & coloc_groups$trait_id == upload_key)
    if ("gwas_upload_id" %in% names(coloc_groups)) {
      own_rows <- own_rows |
        (!is.na(coloc_groups$gwas_upload_id) &
          as.character(coloc_groups$gwas_upload_id) == upload_key)
    }
  }
  if ("study_id" %in% names(coloc_groups) && "existing_study_id" %in% names(coloc_groups)) {
    own_rows <- own_rows |
      (!is.na(coloc_groups$study_id) & is.na(coloc_groups$existing_study_id))
  }
  coloc_groups$trait_id[own_rows] <- as.character(lookup_id)

  if (!is.null(trait_name) && "trait_name" %in% names(coloc_groups)) {
    fill_name <- own_rows & (is.na(coloc_groups$trait_name) | coloc_groups$trait_name == "")
    coloc_groups$trait_name[fill_name] <- as.character(trait_name)[[1]]
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