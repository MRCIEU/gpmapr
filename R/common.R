api_to_package_version <- list(
  "1.0.0" = c("0.0.0.9000"),
  "1.0.1" = c("0.0.1.0", "0.0.2.0", "0.1.0", "0.1.1"),
  "1.1.0" = c("0.1.1")
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
#' @param study_extractions A dataframe of study extractions
#' @return A list of `coloc_groups`, `rare_results` and `study_extractions`, each with the association
#'   of its study at its variant joined on
#' @keywords internal
#' @noRd
merge_associations <- function(coloc_groups, rare_results, associations, study_extractions = NULL) {
  return(list(
    coloc_groups = join_associations(coloc_groups, associations),
    rare_results = join_associations(rare_results, associations),
    study_extractions = join_associations(study_extractions, associations)
  ))
}

# Association columns joined onto coloc groups, rare results and study extractions.
association_value_columns <- c("beta", "se", "p", "eaf", "imputed")

#' Join associations onto a table by study and variant
#'
#' @param table A dataframe with a study id column and `variant_id`
#' @param associations A dataframe of associations
#' @param study_key The study id column of `table`, named by the matching column of `associations`
#' @return `table` with the association columns joined on, or unchanged when either has no rows or
#'   lacks the join columns
#' @keywords internal
#' @noRd
join_associations <- function(table, associations, study_key = c(study_id = "study_id")) {
  table_study <- unname(study_key)
  association_study <- names(study_key)
  if (!is.data.frame(table) || nrow(table) == 0 || !is.data.frame(associations) || nrow(associations) == 0) {
    return(table)
  }
  if (!all(c(table_study, "variant_id") %in% names(table)) ||
        !all(c(association_study, "variant_id") %in% names(associations))) {
    return(table)
  }
  value_columns <- setdiff(intersect(association_value_columns, names(associations)), names(table))
  associations <- associations[!is.na(associations[[association_study]]), , drop = FALSE]
  associations <- associations[, c(association_study, "variant_id", value_columns), drop = FALSE]
  associations <- dplyr::distinct(
    associations, dplyr::across(dplyr::all_of(c(association_study, "variant_id"))), .keep_all = TRUE
  )
  by <- c(stats::setNames(association_study, table_study), variant_id = "variant_id")
  return(dplyr::left_join(table, associations, by = by))
}

#' Fold an api object's associations into its tables
#'
#' Joins `associations` onto `coloc_groups`, `rare_results` and `study_extractions`, then drops it.
#'
#' @param api_object A list of api info
#' @return `api_object` without `associations`
#' @keywords internal
#' @noRd
fold_associations <- function(api_object) {
  merged <- merge_associations(
    api_object$coloc_groups, api_object$rare_results, api_object$associations, api_object$study_extractions
  )
  for (name in names(merged)) {
    if (!is.null(merged[[name]])) api_object[[name]] <- merged[[name]]
  }
  api_object$associations <- NULL
  return(api_object)
}

#' Fold a GWAS upload's associations into its tables
#'
#' Upload associations carry the study either in `study_id` (for the upload's own studies, where it is
#' the upload's id) or in `existing_study_id` (for studies already in the map).
#'
#' @param api_object A list of GWAS upload info, as returned by the api
#' @return `api_object` without `associations`
#' @keywords internal
#' @noRd
fold_upload_associations <- function(api_object) {
  associations <- api_object$associations
  if (is.data.frame(associations) && nrow(associations) > 0) {
    api_object$coloc_groups <- merge_gwas_upload_associations(api_object$coloc_groups, associations)
    api_object$study_extractions <- join_associations(
      api_object$study_extractions, associations, c(existing_study_id = "study_id")
    )
    api_object$upload_study_extractions <- join_associations(
      api_object$upload_study_extractions, associations, c(study_id = "gwas_upload_id")
    )
  }
  api_object$associations <- NULL
  return(api_object)
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
  api_object <- standardise_column_names(api_object)
  if (!is.null(api_object$trait)) {
    api_object$trait <- standardise_record(api_object$trait)
  }
  api_object$coloc_groups <- add_upload_source(api_object$coloc_groups)
  api_object <- add_coloc_pair_details(api_object)
  return(api_object)
}

# Column names the API has used, mapped to the names gpmapr returns.  The package
# always returns the new name, filling it from the old one when an older API sends that.
legacy_column_names <- c(
  # API 1.0.1 -> 1.1.0
  trait = "pipeline_id",
  category = "variable_type",
  study_extraction_id_a = "upload_study_extraction_a_id",
  study_extraction_id_b = "upload_study_extraction_b_id",
  existing_study_extraction_id_a = "existing_study_extraction_a_id",
  existing_study_extraction_id_b = "existing_study_extraction_b_id"
)

# Renames that only apply to GWAS upload tables, which are the tables that have an
# `existing_study_extraction_id`: elsewhere `study_extraction_id` keeps its meaning.
legacy_upload_column_names <- c(
  # API 1.0.1 -> 1.1.0
  study_extraction_id = "upload_study_extraction_id"
)

#' Rename legacy API column names
#'
#' @param x A data frame, or a list of api info containing data frames
#' @return `x` with any legacy column names in its data frames (including data frame columns)
#'   renamed to their current names.  When the API sends both the old and the new name, the new
#'   column is kept.
#' @keywords internal
#' @noRd
standardise_column_names <- function(x) {
  if (is.data.frame(x)) {
    x <- rename_legacy_names(x, legacy_column_names)
    is_legacy_upload_table <- "existing_study_extraction_id" %in% names(x) &&
      !"upload_study_extraction_id" %in% names(x)
    if (is_legacy_upload_table) {
      x <- rename_legacy_names(x, legacy_upload_column_names)
    }
    for (column in names(x)[vapply(x, is.data.frame, logical(1))]) {
      x[[column]] <- standardise_column_names(x[[column]])
    }
    return(x)
  }
  if (is.list(x)) {
    x[] <- lapply(x, standardise_column_names)
  }
  return(x)
}

#' Rename legacy API field names in a record
#'
#' @param x A named list describing one thing, e.g. the `trait` of `trait()` or `get_gwas()`
#' @return `x` with legacy field names renamed, in it and in the records and data frames it contains
#' @keywords internal
#' @noRd
standardise_record <- function(x) {
  if (!is.list(x) || is.data.frame(x)) {
    return(standardise_column_names(x))
  }
  x <- rename_legacy_names(x, legacy_column_names)
  for (field in names(x)[vapply(x, is.list, logical(1))]) {
    x[[field]] <- standardise_record(x[[field]])
  }
  return(x)
}

# Rename the names of `x` found in `legacy_names` (old name -> new name), keeping the new one
# when `x` has both.
rename_legacy_names <- function(x, legacy_names) {
  for (old_name in intersect(names(x), names(legacy_names))) {
    new_name <- legacy_names[[old_name]]
    if (new_name %in% names(x)) {
      x[[old_name]] <- NULL
    } else {
      names(x)[names(x) == old_name] <- new_name
    }
  }
  return(x)
}

#' Say which side of a GWAS upload each coloc group row comes from
#'
#' @param coloc_groups A dataframe of coloc groups
#' @return For GWAS upload coloc groups, `coloc_groups` with a `source` column ("upload" when the row is
#'   one of the upload's study extractions, "existing" when it is a study already in the map) after its
#'   study extraction id columns; other coloc groups unchanged
#' @keywords internal
#' @noRd
add_upload_source <- function(coloc_groups) {
  id_columns <- c("existing_study_extraction_id", "upload_study_extraction_id")
  if (!is.data.frame(coloc_groups) || !all(id_columns %in% names(coloc_groups))) {
    return(coloc_groups)
  }
  is_upload <- !is.na(coloc_groups$upload_study_extraction_id)
  is_existing <- !is.na(coloc_groups$existing_study_extraction_id)
  coloc_groups$source <- ifelse(is_upload, "upload", ifelse(is_existing, "existing", NA_character_))
  last_id_column <- names(coloc_groups)[max(match(id_columns, names(coloc_groups)))]
  coloc_groups <- dplyr::relocate(coloc_groups, "source", .after = dplyr::all_of(last_id_column))
  return(coloc_groups)
}

# Study extraction columns copied onto each side of a coloc pair, as <column>_a and <column>_b.
coloc_pair_detail_columns <- c("trait_name", "data_type", "tissue", "gene")

#' Add the details of each side of a coloc pair
#'
#' Joins the study extractions referred to by each coloc pair onto it, so each pair says which traits
#' it links without a manual join.  For GWAS upload pairs, `source_a` / `source_b` say whether that
#' side is the upload ("upload", id in `upload_study_extraction_*_id`) or a study already in the map
#' ("existing", id in `existing_study_extraction_*_id`).
#'
#' @param api_object A list of api info containing `coloc_pairs` and the `study_extractions` it refers to
#' @return `api_object` with `<column>_a` / `<column>_b` columns added to `coloc_pairs` after its id
#'   columns, `NA` where a study extraction isn't in the response
#' @keywords internal
#' @noRd
add_coloc_pair_details <- function(api_object) {
  pairs <- api_object$coloc_pairs
  if (!is.data.frame(pairs) || nrow(pairs) == 0) {
    return(api_object)
  }
  existing <- as_table(api_object$study_extractions)
  upload_names <- upload_extraction_names(api_object)
  is_upload_pairs <- "upload_study_extraction_a_id" %in% names(pairs)

  id_columns <- grep("study_extraction_[ab]_id$", names(pairs), value = TRUE)
  new_columns <- c()
  for (side in c("a", "b")) {
    existing_id_column <- if (is_upload_pairs) "existing_study_extraction_" else "study_extraction_"
    existing_id <- pairs[[paste0(existing_id_column, side, "_id")]]
    if (is.null(existing_id)) next
    row <- match(existing_id, existing$id)
    details <- lapply(coloc_pair_detail_columns, function(column) {
      if (!column %in% names(existing)) {
        return(rep(NA, nrow(pairs)))
      }
      return(existing[[column]][row])
    })
    names(details) <- paste0(coloc_pair_detail_columns, "_", side)

    if (is_upload_pairs) {
      upload_id <- pairs[[paste0("upload_study_extraction_", side, "_id")]]
      is_upload <- !is.na(upload_id)
      trait_name_column <- paste0("trait_name_", side)
      details[[trait_name_column]][is_upload] <- unname(upload_names[as.character(upload_id[is_upload])])
      source <- ifelse(is_upload, "upload", ifelse(is.na(existing_id), NA_character_, "existing"))
      details <- c(stats::setNames(list(source), paste0("source_", side)), details)
    }
    pairs[names(details)] <- details
    new_columns <- c(new_columns, names(details))
  }

  if (length(id_columns) > 0 && length(new_columns) > 0) {
    pairs <- dplyr::relocate(pairs, dplyr::all_of(new_columns), .after = dplyr::all_of(utils::tail(id_columns, 1)))
  }
  api_object$coloc_pairs <- pairs
  return(api_object)
}

# The name of the GWAS upload each upload study extraction belongs to, keyed by its id: from the
# upload's coloc groups (which cover compared uploads), else the name of the upload being fetched.
upload_extraction_names <- function(api_object) {
  upload <- as_table(api_object$upload_study_extractions)
  if (is.null(upload$id)) {
    return(character(0))
  }
  names <- rep(NA_character_, nrow(upload))
  trait <- api_object$trait
  if (!is.null(trait$id) && !is.null(trait$name) && "gwas_upload_id" %in% names(upload)) {
    names[upload$gwas_upload_id %in% trait$id] <- as.character(trait$name)
  }
  groups <- as_table(api_object$coloc_groups)
  if (all(c("upload_study_extraction_id", "trait_name") %in% names(groups))) {
    group_names <- groups$trait_name[match(upload$id, groups$upload_study_extraction_id)]
    names <- ifelse(is.na(group_names), names, group_names)
  }
  return(stats::setNames(names, upload$id))
}

# A data frame from an api element that may be a data frame, a list of data frames or NULL.
as_table <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.list(x) && length(x) > 0 && all(vapply(x, is.data.frame, logical(1)))) {
    return(dplyr::bind_rows(x))
  }
  return(data.frame())
}
