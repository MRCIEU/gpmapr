api_to_package_version <- list(
  "1.0.0" = c("0.0.0.9000"),
  "1.0.1" = c("0.0.1.0", "0.0.2.0", "0.1.0"),
  "1.1.0" = c("1.0.0")
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
  associations <- mask_padded_associations(associations)
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
#' Joins `associations` onto the upload's `coloc_groups` and `study_extractions`, then drops it.  Rows
#' are keyed by `source`: existing studies by `study_id`, the uploads' own studies by `gwas_upload_id`.
#'
#' @param api_object A list of GWAS upload info, as returned by the api
#' @return `api_object` without `associations`
#' @keywords internal
#' @noRd
fold_upload_associations <- function(api_object) {
  associations <- api_object$associations
  if (is.data.frame(associations) && nrow(associations) > 0) {
    api_object$coloc_groups <- join_upload_associations(api_object$coloc_groups, associations)
    api_object$study_extractions <- join_upload_associations(api_object$study_extractions, associations)
  }
  api_object$associations <- NULL
  return(api_object)
}

# Join GWAS upload associations onto a table with a `source` column, keeping its row order.
join_upload_associations <- function(table, associations) {
  if (!is.data.frame(table) || nrow(table) == 0 || !"source" %in% names(table)) {
    return(table)
  }
  table$.row <- seq_len(nrow(table))
  is_upload <- table$source %in% "upload"
  existing <- join_associations(
    table[!is_upload, , drop = FALSE], associations[associations$source %in% "existing", , drop = FALSE]
  )
  upload <- join_associations(
    table[is_upload, , drop = FALSE], associations[associations$source %in% "upload", , drop = FALSE],
    c(gwas_upload_id = "gwas_upload_id")
  )
  table <- dplyr::bind_rows(existing, upload)
  table <- table[order(table$.row), , drop = FALSE]
  table$.row <- NULL
  rownames(table) <- NULL
  return(table)
}

#' Treat placeholder associations as missing
#'
#' The GPMap pipeline stores `beta = 0, se = 1, p = 1` (with `imputed = TRUE`) for variants a study
#' did not report, with `eaf` taken from the LD reference panel.  These rows are kept, with `beta`,
#' `se`, `p` and `eaf` set to `NA`.
#'
#' @param associations A dataframe of associations
#' @return `associations` with placeholder values set to `NA`
#' @keywords internal
#' @noRd
mask_padded_associations <- function(associations) {
  if (!is.data.frame(associations) || !all(c("beta", "se") %in% names(associations))) {
    return(associations)
  }
  padded <- associations$beta %in% 0 & associations$se %in% 1
  if ("p" %in% names(associations)) padded <- padded & associations$p %in% 1
  if ("imputed" %in% names(associations)) padded <- padded & associations$imputed %in% TRUE
  for (column in intersect(c("beta", "se", "p", "eaf"), names(associations))) {
    associations[[column]][padded] <- NA
  }
  return(associations)
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
#' @return An updated list of api info with null and empty elements removed, and coloc pairs
#'   described (see `add_coloc_pair_details()`)
#' @keywords internal
#' @noRd
cleanup_api_object <- function(api_object) {
  api_object <- api_object[!sapply(api_object, is.null)]
  api_object <- api_object[!sapply(api_object, length) == 0]
  api_object <- add_coloc_pair_details(api_object)
  return(api_object)
}

# Study extraction columns copied onto each side of a coloc pair, as <column>_a and <column>_b.
coloc_pair_detail_columns <- c("trait_name", "data_type", "tissue", "gene")

#' Add the details of each side of a coloc pair
#'
#' Joins the study extraction referred to by each side of a coloc pair onto it, so each pair says which
#' traits it links without a manual join.  A side is identified by `study_extraction_*_id` together with
#' `source_*` (GWAS upload results, where upload and existing study extractions have separate ids), or by
#' `study_extraction_*_id` alone (other results, where every study extraction is an existing one).
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
  extractions <- as_table(api_object$study_extractions)
  extraction_key <- paste(row_sources(extractions), extractions$id)

  new_columns <- c()
  for (side in c("a", "b")) {
    pair_id <- pairs[[paste0("study_extraction_", side, "_id")]]
    if (is.null(pair_id)) next
    pair_source <- pairs[[paste0("source_", side)]]
    if (is.null(pair_source)) pair_source <- rep("existing", nrow(pairs))
    row <- match(paste(pair_source, pair_id), extraction_key)
    details <- lapply(coloc_pair_detail_columns, function(column) {
      if (!column %in% names(extractions)) {
        return(rep(NA, nrow(pairs)))
      }
      return(extractions[[column]][row])
    })
    names(details) <- paste0(coloc_pair_detail_columns, "_", side)
    pairs[names(details)] <- details
    new_columns <- c(new_columns, names(details))
  }

  id_columns <- grep("^(study_extraction_[ab]_id|source_[ab])$", names(pairs), value = TRUE)
  if (length(id_columns) > 0 && length(new_columns) > 0) {
    pairs <- dplyr::relocate(pairs, dplyr::all_of(new_columns), .after = dplyr::all_of(utils::tail(id_columns, 1)))
  }
  api_object$coloc_pairs <- pairs
  return(api_object)
}

# The source of each row of a table: its `source` column, or "existing" for tables without one.
row_sources <- function(table) {
  if ("source" %in% names(table)) {
    return(table$source)
  }
  return(rep("existing", nrow(table)))
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
