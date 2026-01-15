api_to_package_version <- list(
  "1.0.0" = c("0.0.0.9000")
)


#' Merge associations into coloc info
#'
#' @param coloc_info A list of coloc info, as returned by the api
#' @return An updated list of coloc info with associations merged in
merge_associations <- function(coloc_info) {
  if (is.null(coloc_info$associations)) return(coloc_info)

  if (!is.null(coloc_info$coloc_groups) && nrow(coloc_info$coloc_groups) > 0) {
    coloc_info$coloc_groups <- dplyr::left_join(
      coloc_info$coloc_groups,
      coloc_info$associations,
      by = c("study_id", "snp_id")
    )
  }

  if (!is.null(coloc_info$rare_results) && nrow(coloc_info$rare_results) > 0) {
    coloc_info$rare_results <- dplyr::left_join(
      coloc_info$rare_results,
      coloc_info$associations,
      by = c("study_id", "snp_id")
    )
  }

  coloc_info$associations <- NULL
  return(coloc_info)
}

#' Create a source url for a study
#'
#' @param study_name character string specifying the study name
#' @return A character string specifying the source url
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
cleanup_api_object <- function(api_object) {
  api_object <- api_object[!sapply(api_object, is.null)]
  api_object <- api_object[!sapply(api_object, length) == 0]
  return(api_object)
}