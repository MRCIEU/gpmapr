api_to_package_version <- list(
  "1.0.0" = c("0.0.0.9000")
)


#' Merge associations into coloc info
#'
#' @param coloc_groups A dataframe of coloc groups
#' @param rare_results A dataframe of rare results
#' @param associations A dataframe of associations
#' @return An updated list of coloc info with associations merged in
merge_associations <- function(coloc_groups, rare_results, associations) {
  if (is.null(associations)) return(coloc_groups = coloc_groups, rare_results = rare_results)

  if (!is.null(coloc_groups) && nrow(coloc_groups) > 0) {
    coloc_groups <- dplyr::left_join(
      coloc_groups,
      associations,
      by = c("study_id", "snp_id")
    )
  }

  if (!is.null(rare_results) && nrow(rare_results) > 0) {
    rare_results <- dplyr::left_join(
      rare_results,
      associations,
      by = c("study_id", "snp_id")
    )
  }

  return(list(coloc_groups = coloc_groups, rare_results = rare_results))
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