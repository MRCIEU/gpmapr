#' Merge associations into coloc info
#'
#' @param coloc_info A list of coloc info, as returned by the api
#' @return An updated list of coloc info with associations merged in
merge_associations <- function(coloc_info) {
  if (is.null(coloc_info$associations)) return(coloc_info)

  coloc_info$coloc_groups <- dplyr::left_join(
    coloc_info$coloc_groups,
    coloc_info$associations,
    by = c("study_id", "snp_id")
  )

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

#' Clean up an api object
#'
#' @param api_object A list of api info, as returned by the api
#' @return An updated list of api info with null, empty, and NA elements removed
cleanup_api_object <- function(api_object) {
  api_object <- api_object[!sapply(api_object, is.null)]
  api_object <- api_object[!sapply(api_object, length) == 0]
  return(api_object)
}