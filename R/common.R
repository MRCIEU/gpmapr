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
  coloc_info$study_extractions <- dplyr::left_join(
    coloc_info$study_extractions,
    coloc_info$associations,
    by = c("study_id", "snp_id")
  )
  coloc_info$associations <- NULL
  return(coloc_info)
}