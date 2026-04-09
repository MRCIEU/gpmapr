#' @title Get Associations by SNP ID and Study ID
#' @description Get associations from the API by SNP id and study id
#' @param variant_ids A vector of numeric values specifying the SNP IDs
#' @param study_ids A vector of numeric values specifying the Study IDs
#' @return A dataframe containing the associations
#' @inheritSection associations_doc associations_dataframe
#' @export
associations <- function(variant_ids, study_ids) {
  associations <- associations_api(variant_ids, study_ids)
  return(associations$associations)
}