#' @title Get Associations by SNP ID and Study ID
#' @description Get associations from the API by SNP id and study id
#' @param variant_ids A vector of numeric values specifying the SNP IDs
#' @param study_ids A vector of numeric values specifying the Study IDs
#' @return A dataframe containing the associations
#' @inheritSection associations_doc associations_dataframe
#' @examplesIf gpmap_available()
#' associations(variant_ids = 80717, study_ids = c(926, 927))
#' @export
associations <- function(variant_ids, study_ids) {
  associations <- associations_api(variant_ids, study_ids)
  return(mask_padded_associations(associations$associations))
}
