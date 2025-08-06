#' @title Get LD Matrix
#' @description Get LD matrix from the API by Variant ID
#' @param variant_ids A character string specifying the Variant ID.  Variant IDs can be SNP IDs or variant IDs.
#' @return A list containing the LD matrix
#' @inheritSection ld_doc ld dataframe
#' @export
ld_matrix <- function(variant_ids = c()) {
  if (is.numeric(variant_ids)) {
    ld_matrix <- ld_matrix_by_snp_id_api(variant_ids)
  }
  if (is.character(variant_ids)) {
    ld_matrix <- ld_matrix_by_variant_api(variant_ids)
  }

  return(ld_matrix)
}

#' @title Get LD Proxies
#' @description Get LD proxies from the API by Variant ID
#' @param variant_ids A character string specifying the Variant ID.  Variant IDs can be SNP IDs or variant IDs.
#' @return A list containing the LD proxies
#' @inheritSection ld_doc ld dataframe
#' @export
ld_proxies <- function(variant_ids = c()) {
  if (is.numeric(variant_ids)) {
    ld_proxies <- ld_proxies_by_snp_id_api(variant_ids)
  }
  if (is.character(variant_ids)) {
    ld_proxies <- ld_proxies_by_variant_api(variant_ids)
  }

  return(ld_proxies)
}