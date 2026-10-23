#' @title LD Matrix
#' @description Get the LD between the given variants from the API.  The GPMap stores LD only between
#'   lead variants and their proxies, not between every pair of variants, so this returns the pairs among
#'   `variant_ids` where one is a lead variant and the other its proxy.
#' @param variant_ids A character string specifying the Variant ID.  Variant IDs can be SNP IDs or variant IDs.
#' @return A list containing `lds`, a dataframe of the stored LD pairs among `variant_ids`, with no rows
#'   when there are none
#' @inheritSection ld_doc ld_dataframe
#' @examplesIf gpmap_available()
#' ld_matrix(c(5553693, 5553688))
#' @export
ld_matrix <- function(variant_ids = c()) {
  if (is.numeric(variant_ids)) {
    ld_matrix <- ld_matrix_by_variant_id_api(variant_ids)
  }
  if (is.character(variant_ids)) {
    ld_matrix <- ld_matrix_by_variant_api(variant_ids)
  }

  ld_matrix$lds <- add_rsq(ld_matrix$lds)
  return(ld_matrix)
}

#' @title LD Proxies
#' @description Get LD proxies from the API by Variant ID
#' @param variant_ids A character string specifying the Variant ID.  Variant IDs can be SNP IDs or variant IDs.
#' @return A list containing the LD proxies
#' @inheritSection ld_doc ld_dataframe
#' @examplesIf gpmap_available()
#' proxies <- ld_proxies(c(5553693, 5553694))
#' head(proxies$lds)
#' @export
ld_proxies <- function(variant_ids = c()) {
  if (is.numeric(variant_ids)) {
    ld_proxies <- ld_proxies_by_variant_id_api(variant_ids)
  }
  if (is.character(variant_ids)) {
    ld_proxies <- ld_proxies_by_variant_api(variant_ids)
  }

  ld_proxies$lds <- add_rsq(ld_proxies$lds)
  return(ld_proxies)
}

# LD pairs with their r-squared, rounded as in search_gpmap().  No pairs give an empty dataframe
# with the same columns.
add_rsq <- function(lds) {
  if (!is.data.frame(lds) || nrow(lds) == 0) {
    return(data.frame(
      lead_variant_id = integer(0), proxy_variant_id = integer(0), ld_block_id = integer(0),
      r = numeric(0), rsq = numeric(0)
    ))
  }
  lds$rsq <- round(lds$r^2, 2)
  return(lds)
}
