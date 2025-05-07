#' @title Select API
#' @description Select the API to use
#' @param api A character string specifying the API to use. Default is "production".
#' @return NULL
#' @export
select_api <- function(api="production") {
  if (api == "dev") {
    options("gpmap_url" = "http://localhost:8001")
  } else if (api == "production") {
    options("gpmap_url" = "http://gpmap.opengwas.io/api")
  } else {
    stop("Invalid API")
  }
}

#' @title Get API Health
#' @description Get the health status of the API
#' @return A list containing the health status
#' @export
health <- function() {
  url <- paste0(getOption("gpmap_url"), "/health")
  health <- httr::GET(url)
  health <- httr::content(health, "text", encoding = "UTF-8")
  health <- jsonlite::fromJSON(health)
  health
}


studies <- function(study_id=NULL) {
  message("Fetching studies from the API, may take a few minutes...")
  if(is.null(study_id)) {
    url <- paste0(getOption("gpmap_url"), "/v1/studies/")
  } else {
    url <- paste0(getOption("gpmap_url"), "/v1/studies/", study_id)
  }
  studies <- httr::GET(url)
  studies <- httr::content(studies, "text", encoding = "UTF-8")
  studies <- jsonlite::fromJSON(studies)
  studies
}


#' @title Get Gene Information
#' @description Get gene information from the API
#' @param symbol A character string specifying the gene symbol
#' @return A list containing the gene information
#' @export
genes <- function(symbol) {
  url <- paste0(getOption("gpmap_url"), "/v1/genes/", symbol)
  genes <- httr::GET(url)
  genes <- httr::content(genes, "text", encoding = "UTF-8")
  genes <- jsonlite::fromJSON(genes)
  genes
}


# /v1/variants?rsids=rs61770163&rsids=rs11240777


variants_by_rsids <- function(rsids) {
  rsids <- paste(rsids, collapse = "&rsids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?rsids=", rsids)
  variants <- httr::GET(url)
  variants <- httr::content(variants, "text", encoding = "UTF-8")
  variants <- jsonlite::fromJSON(variants)
  variants
}
