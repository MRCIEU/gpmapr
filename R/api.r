#' @title Select API
#' @description Select the API to use
#' @param api A character string specifying the API to use. Default is "production".
#' @return NULL
#' @export
select_api <- function(api = "production") {
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


traits <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/traits/")
  traits <- httr::GET(url)
  traits <- httr::content(traits, "text", encoding = "UTF-8")
  traits <- jsonlite::fromJSON(traits)
  traits
}

trait <- function(trait_id, include_associations = FALSE) {
  url <- paste0(getOption("gpmap_url"), "/v1/traits/", trait_id)
  if (include_associations) {
    url <- paste0(url, "?include_associations=true")
  }
  trait <- httr::GET(url)
  trait <- httr::content(trait, "text", encoding = "UTF-8")
  trait <- jsonlite::fromJSON(trait)
  trait
}


#' @title Get Gene Information
#' @description Get gene information from the API
#' @param symbol A character string specifying the gene symbol
#' @return A list containing the gene information
#' @export
genes <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/genes/")
  genes <- httr::GET(url)
  genes <- httr::content(genes, "text", encoding = "UTF-8")
  genes <- jsonlite::fromJSON(genes)
  genes
}

gene <- function(gene_id, include_associations = FALSE) {
  url <- paste0(getOption("gpmap_url"), "/v1/genes/", gene_id)
  if (include_associations) {
    url <- paste0(url, "?include_associations=true")
  }
  gene <- httr::GET(url)
  gene <- httr::content(gene, "text", encoding = "UTF-8")
  gene <- jsonlite::fromJSON(gene)
  gene
}

variants_by_rsid <- function(rsids) {
  rsids <- paste(rsids, collapse = "&rsids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?rsids=", rsids)
  variants <- httr::GET(url)
  variants <- httr::content(variants, "text", encoding = "UTF-8")
  variants <- jsonlite::fromJSON(variants)
  variants
}

variants_by_snp_id <- function(snp_ids) {
  snp_ids <- paste(snp_ids, collapse = "&snp_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?snp_ids=", snp_ids)
  variants <- httr::GET(url)
  variants <- httr::content(variants, "text", encoding = "UTF-8")
  variants <- jsonlite::fromJSON(variants)
  variants
}

variants_by_variant_id <- function(variant_ids) {
  variant_ids <- paste(variant_ids, collapse = "&variant_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?variant_ids=", variant_ids)
  variants <- httr::GET(url)
  variants <- httr::content(variants, "text", encoding = "UTF-8")
  variants <- jsonlite::fromJSON(variants)
  variants
}

variants_by_grange <- function(chr, start, stop) {
  url <- paste0(getOption("gpmap_url"), "/v1/variants?grange=", chr, ":", start, "-", stop)
  variants <- httr::GET(url)
  variants <- httr::content(variants, "text", encoding = "UTF-8")
  variants <- jsonlite::fromJSON(variants)
  variants
}

variant <- function(snp_id) {
  url <- paste0(getOption("gpmap_url"), "/v1/variants/", snp_id)
  variant <- httr::GET(url)
  variant <- httr::content(variant, "text", encoding = "UTF-8")
  variant <- jsonlite::fromJSON(variant)
  variant
}

ld_proxies_by_variant_id <- function(variant_ids) {
  variant_ids <- paste(variant_ids, collapse = "&variant_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants/ld_proxies?variant_ids=", variant_ids)
  ld_proxies <- httr::GET(url)
  ld_proxies <- httr::content(ld_proxies, "text", encoding = "UTF-8")
  ld_proxies <- jsonlite::fromJSON(ld_proxies)
  ld_proxies
}

ld_proxies_by_snp_id <- function(snp_ids) {
  snp_ids <- paste(snp_ids, collapse = "&snp_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants/ld_proxies?snp_ids=", snp_ids)
  ld_proxies <- httr::GET(url)
  ld_proxies <- httr::content(ld_proxies, "text", encoding = "UTF-8")
  ld_proxies <- jsonlite::fromJSON(ld_proxies)
  ld_proxies
}

ld_matrix_by_variant_id <- function(variant_ids) {
  variant_ids <- paste(variant_ids, collapse = "&variant_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants/ld_matrix?variant_ids=", variant_ids)
  ld_matrix <- httr::GET(url)
  ld_matrix <- httr::content(ld_matrix, "text", encoding = "UTF-8")
  ld_matrix <- jsonlite::fromJSON(ld_matrix)
  ld_matrix
}

ld_matrix_by_snp_id <- function(snp_ids) {
  snp_ids <- paste(snp_ids, collapse = "&snp_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants/ld_matrix?snp_ids=", snp_ids)
  ld_matrix <- httr::GET(url)
  ld_matrix <- httr::content(ld_matrix, "text", encoding = "UTF-8")
  ld_matrix <- jsonlite::fromJSON(ld_matrix)
  ld_matrix
}