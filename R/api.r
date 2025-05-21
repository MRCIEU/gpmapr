#' @title Select API
#' @description Select the API to use
#' @param api A character string specifying the API to use. Default is "production".
#' @return NULL
#' @export
select_api <- function(api = "production") {
  if (api == "dev") {
    options("gpmap_url" = "https://localhost:8000")
  } else if (api == "production") {
    options("gpmap_url" = "https://gpmap.opengwas.io/api")
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

#' @title Get All Traits
#' @description Get all traits from the API
#' @return A list containing all traits
#' @export
traits <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/traits")
  traits <- httr::GET(url)
  traits <- httr::content(traits, "text", encoding = "UTF-8")
  traits <- jsonlite::fromJSON(traits)
  traits
}

#' @title Get a Trait
#' @description Get a trait from the API
#' @param trait_id A character string specifying the trait ID
#' @param include_associations A logical value specifying whether to include associations
#' @return A list containing the trait
#' @export
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
#' @return A list containing the gene information
#' @export
genes <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/genes")
  genes <- httr::GET(url)
  genes <- httr::content(genes, "text", encoding = "UTF-8")
  genes <- jsonlite::fromJSON(genes)
  genes
}

#' @title Get a Gene
#' @description Get a gene from the API
#' @param gene_id A character string specifying the gene ID
#' @param include_associations A logical value specifying whether to include associations
#' @return A list containing the gene information
#' @export
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

#' @title Get Variants by RSID
#' @description Get variants from the API by RSID
#' @param rsids A character string specifying the RSID
#' @return A list containing the variants
#' @export
variants_by_rsid <- function(rsids) {
  rsids <- paste(rsids, collapse = "&rsids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?rsids=", rsids)
  variants <- httr::GET(url)
  variants <- httr::content(variants, "text", encoding = "UTF-8")
  variants <- jsonlite::fromJSON(variants)
  variants
}

#' @title Get Variants by SNP ID
#' @description Get variants from the API by SNP ID
#' @param snp_ids A character string specifying the SNP ID
#' @param include_associations A logical value specifying whether to include associations
#' @param p_value_threshold A numeric value specifying the p-value threshold
#' @return A list containing the variants
#' @export
variants_by_snp_id <- function(snp_ids, include_associations = FALSE, p_value_threshold = NULL) {
  snp_ids <- paste(snp_ids, collapse = "&snp_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?snp_ids=", snp_ids)
  return(get_variants_with_options(url, include_associations, p_value_threshold))
}

#' @title Get Variants by Variant ID
#' @description Get variants from the API by Variant ID
#' @param variants A character string specifying the Variant ID
#' @param include_associations A logical value specifying whether to include associations
#' @param p_value_threshold A numeric value specifying the p-value threshold
#' @return A list containing the variants
#' @export
variants_by_variant <- function(variants, include_associations = FALSE, p_value_threshold = NULL) {
  variants <- paste(variants, collapse = "&variants=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?variants=", variants)
  return(get_variants_with_options(url, include_associations, p_value_threshold))
}

#' @title Get Variants by GRange
#' @description Get variants from the API by GRange
#' @param chr A character string specifying the chromosome
#' @param start A numeric value specifying the start position
#' @param stop A numeric value specifying the stop position
#' @param include_associations A logical value specifying whether to include associations
#' @param p_value_threshold A numeric value specifying the p-value threshold
#' @return A list containing the variants
#' @export
variants_by_grange <- function(chr, start, stop, include_associations = FALSE, p_value_threshold = NULL) {
  url <- paste0(getOption("gpmap_url"), "/v1/variants?grange=", chr, ":", start, "-", stop)
  return(get_variants_with_options(url, include_associations, p_value_threshold))
}

get_variants_with_options <- function(url, include_associations = FALSE, p_value_threshold = NULL) {
  if (include_associations) {
    url <- paste0(url, "&include_associations=true")
  }
  if (!is.null(p_value_threshold)) {
    url <- paste0(url, "&p_value_threshold=", p_value_threshold)
  }
  variants <- httr::GET(url)
  variants <- httr::content(variants, "text", encoding = "UTF-8")
  variants <- jsonlite::fromJSON(variants)
  variants
}

#' @title Get a Variant
#' @description Get a variant from the API
#' @param snp_id A character string specifying the SNP ID
#' @return A list containing the variant
#' @export
variant <- function(snp_id) {
  url <- paste0(getOption("gpmap_url"), "/v1/variants/", snp_id)
  variant <- httr::GET(url)
  variant <- httr::content(variant, "text", encoding = "UTF-8")
  variant <- jsonlite::fromJSON(variant)
  variant
}

#' @title Get LD Proxies by Variant ID
#' @description Get LD proxies from the API by Variant ID
#' @param variants A character string specifying the Variant ID
#' @return A list containing the LD proxies
#' @export
ld_proxies_by_variant <- function(variants) {
  variants <- paste(variants, collapse = "&variants=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants/ld_proxies?variants=", variants)
  ld_proxies <- httr::GET(url)
  ld_proxies <- httr::content(ld_proxies, "text", encoding = "UTF-8")
  ld_proxies <- jsonlite::fromJSON(ld_proxies)
  ld_proxies
}

#' @title Get LD Proxies by SNP ID
#' @description Get LD proxies from the API by SNP ID
#' @param snp_ids A character string specifying the SNP ID
#' @return A list containing the LD proxies
#' @export
ld_proxies_by_snp_id <- function(snp_ids) {
  snp_ids <- paste(snp_ids, collapse = "&snp_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants/ld_proxies?snp_ids=", snp_ids)
  ld_proxies <- httr::GET(url)
  ld_proxies <- httr::content(ld_proxies, "text", encoding = "UTF-8")
  ld_proxies <- jsonlite::fromJSON(ld_proxies)
  ld_proxies
}

#' @title Get LD Matrix by Variant ID
#' @description Get LD matrix from the API by Variant ID
#' @param variants A character string specifying the Variant ID
#' @return A list containing the LD matrix
#' @export
ld_matrix_by_variant <- function(variants) {
  variants <- paste(variants, collapse = "&variants=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants/ld_matrix?variants=", variants)
  ld_matrix <- httr::GET(url)
  ld_matrix <- httr::content(ld_matrix, "text", encoding = "UTF-8")
  ld_matrix <- jsonlite::fromJSON(ld_matrix)
  ld_matrix
}

#' @title Get LD Matrix by SNP ID
#' @description Get LD matrix from the API by SNP ID
#' @param snp_ids A character string specifying the SNP ID
#' @return A list containing the LD matrix
#' @export
ld_matrix_by_snp_id <- function(snp_ids) {
  snp_ids <- paste(snp_ids, collapse = "&snp_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants/ld_matrix?snp_ids=", snp_ids)
  ld_matrix <- httr::GET(url)
  ld_matrix <- httr::content(ld_matrix, "text", encoding = "UTF-8")
  ld_matrix <- jsonlite::fromJSON(ld_matrix)
  ld_matrix
}