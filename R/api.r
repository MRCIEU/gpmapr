timeout_seconds <- 360

#' @title Select API
#' @description Select the API to use
#' @param api A character string specifying the API to use. Default is "production".
#' @return NULL
#' @export
select_api <- function(api = "dev") {
  if (api == "dev") {
    api_option <- options("gpmap_url" = "http://localhost:8000")
  } else if (api == "production") {
    api_option <- options("gpmap_url" = "https://gpmap.opengwas.io/api")
  } else {
    stop("Invalid API")
  }
  return(api_option)
}

#' @title Get API Health
#' @description Get the health status of the API
#' @return A list containing the health status
#' @export
health_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/health")
  health <- httr::GET(url)
  health <- httr::content(health, "text", encoding = "UTF-8")
  health <- jsonlite::fromJSON(health)
  return(health)
}


#' @title Get API Version
#' @description Get the version of the API
#' @return A list containing the version
#' @export 
version_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/info/version")
  version <- httr::GET(url)
  version <- httr::content(version, "text", encoding = "UTF-8")
  version <- jsonlite::fromJSON(version)
  return(version$version)
}

#' @title Search Options API
#' @description Search options from the API
#' @return A list containing the search options
search_options_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/search/options")
  search_options <- httr::GET(url, httr::timeout(timeout_seconds))
  search_options <- httr::content(search_options, "text", encoding = "UTF-8")
  search_options <- jsonlite::fromJSON(search_options)
  return(search_options)
}

#' @title Search Variants API
#' @description Search variants from the API
#' @param query A character string specifying the query
#' @return A list containing the variants
search_variants_api <- function(query) {
  url <- paste0(getOption("gpmap_url"), "/v1/search/variant/", query)
  search_variants <- httr::GET(url, httr::timeout(timeout_seconds))
  search_variants <- httr::content(search_variants, "text", encoding = "UTF-8")
  search_variants <- jsonlite::fromJSON(search_variants)
  return(search_variants)
}

#' @title Get All Traits API
#' @description Get all traits from the API
#' @return A list containing all traits
traits_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/traits")
  traits <- httr::GET(url, httr::timeout(timeout_seconds))
  traits <- httr::content(traits, "text", encoding = "UTF-8")
  traits <- jsonlite::fromJSON(traits)
  return(traits)
}

#' @title Get a Trait API
#' @description Get a trait from the API
#' @param trait_id A character string specifying the trait ID
#' @param include_associations A logical value specifying whether to include associations
#' @return A list containing the trait
trait_api <- function(trait_id, include_associations = FALSE) {
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/traits/", trait_id,
    "?include_associations=", include_associations
  )
  trait <- httr::GET(url, httr::timeout(timeout_seconds))
  trait <- httr::content(trait, "text", encoding = "UTF-8")
  trait <- jsonlite::fromJSON(trait)
  return(trait)
}

trait_coloc_pairs_api <- function(trait_id, h4_threshold = 0.8) {
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/traits/", trait_id, "/coloc-pairs",
    "?h4_threshold=", h4_threshold
  )
  trait_coloc_pairs <- httr::GET(url, httr::timeout(timeout_seconds))
  trait_coloc_pairs <- httr::content(trait_coloc_pairs, "text", encoding = "UTF-8")
  trait_coloc_pairs <- jsonlite::fromJSON(trait_coloc_pairs)
  to_dataframe <- as.data.frame(trait_coloc_pairs$coloc_pair_rows)
  names(to_dataframe) <- trait_coloc_pairs$coloc_pair_column_names
  return(to_dataframe)
}

#' @title Get Gene Information
#' @description Get gene information from the API
#' @return A list containing the gene information
genes_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/genes")
  genes <- httr::GET(url, httr::timeout(timeout_seconds))
  genes <- httr::content(genes, "text", encoding = "UTF-8")
  genes <- jsonlite::fromJSON(genes)
  return(genes)
}

#' @title Get a Gene API
#' @description Get a gene from the API
#' @param gene_id A character string specifying the gene ID
#' @param include_associations A logical value specifying whether to include associations
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs
#' @param include_trans A logical value specifying whether to include trans genetic effects
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs
#' @return A list containing the gene information
gene_api <- function(gene_id,
                     include_associations = FALSE,
                     include_coloc_pairs = FALSE,
                     include_trans = TRUE,
                     h4_threshold = 0.8) {
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/genes/", gene_id,
    "?include_associations=", include_associations,
    "&include_coloc_pairs=", include_coloc_pairs,
    "&include_trans=", include_trans,
    "&h4_threshold=", h4_threshold
  )

  gene <- httr::GET(url, httr::timeout(timeout_seconds))
  gene <- httr::content(gene, "text", encoding = "UTF-8")
  gene <- jsonlite::fromJSON(gene)
  return(gene)
}

#' @title Get a Region API
#' @description Get a region from the API
#' @param region_id A character string specifying the region ID
#' @param include_associations A logical value specifying whether to include associations
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs
#' @return A list containing the region information
region_api <- function(region_id, include_associations = FALSE, include_coloc_pairs = FALSE, h4_threshold = 0.8) {
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/regions/", region_id,
    "?include_associations=", include_associations,
    "&include_coloc_pairs=", include_coloc_pairs,
    "&h4_threshold=", h4_threshold
  )
  region <- httr::GET(url, httr::timeout(timeout_seconds))
  region <- httr::content(region, "text", encoding = "UTF-8")
  region <- jsonlite::fromJSON(region)
  return(region)
}

#' @title Get Variants by RSID
#' @description Get variants from the API by RSID
#' @param rsids A character string specifying the RSID
#' @return A list containing the variants
variants_by_rsid_api <- function(rsids) {
  rsids <- paste(rsids, collapse = "&rsids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?rsids=", rsids)
  variants <- httr::GET(url, httr::timeout(timeout_seconds))
  variants <- httr::content(variants, "text", encoding = "UTF-8")
  variants <- jsonlite::fromJSON(variants)
  return(variants)
}

#' @title Get Variants by SNP ID
#' @description Get variants from the API by SNP ID
#' @param snp_ids A character string specifying the SNP ID
#' @param include_associations A logical value specifying whether to include associations
#' @param p_value_threshold A numeric value specifying the p-value threshold
#' @return A list containing the variants
variants_by_snp_id_api <- function(snp_ids, include_associations = FALSE, p_value_threshold = NULL) {
  snp_ids <- paste(snp_ids, collapse = "&snp_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?snp_ids=", snp_ids)
  return(get_variants_with_options_api(url, include_associations, p_value_threshold))
}

#' @title Get Variants by Variant ID
#' @description Get variants from the API by Variant ID
#' @param variants A character string specifying the Variant ID
#' @param include_associations A logical value specifying whether to include associations
#' @param p_value_threshold A numeric value specifying the p-value threshold
#' @return A list containing the variants
variants_by_variant_api <- function(variants, include_associations = FALSE, p_value_threshold = NULL) {
  variants <- paste(variants, collapse = "&variants=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?variants=", variants)
  return(get_variants_with_options_api(url, include_associations, p_value_threshold))
}

#' @title Get Variants by GRange
#' @description Get variants from the API by GRange
#' @param chr A character string specifying the chromosome
#' @param start A numeric value specifying the start position
#' @param stop A numeric value specifying the stop position
#' @param include_associations A logical value specifying whether to include associations
#' @param p_value_threshold A numeric value specifying the p-value threshold
#' @return A list containing the variants
variants_by_grange_api <- function(chr, start, stop, include_associations = FALSE, p_value_threshold = NULL) {
  url <- paste0(getOption("gpmap_url"), "/v1/variants?grange=", chr, ":", start, "-", stop)
  return(get_variants_with_options_api(url, include_associations, p_value_threshold))
}

get_variants_with_options_api <- function(url, include_associations = FALSE, p_value_threshold = NULL) {
  url <- paste0(url, "&include_associations=", include_associations)
  url <- paste0(url, "&p_value_threshold=", p_value_threshold)

  variants <- httr::GET(url, httr::timeout(timeout_seconds))
  variants <- httr::content(variants, "text", encoding = "UTF-8")
  variants <- jsonlite::fromJSON(variants)
  return(variants)
}

#' @title Get a Variant API
#' @description Get a variant from the API
#' @param snp_id A character string specifying the SNP ID
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs
#' @return A list containing the variant
variant_api <- function(snp_id, include_coloc_pairs = FALSE, h4_threshold = 0.8) {
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/variants/", snp_id,
    "?include_coloc_pairs=", include_coloc_pairs,
    "&h4_threshold=", h4_threshold
  )

  variant <- httr::GET(url, httr::timeout(timeout_seconds))
  variant <- httr::content(variant, "text", encoding = "UTF-8")
  variant <- jsonlite::fromJSON(variant)
  return(variant)
}

#' @title Get Variant Summary Stats API
#' @description Get variant summary statistics from the API
#' @param snp_id A character string specifying the SNP ID
#' @return A list containing dataframes from the TSV files
variant_summary_stats_api <- function(snp_id) {
  temp_zip <- file.path(tempdir(), paste0("summary_stats_", snp_id, ".zip"))

  url <- paste0(getOption("gpmap_url"), "/v1/variants/", snp_id, "/summary-stats")
  response <- httr::GET(url, httr::timeout(timeout_seconds))
  writeBin(httr::content(response, "raw"), temp_zip)
  utils::unzip(temp_zip, exdir = tempdir())

  tsv_files <- list.files(tempdir(), pattern = "\\.tsv\\.gz$", full.names = TRUE)
  dataframes <- list()
  for (file in tsv_files) {
    filename <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(file)))
    filename <- sub("_with_lbfs", "", filename)
    dataframes[[filename]] <- readr::read_tsv(file, show_col_types = FALSE)
  }
  unlink(temp_zip)
  unlink(tsv_files)
  return(dataframes)
}

#' @title Get LD Proxies by Variant ID API
#' @description Get LD proxies from the API by Variant ID
#' @param variants A character string specifying the Variant ID
#' @return A list containing the LD proxies
ld_proxies_by_variant_api <- function(variants) {
  variants <- paste(variants, collapse = "&variants=")
  url <- paste0(getOption("gpmap_url"), "/v1/ld/proxies?variants=", variants)
  ld_proxies <- httr::GET(url, httr::timeout(timeout_seconds))
  ld_proxies <- httr::content(ld_proxies, "text", encoding = "UTF-8")
  ld_proxies <- jsonlite::fromJSON(ld_proxies)
  return(ld_proxies)
}

#' @title Get LD Proxies by SNP ID API
#' @description Get LD proxies from the API by SNP ID
#' @param snp_ids A character string specifying the SNP ID
#' @return A list containing the LD proxies
ld_proxies_by_snp_id_api <- function(snp_ids) {
  snp_ids <- paste(snp_ids, collapse = "&snp_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/ld/proxies?snp_ids=", snp_ids)
  ld_proxies <- httr::GET(url)
  ld_proxies <- httr::content(ld_proxies, "text", encoding = "UTF-8")
  ld_proxies <- jsonlite::fromJSON(ld_proxies)
  return(ld_proxies)
}

#' @title Get LD Matrix by Variant ID API
#' @description Get LD matrix from the API by Variant ID
#' @param variants A character string specifying the Variant ID
#' @return A list containing the LD matrix
ld_matrix_by_variant_api <- function(variants) {
  variants <- paste(variants, collapse = "&variants=")
  url <- paste0(getOption("gpmap_url"), "/v1/ld/matrix?variants=", variants)
  ld_matrix <- httr::GET(url)
  ld_matrix <- httr::content(ld_matrix, "text", encoding = "UTF-8")
  ld_matrix <- jsonlite::fromJSON(ld_matrix)
  return(ld_matrix)
}

#' @title Get LD Matrix by SNP ID API
#' @description Get LD matrix from the API by SNP ID
#' @param snp_ids A character string specifying the SNP ID
#' @return A list containing the LD matrix
ld_matrix_by_snp_id_api <- function(snp_ids) {
  snp_ids <- paste(snp_ids, collapse = "&snp_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/ld/matrix?snp_ids=", snp_ids)
  ld_matrix <- httr::GET(url)
  ld_matrix <- httr::content(ld_matrix, "text", encoding = "UTF-8")
  ld_matrix <- jsonlite::fromJSON(ld_matrix)
  return(ld_matrix)
}

#' @title Get Associations by SNP ID and Study ID API
#' @description Get associations from the API by SNP id and study id
#' @param snp_ids A vector of numeric values specifying the SNP IDs
#' @param study_ids A vector of numeric values specifying the Study IDs
#' @return A list containing the associations
associations_api <- function(snp_ids, study_ids) {
  snp_ids <- paste(snp_ids, collapse = "&snp_ids=")
  study_ids <- paste(study_ids, collapse = "&study_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/associations?snp_ids=", snp_ids, "&study_ids=", study_ids)

  associations <- httr::GET(url, httr::timeout(timeout_seconds))
  associations <- httr::content(associations, "text", encoding = "UTF-8")
  associations <- jsonlite::fromJSON(associations)
  return(associations)
}

#' @title Get All Gene Pleiotropies API
#' @description Get all gene pleiotropies from the API
#' @return A list containing the gene pleiotropies
gene_pleiotropies_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/pleiotropy/genes")
  gene_pleiotropies <- httr::GET(url, httr::timeout(timeout_seconds))
  gene_pleiotropies <- httr::content(gene_pleiotropies, "text", encoding = "UTF-8")
  gene_pleiotropies <- jsonlite::fromJSON(gene_pleiotropies)
  return(gene_pleiotropies)
}

#' @title Get All SNP Pleiotropies API
#' @description Get all SNP pleiotropies from the API
#' @return A list containing the SNP pleiotropies
snp_pleiotropies_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/pleiotropy/snps")
  snp_pleiotropies <- httr::GET(url, httr::timeout(timeout_seconds))
  snp_pleiotropies <- httr::content(snp_pleiotropies, "text", encoding = "UTF-8")
  snp_pleiotropies <- jsonlite::fromJSON(snp_pleiotropies)
  return(snp_pleiotropies)
}

#' @title Get a GWAS from the API
#' @description Get a GWAS from the API
#' @param gwas_id The ID of the GWAS
#' @param include_summary_stats Whether to include summary statistics
#' @return A list containing the GWAS information
#' @export
get_gwas_api <- function(gwas_id, include_summary_stats = FALSE) {
  url <- paste0(getOption("gpmap_url"), "/v1/gwas/", gwas_id)
  gwas <- httr::GET(url, httr::timeout(timeout_seconds))
  gwas <- httr::content(gwas, "text", encoding = "UTF-8")
  gwas <- jsonlite::fromJSON(gwas)
  if (include_summary_stats) {
    summary_stats_df <- get_gwas_summary_stats_api(gwas_id)
    gwas$summary_stats <- summary_stats_df
  }
  return(gwas)
}

get_gwas_summary_stats_api <- function(gwas_id) {
  url <- paste0(getOption("gpmap_url"), "/v1/gwas/", gwas_id, "/summary-stats")
  summary_stats <- httr::GET(url, httr::timeout(timeout_seconds))
  summary_stats <- httr::content(summary_stats, "text", encoding = "UTF-8")
  summary_stats_url <- jsonlite::fromJSON(summary_stats)

  temp_file <- file.path(tempdir(), paste0("gwas_summary_stats_", gwas_id, ".tsv.gz"))
  response <- httr::GET(summary_stats_url, httr::timeout(timeout_seconds))
  writeBin(httr::content(response, "raw"), temp_file)

  summary_stats_df <- readr::read_tsv(temp_file, show_col_types = FALSE)
  unlink(temp_file)

  return(summary_stats_df)
}

#' @title Upload a GWAS to the API
#' @description Upload a GWAS to the API
#' @param file The path to the GWAS file, maximum size is 1GB
#' @param name The name of the GWAS
#' @param p_value_threshold The p-value threshold for the GWAS
#' @param column_names A list of column names in the format of: list(CHR = "chr", BP = "pos"...)
#' @param email The email of the user
#' @param category The category of the GWAS.  Only "continuous" and "categorical" are accepted.
#' @param is_published Whether the GWAS is published
#' @param doi The DOI of the GWAS
#' @param should_be_added Whether the GWAS should be added to the API
#' @param ancestry The ancestry of the GWAS.  Currently only "EUR" is accepted.
#' @param sample_size The sample size of the GWAS
#' @param reference_build The reference build of the GWAS.  Only "GRCh37" and "GRCh38" are accepted.
#' @return A list containing the GWAS information
upload_gwas_api <- function(file,
                            name,
                            p_value_threshold,
                            column_names,
                            email,
                            category,
                            is_published = FALSE,
                            doi = NA,
                            should_be_added = FALSE,
                            ancestry = "EUR",
                            sample_size,
                            reference_build = "GRCh38") {
  url <- paste0(getOption("gpmap_url"), "/v1/gwas")

  gwas_request <- list(
    reference_build = reference_build,
    email = email,
    name = name,
    category = category,
    is_published = is_published,
    doi = doi,
    should_be_added = should_be_added,
    ancestry = ancestry,
    sample_size = sample_size,
    p_value_threshold = p_value_threshold,
    column_names = column_names
  )
  request_json <- jsonlite::toJSON(gwas_request, auto_unbox = TRUE)

  gwas <- httr::POST(url, body = list(
    file = httr::upload_file(file),
    request = request_json
  ), httr::timeout(timeout_seconds))

  gwas <- httr::content(gwas, "text", encoding = "UTF-8")
  gwas <- jsonlite::fromJSON(gwas)
  return(gwas)
}