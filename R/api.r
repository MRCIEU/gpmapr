timeout_seconds <- 360
max_retry_attempts <- 5
retry_base_delay_seconds <- 1
retry_max_delay_seconds <- 32

#' @title GET request to the GPMap API
#' @description Wrapper around `httr::GET()` used by every API call. It checks
#'   once per session that this gpmapr version supports the API (see
#'   `ensure_api_compatible()`), retries with exponential backoff when the API
#'   responds with 429 (Too Many Requests), and turns an unreachable API or an
#'   HTTP error into an informative error.
#' @param ... Passed on to `httr::GET()`
#' @param check_version If `FALSE`, skip the API version check (used by the
#'   version and health endpoints themselves).
#' @return An `httr` response object
#' @keywords internal
#' @noRd
api_get <- function(..., check_version = TRUE) {
  return(api_request(httr::GET, ..., check_version = check_version))
}

#' @title POST request to the GPMap API
#' @description As `api_get()`, for `httr::POST()`.
#' @param ... Passed on to `httr::POST()`
#' @param check_version If `FALSE`, skip the API version check.
#' @return An `httr` response object
#' @keywords internal
#' @noRd
api_post <- function(..., check_version = TRUE) {
  return(api_request(httr::POST, ..., check_version = check_version))
}

api_request <- function(verb, ..., check_version = TRUE) {
  if (check_version) {
    ensure_api_compatible()
  }
  response <- request_with_backoff(verb, ...)
  stop_for_api_status(response)
  return(response)
}

# Session cache of the API URLs whose version this gpmapr version supports.
.api_state <- new.env(parent = emptyenv())

#' @title Check that this gpmapr version supports the selected API
#' @description Fetches the API version once per session per API URL and
#'   stops, with upgrade instructions, when this gpmapr version is not listed
#'   for it in `api_to_package_version`. When the API cannot be reached nothing
#'   is cached and the request itself reports the problem.
#' @return `TRUE` (invisibly) when compatible, `FALSE` when the version could
#'   not be fetched.
#' @keywords internal
#' @noRd
ensure_api_compatible <- function() {
  url <- getOption("gpmap_url")
  if (isTRUE(.api_state$compatible[[url]])) {
    return(invisible(TRUE))
  }
  api_version <- tryCatch(fetch_api_version(url), error = function(e) NULL)
  if (is.null(api_version)) {
    return(invisible(FALSE))
  }
  package_version <- installed_package_version()
  if (!package_version %in% api_to_package_version[[api_version]]) {
    stop(
      "gpmapr ", package_version, " does not support version ", api_version,
      " of the GPMap API at ", url, ". Please update gpmapr with ",
      "install.packages(\"gpmapr\") (or pak::pak(\"MRCIEU/gpmapr\") for the ",
      "development version).",
      call. = FALSE
    )
  }
  .api_state$compatible[[url]] <- TRUE
  return(invisible(TRUE))
}

installed_package_version <- function() {
  return(as.character(utils::packageVersion("gpmapr")))
}

fetch_api_version <- function(url) {
  response <- httr::GET(paste0(url, "/v1/info/version"), httr::timeout(10))
  if (httr::status_code(response) >= 400) {
    return(NULL)
  }
  version <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  return(version$version)
}

#' @title Stop with an informative error for an HTTP error response
#' @param response An `httr` response object
#' @return `response`, invisibly, when its status is below 400
#' @keywords internal
#' @noRd
stop_for_api_status <- function(response) {
  status <- httr::status_code(response)
  if (status < 400) {
    return(invisible(response))
  }
  body <- tryCatch(
    jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8")),
    error = function(e) NULL
  )
  detail <- if (is.list(body) && !is.null(body$detail)) {
    paste(unlist(body$detail), collapse = "; ")
  } else {
    httr::http_status(response)$message
  }
  stop(structure(
    class = c("gpmap_api_error", "error", "condition"),
    list(message = paste0("GPMap API request failed (HTTP ", status, "): ", detail), call = NULL, status = status)
  ))
}

#' @title Retry a request with exponential backoff when rate-limited
#' @description Performs a request via `verb` and, if the API responds with
#'   429, waits (honoring a `Retry-After` header when the API sends one,
#'   otherwise using exponential backoff with jitter) and retries, up to
#'   `max_retries` times. Any other response, including other error
#'   statuses, is returned immediately for the caller to handle as before.
#' @param verb `httr::GET` or `httr::POST`
#' @param ... Passed on to `verb`
#' @param max_retries Maximum number of retries after the initial request
#' @param base_delay Base delay in seconds used for exponential backoff
#' @param max_delay Maximum delay in seconds between retries
#' @param attempt The zero-based retry attempt number, used internally to
#'   track recursion; callers should not need to set this.
#' @return An `httr` response object
#' @keywords internal
#' @noRd
request_with_backoff <- function(verb, ...,
                                 max_retries = max_retry_attempts,
                                 base_delay = retry_base_delay_seconds,
                                 max_delay = retry_max_delay_seconds,
                                 attempt = 0) {
  response <- tryCatch(verb(...), error = function(e) {
    stop(
      "Could not reach the GPMap API at ", getOption("gpmap_url"), ": ",
      conditionMessage(e),
      call. = FALSE
    )
  })
  if (httr::status_code(response) != 429L) {
    return(response)
  }
  if (attempt >= max_retries) {
    stop(
      "GPMap API rate limit exceeded after ", max_retries,
      " retries. Please wait a moment and try again.",
      call. = FALSE
    )
  }
  Sys.sleep(retry_delay_seconds(attempt, response, base_delay, max_delay))
  return(request_with_backoff(
    verb, ...,
    max_retries = max_retries, base_delay = base_delay, max_delay = max_delay,
    attempt = attempt + 1
  ))
}

#' @title Compute the delay before the next retry
#' @description Uses the `Retry-After` header when the API provides one,
#'   otherwise falls back to exponential backoff with full jitter.
#' @param attempt The zero-based retry attempt number
#' @param response The `httr` response that returned status 429
#' @param base_delay Base delay in seconds used for exponential backoff
#' @param max_delay Maximum delay in seconds between retries
#' @return A delay in seconds
#' @keywords internal
#' @noRd
retry_delay_seconds <- function(attempt, response,
                                base_delay = retry_base_delay_seconds,
                                max_delay = retry_max_delay_seconds) {
  retry_after <- httr::headers(response)[["retry-after"]]
  if (!is.null(retry_after)) {
    retry_after_seconds <- suppressWarnings(as.numeric(retry_after))
    if (!is.na(retry_after_seconds)) {
      return(max(retry_after_seconds, 0))
    }
  }
  capped_delay <- min(max_delay, base_delay * 2^attempt)
  return(stats::runif(1, min = 0, max = capped_delay))
}

#' @title Is the GPMap API available?
#' @description Checks whether the GPMap API can be reached and is a version this
#'   version of gpmapr supports, without raising an error. Useful before running
#'   code that needs the API, for example in examples, vignettes or tests.
#' @return `TRUE` if the API's health endpoint responds and its version is
#'   supported, otherwise `FALSE`.
#' @examples
#' gpmap_available()
#' @export
gpmap_available <- function() {
  response <- tryCatch(
    httr::GET(paste0(getOption("gpmap_url"), "/health"), httr::timeout(10)),
    error = function(e) NULL
  )
  if (is.null(response) || httr::status_code(response) >= 400) {
    return(FALSE)
  }
  return(isTRUE(tryCatch(ensure_api_compatible(), error = function(e) FALSE)))
}

#' @title Select API
#' @description Select the API to use
#' @param api A character string specifying the API to use. Default is "production".
#' @return NULL
#' @noRd
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
#' @examplesIf gpmap_available()
#' health_api()$status
#' @export
health_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/health")
  health <- api_get(url, check_version = FALSE)
  health <- httr::content(health, "text", encoding = "UTF-8")
  health <- jsonlite::fromJSON(health)
  return(health)
}


#' @title Get API Version
#' @description Get the version of the API
#' @return A list containing the version
#' @noRd
version_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/info/version")
  version <- api_get(url, check_version = FALSE)
  version <- httr::content(version, "text", encoding = "UTF-8")
  version <- jsonlite::fromJSON(version)
  return(version$version)
}

#' @title Search Options API
#' @description Search options from the API
#' @return A list containing the search options
#' @noRd
search_options_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/search/options")
  search_options <- api_get(url, httr::timeout(timeout_seconds))
  search_options <- httr::content(search_options, "text", encoding = "UTF-8")
  search_options <- jsonlite::fromJSON(search_options)
  return(search_options)
}

#' @title Search Variants API
#' @description Search variants from the API
#' @param query A character string specifying the query
#' @param rsquared_threshold A numeric value specifying the rsquared threshold for the proxy variants
#' @return A list containing the variants
#' @noRd
search_variants_api <- function(query, rsquared_threshold = 0.8) {
  url <- paste0(getOption("gpmap_url"), "/v1/search/variant/", query, "?rsquared_threshold=", rsquared_threshold)
  search_variants <- api_get(url, httr::timeout(timeout_seconds))
  search_variants <- httr::content(search_variants, "text", encoding = "UTF-8")
  search_variants <- jsonlite::fromJSON(search_variants)
  return(search_variants)
}

#' @title Get All Traits API
#' @description Get all traits from the API
#' @return A list containing all traits
#' @noRd
traits_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/traits")
  traits <- api_get(url, httr::timeout(timeout_seconds))
  traits <- httr::content(traits, "text", encoding = "UTF-8")
  traits <- jsonlite::fromJSON(traits)
  return(traits)
}

#' @title Get Trait Duplicates API
#' @description Get all traits marked as duplicates, along with the parent trait they duplicate
#' @return A dataframe containing the trait duplicates
#' @noRd
trait_duplicates_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/traits/duplicates")
  duplicates <- api_get(url, httr::timeout(timeout_seconds))
  duplicates <- httr::content(duplicates, "text", encoding = "UTF-8")
  duplicates <- jsonlite::fromJSON(duplicates)
  return(duplicates$duplicates)
}

#' @title Get Specific Traits API
#' @description Get specific traits from the API
#' @param trait_ids A vector of trait ids
#' @param include_associations A logical value specifying whether to include associations
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs
#' @return A list containing the traits
#' @noRd
specific_traits_api <- function(
  trait_ids,
  include_associations = FALSE,
  include_coloc_pairs = FALSE,
  h4_threshold = 0.8
) {
  trait_ids <- paste(trait_ids, collapse = "&ids=")
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/traits?ids=", trait_ids,
    "&include_associations=", include_associations
  )
  traits <- api_get(url, httr::timeout(timeout_seconds))
  traits <- httr::content(traits, "text", encoding = "UTF-8")
  traits <- jsonlite::fromJSON(traits)
  return(traits)
}

#' @title Get a Trait API
#' @description Get a trait from the API
#' @param trait_id A character string specifying the trait ID
#' @param include_associations A logical value specifying whether to include associations
#' @return A list containing the trait
#' @noRd
trait_api <- function(trait_id, include_associations = FALSE) {
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/traits/", trait_id,
    "?include_associations=", include_associations
  )
  trait <- api_get(url, httr::timeout(timeout_seconds))
  trait <- httr::content(trait, "text", encoding = "UTF-8")
  trait <- jsonlite::fromJSON(trait)
  return(trait)
}

#' @title Get Trait Coloc Pairs API
#' @description Get trait coloc pairs from the API
#' @param trait_id A character string specifying the trait ID
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs
#' @return A list containing the trait coloc pairs
#' @noRd
trait_coloc_pairs_api <- function(trait_id, h4_threshold = 0.8) {
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/traits/", trait_id, "/coloc-pairs",
    "?h4_threshold=", h4_threshold
  )
  trait_coloc_pairs <- api_get(url, httr::timeout(timeout_seconds))
  trait_coloc_pairs <- httr::content(trait_coloc_pairs, "text", encoding = "UTF-8")
  trait_coloc_pairs <- jsonlite::fromJSON(trait_coloc_pairs)
  to_dataframe <- as.data.frame(trait_coloc_pairs$coloc_pair_rows)
  names(to_dataframe) <- trait_coloc_pairs$coloc_pair_column_names
  return(to_dataframe)
}

#' @title Get Trait Full Associations API
#' @description Get full trait associations from the API
#' @param trait_id A character string specifying the trait ID
#' @return A dataframe containing full trait associations
#' @noRd
trait_associations_full_api <- function(trait_id) {
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/traits/", trait_id, "/associations-full"
  )
  response <- api_get(url, httr::timeout(timeout_seconds))
  response <- httr::content(response, "text", encoding = "UTF-8")
  response <- jsonlite::fromJSON(response)
  to_dataframe <- as.data.frame(response$associations_full_rows)
  names(to_dataframe) <- response$associations_full_column_names
  return(to_dataframe)
}

#' @title Get Gene Information
#' @description Get gene information from the API
#' @return A list containing the genes
#' @noRd
genes_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/genes")
  genes <- api_get(url, httr::timeout(timeout_seconds))
  genes <- httr::content(genes, "text", encoding = "UTF-8")
  genes <- jsonlite::fromJSON(genes)
  return(genes)
}

#' @title Get Specific Genes API
#' @description Get specific genes from the API
#' @param gene_ids A vector of gene ids
#' @param include_associations A logical value specifying whether to include associations
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs
#' @param include_trans A logical value specifying whether to include trans genetic effects
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs
#' @return A list containing the genes
#' @noRd
specific_genes_api <- function(
  gene_ids,
  include_associations = FALSE,
  include_coloc_pairs = FALSE,
  include_trans = TRUE,
  h4_threshold = 0.8
) {
  gene_ids <- paste(gene_ids, collapse = "&ids=")
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/genes?ids=", gene_ids, "&include_associations=",
    include_associations, "&include_coloc_pairs=",
    include_coloc_pairs, "&include_trans=", include_trans,
    "&h4_threshold=", h4_threshold
  )
  genes <- api_get(url, httr::timeout(timeout_seconds))
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
#' @noRd
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

  gene <- api_get(url, httr::timeout(timeout_seconds))
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
#' @noRd
region_api <- function(region_id, include_associations = FALSE, include_coloc_pairs = FALSE, h4_threshold = 0.8) {
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/regions/", region_id,
    "?include_associations=", include_associations,
    "&include_coloc_pairs=", include_coloc_pairs,
    "&h4_threshold=", h4_threshold
  )
  region <- api_get(url, httr::timeout(timeout_seconds))
  region <- httr::content(region, "text", encoding = "UTF-8")
  region <- jsonlite::fromJSON(region)
  return(region)
}

#' @title Get Variants API
#' @description Get variants from the API. The API accepts variants (variant_ids, rsids, or strings - auto-detected)
#' or a genomic range (grange).
#' @param variants A vector of variant identifiers (variant_ids, rsids, or strings)
#' @param expand Logical. FALSE (default) returns minimal data. TRUE returns full VariantResponse (max 10, not available with grange)
#' @param include_associations Logical. Whether to include associations for SNPs (only when expand=TRUE)
#' @param include_coloc_pairs Logical. Whether to include coloc pairs for SNPs (only when expand=TRUE)
#' @param h4_threshold Numeric. H4 threshold for coloc pairs, defaults to 0.8
#' @return A list containing the variants
#' @noRd
variants_api <- function(variants, expand = FALSE, include_associations = FALSE, include_coloc_pairs = FALSE, h4_threshold = 0.8) {
  variants <- paste(variants, collapse = "&variants=")
  url <- paste0(getOption("gpmap_url"), "/v1/variants?variants=", variants)
  return(get_variants_with_options_api(url, expand, include_associations, include_coloc_pairs, h4_threshold))
}

#' @title Get Variants by GRange API
#' @description Get variants from the API by genomic range. Expand is not available with grange.
#' @param chr A character string specifying the chromosome
#' @param start A numeric value specifying the start position
#' @param stop A numeric value specifying the stop position
#' @param include_associations Logical. Whether to include associations (only when expand=TRUE, but expand not available with grange)
#' @param include_coloc_pairs Logical. Whether to include coloc pairs
#' @param h4_threshold Numeric. H4 threshold for coloc pairs, defaults to 0.8
#' @return A list containing the variants
#' @noRd
variants_grange_api <- function(chr, start, stop, include_associations = FALSE, include_coloc_pairs = FALSE, h4_threshold = 0.8) {
  url <- paste0(getOption("gpmap_url"), "/v1/variants?grange=", chr, ":", start, "-", stop)
  return(get_variants_with_options_api(url, expand = FALSE, include_associations, include_coloc_pairs, h4_threshold))
}

#' @title Get Variants with Options API
#' @description Get variants from the API with options
#' @param url A character string specifying the URL
#' @param expand A logical value specifying whether to expand the variants
#' @param include_associations A logical value specifying whether to include associations
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs
#' @return A list containing the variants
#' @noRd
get_variants_with_options_api <- function(
  url,
  expand = FALSE,
  include_associations = FALSE,
  include_coloc_pairs = FALSE,
  h4_threshold = 0.8
) {
  sep <- if (grepl("?", url, fixed = TRUE)) "&" else "?"
  url <- paste0(url, sep, "expand=", tolower(expand))
  url <- paste0(url, "&include_associations=", tolower(include_associations))
  url <- paste0(url, "&include_coloc_pairs=", tolower(include_coloc_pairs))
  url <- paste0(url, "&h4_threshold=", h4_threshold)

  variants <- api_get(url, httr::timeout(timeout_seconds))
  variants <- httr::content(variants, "text", encoding = "UTF-8")
  variants <- jsonlite::fromJSON(variants)
  return(variants)
}

#' @title Get a Variant API
#' @description Get a variant from the API
#' @param variant_id A character string specifying the SNP ID
#' @param include_coloc_pairs A logical value specifying whether to include coloc pairs
#' @param h4_threshold A numeric value specifying the h4 threshold for coloc pairs
#' @return A list containing the variant
#' @noRd
variant_api <- function(variant_id, include_coloc_pairs = FALSE, h4_threshold = 0.8) {
  url <- paste0(
    getOption("gpmap_url"),
    "/v1/variants/", variant_id,
    "?include_coloc_pairs=", include_coloc_pairs,
    "&h4_threshold=", h4_threshold
  )

  variant <- api_get(url, httr::timeout(timeout_seconds))
  variant <- httr::content(variant, "text", encoding = "UTF-8")
  variant <- jsonlite::fromJSON(variant)
  return(variant)
}

#' @title Get Variant Summary Stats API
#' @description Get variant summary statistics from the API
#' @param variant_id A character string specifying the SNP ID
#' @return A list containing dataframes from the TSV files
#' @noRd
variant_summary_stats_api <- function(variant_id) {
  temp_zip <- file.path(tempdir(), paste0("summary_stats_", variant_id, ".zip"))

  url <- paste0(getOption("gpmap_url"), "/v1/variants/", variant_id, "/summary-stats")
  response <- api_get(url, httr::timeout(timeout_seconds))
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
#' @noRd
ld_proxies_by_variant_api <- function(variants) {
  variants <- paste(variants, collapse = "&variants=")
  url <- paste0(getOption("gpmap_url"), "/v1/ld/proxies?variants=", variants)
  ld_proxies <- api_get(url, httr::timeout(timeout_seconds))
  ld_proxies <- httr::content(ld_proxies, "text", encoding = "UTF-8")
  ld_proxies <- jsonlite::fromJSON(ld_proxies)
  return(ld_proxies)
}

#' @title Get LD Proxies by SNP ID API
#' @description Get LD proxies from the API by SNP ID
#' @param variant_ids A character string specifying the SNP ID
#' @return A list containing the LD proxies
#' @noRd
ld_proxies_by_variant_id_api <- function(variant_ids) {
  variant_ids <- paste(variant_ids, collapse = "&variant_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/ld/proxies?variant_ids=", variant_ids)
  ld_proxies <- api_get(url)
  ld_proxies <- httr::content(ld_proxies, "text", encoding = "UTF-8")
  ld_proxies <- jsonlite::fromJSON(ld_proxies)
  return(ld_proxies)
}

#' @title Get LD Matrix by Variant ID API
#' @description Get LD matrix from the API by Variant ID
#' @param variants A character string specifying the Variant ID
#' @return A list containing the LD matrix
#' @noRd
ld_matrix_by_variant_api <- function(variants) {
  variants <- paste(variants, collapse = "&variants=")
  url <- paste0(getOption("gpmap_url"), "/v1/ld/matrix?variants=", variants)
  return(get_ld_matrix(url))
}

#' @title Get LD Matrix by SNP ID API
#' @description Get LD matrix from the API by SNP ID
#' @param variant_ids A character string specifying the SNP ID
#' @return A list containing the LD matrix
#' @noRd
ld_matrix_by_variant_id_api <- function(variant_ids) {
  variant_ids <- paste(variant_ids, collapse = "&variant_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/ld/matrix?variant_ids=", variant_ids)
  return(get_ld_matrix(url))
}

# An LD matrix response; `lds` is empty when none of the variants are stored as a lead and its proxy.
get_ld_matrix <- function(url) {
  ld_matrix <- api_get(url)
  ld_matrix <- httr::content(ld_matrix, "text", encoding = "UTF-8")
  ld_matrix <- jsonlite::fromJSON(ld_matrix)
  return(ld_matrix)
}

#' @title Get Associations by SNP ID and Study ID API
#' @description Get associations from the API by SNP id and study id
#' @param variant_ids A vector of numeric values specifying the SNP IDs
#' @param study_ids A vector of numeric values specifying the Study IDs
#' @return A list containing the associations
#' @noRd
associations_api <- function(variant_ids, study_ids) {
  variant_ids <- paste(variant_ids, collapse = "&variant_ids=")
  study_ids <- paste(study_ids, collapse = "&study_ids=")
  url <- paste0(getOption("gpmap_url"), "/v1/associations?variant_ids=", variant_ids, "&study_ids=", study_ids)

  associations <- api_get(url, httr::timeout(timeout_seconds))
  associations <- httr::content(associations, "text", encoding = "UTF-8")
  associations <- jsonlite::fromJSON(associations)
  return(associations)
}

#' @title Get All Gene Pleiotropies API
#' @description Get all gene pleiotropies from the API
#' @return A list containing the gene pleiotropies
#' @noRd
gene_pleiotropies_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/pleiotropy/genes")
  gene_pleiotropies <- api_get(url, httr::timeout(timeout_seconds))
  gene_pleiotropies <- httr::content(gene_pleiotropies, "text", encoding = "UTF-8")
  gene_pleiotropies <- jsonlite::fromJSON(gene_pleiotropies)
  return(gene_pleiotropies)
}

#' @title Get All SNP Pleiotropies API
#' @description Get all SNP pleiotropies from the API
#' @return A list containing the SNP pleiotropies
#' @noRd
variant_pleiotropies_api <- function() {
  url <- paste0(getOption("gpmap_url"), "/v1/pleiotropy/snps")
  variant_pleiotropies <- api_get(url, httr::timeout(timeout_seconds))
  variant_pleiotropies <- httr::content(variant_pleiotropies, "text", encoding = "UTF-8")
  variant_pleiotropies <- jsonlite::fromJSON(variant_pleiotropies)
  return(variant_pleiotropies)
}

#' @title Get a GWAS from the API
#' @description Get a GWAS from the API
#' @param gwas_id The ID of the GWAS
#' @param include_summary_stats Whether to include summary statistics
#' @param include_associations Whether to include associations
#' @return A list containing the GWAS information
#' @noRd
get_gwas_api <- function(gwas_id, include_summary_stats = FALSE, include_associations = FALSE) {
  url <- paste0(getOption("gpmap_url"), "/v1/gwas/", gwas_id)
  if (include_associations) {
    url <- paste0(url, "?include_associations=true")
  }
  gwas <- api_get(url, httr::timeout(timeout_seconds))
  gwas <- httr::content(gwas, "text", encoding = "UTF-8")
  gwas <- jsonlite::fromJSON(gwas)
  if (include_summary_stats) {
    summary_stats_df <- get_gwas_summary_stats_api(gwas_id)
    gwas$summary_stats <- summary_stats_df
  }
  return(gwas)
}

#' @title Get a GWAS Summary Stats API
#' @description Get a GWAS summary stats from the API
#' @param gwas_id The ID of the GWAS
#' @return A list containing the GWAS summary stats
#' @noRd
get_gwas_summary_stats_api <- function(gwas_id) {
  url <- paste0(getOption("gpmap_url"), "/v1/gwas/", gwas_id, "/summary-stats")
  summary_stats <- api_get(url, httr::timeout(timeout_seconds))
  summary_stats <- httr::content(summary_stats, "text", encoding = "UTF-8")
  summary_stats_url <- jsonlite::fromJSON(summary_stats)

  temp_file <- file.path(tempdir(), paste0("gwas_summary_stats_", gwas_id, ".tsv.gz"))
  response <- api_get(summary_stats_url, httr::timeout(timeout_seconds))
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
#' @param column_names A validated named list of column names, e.g. list(CHR = "chr", BP = "pos"...)
#' @param email The email of the user
#' @param variable_type The type of variable the GWAS is of.  Only "continuous" and "categorical" are accepted.
#' @param is_published Whether the GWAS is published
#' @param doi The DOI of the GWAS
#' @param should_be_added Whether the GWAS should be added to the API
#' @param ancestry The ancestry of the GWAS.  Currently only "EUR" is accepted.
#' @param sample_size The sample size of the GWAS
#' @param reference_build The reference build of the GWAS.  Only "GRCh37" and "GRCh38" are accepted.
#' @param compare_with_upload_guids A vector of GUIDs of uploads to compare with
#' @return A list containing the GWAS information
#' @noRd
upload_gwas_api <- function(file,
                            name,
                            p_value_threshold,
                            column_names,
                            email,
                            variable_type,
                            is_published = FALSE,
                            doi = NA,
                            should_be_added = FALSE,
                            ancestry = "EUR",
                            sample_size,
                            reference_build = "GRCh38",
                            compare_with_upload_guids = NA) {
  url <- paste0(getOption("gpmap_url"), "/v1/gwas")

  if (!all(is.na(compare_with_upload_guids))) {
    compare_with_upload_guids <- I(as.character(compare_with_upload_guids))
  }

  gwas_request <- list(
    reference_build = reference_build,
    email = email,
    name = name,
    variable_type = variable_type,
    is_published = is_published,
    doi = doi,
    should_be_added = should_be_added,
    ancestry = ancestry,
    sample_size = sample_size,
    p_value_threshold = p_value_threshold,
    column_names = column_names,
    compare_with_upload_guids = compare_with_upload_guids
  )
  request_json <- jsonlite::toJSON(gwas_request, auto_unbox = TRUE)

  gwas <- api_post(url, body = list(
    file = httr::upload_file(file),
    request = request_json
  ), httr::timeout(timeout_seconds))

  gwas <- httr::content(gwas, "text", encoding = "UTF-8")
  gwas <- jsonlite::fromJSON(gwas)
  return(gwas)
}

#' @title Delete a GWAS API
#' @description Permanently delete a GWAS upload from the API
#' @param guid The GUID of the GWAS upload
#' @param email The email used for the upload
#' @return A list containing the API response message
#' @noRd
delete_gwas_api <- function(guid, email) {
  url <- paste0(getOption("gpmap_url"), "/v1/gwas/", guid)
  response <- httr::DELETE(url, body = list(email = email), encode = "json", httr::timeout(timeout_seconds))
  content <- httr::content(response, "text", encoding = "UTF-8")
  content <- tryCatch(jsonlite::fromJSON(content), error = function(e) NULL)

  if (httr::status_code(response) >= 400) {
    detail <- if (!is.null(content$detail)) content$detail else httr::http_status(response)$message
    stop("Failed to delete GWAS upload ", guid, ": ", detail)
  }
  return(content)
}
