get_gwas <- function(gwas_id) {
  gwas_info <- get_gwas_api(gwas_id)
  return(gwas_info)
}

#' @title Upload a GWAS to the API
#' @description Upload a GWAS to the API
#' @param gwas_file The path to the GWAS file
#' @param gwas_name The name of the GWAS
#' @param p_value_threshold The p-value threshold for the GWAS
#' @param header_map A list of header mappings
#' @param email The email of the user
#' @param category The category of the GWAS.  Only "continuous" and "categorical" are accepted.
#' @param is_published Whether the GWAS is published
#' @param doi The DOI of the GWAS
#' @param should_be_added Whether the GWAS should be added to the API
#' @param ancestry The ancestry of the GWAS.  Currently only "EUR" is accepted.
#' @param sample_size The sample size of the GWAS
#' @param reference_build The reference build of the GWAS.  "GRCh37" and "GRCh38" are accepted.
#' @return A list containing the GWAS information
#' @export
upload_gwas <- function(gwas_file,
                        gwas_name,
                        p_value_threshold = 5e-8,
                        header_map = list(),
                        email = "",
                        category = "continuous",
                        is_published = FALSE,
                        doi = NA,
                        should_be_added = FALSE,
                        ancestry = "EUR",
                        sample_size = NA,
                        reference_build = "GRCh38") {
  if (is.na(gwas_name)) stop("gwas_name is required")
  if (is.na(email)) stop("email is required")
  if (is.na(category) || !category %in% c("continuous", "categorical")) stop("category is required")
  if (is.na(ancestry) || !ancestry %in% c("EUR")) stop("ancestry is required")
  if (is.na(reference_build) || !reference_build %in% c("GRCh37", "GRCh38")) stop("reference_build is required")
  if (is.na(p_value_threshold) || !is.numeric(p_value_threshold) || p_value_threshold > 1.5e-4) {
    stop("p_value_threshold must be a number between 0 and 1.5e-4")
  }
  if (is.na(sample_size) || !is.numeric(sample_size) || sample_size <= 0) {
    stop("sample_size must be a positive number")
  }

  gwas_info <- upload_gwas_api(gwas_file,
    gwas_name,
    p_value_threshold,
    header_map,
    email,
    category,
    is_published,
    doi,
    should_be_added,
    ancestry,
    sample_size,
    reference_build
  )
  return(gwas_info)
}

get_gwas <- function(gwas_id) {
  gwas_info <- get_gwas_api(gwas_id)
  return(gwas_info)
}