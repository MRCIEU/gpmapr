#' @title Get a GWAS from the API
#' @description Get a GWAS from the API
#' @param gwas_id The ID of the GWAS
#' @param include_summary_stats Whether to include summary statistics
#' @return A list containing the GWAS information
#' @export
get_gwas <- function(gwas_id, include_summary_stats = FALSE) {
  gwas_info <- get_gwas_api(gwas_id, include_summary_stats = include_summary_stats)
  gwas_info <- cleanup_api_object(gwas_info)
  return(gwas_info)
}

#' @title Upload a GWAS to the API
#' @description Upload a GWAS to the API
#' @param file The path to the GWAS file, maximum size is 1GB
#' @param name The name of the GWAS
#' @param p_value_threshold The p-value threshold for the GWAS
#' @param column_names A list of column names in the format of: list(CHR = "chr", BP = "pos"...)
#' * CHR: chromosome
#' * BP: base pair position
#' * P: p-value
#' * EA: allele 1
#' * OA: allele 2
#' * EAF: allele frequency
#' And either BETA and SE, or OR, LB, and UB
#' * BETA: beta
#' * SE: standard error
#' * OR: odds ratio
#' * LB: lower bound of the confidence interval
#' * UB: upper bound of the confidence interval
#' @param email The email of the user
#' @param category The category of the GWAS.  Only "continuous" and "categorical" are accepted.
#' @param is_published Whether the GWAS is published
#' @param doi The DOI of the GWAS
#' @param should_be_added Whether the GWAS should be added to the API
#' @param ancestry The ancestry of the GWAS.  Currently only "EUR" is accepted.
#' @param sample_size The sample size of the GWAS
#' @param reference_build The reference build of the GWAS.  Only "GRCh37" and "GRCh38" are accepted.
#' @return A list containing the GWAS information
#' @export
upload_gwas <- function(file,
                        name,
                        p_value_threshold = 5e-8,
                        column_names = list(),
                        email = NA,
                        category = "continuous",
                        is_published = FALSE,
                        doi = NA,
                        should_be_added = FALSE,
                        ancestry = "EUR",
                        sample_size = NA,
                        reference_build = "GRCh38") {
  if (is.na(file) || !file.exists(file) || file.info(file)$size > 1024^3) {
    stop("file must be a valid file and less than 1GB")
  }
  if (is.na(name)) stop("name is required")
  if (is.na(email)) stop("email is required")
  if (is.na(category) || !category %in% c("continuous", "categorical")) stop("category is required")
  if (is.na(ancestry) || !ancestry %in% c("EUR")) stop("ancestry is required")
  if (is.na(reference_build) || !reference_build %in% c("GRCh37", "GRCh38")) stop("reference_build is required")
  if (is.na(p_value_threshold) || !is.numeric(p_value_threshold) || p_value_threshold >= 1e-5) {
    stop("p_value_threshold must be a number between 0 and 1e-5")
  }
  if (is.na(sample_size) || !is.numeric(sample_size) || sample_size <= 0) {
    stop("sample_size must be a positive number")
  }

  gwas_info <- upload_gwas_api(file,
    name,
    p_value_threshold,
    column_names,
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