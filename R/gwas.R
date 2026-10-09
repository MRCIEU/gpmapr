#' @title Get a GWAS from the API
#' @description Get a GWAS from the API
#' @param gwas_id The ID of the GWAS
#' @param include_summary_stats Whether to include summary statistics
#' @param include_associations Whether to include associations
#' @return A list containing the GWAS information
#' @export
get_gwas <- function(gwas_id, include_associations = FALSE, include_summary_stats = FALSE) {
  gwas_info <- get_gwas_api(
    gwas_id,
    include_summary_stats = include_summary_stats,
    include_associations = include_associations
  )
  gwas_info <- cleanup_api_object(gwas_info)
  return(gwas_info)
}

#' @title Upload a GWAS to the API
#' @description Upload a GWAS to the API.  Column names are specified with the `*_col` arguments, e.g.
#' `chr_col = "chr"`.  `chr_col`, `bp_col`, `ea_col`, `oa_col` and `p_col` are required, along with either
#' `beta_col` and `se_col`, or `or_col`, `or_lb_col` and `or_ub_col`.
#' @param file The path to the GWAS file, maximum size is 1GB
#' @param name The name of the GWAS
#' @param p_value_threshold The p-value threshold for the GWAS
#' @param column_names Deprecated, use the `*_col` arguments instead.
#' A named list of column names in the format of: list(CHR = "chr", BP = "pos"...).
#' Accepted names are SNP, RSID, CHR, BP, EA, OA, P, EAF, BETA, SE, OR, OR_LB (or LB) and OR_UB (or UB).
#' @param email The email of the user
#' @param category The category of the GWAS.  Only "continuous" and "categorical" are accepted.
#' @param is_published Whether the GWAS is published
#' @param doi The DOI of the GWAS
#' @param should_be_added Whether the GWAS should be added to the API
#' @param ancestry The ancestry of the GWAS.  Currently only "EUR" is accepted.
#' @param sample_size The sample size of the GWAS
#' @param reference_build The reference build of the GWAS.  Only "GRCh37" and "GRCh38" are accepted.
#' @param compare_with_upload_guids A vector of GUIDs of uploads to compare with
#' @param chr_col The name of the chromosome column (required)
#' @param bp_col The name of the base pair position column (required)
#' @param ea_col The name of the effect allele column (required)
#' @param oa_col The name of the other allele column (required)
#' @param p_col The name of the p-value column (required)
#' @param eaf_col The name of the effect allele frequency column (optional)
#' @param beta_col The name of the beta column (required with `se_col`, unless using odds ratios)
#' @param se_col The name of the standard error column (required with `beta_col`, unless using odds ratios)
#' @param or_col The name of the odds ratio column (required with `or_lb_col` and `or_ub_col`, unless using beta)
#' @param or_lb_col The name of the odds ratio lower confidence bound column
#' @param or_ub_col The name of the odds ratio upper confidence bound column
#' @param snp_col The name of the SNP identifier column (optional)
#' @param rsid_col The name of the rsID column (optional)
#' @return A list containing the GWAS information, including the `guid` used to fetch results
#' @export
upload_gwas <- function(file,
                        name,
                        p_value_threshold = 5e-8,
                        column_names = NULL,
                        email = NA,
                        category = "continuous",
                        is_published = FALSE,
                        doi = NA,
                        should_be_added = FALSE,
                        ancestry = "EUR",
                        sample_size = NA,
                        reference_build = "GRCh38",
                        compare_with_upload_guids = NA,
                        chr_col = NULL,
                        bp_col = NULL,
                        ea_col = NULL,
                        oa_col = NULL,
                        p_col = NULL,
                        eaf_col = NULL,
                        beta_col = NULL,
                        se_col = NULL,
                        or_col = NULL,
                        or_lb_col = NULL,
                        or_ub_col = NULL,
                        snp_col = NULL,
                        rsid_col = NULL) {

  if (is.na(file) || !file.exists(file) || file.info(file)$size > 1024^3) {
    stop("file must be a valid file and less than 1GB")
  }
  if (is.na(name)) stop("name is required")
  if (is.na(email)) stop("email is required")
  if (is.na(category) || !category %in% c("continuous", "categorical")) stop("category is required")
  if (is.na(ancestry) || !ancestry %in% c("EUR")) stop("ancestry is required")
  if (is.na(reference_build) || !reference_build %in% c("GRCh37", "GRCh38")) stop("reference_build is required")
  if (is.na(p_value_threshold) || !is.numeric(p_value_threshold) || p_value_threshold > 1e-5) {
    stop("p_value_threshold must be a number between 0 and 1e-5")
  }
  if (is.na(sample_size) || !is.numeric(sample_size) || sample_size <= 0) {
    stop("sample_size must be a positive number")
  }
  if (!all(is.na(compare_with_upload_guids))) {
    if (!all(sapply(compare_with_upload_guids, is_guid))) {
      stop("compare_with_upload_guids must be a vector of GUIDs")
    }
  }

  col_args <- list(
    chr_col = chr_col, bp_col = bp_col, ea_col = ea_col, oa_col = oa_col, p_col = p_col,
    eaf_col = eaf_col, beta_col = beta_col, se_col = se_col, or_col = or_col,
    or_lb_col = or_lb_col, or_ub_col = or_ub_col, snp_col = snp_col, rsid_col = rsid_col
  )
  column_names <- build_gwas_column_names(column_names, col_args)
  check_gwas_file_columns(file, column_names)

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
    reference_build,
    compare_with_upload_guids
  )
  return(gwas_info)
}

#' @title Delete a GWAS upload
#' @description Permanently delete a GWAS upload and all of its results and files.
#' The email must match the one used for the upload, and uploads that are still processing cannot be deleted.
#' @param gwas_id The GUID of the GWAS upload
#' @param email The email address used for the upload
#' @return Invisibly, a list containing the API response message
#' @export
delete_gwas <- function(gwas_id, email) {
  if (missing(gwas_id) || length(gwas_id) != 1 || is.na(gwas_id) || !is_guid(gwas_id)) {
    stop("gwas_id must be the GUID of a GWAS upload")
  }
  if (missing(email) || length(email) != 1 || is.na(email) || !nzchar(trimws(email))) {
    stop("email is required")
  }

  response <- delete_gwas_api(gwas_id, email)
  return(invisible(response))
}

# API column name -> upload_gwas() argument name
gwas_column_args <- c(
  SNP = "snp_col", RSID = "rsid_col", CHR = "chr_col", BP = "bp_col", EA = "ea_col", OA = "oa_col",
  P = "p_col", EAF = "eaf_col", BETA = "beta_col", SE = "se_col",
  OR = "or_col", OR_LB = "or_lb_col", OR_UB = "or_ub_col"
)

#' Build and validate the column names list sent to the API
#' @param column_names The deprecated column_names list, or NULL
#' @param col_args A named list of the *_col arguments
#' @return A named list of column names keyed by API column name (CHR, BP, ...)
#' @noRd
build_gwas_column_names <- function(column_names, col_args) {
  col_args <- col_args[!sapply(col_args, is.null)]

  if (length(column_names) > 0) {
    if (length(col_args) > 0) {
      stop("Specify columns with either the *_col arguments or column_names, not both")
    }
    warning("column_names is deprecated, please use the *_col arguments instead (e.g. chr_col = \"chr\")",
      call. = FALSE
    )
    if (is.null(names(column_names)) || any(names(column_names) == "")) {
      stop("column_names must be a named list, e.g. list(CHR = \"chr\")")
    }
    keys <- toupper(names(column_names))
    keys[keys == "LB"] <- "OR_LB"
    keys[keys == "UB"] <- "OR_UB"
    unknown <- setdiff(keys, names(gwas_column_args))
    if (length(unknown) > 0) {
      stop("Unknown column_names: ", paste(unknown, collapse = ", "),
        ". Accepted names are: ", paste(names(gwas_column_args), collapse = ", ")
      )
    }
    col_args <- as.list(column_names)
    names(col_args) <- gwas_column_args[keys]
  }

  is_column_name <- vapply(col_args, function(x) {
    is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
  }, logical(1))
  if (!all(is_column_name)) {
    stop(paste(names(col_args)[!is_column_name], collapse = ", "), " must be a single column name")
  }

  missing_args <- setdiff(c("chr_col", "bp_col", "ea_col", "oa_col", "p_col"), names(col_args))
  if (length(missing_args) > 0) {
    stop("Missing required column arguments: ", paste(missing_args, collapse = ", "))
  }
  has_beta <- all(c("beta_col", "se_col") %in% names(col_args))
  has_or <- all(c("or_col", "or_lb_col", "or_ub_col") %in% names(col_args))
  if (!has_beta && !has_or) {
    stop("Either beta_col and se_col, or or_col, or_lb_col and or_ub_col must be specified")
  }

  column_names <- col_args
  names(column_names) <- names(gwas_column_args)[match(names(col_args), gwas_column_args)]
  return(column_names)
}

#' Check that the specified columns exist in the GWAS file header
#' @param file The path to the GWAS file
#' @param column_names A named list of column names
#' @return NULL, errors if any columns are missing from the file.  If the header can't be read, the check is skipped.
#' @noRd
check_gwas_file_columns <- function(file, column_names) {
  header <- tryCatch(
    names(readr::read_delim(file, n_max = 0, show_col_types = FALSE, progress = FALSE)),
    error = function(e) NULL
  )
  if (length(header) <= 1) return(invisible(NULL))

  missing_columns <- setdiff(unlist(column_names), header)
  if (length(missing_columns) > 0) {
    stop("Columns not found in file: ", paste(missing_columns, collapse = ", "))
  }
  return(invisible(NULL))
}
