#' @title Pathway Enrichment
#' @description Test a set of genes for pathway enrichment using the GPMap pathway database.
#' @param genes A vector of numeric gene IDs (from `all_genes()$id` or coloc group `gene_id`
#'   values) or gene names (e.g. `"APOE"`).
#' @param source Optional pathway source to filter by: `"Reactome"`, `"KEGG"`, or `"HP"`.
#' @param p_value_threshold FDR-adjusted p-value threshold for filtering results. Defaults to 0.05.
#' @param minimum_count_in_network Optional minimum number of input genes that must overlap a
#'   pathway for it to be tested. If `NULL`, the API default is used.
#' @return A list with:
#'   \itemize{
#'     \item results: a dataframe of enriched pathways with columns term_id, source,
#'       description, pathway_size, background_size, overlap, p_value, fdr, gene_ids
#'       (input genes overlapping the pathway), and pathway_gene_ids (all genes in the pathway)
#'     \item input_gene_count: number of genes submitted
#'     \item matched_gene_count: number of genes matched in the pathway database
#'     \item source: pathway source filter applied (if any)
#'     \item p_value_threshold: FDR threshold used
#'     \item minimum_count_in_network: minimum overlap threshold used (if returned by API)
#'     \item total_terms_tested: total pathway terms tested
#'   }
#' @export
pathway_enrichment <- function(genes,
                               source = NULL,
                               p_value_threshold = 0.05,
                               minimum_count_in_network = NULL) {
  if (is.null(genes) || length(genes) == 0) {
    stop("genes is required")
  }
  if (any(is.na(genes))) {
    stop("genes must not contain NA values")
  }
  if (!is.numeric(genes) && !is.character(genes)) {
    stop("genes must be numeric gene IDs or character gene names")
  }
  if (is.character(genes) && any(genes == "")) {
    stop("genes must not contain empty strings")
  }

  genes <- unique(genes)

  if (!is.null(source)) {
    valid_sources <- c("Reactome", "KEGG", "HP")
    if (!source %in% valid_sources) {
      stop("source must be one of: ", paste(valid_sources, collapse = ", "))
    }
  }

  if (!is.numeric(p_value_threshold) || p_value_threshold <= 0 || p_value_threshold > 1) {
    stop("p_value_threshold must be a number between 0 and 1")
  }

  if (!is.null(minimum_count_in_network)) {
    if (!is.numeric(minimum_count_in_network) ||
        minimum_count_in_network < 1 ||
        minimum_count_in_network != as.integer(minimum_count_in_network)) {
      stop("minimum_count_in_network must be a positive integer")
    }
  }

  return(pathway_enrichment_api(
    genes = genes,
    source = source,
    p_value_threshold = p_value_threshold,
    minimum_count_in_network = minimum_count_in_network
  ))
}
