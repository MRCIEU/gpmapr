#' @title Pathway Enrichment
#' @description Test a set of genes for pathway enrichment using the GPMap pathway database.
#' @param gene_ids A vector of numeric gene IDs (from `all_genes()$id` or coloc group `gene_id` values).
#' @param source Optional pathway source to filter by: `"Reactome"`, `"KEGG"`, or `"HP"`.
#' @param p_value_threshold FDR-adjusted p-value threshold for filtering results. Defaults to 0.05.
#' @return A list with:
#'   \itemize{
#'     \item results: a dataframe of enriched pathways with columns term_id, source,
#'       description, pathway_size, background_size, overlap, p_value, fdr, gene_ids
#'     \item input_gene_count: number of gene IDs submitted
#'     \item matched_gene_count: number of gene IDs matched in the pathway database
#'     \item source: pathway source filter applied (if any)
#'     \item p_value_threshold: FDR threshold used
#'     \item total_terms_tested: total pathway terms tested
#'   }
#' @export
pathway_enrichment <- function(gene_ids,
                               source = NULL,
                               p_value_threshold = 0.05) {
  if (is.null(gene_ids) || length(gene_ids) == 0) {
    stop("gene_ids is required")
  }
  if (any(is.na(gene_ids))) {
    stop("gene_ids must not contain NA values")
  }
  if (!is.numeric(gene_ids)) {
    stop("gene_ids must be numeric gene IDs")
  }

  gene_ids <- unique(gene_ids)

  if (!is.null(source)) {
    valid_sources <- c("Reactome", "KEGG", "HP")
    if (!source %in% valid_sources) {
      stop("source must be one of: ", paste(valid_sources, collapse = ", "))
    }
  }

  if (!is.numeric(p_value_threshold) || p_value_threshold <= 0 || p_value_threshold > 1) {
    stop("p_value_threshold must be a number between 0 and 1")
  }

  pathway_enrichment_api(
    gene_ids = gene_ids,
    source = source,
    p_value_threshold = p_value_threshold
  )
}
