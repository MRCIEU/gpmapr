#' @title Get All Gene Pleiotropies
#' @description Get gene pleiotropy from the API by gene id
#' @return A list containing the gene pleiotropy
#'   \itemize{
#'     \item gene_id: the id of the gene
#'     \item gene: the name of the gene
#'     \item distinct_trait_categories: the number of trait categories that the gene is associated with via coloc groups
#'     \item distinct_protein_coding_genes: the number of genes that the gene is associated with via coloc groups
#'   }
#' @export
get_all_gene_pleiotropies <- function() {
  gene_pleiotropies <- gene_pleiotropies_api()
  return(gene_pleiotropies$genes)
}

#' @title Get All SNP Pleiotropies
#' @description Get all SNP pleiotropies from the API
#' @return A list containing the SNP pleiotropies
#'   \itemize{
#'     \item variant_id: the id of the SNP
#'     \item display_snp: the name of the SNP
#'     \item distinct_trait_categories: the number of trait categories that the SNP is associated with via coloc groups
#'     \item distinct_protein_coding_genes: the number of genes that the SNP is associated with via coloc groups
#'   }
#' @export
get_all_variant_pleiotropies <- function() {
  variant_pleiotropies <- variant_pleiotropies_api()
  return(variant_pleiotropies$snps)
}

#' @title Module-Level SNP Pleiotropy Summary
#' @description Summarise per-SNP pleiotropy within SNP groups/modules: the
#'   number of distinct trait categories and the number of distinct
#'   protein-coding genes each SNP is associated with via coloc groups. A module
#'   of broadly pleiotropic SNPs (high category count) is more hub-like or
#'   generic; a module of narrowly pleiotropic SNPs is more specific. This is a
#'   descriptive annotation, not a reliability gate — combine with
#'   `summarise_snp_module_quality()` to interpret modules.
#' @param groups Named vector mapping SNP ids to group/module ids, or a dataframe
#'   with snp/variant and group/cluster columns (see `.normalize_snp_groups()`).
#' @param variant_pleiotropies Dataframe of per-SNP pleiotropy, e.g. the output
#'   of `get_all_variant_pleiotropies()`.
#' @param snp_key Column in `variant_pleiotropies` matching the SNP ids.
#'   Defaults to `"variant_id"`.
#' @param min_group_size Only summarise groups with at least this many SNPs.
#'   Defaults to 0.
#' @return A dataframe with one row per group:
#'   \itemize{
#'     \item group
#'     \item n_snps
#'     \item n_snps_with_pleiotropy: SNPs matched in `variant_pleiotropies`
#'     \item mean_trait_category_pleiotropy / median_trait_category_pleiotropy
#'     \item mean_gene_pleiotropy / median_gene_pleiotropy
#'   }
#' @export
summarise_snp_group_pleiotropy <- function(groups,
                                           variant_pleiotropies,
                                           snp_key = "variant_id",
                                           min_group_size = 0) {
  if (is.null(variant_pleiotropies) || nrow(variant_pleiotropies) == 0) {
    stop("variant_pleiotropies is required")
  }
  if (!snp_key %in% names(variant_pleiotropies)) {
    stop("variant_pleiotropies must include column: ", snp_key)
  }
  if (!all(c("distinct_trait_categories", "distinct_protein_coding_genes") %in% names(variant_pleiotropies))) {
    stop(
      "variant_pleiotropies must include distinct_trait_categories and ",
      "distinct_protein_coding_genes"
    )
  }
  if (!is.numeric(min_group_size) || min_group_size < 0) {
    stop("min_group_size must be a non-negative number")
  }

  group_df <- .normalize_snp_groups(groups)
  pleio <- variant_pleiotropies |>
    dplyr::mutate(
      .join_key = as.character(.data[[snp_key]])
    ) |>
    dplyr::select(
      .join_key,
      distinct_trait_categories,
      distinct_protein_coding_genes
    ) |>
    dplyr::distinct(.join_key, .keep_all = TRUE)

  joined <- group_df |>
    dplyr::left_join(
      pleio,
      by = c("snp_id" = ".join_key")
    )

  out <- joined |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      n_snps = dplyr::n(),
      n_snps_with_pleiotropy = sum(!is.na(distinct_trait_categories)),
      mean_trait_category_pleiotropy = mean(distinct_trait_categories, na.rm = TRUE),
      median_trait_category_pleiotropy = stats::median(distinct_trait_categories, na.rm = TRUE),
      mean_gene_pleiotropy = mean(distinct_protein_coding_genes, na.rm = TRUE),
      median_gene_pleiotropy = stats::median(distinct_protein_coding_genes, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(n_snps >= min_group_size) |>
    dplyr::arrange(as.integer(group))

  return(out)
}