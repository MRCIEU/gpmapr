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
#'     \item snp_id: the id of the SNP
#'     \item display_snp: the name of the SNP
#'     \item distinct_trait_categories: the number of trait categories that the SNP is associated with via coloc groups
#'     \item distinct_protein_coding_genes: the number of genes that the SNP is associated with via coloc groups
#'   }
#' @export
get_all_snp_pleiotropies <- function() {
  snp_pleiotropies <- snp_pleiotropies_api()
  return(snp_pleiotropies$snps)
}