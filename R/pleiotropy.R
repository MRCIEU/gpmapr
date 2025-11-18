#' @title Get All Gene Pleiotropies
#' @description Get gene pleiotropy from the API by gene id
#' @return A list containing the gene pleiotropy
#' @inheritSection gene_pleiotropy_doc gene_pleiotropy_dataframe
#' @export
get_all_gene_pleiotropies <- function() {
  gene_pleiotropies <- gene_pleiotropies_api()
  return(gene_pleiotropies)
}

#' @title Get All SNP Pleiotropies
#' @description Get all SNP pleiotropies from the API
#' @return A list containing the SNP pleiotropies
#' @inheritSection snp_pleiotropy_doc snp_pleiotropy_dataframe
#' @export
get_all_snp_pleiotropies <- function() {
  snp_pleiotropies <- snp_pleiotropies_api()
  return(snp_pleiotropies)
}