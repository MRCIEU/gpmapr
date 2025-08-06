#' @title Coloc Groups Dataframe Documentation
#' @description Shared documentation for coloc_groups dataframe columns
#' @section coloc_groups dataframe:
#' The coloc_groups dataframe contains information about which studies have coloc results.
#' It has the following columns:
#' \itemize{
#'   \item coloc_group_id: the unique id for this group of colocalised results
#'   \item study_id: the id of the study
#'   \item study_extraction_id: the id of the study extraction
#'   \item snp_id: the id of the SNP
#'   \item ld_block_id: the id of the LD block
#'   \item group_threshold: the threshold for this group of colocalised results.
#'     'strong' signifies h4 > 0.8, 'moderate' signifies h4 > 0.6
#'   \item chr: the chromosome of the SNP
#'   \item bp: the base pair position of the SNP
#'   \item min_p: the minimum p-value related to the study_extraction_id
#'   \item cis_trans: the cis/trans status of the SNP
#'   \item ld_block: the LD block of the SNP
#'   \item display_snp: the display SNP name
#'   \item gene: the gene associated with the SNP
#'   \item gene_id: the id of the gene
#'   \item trait_id: the id of the trait
#'   \item trait_name: the name of the trait
#'   \item trait_category: the category of the trait
#'   \item data_type: the data type of the trait
#'   \item tissue: the tissue of the trait
#' }
#' @keywords internal
coloc_groups_doc <- function() {}

#' @title Study Extractions Dataframe Documentation
#' @description Shared documentation for study_extractions dataframe columns
#' @section study_extractions dataframe:
#' The study_extractions dataframe contains information about which studies have coloc results.
#' It has the following columns:
#' \itemize{
#'   \item id: the unique id for this study extraction
#'   \item study_id: the id of the study associated with this study extraction
#'   \item snp_id: the id of the SNP
#'   \item snp: the SNP name
#'   \item ld_block_id: the id of the LD block
#'   \item unique_study_id: the unique id for this study
#'   \item study: the study name
#'   \item file: the file name
#'   \item svg_file: the SVG file name
#'   \item file_with_lbfs: the file name with lbfs
#'   \item chr: the chromosome of the SNP
#'   \item bp: the base pair position of the SNP
#'   \item min_p: the minimum p-value related to the study_extraction_id
#'   \item cis_trans: the cis/trans status of the SNP
#'   \item ld_block: the LD block of the SNP
#'   \item gene: the gene associated with the SNP
#'   \item gene_id: the id of the gene
#'   \item trait_id: the id of the trait
#'   \item trait_name: the name of the trait
#'   \item trait_category: the category of the trait
#'   \item data_type: the data type of the trait
#'   \item tissue: the tissue of the trait
#' }
#' @keywords internal
study_extractions_doc <- function() {}

#' @title Rare Results Dataframe Documentation
#' @description Shared documentation for rare_results dataframe columns
#' @section rare_results dataframe:
#' The rare_results dataframe contains information about which studies have coloc results.
#' It has the following columns:
#' \itemize{
#'   \item rare_result_group_id: the unique id for this rare result group
#'   \item study_id: the id of the study associated with this rare result
#'   \item study_extraction_id: the id of the study extraction associated with this rare result
#'   \item snp_id: the id of the SNP
#'   \item ld_block_id: the id of the LD block
#'   \item chr: the chromosome of the SNP
#'   \item bp: the base pair position of the SNP
#'   \item min_p: the minimum p-value related to the study_extraction_id
#'   \item display_snp: the display SNP name
#'   \item gene: the gene associated with the SNP
#'   \item gene_id: the id of the gene
#'   \item trait_id: the id of the trait
#'   \item trait_name: the name of the trait
#'   \item trait_category: the category of the trait
#'   \item data_type: the data type of the trait
#'   \item tissue: the tissue of the trait
#'   \item ld_block: the LD block of the SNP
#' }
#' @keywords internal
rare_results_doc <- function() {}

#' @title Coloc Pairs Dataframe Documentation
#' @description Shared documentation for coloc_pairs dataframe columns
#' @section coloc_pairs dataframe:
#' The coloc_pairs dataframe contains information about which studies have coloc pairs.
#' It has the following columns:
#' \itemize{
#'   \item study_extraction_a_id: the id of the study extraction associated with this coloc pair
#'   \item study_extraction_b_id: the id of the study extraction associated with this coloc pair
#'   \item ld_block_id: the id of the LD block
#'   \item h3: the h3 value for this coloc pair
#'   \item h4: the h4 value for this coloc pair
#'   \item spurious: whether this coloc pair is spurious
#' }
#' @keywords internal
coloc_pairs_doc <- function() {}

#' @title Associations Dataframe Documentation
#' @description Shared documentation for associations dataframe columns
#' @section associations dataframe:
#' The associations dataframe contains information about which studies have association results.
#' It has the following columns:
#' \itemize{
#'   \item snp_id: the id of the SNP associated with this association
#'   \item study_id: the id of the study associated with this association
#'   \item beta: the beta value of the association
#'   \item se: the standard error of the association
#'   \item p: the p-value of the association
#'   \item eaf: the estimated allele frequency of the association
#'   \item imputed: whether the association is imputed
#' }
#' @keywords internal
associations_doc <- function() {}

#' @title Summary Statistics Dataframe Documentation
#' @description Shared documentation for summary_statistics dataframe columns
#' @section summary_statistics dataframe:
#' The summary_statistics dataframe contains information about which studies have summary statistics.
#' It has the following columns:
#' \itemize{
#'   \item snp_id: the id of the SNP
#'   \item chr: the chromosome of the SNP
#'   \item bp: the base pair position of the SNP
#'   \item ea: the effect allele
#'   \item oa: the other allele
#'   \item eaf: the estimated allele frequency
#'   \item z: the z-score
#'   \item beta: the beta value
#'   \item se: the standard error
#'   \item p: the p-value
#'   \item imputed: whether the summary statistics are imputed
#'   \item lbf_*: all different finemapped log-bayes factors for each credible set.
#'     Each credible set is numbered from 1 to 10.  If finemapped failed or only returned 1 credible set, the lbf_1
#'     column is just converted directly from the z-score.
#' }
#' @keywords internal
summary_statistics_doc <- function() {}

#' @title Variants Dataframe Documentation
#' @description Shared documentation for variants dataframe columns
#' @section variants dataframe:
#' The variants dataframe contains variant information that is pulled from the Variant Effect Predictor (VEP) database.
#' It has the following columns, along side many more columns from VEP:
#' \itemize{
#'   \item id: the id of the SNP
#'   \item gene_id: the id of the gene as predicted by VEP
#'   \item gene: the gene name as predicted by VEP
#' }
#' @keywords internal
variants_doc <- function() {}


#' @title Genes in Region Dataframe Documentation
#' @description Shared documentation for genes_in_region dataframe columns
#' @section genes_in_region dataframe:
#' The genes_in_region dataframe contains information about which genes are in a region.
#' It has the following columns:
#' \itemize{
#'   \item id: the id of the gene
#'   \item ensembl_id: the ensembl id of the gene
#'   \item gene: the name of the gene
#'   \item description: the description of the gene
#'   \item gene_biotype: the gene biotype
#'   \item chr: the chromosome of the gene
#'   \item start: the start position of the gene
#'   \item stop: the stop position of the gene
#'   \item strand: the strand of the gene
#'   \item source: the source of the gene
#' }
#' @keywords internal
genes_in_region_doc <- function() {}

#' @title LD Dataframe Documentation
#' @description Shared documentation for ld dataframe columns
#' @section ld dataframe:
#' The ld dataframe contains information about the LD matrix.
#' It has the following columns:
#' \itemize{
#'   \item lead_snp_id: the id of the lead SNP
#'   \item variant_snp_id: the id of the variant SNP
#'   \item ld_block_id: the id of the LD block
#'   \item r: the r value between the lead and variant SNPs
#' }
#' @keywords internal
ld_doc <- function() {}