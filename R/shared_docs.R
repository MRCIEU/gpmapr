#' @title Coloc Groups Dataframe Documentation
#' @description Shared documentation for coloc_groups dataframe columns
#' @section coloc_groups_dataframe:
#' The coloc_groups dataframe contains information about which studies have coloc results.
#' It has the following columns:
#' \itemize{
#'   \item coloc_group_id: the unique id for this group of colocalised results
#'   \item study_id: the id of the study
#'   \item study_extraction_id: the id of the study extraction
#'   \item variant_id: the id of the SNP
#'   \item ld_block_id: the id of the LD block
#'   \item h4_connectedness: the proportion of the pairs of studies in the group that colocalise
#'     (posterior probability of H4, one shared causal variant, of at least 0.8)
#'   \item h3_connectedness: the proportion of the pairs of studies in the group with a posterior probability
#'     of H3 (two distinct causal variants) of at least 0.8
#'   \item chr: the chromosome of the SNP
#'   \item bp: the base pair position of the SNP
#'   \item min_p: the minimum p-value related to the study_extraction_id
#'   \item cis_trans: the cis/trans status of the SNP
#'   \item ld_block: the LD block of the SNP
#'   \item display_snp: the display SNP name
#'   \item rsid: the rsID of the SNP
#'   \item gene: the gene whose expression, splicing or protein level the study measures (for molecular QTLs)
#'   \item gene_id: the id of the gene
#'   \item situated_gene: the gene the measured site is in, for methylation studies (which have no `gene`)
#'   \item situated_gene_id: the id of the situated gene
#'   \item trait_id: the id of the trait
#'   \item trait_name: the name of the trait
#'   \item trait_category: the category of the trait
#'   \item data_type: the data type of the trait
#'   \item tissue: the tissue of the trait
#'   \item cell_type: the cell type of the trait
#'   \item source_id: the id of the source the study comes from
#'   \item source_name: the name of the source the study comes from
#'   \item source_url: a link to the study at its source
#'   \item beta, se, p, eaf, imputed: (with `include_associations = TRUE`) the association of the study at the SNP
#'   \item source, gwas_upload_id: (GWAS upload results only) whether the row is from an uploaded GWAS
#'     ("upload", with the upload's id in `gwas_upload_id`) or already in the GPMap ("existing").  Identify a
#'     study extraction by `source` and `study_extraction_id` together.
#' }
#' @keywords internal
coloc_groups_doc <- function() {}

#' @title Study Extractions Dataframe Documentation
#' @description Shared documentation for study_extractions dataframe columns
#' @section study_extractions_dataframe:
#' The study_extractions dataframe contains information about which studies have coloc results.
#' It has the following columns:
#' \itemize{
#'   \item id: the unique id for this study extraction
#'   \item study_id: the id of the study associated with this study extraction
#'   \item variant_id: the id of the SNP
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
#'   \item source, gwas_upload_id: (GWAS upload results only) whether the study extraction is from an
#'     uploaded GWAS ("upload", with the upload's id in `gwas_upload_id`) or already in the GPMap
#'     ("existing").  Upload and existing study extractions have separate ids, so identify a study
#'     extraction by `source` and `id` together.
#' }
#' @keywords internal
study_extractions_doc <- function() {}

#' @title Rare Results Dataframe Documentation
#' @description Shared documentation for rare_results dataframe columns
#' @section rare_results_dataframe:
#' The rare_results dataframe contains information about which studies have coloc results.
#' It has the following columns:
#' \itemize{
#'   \item rare_result_group_id: the unique id for this rare result group
#'   \item study_id: the id of the study associated with this rare result
#'   \item study_extraction_id: the id of the study extraction associated with this rare result
#'   \item variant_id: the id of the SNP
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
#' @section coloc_pairs_dataframe:
#' The coloc_pairs dataframe contains information about which studies have coloc pairs.
#' It has the following columns:
#' \itemize{
#'   \item study_extraction_a_id: the id of the study extraction associated with this coloc pair
#'   \item study_extraction_b_id: the id of the study extraction associated with this coloc pair
#'   \item source_a, source_b: (GWAS upload results only) whether each study extraction is from an uploaded
#'     GWAS ("upload") or already in the GPMap ("existing")
#'   \item trait_name_a, data_type_a, tissue_a, gene_a: the trait name, data type, tissue and gene of
#'     study extraction a (and the same with `_b` for study extraction b), joined from `study_extractions`.
#'     Not added to the coloc pairs of `trait()`, whose endpoint doesn't return the study extractions
#'     its pairs refer to.
#'   \item ld_block_id: the id of the LD block
#'   \item h3: the h3 value for this coloc pair
#'   \item h4: the h4 value for this coloc pair
#'   \item spurious: whether this coloc pair is spurious
#' }
#' @keywords internal
coloc_pairs_doc <- function() {}

#' @title Associations Dataframe Documentation
#' @description Shared documentation for associations dataframe columns
#' @section associations_dataframe:
#' The associations dataframe contains information about which studies have association results.
#' It has the following columns:
#' \itemize{
#'   \item variant_id: the id of the SNP associated with this association
#'   \item study_id: the id of the study associated with this association
#'   \item beta: the beta value of the association
#'   \item se: the standard error of the association
#'   \item p: the p-value of the association
#'   \item eaf: the estimated allele frequency of the association
#'   \item imputed: whether the association is imputed
#' }
#' Associations of variants a study did not report are placeholders in the GPMap (beta 0, se 1, p 1),
#' and are returned with `beta`, `se`, `p` and `eaf` set to `NA` (and `imputed` TRUE).
#' @keywords internal
associations_doc <- function() {}

#' @title Summary Statistics Dataframe Documentation
#' @description Shared documentation for summary_statistics dataframe columns
#' @section summary_statistics_dataframe:
#' The summary_statistics dataframe contains information about which studies have summary statistics.
#' From the API, column names are typically uppercase (SNP, CHR, BP, EA, OA, EAF, Z, BETA, SE, P, LBF_1, etc.).
#' It has the following columns (names may be upper or lower case depending on source):
#' \itemize{
#'   \item SNP / variant_id: the id of the SNP
#'   \item CHR / chr: the chromosome of the SNP
#'   \item BP / bp: the base pair position of the SNP
#'   \item EA / ea: the effect allele
#'   \item OA / oa: the other allele
#'   \item EAF / eaf: the estimated allele frequency
#'   \item Z / z: the z-score
#'   \item BETA / beta: the beta value
#'   \item SE / se: the standard error
#'   \item P / p: the p-value
#'   \item imputed: whether the summary statistics are imputed
#'   \item LBF_* / lbf_*: all different finemapped log-bayes factors for each credible set.
#'     Each credible set is numbered from 1 to 10.  If finemapped failed or only returned 1 credible set, the LBF_1
#'     column is just converted directly from the z-score.
#' }
#' @keywords internal
summary_statistics_doc <- function() {}

#' @title Variants Dataframe Documentation
#' @description Shared documentation for variants dataframe columns
#' @section variants_dataframe:
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
#' @section genes_in_region_dataframe:
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
#' @section ld_dataframe:
#' The ld dataframe contains information about the LD matrix.
#' It has the following columns:
#' \itemize{
#'   \item lead_variant_id: the id of the lead SNP
#'   \item proxy_variant_id: the id of the variant SNP
#'   \item ld_block_id: the id of the LD block
#'   \item r: the r value between the lead and variant SNPs
#'   \item rsq: the r-squared between the lead and variant SNPs, rounded to 2 decimal places
#' }
#' @keywords internal
ld_doc <- function() {}

#' @title Gene Pleiotropies Dataframe Documentation
#' @description Shared documentation for gene_pleiotropies dataframe columns
#' @section gene_pleiotropies_dataframe:
#' The gene_pleiotropies dataframe contains information about which genes are pleiotropic.
#' It has the following columns:
#' \itemize{
#'   \item gene_id: the id of the gene
#'   \item gene: the name of the gene
#'   \item distinct_trait_categories: the number of distinct trait categories in the coloc groups the gene's own
#'     QTLs are in
#'     (the same as distinct_trait_categories_cis_trans_all)
#'   \item distinct_protein_coding_genes: the number of distinct protein coding genes in the coloc groups the
#'     gene's own QTLs are in
#'     (the same as distinct_protein_coding_genes_cis_trans_all)
#'   \item distinct_trait_categories_cis_gws, distinct_trait_categories_cis_all,
#'     distinct_trait_categories_cis_trans_gws, distinct_trait_categories_cis_trans_all: the number of distinct
#'     trait categories, counting only cis results (\code{_cis_}) or cis and trans results (\code{_cis_trans_}),
#'     and only genome-wide significant study extractions, with p <= 5e-8 (\code{_gws}), or all of them
#'     (\code{_all})
#'   \item distinct_protein_coding_genes_cis_gws, distinct_protein_coding_genes_cis_all,
#'     distinct_protein_coding_genes_cis_trans_gws, distinct_protein_coding_genes_cis_trans_all: the number of
#'     distinct protein coding genes, with the same cis / cis_trans and gws / all variants
#' }
#' @keywords internal
gene_pleiotropies_doc <- function() {}

#' @title SNP Pleiotropies Dataframe Documentation
#' @description Shared documentation for variant_pleiotropies dataframe columns
#' @section variant_pleiotropies_dataframe:
#' The variant_pleiotropies dataframe contains information about which SNPs are pleiotropic.
#' It has the following columns:
#' \itemize{
#'   \item variant_id: the id of the SNP
#'   \item display_snp: the display name of the SNP
#'   \item distinct_trait_categories: the number of distinct trait categories in the coloc groups at the SNP
#'     (the same as distinct_trait_categories_cis_trans_all)
#'   \item distinct_protein_coding_genes: the number of distinct protein coding genes in the coloc groups at the SNP
#'     (the same as distinct_protein_coding_genes_cis_trans_all)
#'   \item distinct_trait_categories_cis_gws, distinct_trait_categories_cis_all,
#'     distinct_trait_categories_cis_trans_gws, distinct_trait_categories_cis_trans_all: the number of distinct
#'     trait categories, counting only cis results (\code{_cis_}) or cis and trans results (\code{_cis_trans_}),
#'     and only genome-wide significant study extractions, with p <= 5e-8 (\code{_gws}), or all of them
#'     (\code{_all})
#'   \item distinct_protein_coding_genes_cis_gws, distinct_protein_coding_genes_cis_all,
#'     distinct_protein_coding_genes_cis_trans_gws, distinct_protein_coding_genes_cis_trans_all: the number of
#'     distinct protein coding genes, with the same cis / cis_trans and gws / all variants
#' }
#' @keywords internal
variant_pleiotropies_doc <- function() {}
