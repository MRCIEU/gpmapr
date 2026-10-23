#' @section Data model:
#' Most functions return a list of dataframes describing the same few kinds of record:
#' \itemize{
#'   \item `study_extractions`: a study's association signal in one LD block (a "study extraction"),
#'     identified by `id`.
#'   \item `coloc_groups`: groups of study extractions that colocalise, one row per study extraction in
#'     each group, identified by `coloc_group_id` and pointing at their study extraction with
#'     `study_extraction_id`.
#'   \item `coloc_pairs`: pairwise colocalisation results between two study extractions,
#'     `study_extraction_a_id` and `study_extraction_b_id`, with their trait names alongside.
#'   \item `rare_results`: rare-variant results, grouped by `rare_result_group_id`.
#' }
#' Join `coloc_groups$study_extraction_id` (and the pair ids) to `study_extractions$id`.  `trait_id`,
#' `study_id`, `variant_id`, `gene_id` and `ld_block_id` refer to traits, studies, variants, genes and LD
#' blocks.  With `include_associations = TRUE`, each study's association at the row's variant (`beta`, `se`,
#' `p`, `eaf`, `imputed`) is added to `coloc_groups`, `rare_results` and `study_extractions`.
#'
#' Associations are not stored for every study at every variant. Every study with a study extraction in an
#' LD block has an association at the variants of that block's coloc groups and rare results, but at a study
#' extraction's own lead variant only that study does. Compare studies at the variant of their coloc group;
#' see "Which associations are stored" in [associations()].
#'
#' GWAS upload results ([get_gwas()], or [trait()] with a GUID) have the same tables, mixing the upload's
#' own rows with rows for studies already in the GPMap.  Their `source` column ("upload" or "existing")
#' says which: identify a study extraction by `source` and its id together, and find the upload's own rows
#' with `source == "upload" & gwas_upload_id == result$trait$id`.  [compare_upload()] lines an upload up
#' with another trait.
#' @keywords internal
"_PACKAGE"

#' @importFrom dplyr .data
NULL
