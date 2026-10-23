# Build the hand-written GWAS upload fixtures used by the mocked tests.
#
# The fixture recorder (record_api_fixtures.R) records GET /v1/gwas/{guid} only
# from a server that has the example upload. Until production serves API 1.1.0,
# that upload exists nowhere in the 1.1.0 shape, so these small responses in
# that shape stand in for it. They cover the cases upload code gets wrong:
#   * the upload's own rows, an existing trait's rows and a compared upload's
#     rows in the same coloc group, with an unrelated existing study that has a
#     smaller p-value than the upload;
#   * an upload and an existing study extraction with the same id (7223);
#   * a placeholder association (beta 0, se 1, p 1, imputed);
#   * a coloc pair at a locus with no shared coloc group.
#
# Run from the package root:
#   Rscript data-raw/build_upload_fixture.R

fixture_dir <- file.path("tests", "testthat", "fixtures")
guid <- "7a289615-c1b4-91f3-3d97-887f60de9155"
compared_guid <- "11111111-2222-3333-4444-555555555555"

trait <- list(
  id = 90, guid = guid, email = NULL, name = "AD GWAS", sample_size = 50000, ancestry = "EUR",
  variable_type = "categorical", is_published = FALSE, doi = NULL, should_be_added = FALSE,
  upload_metadata = NULL, status = "completed", failure_reason = NULL,
  created_at = "2026-01-01T00:00:00", updated_at = "2026-01-01T00:00:00", message = NULL
)

extraction <- function(id, source, variant_id, ld_block_id, min_p, trait_name, gwas_upload_id = NULL,
                       study_id = NULL, trait_id = NULL, study = NULL, data_type = "Phenotype") {
  return(list(
    id = id, study_id = study_id, variant_id = variant_id, display_snp = paste0("1:", variant_id, " A/G"),
    rsid = paste0("rs", variant_id), ld_block_id = ld_block_id, unique_study_id = paste0(study, "_", id),
    study = study, file = paste0(id, ".tsv.gz"), svg_file = NULL, file_with_lbfs = NULL, chr = 1,
    bp = variant_id * 1000, min_p = min_p, cis_trans = NULL, ld_block = paste0("EUR/1/", ld_block_id),
    credible_set = NULL, gene_id = NULL, situated_gene_id = NULL, gene = NULL, trait_id = trait_id,
    trait_name = trait_name, trait_category = NULL, data_type = data_type, tissue = NULL, cell_type = NULL,
    source = source, gwas_upload_id = gwas_upload_id
  ))
}

study_extractions <- list(
  extraction(32958, "existing", 100, 9, 1e-12, "Allergic rhinitis", study_id = 3001, trait_id = 2527,
             study = "ebi-a-GCST000001"),
  extraction(50000, "existing", 100, 9, 1e-50, "Type 1 diabetes", study_id = 3500, trait_id = 3500,
             study = "ebi-a-GCST000002"),
  extraction(7223, "existing", 200, 58, 1e-10, "Granulocyte percentage", study_id = 21, trait_id = 21,
             study = "ebi-a-GCST000003"),
  extraction(40000, "existing", 300, 60, 1e-9, "Allergic rhinitis", study_id = 3001, trait_id = 2527,
             study = "ebi-a-GCST000001"),
  extraction(7223, "upload", 100, 9, 1e-9, "AD GWAS", gwas_upload_id = 90, study = guid),
  extraction(7224, "upload", 200, 58, 1e-8, "AD GWAS", gwas_upload_id = 90, study = guid),
  extraction(7225, "upload", 300, 60, 1e-8, "AD GWAS", gwas_upload_id = 90, study = guid),
  extraction(8001, "upload", 100, 9, 1e-9, "Compared GWAS", gwas_upload_id = 144, study = compared_guid)
)

# A coloc group row for a study extraction: its columns plus the group's.
group_row <- function(coloc_group_id, e, h4_connectedness) {
  return(list(
    coloc_group_id = coloc_group_id, study_id = e$study_id, study_extraction_id = e$id, variant_id = e$variant_id,
    ld_block_id = e$ld_block_id, situated_gene_id = NULL, h4_connectedness = h4_connectedness,
    h3_connectedness = 0, chr = e$chr, bp = e$bp, min_p = e$min_p, cis_trans = NULL, ld_block = e$ld_block,
    display_snp = e$display_snp, rsid = e$rsid, gene = NULL, gene_id = NULL, situated_gene = NULL,
    trait_id = e$trait_id, trait_name = e$trait_name, trait_category = NULL, data_type = e$data_type,
    tissue = NULL, cell_type = NULL, source_id = NULL, source_name = NULL, source_url = NULL,
    unique_study_id = e$unique_study_id, study = e$study, file = e$file, svg_file = NULL, file_with_lbfs = NULL,
    sample_size = NULL, variable_type = NULL, ancestry = "EUR", heritability = NULL, heritability_se = NULL,
    source = e$source, gwas_upload_id = e$gwas_upload_id
  ))
}
by_key <- function(source, id) {
  return(Filter(function(e) e$source == source && e$id == id, study_extractions)[[1]])
}
coloc_groups <- list(
  # Variant 100: the upload, allergic rhinitis, the compared upload, and type 1 diabetes (smallest p).
  group_row(1, by_key("upload", 7223), 0.9),
  group_row(1, by_key("existing", 32958), 0.9),
  group_row(1, by_key("upload", 8001), 0.9),
  group_row(1, by_key("existing", 50000), 0.9),
  # Variant 200: the upload and an existing study extraction with the same id as one of the upload's.
  group_row(2, by_key("upload", 7224), 1),
  group_row(2, by_key("existing", 7223), 1)
)

pair <- function(source_a, id_a, source_b, id_b, ld_block_id, h3, h4) {
  return(list(
    study_extraction_a_id = id_a, source_a = source_a, study_extraction_b_id = id_b, source_b = source_b,
    ld_block_id = ld_block_id, h3 = h3, h4 = h4, false_positive = FALSE, false_negative = FALSE, ignore = FALSE
  ))
}
coloc_pairs <- list(
  pair("upload", 7223, "existing", 32958, 9, 0.05, 0.9),
  pair("existing", 40000, "upload", 7225, 60, 0.02, 0.95),
  pair("upload", 7223, "upload", 8001, 9, 0.01, 0.99),
  pair("existing", 7223, "upload", 7224, 58, 0.1, 0.85)
)

association <- function(variant_id, source, beta, se, p, study_id = NULL, gwas_upload_id = NULL, imputed = FALSE) {
  return(list(
    variant_id = variant_id, source = source, study_id = study_id, gwas_upload_id = gwas_upload_id,
    beta = beta, se = se, p = p, eaf = 0.3, imputed = imputed
  ))
}
associations <- list(
  association(100, "upload", 0.5, 0.1, 1e-9, gwas_upload_id = 90),
  association(200, "upload", 0.3, 0.05, 1e-8, gwas_upload_id = 90),
  association(300, "upload", 0.25, 0.04, 1e-8, gwas_upload_id = 90),
  association(100, "upload", 0.45, 0.1, 1e-9, gwas_upload_id = 144),
  association(100, "existing", 0.2, 0.02, 1e-12, study_id = 3001),
  association(100, "existing", 1.2, 0.08, 1e-50, study_id = 3500),
  association(200, "existing", 0, 1, 1, study_id = 21, imputed = TRUE),
  association(300, "existing", 0.1, 0.01, 1e-9, study_id = 3001)
)

response <- list(
  trait = trait, coloc_groups = coloc_groups, coloc_pairs = coloc_pairs, study_extractions = study_extractions,
  rare_results = list(), associations = NULL, queue_position = NULL, queue_status = NULL
)
write_fixture <- function(body, file) {
  json <- jsonlite::toJSON(body, auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = TRUE)
  writeLines(json, file.path(fixture_dir, file))
  return(invisible(file))
}
write_fixture(response, "gwas-upload.json")
response$associations <- associations
write_fixture(response, "gwas-upload-associations.json")
message("Wrote gwas-upload.json and gwas-upload-associations.json to ", fixture_dir)
