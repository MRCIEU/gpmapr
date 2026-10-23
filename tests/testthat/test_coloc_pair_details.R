library(testthat)

# add_coloc_pair_details() joins each side of a coloc pair to its study extraction.

existing_extractions <- data.frame(
  id = c(10, 20),
  trait_name = c("LDL cholesterol", "PCSK9 Liver eQTL"),
  data_type = c("Phenotype", "Gene Expression"),
  tissue = c(NA, "Liver"),
  gene = c(NA, "PCSK9")
)

test_that("both sides of a coloc pair get their study extraction details, after the id columns", {
  api_object <- list(
    coloc_pairs = data.frame(variant_id = 1, study_extraction_a_id = 10, study_extraction_b_id = 20, h4 = 0.9),
    study_extractions = existing_extractions
  )
  pairs <- add_coloc_pair_details(api_object)$coloc_pairs
  expect_equal(names(pairs), c(
    "variant_id", "study_extraction_a_id", "study_extraction_b_id",
    "trait_name_a", "data_type_a", "tissue_a", "gene_a",
    "trait_name_b", "data_type_b", "tissue_b", "gene_b", "h4"
  ))
  expect_equal(pairs$trait_name_a, "LDL cholesterol")
  expect_equal(pairs$trait_name_b, "PCSK9 Liver eQTL")
  expect_equal(pairs$tissue_b, "Liver")
  expect_false("source_a" %in% names(pairs))
})

test_that("a study extraction missing from the response gives NA details", {
  api_object <- list(
    coloc_pairs = data.frame(study_extraction_a_id = c(10, 99), study_extraction_b_id = c(20, 20)),
    study_extractions = existing_extractions
  )
  pairs <- add_coloc_pair_details(api_object)$coloc_pairs
  expect_equal(pairs$trait_name_a, c("LDL cholesterol", NA))
  expect_equal(pairs$trait_name_b, c("PCSK9 Liver eQTL", "PCSK9 Liver eQTL"))
})

test_that("study extractions given as a list of data frames are joined too", {
  api_object <- list(
    coloc_pairs = data.frame(study_extraction_a_id = 10, study_extraction_b_id = 20),
    study_extractions = list(existing_extractions[1, ], existing_extractions[2, ])
  )
  pairs <- add_coloc_pair_details(api_object)$coloc_pairs
  expect_equal(pairs$trait_name_b, "PCSK9 Liver eQTL")
})

test_that("upload pairs are matched to their study extractions by source and id", {
  # Upload and existing study extractions have separate ids: 7223 is both here.
  api_object <- list(
    coloc_pairs = data.frame(
      study_extraction_a_id = c(7223, 7223, 10), source_a = c("upload", "upload", "existing"),
      study_extraction_b_id = c(10, 7223, 7223), source_b = c("existing", "existing", "upload"),
      h4 = c(0.9, 0.95, 0.85)
    ),
    study_extractions = data.frame(
      id = c(10, 7223, 7223), source = c("existing", "existing", "upload"),
      gwas_upload_id = c(NA, NA, 90), trait_name = c("LDL cholesterol", "Granulocyte percentage", "My GWAS"),
      data_type = "Phenotype"
    )
  )
  pairs <- add_coloc_pair_details(api_object)$coloc_pairs
  expect_equal(pairs$trait_name_a, c("My GWAS", "My GWAS", "LDL cholesterol"))
  expect_equal(pairs$trait_name_b, c("LDL cholesterol", "Granulocyte percentage", "My GWAS"))
  expect_equal(which(names(pairs) == "trait_name_a"), which(names(pairs) == "source_b") + 1)
})

test_that("objects without coloc pairs are returned unchanged", {
  api_object <- list(coloc_groups = data.frame(id = 1), coloc_pairs = data.frame())
  expect_equal(add_coloc_pair_details(api_object), api_object)
  expect_equal(add_coloc_pair_details(list(trait = list(id = 1))), list(trait = list(id = 1)))
})

test_that("get_gwas() describes both sides of its coloc pairs", {
  local_mock_api()
  pairs <- get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155")$coloc_pairs
  expect_equal(pairs$trait_name_a, c("AD GWAS", "Allergic rhinitis", "AD GWAS", "Granulocyte percentage"))
  expect_equal(pairs$trait_name_b, c("Allergic rhinitis", "AD GWAS", "Compared GWAS", "AD GWAS"))
})

test_that("trait() coloc pairs are left without details", {
  local_mocked_bindings(
    trait_api = function(...) list(trait = list(id = 1), study_extractions = existing_extractions),
    trait_coloc_pairs_api = function(...) data.frame(study_extraction_a_id = 10, study_extraction_b_id = 20)
  )
  pairs <- trait(1, include_coloc_pairs = TRUE)$coloc_pairs
  expect_equal(names(pairs), c("study_extraction_a_id", "study_extraction_b_id"))
})
