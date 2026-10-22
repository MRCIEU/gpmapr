library(testthat)

# Live tests: these call the real GPMap API. Each starts with skip_if_no_api()
# (helper-live-api.R), so they are skipped on CRAN and when the API is down.
# The same functions are tested offline in test_mocked_api.R.

test_that("health() returns expected output", {
  skip_if_no_api()
  result <- health_api()
  expect_type(result, "list")
  expect_true(result$status %in% c("healthy", "unhealthy"))
})

test_that("search_gpmap() returns expected output", {
  skip_if_no_api()
  result <- search_gpmap("haemoglobin")
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
  expected_names <- c(
    "call", "name", "type", "type_id", "num_coloc_groups",
    "num_coloc_studies", "num_rare_results", "num_study_extractions"
  )
  expect_true(all(expected_names %in% names(result)))
})

test_that("trait() returns expected output", {
  skip_if_no_api()
  trait_id <- 5020
  result <- trait(trait_id)
  expect_type(result, "list")
  expect_true(result$trait$id == trait_id)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("trait() returns full_associations when requested", {
  skip_if_no_api()
  trait_id <- 5020
  result <- trait(trait_id, include_full_associations = TRUE)
  expect_type(result, "list")
  expect_true(is.data.frame(result$full_associations))
  expect_true(nrow(result$full_associations) > 0)
  expected_names <- c("variant_id", "study_id", "beta", "se", "p", "eaf", "imputed")
  expect_true(all(expected_names %in% names(result$full_associations)))
})

test_that("variant() returns expected output", {
  skip_if_no_api()
  variant_id <- 5553693
  result <- variant(variant_id)
  expect_type(result, "list")
  expect_true(result$variant$id == variant_id)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("gene() returns expected output", {
  skip_if_no_api()
  gene_id <- "WNT7B"
  result <- gene(gene_id)
  expect_type(result, "list")
  expect_true(result$gene$gene == gene_id)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("associations() returns expected output", {
  skip_if_no_api()
  variant_ids <- c(80717)
  study_ids <- c(926, 927)
  result <- associations(variant_ids, study_ids)
  expect_type(result, "list")
  expect_gt(nrow(result), 0)
  expect_true(all(result$variant_id %in% variant_ids))
  expect_true(all(result$study_id %in% study_ids))
})

test_that("get_all_gene_pleiotropies() returns expected output", {
  skip_if_no_api()
  result <- get_all_gene_pleiotropies()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("get_all_variant_pleiotropies() returns expected output", {
  skip_if_no_api()
  result <- get_all_variant_pleiotropies()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("all_genes() returns expected output", {
  skip_if_no_api()
  result <- all_genes()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("all_traits() returns expected output", {
  skip_if_no_api()
  result <- all_traits()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("traits(trait_ids) returns expected output", {
  skip_if_no_api()
  trait_ids <- c(4405, 4872)
  result <- traits(trait_ids = trait_ids, include_associations = TRUE)
  expect_type(result, "list")
  expect_true(all(result$traits$id %in% trait_ids))
  expect_true(nrow(result$coloc_groups) > 0)
  expect_true(!is.null(result$study_extractions))
  expect_true(nrow(result$study_extractions) > 0)
})

test_that("genes(gene_ids) returns expected output", {
  skip_if_no_api()
  gene_ids <- c("WNT7B", "WNT7A")
  result <- genes(gene_ids = gene_ids, include_associations = TRUE)
  expect_type(result, "list")
  expect_true(all(result$genes$gene %in% gene_ids))
  expect_true(nrow(result$coloc_groups) > 0)
  expect_true(!is.null(result$study_extractions))
  expect_true(nrow(result$study_extractions) > 0)
})

test_that("variants() returns expected output", {
  skip_if_no_api()
  variant_ids <- c(5553693, 5553694)
  result <- variants(variants = variant_ids, include_associations = TRUE, expand = TRUE)
  expect_type(result, "list")
  expect_true(all(result$variants$variant_id %in% variant_ids))
  expect_true(nrow(result$coloc_groups) > 0)
  expect_true(!is.null(result$study_extractions))
  expect_true(nrow(result$study_extractions) > 0)
})

test_that("search_gpmap() orders results by coloc groups and rare results", {
  skip_if_no_api()
  result <- search_gpmap("haemoglobin")
  importance <- result$num_coloc_groups + result$num_rare_results
  expect_false(is.unsorted(rev(importance)))
})

test_that("trait_duplicates() returns expected output", {
  skip_if_no_api()
  result <- trait_duplicates()
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  # Older API versions name each duplicate's parent; newer ones list duplicate
  # groups with a representative. remove_duplicate_traits() reads both.
  parent_shape <- c("trait_id", "trait_name", "parent_trait_id", "parent_trait_name")
  group_shape <- c("trait_id", "trait_name", "duplicate_group", "is_representative")
  expect_true(all(parent_shape %in% names(result)) || all(group_shape %in% names(result)))
})

test_that("delete_gwas() surfaces API errors", {
  skip_if_no_api()
  expect_error(
    delete_gwas("00000000-0000-0000-0000-000000000000", "user@example.com"),
    "not found"
  )
})

test_that("region() with include_associations works against the live API", {
  skip_if_no_api()
  result <- region(459, include_associations = TRUE)
  expect_gt(nrow(result$coloc_groups), 0)
})

test_that("search_gpmap() returns a well-formed table for a variant", {
  skip_if_no_api()
  result <- search_gpmap("rs9939609")
  expect_equal(result$type[1], "original_variant")
  expect_true(all(c("type", "name", "type_id", "rsq", "num_coloc_studies", "call") %in% names(result)))
})

test_that("trait tables and records use pipeline_id and variable_type", {
  skip_if_no_api()
  traits_table <- all_traits()
  expect_true(all(c("pipeline_id", "variable_type") %in% names(traits_table)))
  result <- trait(5020)
  expect_true(is.character(result$trait$pipeline_id))
  expect_match(result$trait$source_url, "^https://opengwas.io/datasets/")
})

test_that("include_associations folds associations into the tables", {
  skip_if_no_api()
  result <- gene("WNT7B", include_associations = TRUE)
  expect_null(result$associations)
  expect_true(all(c("beta", "se", "p") %in% names(result$coloc_groups)))
  expect_true(all(c("beta", "se", "p") %in% names(result$study_extractions)))
})

test_that("ld_matrix() returns stored pairs, and an empty table when there are none", {
  skip_if_no_api()
  pair <- ld_matrix(c(5553693, 5553688))$lds
  expect_gt(nrow(pair), 0)
  expect_true(all(c("r", "rsq") %in% names(pair)))
  expect_equal(nrow(ld_matrix(c(5553693, 80717))$lds), 0)
})

test_that("get_gwas() returns an upload in the same shape as other results", {
  skip_if_no_api()
  guid <- "7a289615-c1b4-91f3-3d97-887f60de9155"
  result <- tryCatch(get_gwas(guid, include_associations = TRUE), gpmap_api_error = function(e) NULL)
  skip_if(is.null(result), "The example GWAS upload is not on this API")
  expect_null(result$associations)
  expect_null(result$upload_study_extractions)
  groups <- result$coloc_groups
  expect_true(all(c("study_extraction_id", "source", "gwas_upload_id", "beta") %in% names(groups)))
  expect_true(all(groups$source %in% c("upload", "existing")))
  expect_true(all(is.na(groups$gwas_upload_id[groups$source == "existing"])))
  expect_true(any(groups$source == "upload" & groups$gwas_upload_id == result$trait$id))
  extraction_keys <- paste(result$study_extractions$source, result$study_extractions$id)
  expect_true(all(paste(groups$source, groups$study_extraction_id) %in% extraction_keys))
  expect_true(all(c("study_extraction_a_id", "source_a", "trait_name_a") %in% names(result$coloc_pairs)))

  comparison <- compare_upload(result, 2527)
  expect_true(all(comparison$coloc_groups$coloc_group_id %in% groups$coloc_group_id))
})
