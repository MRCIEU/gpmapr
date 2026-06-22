library(testthat)

test_that("health() returns expected output", {
  result <- health_api()
  expect_type(result, "list")
  expect_true(result$status %in% c("healthy", "unhealthy"))
})

test_that("search_gpmap() returns expected output", {
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
  trait_id <- 5020
  result <- trait(trait_id)
  expect_type(result, "list")
  expect_true(result$trait$id == trait_id)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("trait() returns full_associations when requested", {
  trait_id <- 5020
  result <- trait(trait_id, include_full_associations = TRUE)
  expect_type(result, "list")
  expect_true(is.data.frame(result$full_associations))
  expect_true(nrow(result$full_associations) > 0)
  expected_names <- c("variant_id", "study_id", "beta", "se", "p", "eaf", "imputed")
  expect_true(all(expected_names %in% names(result$full_associations)))
})

test_that("variant() returns expected output", {
  variant_id <- 5553693
  result <- variant(variant_id)
  expect_type(result, "list")
  expect_true(result$variant$id == variant_id)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("gene() returns expected output", {
  gene_id <- "WNT7B"
  result <- gene(gene_id)
  expect_type(result, "list")
  expect_true(result$gene$gene == gene_id)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("associations() returns expected output", {
  variant_ids <- c(80750)
  study_ids <- c(5020, 4870)
  result <- associations(variant_ids, study_ids)
  expect_type(result, "list")
  expect_true(all(result$variant_id %in% variant_ids))
  expect_true(all(result$study_id %in% study_ids))
})

test_that("get_all_gene_pleiotropies() returns expected output", {
  result <- get_all_gene_pleiotropies()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("get_all_variant_pleiotropies() returns expected output", {
  result <- get_all_variant_pleiotropies()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("all_genes() returns expected output", {
  result <- all_genes()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("all_traits() returns expected output", {
  result <- all_traits()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("traits(trait_ids) returns expected output", {
  trait_ids <- c(4405, 4872)
  result <- traits(trait_ids = trait_ids, include_associations = TRUE)
  expect_type(result, "list")
  expect_true(all(result$traits$id %in% trait_ids))
  expect_true(nrow(result$coloc_groups) > 0)
  expect_true(!is.null(result$study_extractions))
  expect_true(nrow(result$study_extractions) > 0)
})

test_that("genes(gene_ids) returns expected output", {
  gene_ids <- c("WNT7B", "WNT7A")
  result <- genes(gene_ids = gene_ids, include_associations = TRUE)
  expect_type(result, "list")
  expect_true(all(result$genes$gene %in% gene_ids))
  expect_true(nrow(result$coloc_groups) > 0)
  expect_true(!is.null(result$study_extractions))
  expect_true(nrow(result$study_extractions) > 0)
})

test_that("variants() returns expected output", {
  variant_ids <- c(5553693, 5553694)
  result <- variants(variants = variant_ids, include_associations = TRUE, expand = TRUE)
  expect_type(result, "list")
  expect_true(all(result$variants$variant_id %in% variant_ids))
  expect_true(nrow(result$coloc_groups) > 0)
  expect_true(!is.null(result$study_extractions))
  expect_true(nrow(result$study_extractions) > 0)
})

test_that("pathway_enrichment() returns expected output with gene IDs", {
  genes <- all_genes()
  gene_ids <- genes[genes$gene %in% c("TREM2", "APOE"), "id"]

  result <- pathway_enrichment(gene_ids, minimum_count_in_network = 1)
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_true("input_gene_count" %in% names(result))
  expect_true("matched_gene_count" %in% names(result))
  expect_true("p_value_threshold" %in% names(result))
  expect_true("total_terms_tested" %in% names(result))
  expect_equal(result$input_gene_count, length(gene_ids))
  expect_true(is.data.frame(result$results) || is.null(result$results))
  if (is.data.frame(result$results) && nrow(result$results) > 0) {
    expected_cols <- c(
      "term_id", "source", "description", "pathway_size",
      "background_size", "overlap", "p_value", "fdr", "gene_ids",
      "pathway_gene_ids"
    )
    expect_true(all(expected_cols %in% names(result$results)))
    expect_true(all(result$results$fdr <= 0.05))
    expect_true(is.list(result$results$gene_ids))
    expect_true(is.list(result$results$pathway_gene_ids))
    expect_true(length(result$results$pathway_gene_ids[[1]]) >=
                  length(result$results$gene_ids[[1]]))
  }
})

test_that("pathway_enrichment() accepts gene names", {
  result <- pathway_enrichment(
    c("TREM2", "APOE"),
    minimum_count_in_network = 1
  )
  expect_type(result, "list")
  expect_equal(result$input_gene_count, 2L)
})

test_that("pathway_enrichment() validates inputs", {
  expect_error(pathway_enrichment(NULL), "genes is required")
  expect_error(pathway_enrichment(c(1, NA)), "genes must not contain NA")
  expect_error(pathway_enrichment(c(TRUE)), "genes must be numeric gene IDs or character gene names")
  expect_error(pathway_enrichment(1, source = "Invalid"), "source must be one of")
  expect_error(pathway_enrichment(1, p_value_threshold = 2), "p_value_threshold must be")
  expect_error(pathway_enrichment(1, minimum_count_in_network = 0),
               "minimum_count_in_network must be a positive integer")
})