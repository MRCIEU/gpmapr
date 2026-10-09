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
test_that("search_gpmap() orders results by coloc groups and rare results", {
  result <- search_gpmap("haemoglobin")
  importance <- result$num_coloc_groups + result$num_rare_results
  expect_false(is.unsorted(rev(importance)))
})

test_that("trait_duplicates() returns expected output", {
  result <- trait_duplicates()
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  expected_names <- c("trait_id", "trait_name", "parent_trait_id", "parent_trait_name")
  expect_true(all(expected_names %in% names(result)))
})

test_that("delete_gwas() validates its inputs", {
  guid <- "00000000-0000-0000-0000-000000000000"
  expect_error(delete_gwas("not-a-guid", "user@example.com"), "GUID")
  expect_error(delete_gwas(guid), "email is required")
  expect_error(delete_gwas(guid, ""), "email is required")
})

test_that("delete_gwas() surfaces API errors", {
  expect_error(
    delete_gwas("00000000-0000-0000-0000-000000000000", "user@example.com"),
    "not found"
  )
})

test_that("upload_gwas() validates column arguments", {
  file <- tempfile(fileext = ".tsv")
  writeLines(c("chr\tpos\tea\toa\tp\tbeta\tse", "1\t100\tA\tG\t1e-9\t0.1\t0.01"), file)
  on.exit(unlink(file))
  upload <- function(...) {
    upload_gwas(file, name = "test", email = "user@example.com", sample_size = 1000, ...)
  }

  expect_error(
    upload(chr_col = "chr", bp_col = "pos", ea_col = "ea", oa_col = "oa", beta_col = "beta", se_col = "se"),
    "Missing required column arguments: p_col"
  )
  expect_error(
    upload(chr_col = "chr", bp_col = "pos", ea_col = "ea", oa_col = "oa", p_col = "p", beta_col = "beta"),
    "Either beta_col and se_col"
  )
  expect_error(
    upload(chr_col = "chr", bp_col = "pos", ea_col = "ea", oa_col = "oa", p_col = "pval",
      beta_col = "beta", se_col = "se"),
    "Columns not found in file: pval"
  )
  expect_error(
    upload(column_names = list(CHR = "chr"), chr_col = "chr"),
    "not both"
  )
  expect_error(
    suppressWarnings(upload(column_names = list(CHR = "chr", POS = "pos"))),
    "Unknown column_names: POS"
  )
  expect_warning(
    expect_error(upload(column_names = list(CHR = "chr", BP = "pos")), "Missing required"),
    "deprecated"
  )
})

test_that("legacy column_names are mapped to API column names", {
  result <- suppressWarnings(build_gwas_column_names(
    list(chr = "chr", BP = "pos", EA = "ea", OA = "oa", P = "p", OR = "or", LB = "lb", UB = "ub"),
    list()
  ))
  expect_equal(names(result), c("CHR", "BP", "EA", "OA", "P", "OR", "OR_LB", "OR_UB"))
  expect_equal(result$OR_LB, "lb")
})
