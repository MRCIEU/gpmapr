library(testthat)

# Offline tests: every request is answered from the recorded fixtures in
# fixtures/ (see helper-mock-api.R), so these run without network access,
# including on CRAN. test_live_api.R checks the same functions against the
# real API.

test_that("health_api() parses the health response", {
  local_mock_api()
  result <- health_api()
  expect_type(result, "list")
  expect_true(result$status %in% c("healthy", "unhealthy"))
})

test_that("all_traits() returns a trait table", {
  local_mock_api()
  result <- all_traits()
  expect_s3_class(result, "data.frame")
  expect_true(all(c("id", "trait_name", "trait_category") %in% names(result)))
  expect_gt(nrow(result), 0)
})

test_that("trait_duplicates() returns the duplicate groups", {
  local_mock_api()
  result <- trait_duplicates()
  expect_s3_class(result, "data.frame")
  expect_true(all(c("trait_id", "duplicate_group", "is_representative") %in% names(result)))
})

test_that("trait() returns the trait and its coloc groups", {
  local_mock_api()
  result <- trait(5020)
  expect_equal(result$trait$id, 5020)
  expect_s3_class(result$coloc_groups, "data.frame")
  expect_gt(nrow(result$coloc_groups), 0)
  expect_false("beta" %in% names(result$coloc_groups))
})

test_that("trait(include_associations = TRUE) joins the associations onto the coloc groups", {
  local_mock_api()
  result <- trait(5020, include_associations = TRUE)
  expect_true(all(c("beta", "se") %in% names(result$coloc_groups)))
})

test_that("trait(include_full_associations = TRUE) adds the full associations", {
  local_mock_api()
  result <- trait(5020, include_full_associations = TRUE)
  expect_s3_class(result$full_associations, "data.frame")
  expect_true(all(c("variant_id", "study_id", "beta", "se", "p") %in% names(result$full_associations)))
  expect_equal(nrow(result$full_associations), 3)
  expect_type(result$full_associations$beta, "double")
  expect_type(result$full_associations$imputed, "logical")
})

test_that("trait(include_coloc_pairs = TRUE) adds the coloc pairs from the columnar response", {
  local_mock_api()
  result <- trait(5020, include_coloc_pairs = TRUE)
  expect_s3_class(result$coloc_pairs, "data.frame")
  expect_equal(
    names(result$coloc_pairs),
    c("variant_id", "study_extraction_a_id", "study_extraction_b_id", "ld_block_id", "h3", "h4",
      "false_positive", "false_negative", "in_coloc_group")
  )
  expect_equal(nrow(result$coloc_pairs), 3)
  expect_type(result$coloc_pairs$h4, "double")
  expect_type(result$coloc_pairs$in_coloc_group, "logical")
})

test_that("tabular_response_to_dataframe() reads columnar and row responses alike", {
  columnar <- jsonlite::fromJSON(
    '{"coloc_pair_column_names": ["variant_id", "h4", "in_coloc_group"],
      "coloc_pair_columns": {"in_coloc_group": [true, false], "variant_id": [1, 2], "h4": [0.9, null]}}'
  )
  rows <- jsonlite::fromJSON(
    '{"coloc_pair_column_names": ["variant_id", "h4", "in_coloc_group"],
      "coloc_pair_rows": [[1, 0.9, true], [2, null, false]]}'
  )
  from_columns <- tabular_response_to_dataframe(columnar, "coloc_pair")
  from_rows <- tabular_response_to_dataframe(rows, "coloc_pair")

  expect_equal(names(from_columns), c("variant_id", "h4", "in_coloc_group"))
  expect_equal(from_columns$variant_id, c(1, 2))
  expect_equal(from_columns$h4, c(0.9, NA))
  expect_equal(from_columns$in_coloc_group, c(TRUE, FALSE))
  expect_equal(names(from_rows), names(from_columns))
  expect_equal(from_rows$variant_id, from_columns$variant_id)
  expect_equal(from_rows$h4, from_columns$h4)
})

test_that("tabular_response_to_dataframe() returns an empty data frame for an empty response", {
  columnar <- jsonlite::fromJSON('{"associations_full_column_names": [], "associations_full_columns": {}}')
  rows <- jsonlite::fromJSON('{"associations_full_column_names": [], "associations_full_rows": []}')
  expect_equal(dim(tabular_response_to_dataframe(columnar, "associations_full")), c(0, 0))
  expect_equal(dim(tabular_response_to_dataframe(rows, "associations_full")), c(0, 0))

  no_rows <- jsonlite::fromJSON(
    '{"associations_full_column_names": ["variant_id", "p"], "associations_full_columns": {"variant_id": [], "p": []}}'
  )
  result <- tabular_response_to_dataframe(no_rows, "associations_full")
  expect_equal(names(result), c("variant_id", "p"))
  expect_equal(nrow(result), 0)
})

test_that("traits() returns the requested traits with their coloc groups", {
  local_mock_api()
  result <- traits(c(4405, 4872), include_associations = TRUE)
  expect_true(all(result$traits$id %in% c(4405, 4872)))
  expect_gt(nrow(result$coloc_groups), 0)
  expect_gt(nrow(result$study_extractions), 0)
})

test_that("all_genes() returns a gene table", {
  local_mock_api()
  result <- all_genes()
  expect_s3_class(result, "data.frame")
  expect_true(all(c("id", "gene", "chr") %in% names(result)))
})

test_that("gene() returns the gene and its coloc groups", {
  local_mock_api()
  result <- gene("WNT7B")
  expect_equal(result$gene$gene, "WNT7B")
  expect_gt(nrow(result$coloc_groups), 0)
  expect_null(result$tissues)
})

test_that("genes() returns the requested genes", {
  local_mock_api()
  result <- genes(c("WNT7B", "WNT7A"), include_associations = TRUE)
  expect_true(all(result$genes$gene %in% c("WNT7B", "WNT7A")))
  expect_gt(nrow(result$coloc_groups), 0)
})

test_that("variant() returns the variant and its coloc groups", {
  local_mock_api()
  result <- variant(5553693)
  expect_equal(result$variant$id, 5553693)
  expect_gt(nrow(result$coloc_groups), 0)
})

test_that("variants() returns the requested variants", {
  local_mock_api()
  result <- variants(c(5553693, 5553694), include_associations = TRUE, expand = TRUE)
  expect_true(all(result$variants$variant_id %in% c(5553693, 5553694)))
  expect_gt(nrow(result$coloc_groups), 0)
})

test_that("region() returns the region with and without associations", {
  local_mock_api()
  plain <- region(459)
  expect_gt(nrow(plain$coloc_groups), 0)
  # The region endpoint sends no associations, which used to make this error.
  with_associations <- region(459, include_associations = TRUE)
  expect_gt(nrow(with_associations$coloc_groups), 0)
  expect_equal(names(with_associations$coloc_groups), names(plain$coloc_groups))
})

test_that("associations() returns the requested variant-study associations", {
  local_mock_api()
  result <- associations(80717, c(926, 927))
  expect_s3_class(result, "data.frame")
  expect_true(all(result$variant_id == 80717))
  expect_true(all(result$study_id %in% c(926, 927)))
  expect_true(all(c("beta", "se", "p") %in% names(result)))
})

test_that("ld_proxies() returns proxies of the requested variants", {
  local_mock_api()
  result <- ld_proxies(c(5553693, 5553694))
  expect_true(all(c("lead_variant_id", "proxy_variant_id", "r") %in% names(result$lds)))
})

test_that("pleiotropy summaries are returned for genes and variants", {
  local_mock_api()
  genes_result <- get_all_gene_pleiotropies()
  variants_result <- get_all_variant_pleiotropies()
  expect_true(all(c("gene_id", "distinct_trait_categories") %in% names(genes_result)))
  expect_true(all(c("variant_id", "distinct_trait_categories") %in% names(variants_result)))
})

test_that("all_traits() and trait() return pipeline_id and variable_type", {
  local_mock_api()
  expect_true(all(c("pipeline_id", "variable_type") %in% names(all_traits())))
  result <- trait(5020)
  expect_true(is.character(result$trait$pipeline_id))
  expect_equal(result$trait$source_url, create_source_url(result$trait$pipeline_id))
})

test_that("get_gwas() returns an upload in the same shape as other results", {
  local_mock_api()
  result <- get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155")
  expect_null(result$upload_study_extractions)
  expect_true(all(c("id", "source", "gwas_upload_id", "trait_name") %in% names(result$study_extractions)))
  groups <- result$coloc_groups
  extraction_keys <- paste(result$study_extractions$source, result$study_extractions$id)
  expect_true(all(paste(groups$source, groups$study_extraction_id) %in% extraction_keys))
  expect_equal(trait("7a289615-c1b4-91f3-3d97-887f60de9155")$coloc_groups, groups)
})

test_that("traits() combines trait ids and GWAS upload GUIDs", {
  local_mock_api()
  result <- traits(c(4405, 4872, "7a289615-c1b4-91f3-3d97-887f60de9155"), include_associations = TRUE)
  expect_true(all(c(4405, 4872, 90) %in% result$traits$id))
  expect_equal(result$traits$trait_name[result$traits$id == 90], "AD GWAS")
  expect_true(all(c("existing", "upload") %in% result$coloc_groups$source))
})

test_that("get_gwas() returns an uploaded GWAS", {
  local_mock_api()
  result <- get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155")
  expect_type(result$trait, "list")
  expect_s3_class(result$coloc_groups, "data.frame")
})

test_that("search_gpmap() finds traits by name, most colocalised first", {
  local_mock_api()
  rm(list = ls(.search_cache), envir = .search_cache)
  result <- search_gpmap("haemoglobin")
  rm(list = ls(.search_cache), envir = .search_cache)
  expect_gt(nrow(result), 0)
  expect_true(all(grepl("haemoglobin", tolower(result$name))))
  expect_true(all(c("type", "name", "type_id", "call", "num_coloc_groups") %in% names(result)))
  importance <- result$num_coloc_groups + result$num_rare_results
  expect_false(is.unsorted(rev(importance)))
})

test_that("search_gpmap() finds a variant and its proxies", {
  local_mock_api()
  result <- search_gpmap("rs9939609")
  expect_equal(result$type[1], "original_variant")
  expect_true(all(c("rsq", "num_coloc_studies", "call") %in% names(result)))
  expect_true(all(grepl("^variant\\(", result$call)))
})

test_that("an HTTP error from the API gives an informative error", {
  response <- fake_api_response('{"detail": "Trait not found"}', status = 404L)
  expect_error(stop_for_api_status(response), "HTTP 404.*Trait not found")
  expect_invisible(stop_for_api_status(fake_api_response("{}")))
})

test_that("api_get() reports an HTTP error instead of returning it", {
  local_mocked_bindings(ensure_api_compatible = function() invisible(TRUE))
  local_mocked_bindings(GET = function(...) fake_api_response('{"detail": "boom"}', status = 500L), .package = "httr")
  expect_error(api_get("https://gpmap.opengwas.io/api/v1/traits"), "HTTP 500.*boom")
})

test_that("an unreachable API gives an informative error", {
  unreachable <- function(...) stop("Could not resolve host: gpmap.invalid")
  expect_error(request_with_backoff(unreachable), "Could not reach the GPMap API.*Could not resolve host")
})
