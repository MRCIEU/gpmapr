library(testthat)

coloc_groups_fixture <- function() {
  return(data.frame(
    coloc_group_id = c(1L, 1L, 2L),
    study_id = c(10L, 11L, 12L),
    variant_id = c(100L, 100L, 200L),
    stringsAsFactors = FALSE
  ))
}

rare_results_fixture <- function() {
  return(data.frame(
    rare_result_group_id = 1L,
    study_id = 13L,
    variant_id = 300L,
    stringsAsFactors = FALSE
  ))
}

test_that("merge_associations() returns its inputs unchanged when there are no associations", {
  result <- merge_associations(coloc_groups_fixture(), rare_results_fixture(), NULL)
  expect_equal(
    result,
    list(coloc_groups = coloc_groups_fixture(), rare_results = rare_results_fixture(), study_extractions = NULL)
  )
})

test_that("merge_associations() joins associations on study and variant", {
  associations <- data.frame(
    study_id = c(10L, 12L, 13L),
    variant_id = c(100L, 200L, 300L),
    beta = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )
  study_extractions <- data.frame(id = c(1L, 2L), study_id = c(10L, 14L), variant_id = c(100L, 400L))
  result <- merge_associations(coloc_groups_fixture(), rare_results_fixture(), associations, study_extractions)
  expect_equal(result$coloc_groups$beta, c(0.1, NA, 0.2))
  expect_equal(result$rare_results$beta, 0.3)
  expect_equal(result$study_extractions$beta, c(0.1, NA))
})

test_that("functions fold the associations into their tables and drop them", {
  local_mocked_bindings(gene_api = function(...) {
    return(list(
      coloc_groups = coloc_groups_fixture(),
      study_extractions = data.frame(id = 1L, study_id = 10L, variant_id = 100L),
      associations = data.frame(study_id = 10L, variant_id = 100L, beta = 0.1, se = 0.01, p = 1e-9)
    ))
  })
  result <- gene("WNT7B", include_associations = TRUE)
  expect_null(result$associations)
  expect_equal(result$coloc_groups$beta, c(0.1, NA, NA))
  expect_equal(result$study_extractions$se, 0.01)
})

test_that("GWAS upload associations are folded in by source", {
  local_mock_api()
  result <- get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155", include_associations = TRUE)
  expect_null(result$associations)
  groups <- result$coloc_groups
  # Existing rows join on study_id, upload rows on gwas_upload_id; the same id in each must not mix.
  expect_equal(groups$beta[groups$source == "upload" & groups$gwas_upload_id == 90], c(0.5, 0.3))
  expect_equal(groups$beta[groups$source == "upload" & groups$gwas_upload_id == 144], 0.45)
  expect_equal(groups$beta[groups$source == "existing" & groups$trait_id == 2527], 0.2)
  extractions <- result$study_extractions
  expect_equal(extractions$beta[extractions$source == "upload" & extractions$id == 7223], 0.5)
  expect_true(is.na(extractions$beta[extractions$source == "existing" & extractions$id == 7223]))
})

test_that("placeholder associations are returned as missing", {
  associations <- data.frame(
    variant_id = 1:3, study_id = 1, beta = c(0, 0, 0.2), se = c(1, 0.5, 1), p = c(1, 0.9, 1e-3),
    eaf = 0.3, imputed = c(TRUE, TRUE, FALSE)
  )
  masked <- mask_padded_associations(associations)
  expect_equal(masked$beta, c(NA, 0, 0.2))
  expect_equal(masked$eaf, c(NA, 0.3, 0.3))
  expect_equal(masked$imputed, associations$imputed)
  expect_equal(nrow(masked), 3)

  local_mock_api()
  groups <- get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155", include_associations = TRUE)$coloc_groups
  expect_true(is.na(groups$beta[groups$source == "existing" & groups$study_extraction_id == 7223]))
})

# Each of these reaches merge_associations() even when the API sends no
# associations, which the live region endpoint never does.
api_object_without_associations <- function(...) {
  return(list(coloc_groups = coloc_groups_fixture(), rare_results = rare_results_fixture()))
}

test_that("region() with include_associations works when the API returns no associations", {
  local_mocked_bindings(region_api = api_object_without_associations)
  result <- region(459, include_associations = TRUE)
  expect_equal(result$coloc_groups, coloc_groups_fixture())
  expect_equal(result$rare_results, rare_results_fixture())
})

test_that("gene() with include_associations works when the API returns no associations", {
  local_mocked_bindings(gene_api = api_object_without_associations)
  result <- gene("ENSG00000140718", include_associations = TRUE)
  expect_equal(result$coloc_groups, coloc_groups_fixture())
})

test_that("trait() with include_associations works when the API returns no associations", {
  local_mocked_bindings(trait_api = api_object_without_associations)
  result <- trait(1992L, include_associations = TRUE)
  expect_equal(result$coloc_groups, coloc_groups_fixture())
})

test_that("variant() works when the API returns no associations", {
  local_mocked_bindings(variant_api = api_object_without_associations)
  result <- variant(100L)
  expect_equal(result$coloc_groups, coloc_groups_fixture())
})
