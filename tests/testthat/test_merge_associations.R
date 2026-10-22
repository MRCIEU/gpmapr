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

test_that("GWAS upload associations are folded into the upload's tables", {
  api_object <- list(
    trait = list(id = 90L, name = "My GWAS"),
    coloc_groups = data.frame(
      study_id = c(90L, NA), existing_study_id = c(NA, 21L), variant_id = c(100L, 100L),
      existing_study_extraction_id = c(NA, 5L), study_extraction_id = c(7L, NA)
    ),
    study_extractions = data.frame(id = 5L, study_id = 21L, variant_id = 100L),
    upload_study_extractions = data.frame(id = 7L, gwas_upload_id = 90L, variant_id = 100L),
    associations = data.frame(
      gwas_upload_id = 90L, variant_id = 100L, study_id = c(90L, NA), existing_study_id = c(NA, 21L),
      beta = c(0.5, 0.2), se = 0.1, p = 1e-9, eaf = 0.3, imputed = FALSE
    )
  )
  local_mocked_bindings(get_gwas_api = function(...) api_object)
  result <- get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155", include_associations = TRUE)
  expect_null(result$associations)
  expect_equal(result$coloc_groups$beta[result$coloc_groups$source == "upload"], 0.5)
  expect_equal(result$coloc_groups$beta[result$coloc_groups$source == "existing"], 0.2)
  expect_equal(result$study_extractions$beta, 0.2)
  expect_equal(result$upload_study_extractions$beta, 0.5)
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
