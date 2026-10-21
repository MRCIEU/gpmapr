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
  expect_equal(result, list(coloc_groups = coloc_groups_fixture(), rare_results = rare_results_fixture()))
})

test_that("merge_associations() joins associations on study and variant", {
  associations <- data.frame(
    study_id = c(10L, 12L, 13L),
    variant_id = c(100L, 200L, 300L),
    beta = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )
  result <- merge_associations(coloc_groups_fixture(), rare_results_fixture(), associations)
  expect_equal(result$coloc_groups$beta, c(0.1, NA, 0.2))
  expect_equal(result$rare_results$beta, 0.3)
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
