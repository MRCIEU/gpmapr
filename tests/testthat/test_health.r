library(testthat)

test_that("health() returns expected output", {
  # Example: Assuming health() returns a status message
  result <- health()
  expect_type(result, "list")
  expect_true(result$status %in% c("healthy", "unhealthy"))
})

test_that("search_gpmap() returns expected output", {
  result <- search_gpmap("hemo")
  expect_type(result, "list")
  expect_true(result$gpmap$id == 1)
})

test_that("trait() returns expected output", {
  result <- trait(1)
  expect_type(result, "list")
  expect_true(result$trait$id == 1)
})

test_that("variant() returns expected output", {
  result <- variant(1)
  expect_type(result, "list")
  expect_true(result$variant$id == 1)
})

test_that("study_extractions() returns expected output", {
  result <- study_extractions(1)
  expect_type(result, "list")
  expect_true(result$study_extractions$id == 1)
})

test_that("gene() returns expected output", {
  result <- gene(1)
  expect_type(result, "list")
  expect_true(result$gene$id == 1)
})
