library(testthat)

test_that("health() returns expected output", {
  result <- health_api()
  expect_type(result, "list")
  expect_true(result$status %in% c("healthy", "unhealthy"))
})

test_that("search_gpmap() returns expected output", {
  result <- search_gpmap("hemoglobin")
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
  expect_true(all(names(result) == c("call", "name", "type", "type_id", "info")))
})

test_that("trait() returns expected output", {
  result <- trait(1)
  expect_type(result, "list")
  expect_true(result$trait$id == 1)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("variant() returns expected output", {
  result <- variant(8466253)
  expect_type(result, "list")
  expect_true(result$variant$id == 8466253)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("gene() returns expected output", {
  result <- gene("SLC44A1")
  expect_type(result, "list")
  expect_true(result$gene$gene == "SLC44A1")
  expect_true(nrow(result$coloc_groups) > 0)
})
