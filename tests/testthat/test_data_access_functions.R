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
  trait_id <- 5020
  result <- trait(trait_id)
  expect_type(result, "list")
  expect_true(result$trait$id == trait_id)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("variant() returns expected output", {
  variant_id <- 8466253
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
