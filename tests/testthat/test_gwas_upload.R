library(testthat)

# Offline checks of the GWAS upload / delete helpers: argument validation only,
# so no API access is needed.

test_that("delete_gwas() validates its inputs", {
  guid <- "00000000-0000-0000-0000-000000000000"
  expect_error(delete_gwas("not-a-guid", "user@example.com"), "GUID")
  expect_error(delete_gwas(guid), "email is required")
  expect_error(delete_gwas(guid, ""), "email is required")
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
