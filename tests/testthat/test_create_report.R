library(testthat)

test_that("create_report has no method for other objects", {
  expect_error(create_report(1), "gpmap_correspondence, gpmap_decomposition")
})

test_that("create_report renders decomposition and correspondence reports", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
  wd <- getwd()
  before <- list.files(wd)
  out_dir <- tempfile("report-test-")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  d <- simulated_decompositions()[[1]]
  f <- create_report(d, output_file = file.path(out_dir, "decomposition.html"), open = FALSE)
  expect_true(file.exists(f))
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  for (heading in c("Traits removed", "Programs and their confidence", "What the programs are about")) {
    expect_match(html, heading, fixed = TRUE)
  }

  x <- simulated_correspondence()
  f <- create_report(x, output_file = file.path(out_dir, "correspondence.html"), open = FALSE)
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  for (heading in c("How many program pairs link", "Families", "Do families recover the same biology?")) {
    expect_match(html, heading, fixed = TRUE)
  }
  expect_false(grepl("Trait: Simulated trait 9001", html, fixed = TRUE))

  f <- create_report(x, output_file = file.path(out_dir, "with-traits.html"), open = FALSE, include_traits = TRUE)
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, "Trait: Simulated trait 9001", fixed = TRUE)
  expect_match(html, "Trait: Simulated trait 9002", fixed = TRUE)

  expect_identical(getwd(), wd)
  expect_setequal(list.files(wd), before)
})

test_that("create_report writes to a temporary file by default", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
  f <- create_report(simulated_decompositions()[[2]], open = FALSE)
  expect_true(file.exists(f))
  expect_true(startsWith(normalizePath(f), normalizePath(tempdir())))
  unlink(f)
})
