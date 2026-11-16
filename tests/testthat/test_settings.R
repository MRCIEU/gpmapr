library(testthat)

# _ebmf_settings.yml lives in vignettes/, which R CMD check does not keep.
skip_if_no_settings_file <- function() {
  found <- tryCatch(nzchar(ebmf_settings_path()), error = function(e) return(FALSE))
  return(skip_if_not(found, "vignettes/_ebmf_settings.yml is not available"))
}

test_that("gpmap_settings defaults match the shared settings file", {
  skip_if_no_settings_file()
  s <- gpmap_settings()
  shared <- read_ebmf_settings()
  for (key in names(shared)) {
    expect_equal(s[[key]], shared[[key]], ignore_attr = TRUE, label = paste("gpmap_settings()$", key))
  }
})

test_that("gpmap_settings defaults match the function formals", {
  s <- gpmap_settings()
  expect_s3_class(s, "gpmap_settings")
  expect_equal(s$n_rep, formals(summarise_ebmf_programs)$n_rep)
  expect_equal(s$seed, formals(summarise_ebmf_programs)$seed)
  link <- formals(compare_program_pairs_loadings)
  for (key in c("n_perm", "phi_high", "phi_low", "min_shared_loci", "min_eff_loci",
                "min_alignment", "fdr_threshold")) {
    expect_equal(s[[key]], eval(link[[key]]), label = paste("gpmap_settings()$", key))
  }
  expect_false(s$mutual_best)
})

test_that("gpmap_settings takes overrides as arguments or lists and rejects unknown keys", {
  s <- gpmap_settings(list(n_perm = 3000, phi_high = 0.9), phi_high = 0.95, seed = NULL)
  expect_equal(s$n_perm, 3000)
  expect_equal(s$phi_high, 0.95)
  expect_equal(s$seed, 1)
  expect_equal(gpmap_settings(list(ebmf_prior = "point_normal"))$ebmf_prior, "point_normal")
  expect_error(gpmap_settings(n_null = 5), "Unknown setting")
  expect_error(gpmap_settings(n_perm = "a lot"), "numeric")
  expect_error(gpmap_settings(mutual_best = 1), "TRUE or FALSE")
  expect_error(gpmap_settings(phi_high = c(0.8, 0.9)), "single")
  expect_output(print(gpmap_settings()), "Cross-trait links")
})

test_that("gpmap_settings accepts a resolved settings file", {
  skip_if_no_settings_file()
  expect_equal(gpmap_settings(resolve_ebmf_settings())$ebmf_prior, "point_laplace")
})

test_that("settings split into the pipeline and link parts", {
  s <- gpmap_settings(similarity_threshold = 0.3, phi_low = 0.6)
  pipeline <- .settings_for(s, "pipeline")
  expect_equal(pipeline$similarity_threshold, 0.3)
  expect_false(any(c("n_perm", "phi_low", "n_rep", "seed") %in% names(pipeline)))
  expect_silent(.ebmf_pipeline_args(pipeline))
  expect_equal(.settings_for(s, "link")$phi_low, 0.6)
  expect_setequal(names(.settings_for(s, "link")), setdiff(names(formals(assign_link_tiers)), "out"))
})

test_that(".default_cores is a positive integer", {
  n <- .default_cores()
  expect_true(is.integer(n) && length(n) == 1 && n >= 1L)
  expect_equal(.default_cores(NA_integer_), 1L)
  expect_equal(.default_cores(1L), 1L)
  expect_equal(.default_cores(10L), 5L)
})
