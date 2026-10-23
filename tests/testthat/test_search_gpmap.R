library(testthat)

# search_gpmap() caches the search options in .search_cache for the session;
# clear it around each test so mocked options never leak into live tests.
with_clean_search_cache <- function(code) {
  clear <- function() rm(list = ls(.search_cache), envir = .search_cache)
  clear()
  on.exit(clear())
  return(code)
}

coloc_groups_of <- function(n) {
  return(data.frame(coloc_group_id = seq_len(n)))
}

# One searched variant and three proxies whose coloc groups, rare results and
# LD to the searched variant all differ, so any per-row mix-up shows.
mock_variant_search <- function() {
  original <- data.frame(id = 10L, display_snp = "rsQ", num_colocs = 5L)
  original$coloc_groups <- list(coloc_groups_of(5))
  original$rare_results <- list(NULL)
  proxies <- data.frame(id = c(11L, 12L, 13L), display_snp = c("rsA", "rsB", "rsC"), num_colocs = c(1L, 4L, 2L))
  proxies$coloc_groups <- list(coloc_groups_of(1), coloc_groups_of(4), coloc_groups_of(2))
  proxies$rare_results <- list(coloc_groups_of(0), coloc_groups_of(3), NULL)
  proxies$ld_proxies <- list(data.frame(r = 0.95), data.frame(r = -0.90), data.frame(r = 0.85))
  return(list(original_variants = original, proxy_variants = proxies))
}

test_that("search_gpmap() reports each proxy variant's own counts and r-squared", {
  local_mocked_bindings(
    search_variants_api = function(...) mock_variant_search(),
    search_options_api = function() stop("search options are not needed for a variant search")
  )
  result <- with_clean_search_cache(search_gpmap("rsq"))

  expect_equal(result$name[1], "rsQ")
  expect_equal(result$type[1], "original_variant")
  row <- function(name) result[result$name == name, ]
  expect_equal(row("rsQ")$num_coloc_studies, 5)
  expect_equal(
    vapply(c("rsA", "rsB", "rsC"), function(n) row(n)$num_coloc_studies, numeric(1)),
    c(rsA = 1, rsB = 4, rsC = 2)
  )
  expect_equal(
    vapply(c("rsA", "rsB", "rsC"), function(n) row(n)$num_rare_results, numeric(1)),
    c(rsA = 0, rsB = 3, rsC = 0)
  )
  expect_true(is.na(row("rsQ")$rsq))
  expect_equal(
    vapply(c("rsA", "rsB", "rsC"), function(n) row(n)$rsq, numeric(1)),
    c(rsA = 0.9, rsB = 0.81, rsC = 0.72)
  )
})

test_that("search_gpmap() does not download the search options for a variant search", {
  local_mocked_bindings(
    search_variants_api = function(...) mock_variant_search(),
    search_options_api = function() stop("search options are not needed for a variant search")
  )
  expect_no_error(with_clean_search_cache(search_gpmap("16:53786615")))
})

test_that("search_gpmap() still searches the options for text", {
  options_calls <- 0
  local_mocked_bindings(
    search_options_api = function() {
      options_calls <<- options_calls + 1
      return(list(search_terms = data.frame(
        type = c("trait", "gene"),
        name = c("Body mass index", "FTO"),
        alt_name = c(NA, "fat mass and obesity associated"),
        type_id = c("1992", "ENSG00000140718"),
        num_coloc_groups = c(10L, 3L),
        num_rare_results = c(0L, 1L),
        stringsAsFactors = FALSE
      )))
    },
    search_variants_api = function(...) stop("not a variant search")
  )
  result <- with_clean_search_cache({
    first <- search_gpmap("body mass")
    second <- search_gpmap("fto")
    list(first = first, second = second)
  })
  expect_equal(result$first$call, "trait(1992)")
  expect_equal(result$second$call, "gene('ENSG00000140718')")
  expect_false("rsq" %in% names(result$first))
  # Cached after the first text search.
  expect_equal(options_calls, 1)
})
