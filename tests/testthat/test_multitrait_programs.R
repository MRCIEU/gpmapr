library(testthat)

make_program_loadings <- function(trait_id, program_id, loadings, high_confidence,
                                  snp_ids = paste0("rs", seq_along(loadings))) {
  return(data.frame(
    program_id = program_id,
    trait_id = trait_id,
    trait_name = trait_id,
    program = 1L,
    snp_id = snp_ids,
    loading = loadings,
    abs_loading = abs(loadings),
    lfsr = 0.001,
    high_confidence = high_confidence,
    stringsAsFactors = FALSE
  ))
}

test_that("compare_program_pairs detects concordant shared programs", {
  set.seed(1)
  n <- 200
  signal <- stats::rnorm(n)
  a <- make_program_loadings("t1", "t1:1", signal, abs(signal) > 0.5)
  b <- make_program_loadings("t2", "t2:1", signal + stats::rnorm(n, sd = 0.05),
                             abs(signal) > 0.5)
  res <- compare_program_pairs(dplyr::bind_rows(a, b), n_perm = 200, seed = 1)
  expect_equal(nrow(res$pairs), 1L)
  expect_gt(res$pairs$r, 0.9)
  expect_true(res$pairs$shared)
  expect_equal(res$pairs$direction, "concordant")
  expect_gt(res$pairs$p_emp, 0)
  expect_lte(res$pairs$p_emp, 1)
})

test_that("compare_program_pairs flags antagonistic programs from signed r", {
  n <- 150
  signal <- rep(c(1, -1), length.out = n)
  a <- make_program_loadings("t1", "t1:1", signal, rep(TRUE, n))
  b <- make_program_loadings("t2", "t2:1", -signal, rep(TRUE, n))
  res <- compare_program_pairs(dplyr::bind_rows(a, b), n_perm = 200, seed = 1)
  expect_equal(res$pairs$r, -1)
  expect_equal(res$pairs$direction, "antagonistic")
  expect_true(res$pairs$shared)
})

test_that("compare_program_pairs does not call orthogonal programs shared", {
  n <- 200
  a_load <- rep(c(1, -1), n / 2)
  b_load <- rep(c(1, 1, -1, -1), n / 4)
  a <- make_program_loadings("t1", "t1:1", a_load, rep(TRUE, n))
  b <- make_program_loadings("t2", "t2:1", b_load, rep(TRUE, n))
  res <- compare_program_pairs(dplyr::bind_rows(a, b), n_perm = 200, seed = 1)
  expect_equal(res$pairs$r, 0)
  expect_false(res$pairs$shared)
})

test_that("compare_program_pairs only forms cross-trait pairs", {
  n <- 60
  a <- make_program_loadings("t1", "t1:1", stats::rnorm(n), rep(TRUE, n))
  b <- make_program_loadings("t1", "t1:2", stats::rnorm(n), rep(TRUE, n))
  res <- compare_program_pairs(dplyr::bind_rows(a, b), n_perm = 20, seed = 1)
  expect_equal(nrow(res$pairs), 0L)
  expect_setequal(
    names(res$pairs),
    c("program_id_a", "program_id_b", "trait_id_a", "trait_id_b", "program_a",
      "program_b", "n_snps_union", "n_snps_used", "r", "jaccard",
      "n_high_conf_shared", "n_high_conf_union", "p_emp", "fdr", "shared",
      "direction")
  )
})

test_that("compare_program_pairs accepts a list of extraction results", {
  n <- 80
  signal <- stats::rnorm(n)
  a <- make_program_loadings("t1", "t1:1", signal, abs(signal) > 0.5)
  b <- make_program_loadings("t2", "t2:1", signal, abs(signal) > 0.5)
  res <- compare_program_pairs(
    list(list(loadings = a), list(loadings = b)),
    n_perm = 50,
    seed = 1
  )
  expect_equal(nrow(res$pairs), 1L)
  expect_gt(res$pairs$r, 0.99)
  expect_equal(res$pairs$jaccard, 1)
})

test_that("build_program_families groups chained programs and keeps singletons", {
  pairs <- data.frame(
    program_id_a = c("t1:1", "t2:1"),
    program_id_b = c("t2:1", "t3:1"),
    trait_id_a = c("t1", "t2"),
    trait_id_b = c("t2", "t3"),
    program_a = c(1L, 1L),
    program_b = c(1L, 1L),
    n_snps_union = c(10L, 10L),
    n_snps_used = c(10L, 10L),
    r = c(0.8, -0.6),
    jaccard = c(0.5, 0.3),
    n_high_conf_shared = c(2L, 1L),
    n_high_conf_union = c(4L, 4L),
    p_emp = c(0.01, 0.02),
    fdr = c(0.01, 0.02),
    shared = c(TRUE, TRUE),
    direction = c("concordant", "antagonistic"),
    stringsAsFactors = FALSE
  )
  program_info <- data.frame(
    program_id = c("t1:1", "t2:1", "t3:1", "t4:1"),
    trait_id = c("t1", "t2", "t3", "t4"),
    stringsAsFactors = FALSE
  )
  fam <- build_program_families(pairs, program_info = program_info)
  expect_equal(nrow(fam$families), 2L)
  big <- fam$families[fam$families$n_programs == 3L, ]
  expect_equal(big$n_traits, 3L)
  expect_equal(big$n_edges, 2L)
  expect_equal(big$n_concordant, 1L)
  expect_equal(big$n_antagonistic, 1L)
  expect_equal(nrow(fam$families[fam$families$n_programs == 1L, ]), 1L)
  expect_false(is.null(fam$graph))
  expect_true(inherits(fam$graph, "igraph"))
})

test_that("build_program_families handles no significant edges", {
  pairs <- compare_program_pairs(
    dplyr::bind_rows(
      make_program_loadings("t1", "t1:1", rep(c(1, -1), 50), rep(TRUE, 100)),
      make_program_loadings("t2", "t2:1", rep(c(1, 1, -1, -1), 25), rep(TRUE, 100))
    ),
    n_perm = 20,
    seed = 1
  )
  fam <- build_program_families(pairs)
  expect_equal(nrow(fam$families), 0L)
  expect_null(fam$graph)
  expect_true(nrow(fam$nodes) >= 2)
})

test_that("extract_program_loadings pulls valid programs from an EBMF fit", {
  sim <- simulate_trait(
    n_coloc_groups = 40,
    K = 2,
    module_sizes = c(10, 10),
    n_background_snps = 20,
    n_traits_per_module = 4,
    n_background_traits = 12,
    log_se_sd = 0.5,
    seed = 11
  )
  res <- run_univariate_clustering(
    sim$trait_object,
    min_snp_signals = 2,
    min_module_size = 3
  )
  summary <- summarise_ebmf_programs(res, n_null = 3, n_rep = 0, verbose = FALSE)
  ex <- extract_program_loadings(res, summary)
  expect_setequal(names(ex), c("loadings", "programs"))
  expect_setequal(
    names(ex$loadings),
    c("program_id", "trait_id", "trait_name", "program", "snp_id", "loading",
      "abs_loading", "lfsr", "high_confidence")
  )
  valid_ids <- summary$programs$program[summary$programs$status == "valid"]
  expect_setequal(unique(ex$loadings$program), valid_ids)
  expect_true(all(ex$loadings$program_id ==
                    paste0(ex$loadings$trait_id, ":", ex$loadings$program)))
  expect_true(all(ex$programs$n_high_confidence >= 0))
  expect_true(all(abs(ex$loadings$loading) == ex$loadings$abs_loading))
})
