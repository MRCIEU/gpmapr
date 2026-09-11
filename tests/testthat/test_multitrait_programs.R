library(testthat)

make_program <- function(trait_id, program, loci, locus_loadings, high_confidence,
                         features, profile_loadings, lfsr = 0.001) {
  return(list(
    loadings = data.frame(
      program_id = paste0(trait_id, ":", program),
      trait_id = trait_id,
      trait_name = trait_id,
      program = as.integer(program),
      snp_id = paste0(trait_id, "_", program, "_", loci),
      coloc_group_id = loci,
      chr = NA_character_,
      bp = NA_real_,
      loading = locus_loadings,
      abs_loading = abs(locus_loadings),
      lfsr = lfsr,
      high_confidence = high_confidence,
      stringsAsFactors = FALSE
    ),
    profiles = data.frame(
      program_id = paste0(trait_id, ":", program),
      trait_id = trait_id,
      trait_name = trait_id,
      program = as.integer(program),
      feature_trait_id = features,
      feature_trait_name = features,
      loading = profile_loadings,
      lfsr = lfsr,
      stringsAsFactors = FALSE
    )
  ))
}

pair_columns <- c(
  "program_id_a", "program_id_b", "trait_id_a", "trait_id_b", "program_a",
  "program_b", "n_loci_a", "n_loci_b", "n_loci_shared", "n_loci_union",
  "jaccard_loci", "p_locus", "candidate", "n_features_shared",
  "n_features_union", "profile_r", "profile_r_resid", "p_profile",
  "q_profile", "direction", "shared"
)

test_that("locus membership uses the high-confidence coloc groups", {
  loci <- paste0("g", 1:20)
  hc <- c(rep(TRUE, 4), rep(FALSE, 16))
  a <- make_program("t1", 1, loci, stats::rnorm(20), hc, paste0("f", 1:8), stats::rnorm(8))
  b <- make_program("t2", 1, loci, stats::rnorm(20), hc, paste0("f", 1:8), stats::rnorm(8))
  res <- compare_program_pairs(list(a, b), n_perm = 50, seed = 1)
  expect_equal(res$pairs$n_loci_a, 4L)
  expect_equal(res$pairs$n_loci_b, 4L)
  expect_equal(res$pairs$n_loci_shared, 4L)
  expect_true(res$pairs$candidate)
})

test_that("a planted shared program is significant, shared and concordant", {
  set.seed(1)
  loci <- paste0("g", 1:12)
  feats <- paste0("f", 1:15)
  hc <- rep(TRUE, 12)
  profile <- stats::rnorm(15)
  a <- make_program("t1", 1, loci, stats::rnorm(12), hc, feats, profile)
  b <- make_program("t2", 1, loci, stats::rnorm(12), hc, feats, profile)
  res <- compare_program_pairs(list(a, b), n_perm = 200, seed = 1)
  expect_equal(res$pairs$n_loci_shared, 12L)
  expect_true(res$pairs$shared)
  expect_equal(res$pairs$direction, "concordant")
  expect_true(all(res$matches$significant))
  expect_equal(res$matches$best_partner[1], "t2:1")
  expect_equal(res$matches$best_locus_jaccard[1], 1)
})

test_that("programs with disjoint loci are not locus-matched candidates", {
  set.seed(2)
  feats <- paste0("f", 1:10)
  profile <- stats::rnorm(10)
  a <- make_program(
    "t1", 1, paste0("a", 1:10), stats::rnorm(10), rep(TRUE, 10), feats, profile
  )
  b <- make_program(
    "t2", 1, paste0("b", 1:10), stats::rnorm(10), rep(TRUE, 10), feats, profile
  )
  res <- compare_program_pairs(list(a, b), n_perm = 50, seed = 1)
  expect_equal(res$pairs$n_loci_shared, 0L)
  expect_false(res$pairs$candidate)
  expect_false(res$pairs$shared)
  expect_true(all(is.na(res$matches$q)))
})

test_that("partialling out the common profile removes a shared global component", {
  set.seed(3)
  loci <- paste0("g", 1:12)
  feats <- paste0("f", 1:20)
  hc <- rep(TRUE, 12)
  common <- stats::rnorm(20)
  programs <- list()
  for (tr in c("t1", "t2")) {
    for (p in 1:3) {
      programs[[length(programs) + 1]] <- make_program(
        tr, p, loci, stats::rnorm(12), hc, feats,
        common + stats::rnorm(20, sd = 0.01)
      )
    }
  }
  raw <- compare_program_pairs(programs, n_perm = 20, seed = 1,
                               partial_out_common = FALSE)
  resid <- compare_program_pairs(programs, n_perm = 20, seed = 1,
                                 partial_out_common = TRUE)
  expect_gt(mean(abs(raw$pairs$profile_r), na.rm = TRUE), 0.9)
  expect_lt(mean(abs(resid$pairs$profile_r_resid), na.rm = TRUE), 0.2)
})

test_that("top_k locus membership is honoured", {
  set.seed(4)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:8)
  hc <- rep(FALSE, 20)  # force top_k
  locus_load <- stats::rnorm(20)
  a <- make_program("t1", 1, loci, locus_load, hc, feats, stats::rnorm(8))
  b <- make_program("t2", 1, loci, locus_load, hc, feats, stats::rnorm(8))
  res <- compare_program_pairs(
    list(a, b), n_perm = 20, seed = 1,
    locus_membership = "top_k", top_k = 5L
  )
  expect_equal(res$pairs$n_loci_a, 5L)
  expect_equal(res$pairs$n_loci_shared, 5L)
})

test_that("compare_program_pairs only forms cross-trait pairs", {
  set.seed(5)
  loci <- paste0("g", 1:10)
  feats <- paste0("f", 1:8)
  a <- make_program("t1", 1, loci, stats::rnorm(10), rep(TRUE, 10), feats, stats::rnorm(8))
  b <- make_program("t1", 2, loci, stats::rnorm(10), rep(TRUE, 10), feats, stats::rnorm(8))
  res <- compare_program_pairs(list(a, b), n_perm = 20, seed = 1)
  expect_equal(nrow(res$pairs), 0L)
  expect_setequal(names(res$pairs), pair_columns)
  expect_setequal(
    names(res$matches),
    c("program_id", "trait_id", "trait_name", "program", "best_partner",
      "best_profile_r", "best_locus_jaccard", "best_n_loci_shared", "emp_p",
      "q", "significant")
  )
})

test_that("build_program_families builds reciprocal best-match families", {
  set.seed(6)
  loci <- paste0("g", 1:12)
  feats <- paste0("f", 1:15)
  hc <- rep(TRUE, 12)
  profile <- stats::rnorm(15)
  programs <- list(
    make_program("t1", 1, loci, stats::rnorm(12), hc, feats, profile),
    make_program("t1", 2, loci, stats::rnorm(12), hc, feats, stats::rnorm(15)),
    make_program("t2", 1, loci, stats::rnorm(12), hc, feats, profile),
    make_program("t2", 2, loci, stats::rnorm(12), hc, feats, stats::rnorm(15))
  )
  res <- compare_program_pairs(programs, n_perm = 200, seed = 1,
                               partial_out_common = FALSE)
  fam <- build_program_families(res, program_info = res$matches)
  fam_of <- function(pid) fam$nodes$family[fam$nodes$program_id == pid]
  expect_equal(nrow(fam$nodes), 4L)
  expect_equal(fam_of("t1:1"), fam_of("t2:1"))
  expect_false(fam_of("t1:1") == fam_of("t1:2"))
})

test_that("build_program_families handles no shared edges", {
  set.seed(7)
  loci <- paste0("g", 1:10)
  feats <- paste0("f", 1:10)
  a <- make_program("t1", 1, loci, stats::rnorm(10), rep(TRUE, 10), feats, stats::rnorm(10))
  b <- make_program("t2", 1, loci, stats::rnorm(10), rep(TRUE, 10), feats, stats::rnorm(10))
  res <- compare_program_pairs(list(a, b), n_perm = 20, seed = 1)
  fam <- build_program_families(res)
  expect_null(fam$graph)
  expect_equal(nrow(fam$families), 0L)
  expect_true(nrow(fam$nodes) >= 2)
})

test_that("extract_program_loadings returns loci, profiles and meta", {
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
  expect_setequal(names(ex), c("loadings", "profiles", "programs"))
  expect_true(all(c(
    "program_id", "trait_id", "trait_name", "program", "snp_id",
    "coloc_group_id", "chr", "bp", "loading", "abs_loading", "lfsr",
    "high_confidence"
  ) %in% names(ex$loadings)))
  expect_true(all(c(
    "program_id", "trait_id", "trait_name", "program", "feature_trait_id",
    "feature_trait_name", "loading", "lfsr"
  ) %in% names(ex$profiles)))
  valid_ids <- summary$programs$program[summary$programs$status == "valid"]
  expect_setequal(unique(ex$loadings$program), valid_ids)
  expect_setequal(unique(ex$profiles$program), valid_ids)
  expect_true(all(abs(ex$loadings$loading) == ex$loadings$abs_loading))
  expect_true(all(ex$programs$n_high_confidence >= 0))
  expect_true(all(ex$programs$n_features > 0))
})
