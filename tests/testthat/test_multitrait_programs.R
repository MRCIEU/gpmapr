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

test_that("build_program_families builds reciprocal best-match families", {
  set.seed(6)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:15)
  hc <- rep(TRUE, 20)
  # The signal is planted on the locus axis, which is what the concordance
  # statistic reads: t1:1 and t2:1 load on the same loci, the others do not.
  shared <- c(rep(1.5, 10), rep(0, 10))
  programs <- list(
    make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), hc,
                 feats, stats::rnorm(15)),
    make_program("t1", 2, loci, stats::rnorm(20), hc, feats, stats::rnorm(15)),
    make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), hc,
                 feats, stats::rnorm(15)),
    make_program("t2", 2, loci, stats::rnorm(20), hc, feats, stats::rnorm(15))
  )
  res <- compare_program_pairs_loadings(programs, n_perm = 200, seed = 1)
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
  res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)
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
