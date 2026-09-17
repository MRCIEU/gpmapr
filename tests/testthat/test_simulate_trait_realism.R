library(testthat)

sim_graph <- function(..., seed = 1) {
  sim <- simulate_trait(
    K = 0, n_traits = 600, n_coloc_groups = 300, noise_sd = 1,
    effect_size = 6, background_effect_scale = 0.5,
    target_pattern = "dense", seed = seed, ...
  )
  res <- run_univariate_clustering(
    sim$trait_object, trait_subset = "phenotypic", min_snp_signals = 2,
    compress_method = "asinh", compress_scale = 5, ebmf_greedy_Kmax = 2L
  )
  simulated_graph_diagnostics(
    res$x_matrix,
    target_trait_id = res$parameters$target_trait_id,
    n_probe = 150
  )
}

test_that("background_corr = 0 reproduces the independent-background behaviour", {
  # Regression guard: existing simulation results were produced with
  # independent, random-signed background, which contributes nothing to
  # SNP-SNP similarity. That must remain the default.
  expect_equal(formals(simulate_trait)$background_corr, 0)
  expect_equal(formals(simulate_trait)$snp_pleiotropy_sd, 0)

  d <- sim_graph(p_active_background = 0.10, background_sparsity_sd = 1.2)
  expect_lt(abs(d$mean_similarity_no_target), 0.02)
})

test_that("background_corr raises similarity that survives removing the target", {
  low <- sim_graph(p_active_background = 0.02, background_sparsity_sd = 1.8,
                   background_corr = 0)
  high <- sim_graph(p_active_background = 0.02, background_sparsity_sd = 1.8,
                    background_corr = 0.9)
  expect_gt(high$mean_similarity_no_target, low$mean_similarity_no_target)
  expect_gt(high$mean_similarity, low$mean_similarity)
})

test_that("snp_pleiotropy_sd lowers pairwise overlap and lifts similarity", {
  flat <- sim_graph(p_active_background = 0.02, background_sparsity_sd = 1.8,
                    background_corr = 0.7, snp_pleiotropy_sd = 0)
  varied <- sim_graph(p_active_background = 0.02, background_sparsity_sd = 1.8,
                      background_corr = 0.7, snp_pleiotropy_sd = 2)
  # Heterogeneous per-SNP pleiotropy makes a typical pair share fewer traits,
  # and the sparse SNPs that result give high cosines on what little they share.
  expect_lt(varied$median_pair_overlap, flat$median_pair_overlap)
  expect_gt(varied$mean_similarity, flat$mean_similarity)
  # The point of the calibration: a random SNP set must be able to clear the
  # old absolute 0.3 gate, as it does ~94% of the time on real data. Where it
  # cannot, a null simulation study cannot detect that the gate is
  # uninformative.
  expect_gt(varied$random_set_pass_rate, 0.5)
  expect_equal(flat$random_set_pass_rate, 0, tolerance = 0.05)
})

test_that("background_sparsity_sd varies shape at roughly constant density", {
  # The lognormal multiplier is normalised to mean 1. Unnormalised it also
  # raised mean density (about 2x at sd 1.2, 5x at sd 1.8), which confounded
  # the null study's background_sparsity_sd axis with a density change.
  low <- sim_graph(p_active_background = 0.05, background_sparsity_sd = 0)
  high <- sim_graph(p_active_background = 0.05, background_sparsity_sd = 1.8)
  expect_equal(high$density, low$density, tolerance = 0.5)
})

test_that("simulated_graph_diagnostics reports the calibration columns", {
  d <- sim_graph(p_active_background = 0.02, background_sparsity_sd = 1.2)
  targets <- real_bmi_graph_targets()
  expect_true(all(names(targets) %in% names(d)))
  expect_true(is.finite(d$mean_similarity))
  expect_true(d$density > 0 && d$density < 1)
  expect_true(d$random_set_pass_rate >= 0 && d$random_set_pass_rate <= 1)
})
