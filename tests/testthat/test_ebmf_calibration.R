library(testthat)

make_ebmf_result <- function() {
  sim <- simulate_trait(
    n_coloc_groups = 40,
    K = 2,
    module_sizes = c(10, 10),
    n_background_snps = 20,
    drivers_per_module = 4,
    n_background_traits = 12,
    log_se_sd = 0.5,
    seed = 11
  )
  res <- run_univariate_clustering(
    sim$trait_object,
    cluster_type = "ebmf",
    min_snp_signals = 2,
    min_module_size = 3
  )
  list(sim = sim, res = res)
}

test_that("ebmf_posterior_table returns one row per SNP x program", {
  r <- make_ebmf_result()
  expect_length(r$res$clusters, 0)
  expect_true(is.matrix(r$res$cluster_membership))
  tab <- ebmf_posterior_table(r$res)
  fit <- r$res$cluster_details$flash_fit
  expect_equal(nrow(tab), nrow(fit$F_pm) * fit$n_factors)
  expect_setequal(
    names(tab),
    c("snp_id", "program", "loading", "abs_loading", "lfsr")
  )
  expect_true(all(tab$abs_loading >= 0))
})

test_that("ebmf_posterior_table rejects non-EBMF results", {
  sim <- simulate_trait(n_coloc_groups = 30, K = 0, seed = 3)
  res <- run_univariate_clustering(sim$trait_object)
  expect_error(ebmf_posterior_table(res), "cluster_type")
})

test_that("calibrate_ebmf_programs returns calibrated gates", {
  r <- make_ebmf_result()
  cal <- calibrate_ebmf_programs(
    r$res,
    n_null = 3,
    verbose = FALSE
  )
  expect_setequal(names(cal), c(
    "programs", "memberships", "null_summary", "settings"
  ))
  expect_true(all(is.finite(unlist(cal$null_summary$max_masses))))
  expect_setequal(
    names(cal$programs),
    c("program", "loading_mass", "n_core", "n_candidate", "detected")
  )
})

test_that("stability scores EBMF programs across trait subsamples", {
  r <- make_ebmf_result()
  out <- stability_ebmf_programs(
    r$res,
    n_rep = 2,
    frac_traits = 0.8,
    top_n = 5,
    verbose = FALSE
  )
  expect_true(all(c("program", "n_ref", "replication", "sd_replication") %in%
                    names(out)))
  expect_true(all(out$replication >= 0 & out$replication <= 1))
})

test_that("louvain path is untouched by calibration helpers", {
  sim <- simulate_trait(n_coloc_groups = 40, K = 2, seed = 5)
  res <- run_univariate_clustering(sim$trait_object, louvain_gamma = 1.5,
                                   similarity_threshold = 0.2)
  expect_identical(res$parameters$cluster_type, "louvain")
  expect_error(ebmf_posterior_table(res), "cluster_type")
  expect_error(
    calibrate_ebmf_programs(res, n_null = 1, verbose = FALSE),
    "cluster_type"
  )
})
