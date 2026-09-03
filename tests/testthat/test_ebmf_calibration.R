library(testthat)

make_ebmf_result <- function() {
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

test_that("calibrate_ebmf_programs returns factor strength and null calibration", {
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
    c("program", "loading_mass", "raw_factor_signal", "factor_strength",
      "n_candidate")
  )
  expect_true(all(is.finite(cal$programs$raw_factor_signal)))
  expect_true(all(cal$programs$factor_strength >= 0))
  expect_true(all(cal$programs$factor_strength <= 1))
  expect_false("factor_strengths" %in% names(cal$null_summary))
  expect_false("factor_strength_quantiles" %in% names(cal$null_summary))
  expect_false("n_core" %in% names(cal$programs))
  expect_false("detected" %in% names(cal$programs))
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
  posterior <- ebmf_posterior_table(r$res)
  filtered_counts <- posterior |>
    dplyr::filter(
      is.finite(abs_loading), abs_loading > 0,
      !is.na(lfsr), lfsr < r$res$parameters$ebmf_lfsr_threshold,
      abs_loading > r$res$parameters$ebmf_magnitude_threshold
    ) |>
    dplyr::count(program, name = "n_filtered")
  expect_equal(
    out$n_ref,
    pmin(5L, filtered_counts$n_filtered[
      match(out$program, filtered_counts$program)
    ])
  )
})

test_that("parallel stability replicates match serial results", {
  r <- make_ebmf_result()
  s1 <- stability_ebmf_programs(r$res, n_rep = 3, top_n = 5, cores = 1,
                                verbose = FALSE)
  s5 <- stability_ebmf_programs(r$res, n_rep = 3, top_n = 5, cores = 5,
                                verbose = FALSE)
  expect_identical(s1, s5)
})

test_that("summarise_ebmf_programs folds filters and additional scores", {
  r <- make_ebmf_result()
  ps <- summarise_ebmf_programs(
    r$res,
    n_null = 3,
    n_rep = 2,
    verbose = FALSE
  )
  expect_true(all(c(
    "program", "n_snps", "n_snps_filtered", "mean_internal_similarity",
    "connectedness",
    "raw_factor_signal", "factor_strength", "factor_strength_per_snp",
    "replication",
    "posterior_evidence", "max_abs_factor_corr", "redundant",
    "size_pass", "internal_pass", "stability_pass",
    "status"
  ) %in% names(ps$programs)))
  expect_true(all(ps$programs$n_snps >= ps$programs$n_snps_filtered))
  filtered_counts <- ps$assigned |>
    dplyr::count(program, name = "n")
  expect_equal(
    ps$programs$n_snps_filtered,
    filtered_counts$n[match(ps$programs$program, filtered_counts$program)]
  )
  expect_true(all(is.na(ps$programs$factor_strength) |
                    ps$programs$factor_strength >= 0))
  expect_true(all(is.na(ps$programs$factor_strength_per_snp) |
                    ps$programs$factor_strength_per_snp >= 0))
  expect_true(all(ps$programs$status == "valid" |
                    grepl("^(size|internal_similarity|connectedness|stability|redundancy)",
                          ps$programs$status)))
  expect_true(is.null(ps$factor_correlation) || is.matrix(ps$factor_correlation))
  expect_true(all(c(
    "max_pair_redundancy", "most_redundant_program", "redundancy_pass"
  ) %in% names(ps$programs)))
  expect_true(all(ps$programs$redundancy_pass %in% c(TRUE, FALSE)))
  expect_true(all(is.na(ps$programs$max_pair_redundancy) |
                    (ps$programs$max_pair_redundancy >= 0 &
                       ps$programs$max_pair_redundancy <= 1)))
  expect_setequal(
    names(ps),
    c("programs", "memberships", "assigned", "factor_correlation",
      "null_summary", "settings")
  )
})

test_that("factor strength is reported but does not gate status", {
  r <- make_ebmf_result()
  ps <- summarise_ebmf_programs(
    r$res,
    n_null = 3,
    n_rep = 0,
    verbose = FALSE
  )
  expect_false("null_factor_strength_q" %in% names(ps$programs))
  expect_false("factor_strength_null_q" %in% names(ps$settings))
  expect_false("factor_strength_pass" %in% names(ps$programs))
  expect_false("factor_strength_min" %in% names(ps$settings))
  expect_true(all(c("n_snps", "factor_strength", "factor_strength_per_snp") %in%
                    names(ps$programs)))
  tol <- 1e-12
  expect_true(all(is.na(ps$programs$factor_strength_per_snp) |
                    abs(ps$programs$factor_strength_per_snp -
                          ps$programs$factor_strength /
                            sqrt(ps$programs$n_snps_filtered)) < tol))
  expect_true(all(ps$programs$status == "valid" |
                    grepl("^(size|internal_similarity|connectedness|stability|redundancy)",
                          ps$programs$status)))
})

test_that("reciprocal SNP-containment redundancy flags near-duplicate programs", {
  # A = {s1..s10}, B = {s1..s9}: shared = 9, pair = min(9/10, 9/9) = 0.9
  assigned <- data.frame(
    snp_id = c(paste0("s", 1:10), paste0("s", 1:9)),
    program = c(rep(1L, 10), rep(2L, 9))
  )
  red <- .program_membership_redundancy(assigned)
  expect_equal(red$program, 1:2)
  expect_equal(round(red$max_pair_redundancy, 3), c(0.9, 0.9))
  expect_equal(red$most_redundant_program, c(2L, 1L))

  # Disjoint programs -> pair redundancy 0
  assigned2 <- data.frame(
    snp_id = c(paste0("s", 1:10), paste0("t", 1:10)),
    program = c(rep(1L, 10), rep(2L, 10))
  )
  red2 <- .program_membership_redundancy(assigned2)
  expect_equal(red2$max_pair_redundancy, c(0, 0))
  expect_true(all(red2$max_pair_redundancy < 0.9))

  # Disjoint programs -> pair redundancy 0 (not NA), partner still reported
  assigned3 <- data.frame(
    snp_id = c("s1", "s2", "s3"),
    program = c(1L, 1L, 2L)
  )
  red3 <- .program_membership_redundancy(assigned3)
  expect_equal(red3$program, 1:2)
  expect_equal(red3$max_pair_redundancy, c(0, 0))
  expect_equal(red3$most_redundant_program, c(2L, 1L))

  # A single program (no partner) reports NA redundancy
  red4 <- .program_membership_redundancy(
    data.frame(snp_id = c("s1", "s2"), program = c(1L, 1L))
  )
  expect_equal(red4$program, 1L)
  expect_true(is.na(red4$max_pair_redundancy[1]))
  expect_true(is.na(red4$most_redundant_program[1]))
})

test_that("summarise_ebmf_programs checks stability only for programs passing coherence", {
  r <- make_ebmf_result()
  ps <- summarise_ebmf_programs(
    r$res,
    n_null = 3,
    n_rep = 2,
    verbose = FALSE
  )
  coherent <- ps$programs$internal_pass %in% TRUE
  expect_true(all(ps$programs$stability_checked[coherent]))
  expect_true(all(!ps$programs$stability_checked[!coherent]))
  expect_true(all(is.finite(ps$programs$replication[coherent])))
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
