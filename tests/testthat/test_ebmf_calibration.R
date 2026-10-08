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
    min_snp_signals = 2,
    min_module_size = 3
  )
  list(sim = sim, res = res)
}

test_that("ebmf_posterior_table returns one row per SNP x program", {
  r <- make_ebmf_result()
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
  # The reference set is size-relative: capped at top_n, but never more than
  # half the program's filtered SNPs, so programs are matched on comparable
  # fractions of themselves rather than on a flat top_n that makes replication
  # track program size.
  n_filtered <- filtered_counts$n_filtered[
    match(out$program, filtered_counts$program)
  ]
  expect_equal(
    out$n_ref,
    pmax(1L, pmin(5L, as.integer(ceiling(n_filtered / 2))))
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

test_that("summarise_ebmf_programs reports metrics and a confidence tier", {
  r <- make_ebmf_result()
  ps <- summarise_ebmf_programs(
    r$res,
    n_null = 3,
    n_rep = 2,
    verbose = FALSE
  )
  expect_true(all(c(
    "program", "n_snps", "n_snps_filtered", "mean_internal_similarity",
    "connectedness", "excess_similarity", "excess_connectedness",
    "similarity_lift", "similarity_z", "similarity_emp_p",
    "similarity_q", "factor_pve", "replication", "sd_replication",
    "parent_program", "parent_containment", "posterior_evidence",
    "max_abs_factor_corr", "redundant", "max_pair_redundancy",
    "most_redundant_program", "effective_n_traits", "confidence_tier",
    "depth", "program_label"
  ) %in% names(ps$programs)))
  expect_false(any(c(
    "weighted_internal", "weighted_quorum", "n_eff", "coherence_q",
    "quorum_q", "factor_strength", "factor_strength_per_snp",
    "raw_factor_signal", "strength_q", "strength_pass", "internal_pass",
    "stability_pass", "stability_checked", "redundancy_pass",
    "internal_similarity_pass", "connectedness_pass", "fail_reason", "status",
    "size_pass", "similarity_pass", "mean_external_similarity", "separation"
  ) %in% names(ps$programs)))
  expect_true(all(ps$programs$n_snps >= ps$programs$n_snps_filtered))
  filtered_counts <- ps$assigned |>
    dplyr::count(program, name = "n")
  expect_equal(
    ps$programs$n_snps_filtered,
    filtered_counts$n[match(ps$programs$program, filtered_counts$program)]
  )
  expect_true(is.ordered(ps$programs$confidence_tier))
  expect_equal(
    as.character(ps$programs$confidence_tier),
    as.character(.program_confidence_tier(
      ps$programs$similarity_q, ps$programs$connectedness,
      ps$settings$similarity_q_threshold, ps$settings$high_similarity_q,
      ps$settings$high_min_connectedness,
      effective_n_traits = ps$programs$effective_n_traits,
      min_effective_traits = ps$settings$min_effective_traits
    ))
  )
  # Programs below the size floor are dropped before any check.
  expect_true(all(ps$programs$n_snps_filtered >= ps$settings$min_module_size))
  expect_true(all(ps$dropped_programs$n_snps_filtered < ps$settings$min_module_size))
  expect_false(any(ps$dropped_programs$program %in% ps$programs$program))
  expect_true(all(ps$assigned$program %in% ps$programs$program))
  expect_true(all(c("mean_similarity", "frac_pairs_edge", "frac_pairs_abs_edge",
                    "median_shared_traits") %in% names(ps$background)))
  fit <- r$res$cluster_details$flash_fit
  expect_equal(ps$programs$factor_pve, as.numeric(fit$pve[ps$programs$program]))
  # Effective number of traits from the trait loadings, target row excluded.
  target_id <- as.character(r$res$parameters$target_trait_id)
  L <- fit$L_pm[rownames(fit$L_pm) != target_id, , drop = FALSE]
  expected_traits <- vapply(ps$programs$program, function(k) {
    p <- L[, k]^2 / sum(L[, k]^2)
    return(1 / sum(p^2))
  }, numeric(1))
  expect_equal(ps$programs$effective_n_traits, expected_traits)
  # Replication is reported for every program with members.
  has_members <- ps$programs$n_snps_filtered > 0
  expect_true(all(is.finite(ps$programs$replication[has_members])))
  expect_true(is.null(ps$factor_correlation) || is.matrix(ps$factor_correlation))
  expect_true(all(is.na(ps$programs$max_pair_redundancy) |
                    (ps$programs$max_pair_redundancy >= 0 &
                       ps$programs$max_pair_redundancy <= 1)))
  expect_true(all(c("similarity_n_perm", "similarity_q_threshold",
                    "high_similarity_q", "high_min_connectedness",
                    "min_effective_traits") %in%
                    names(ps$settings)))
  expect_false(any(c("min_excess_similarity", "min_excess_connectedness") %in%
                     names(ps$settings)))
  expect_false(any(c("stability_threshold", "strength_q_threshold",
                     "snp_redundancy_threshold") %in% names(ps$settings)))
  expect_false("factor_pve" %in% names(ps$null_summary))
  expect_setequal(
    names(ps),
    c("programs", "dropped_programs", "background", "memberships", "assigned",
      "factor_correlation", "null_summary", "settings")
  )
})

test_that("summarise_ebmf_programs rows follow order_programs_by_confidence()", {
  r <- make_ebmf_result()
  ps <- summarise_ebmf_programs(r$res, n_null = 3, n_rep = 0, verbose = FALSE)
  reordered <- order_programs_by_confidence(
    ps$programs[rev(seq_len(nrow(ps$programs))), , drop = FALSE]
  )
  expect_equal(reordered$program, ps$programs$program)
})

test_that("the similarity matrix excludes the target row", {
  r <- make_ebmf_result()
  target_id <- as.character(r$sim$trait_object$trait$id)
  scaled <- r$sim$trait_object
  on_target <- as.character(scaled$coloc_groups$trait_id) == target_id
  # Rescaling the target keeps its signs, so orientation and every other row
  # are unchanged; only the target row of X_star moves.
  scaled$coloc_groups$beta[on_target] <- scaled$coloc_groups$beta[on_target] * 3
  res2 <- run_univariate_clustering(scaled, min_snp_signals = 2, min_module_size = 3)
  expect_false(isTRUE(all.equal(
    r$res$x_star[target_id, ], res2$x_star[target_id, ]
  )))
  expect_equal(r$res$s_matrix, res2$s_matrix)
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
