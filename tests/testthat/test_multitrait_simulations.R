library(testthat)

small_pair <- function(...) {
  # modifyList, not `...` straight through: a test that overrides one of the
  # defaults below would otherwise pass it twice and error.
  defaults <- list(
    n_loci_per_trait = 120L, n_shared_loci = 50L, K_shared = 3L,
    K_specific = 1L, module_size = 12L, n_drivers_per_program = 8L,
    n_background_traits = 150L, p_active_background = 0.02, seed = 1
  )
  return(do.call(
    simulate_trait_pair, utils::modifyList(defaults, list(...))
  ))
}

test_that("locus_overlap controls how many loci a shared program has in common", {
  for (lo in c(0, 0.25, 0.5, 1)) {
    sp <- small_pair(locus_overlap = lo, profile_overlap = 0.5)
    expect_equal(
      unique(sp$ground_truth$correspondence$n_loci_shared),
      round(12 * lo)
    )
  }
})

test_that("profile_overlap controls how many drivers a shared program has in common", {
  for (po in c(0, 0.5, 1)) {
    sp <- small_pair(locus_overlap = 0.5, profile_overlap = po)
    expect_equal(
      unique(sp$ground_truth$correspondence$n_drivers_shared),
      round(8 * po)
    )
  }
})

test_that("ground truth carries each program's per-locus effect, shared at shared loci", {
  sp <- small_pair(locus_overlap = 0.5, profile_overlap = 1)
  gt <- sp$ground_truth
  expect_setequal(names(gt$program_effects), names(gt$program_snps))
  for (k in names(gt$program_snps)) {
    expect_equal(names(gt$program_effects[[k]]), gt$program_snps[[k]])
    expect_true(all(gt$program_effects[[k]] > 0))
  }
  # A shared locus draws its effect once, so both traits see the same a_s.
  a <- gt$program_effects[["9001:shared1"]]
  b <- gt$program_effects[["9002:shared1"]]
  common <- intersect(names(a), names(b))
  expect_length(common, 6)
  expect_equal(a[common], b[common])
})

test_that("the two axes are independent", {
  # Disjoint loci with fully shared drivers: the case the profile axis exists
  # for, and the case the locus axis cannot see at all.
  sp <- small_pair(locus_overlap = 0, profile_overlap = 1)
  a <- unique(sp$traits[[1]]$coloc_groups$coloc_group_id)
  b <- unique(sp$traits[[2]]$coloc_groups$coloc_group_id)

  shared_snp_names <- lapply(names(sp$ground_truth$program_snps), function(k) {
    sp$ground_truth$program_snps[[k]]
  })
  names(shared_snp_names) <- names(sp$ground_truth$program_snps)
  a_shared1 <- shared_snp_names[["9001:shared1"]]
  b_shared1 <- shared_snp_names[["9002:shared1"]]
  expect_length(intersect(a_shared1, b_shared1), 0)

  # but the drivers are identical
  expect_equal(
    sp$ground_truth$program_drivers[["9001:shared1"]],
    sp$ground_truth$program_drivers[["9002:shared1"]]
  )
  # and the traits still share a background-study universe
  fa <- unique(sp$traits[[1]]$coloc_groups$trait_id)
  fb <- unique(sp$traits[[2]]$coloc_groups$trait_id)
  expect_gt(length(intersect(fa, fb)), 50)
})

test_that("directions alternate so sign is scoreable both ways", {
  sp <- small_pair(locus_overlap = 0.5, profile_overlap = 0.8)
  expect_setequal(
    sp$ground_truth$correspondence$direction,
    c("concordant", "antagonistic")
  )
})

test_that("K_shared = 0 gives a usable multi-trait null", {
  sp <- small_pair(locus_overlap = 0, profile_overlap = 0, K_shared = 0L)
  expect_equal(nrow(sp$ground_truth$correspondence), 0)
  # Both traits still carry programs of their own, so there is something to
  # compare -- it just must not link.
  expect_gt(nrow(sp$ground_truth$programs), 0)
  expect_equal(length(sp$traits), 2)
})

test_that("shared loci carry the same coloc_group_id in both traits", {
  # This is what makes a locus "the same" locus across traits; if the ids were
  # per-trait column positions the locus axis could never see any overlap.
  sp <- small_pair(locus_overlap = 1, profile_overlap = 0.5)
  a <- unique(sp$traits[[1]]$coloc_groups$coloc_group_id)
  b <- unique(sp$traits[[2]]$coloc_groups$coloc_group_id)
  expect_gte(length(intersect(a, b)), 36)
})

# Target effect of one trait at a planted program's loci, from its coloc groups.
target_beta <- function(sp, tid, snps) {
  cg <- sp$traits[[tid]]$coloc_groups
  cg <- cg[as.character(cg$trait_id) == tid, , drop = FALSE]
  return(cg$beta[match(snps, cg$variant_id)])
}

test_that("an antagonistic program has opposite-sign target effects in a common allele frame", {
  sp <- small_pair(
    K_shared = 2L, locus_overlap = 1, profile_overlap = 1, directions = c(1, -1)
  )
  snps_conc <- sp$ground_truth$program_snps[["9002:shared1"]]
  snps_anta <- sp$ground_truth$program_snps[["9002:shared2"]]

  expect_gt(mean(target_beta(sp, "9001", snps_anta) > 0), 0.9)
  expect_gt(mean(target_beta(sp, "9002", snps_anta) < 0), 0.9)
  expect_gt(mean(target_beta(sp, "9002", snps_conc) > 0), 0.9)

  # The same GWAS: driver effects carry the same sign in both traits' matrices.
  drivers <- sp$ground_truth$program_drivers[["9002:shared2"]]
  for (tid in c("9001", "9002")) {
    cells <- sp$x_matrices[[tid]][drivers, snps_anta]
    expect_gt(mean(cells > 0, na.rm = TRUE), 0.9)
  }
})

test_that("orienting to the target's risk allele turns antagonism into negative drivers", {
  sp <- small_pair(
    K_shared = 2L, locus_overlap = 1, profile_overlap = 1, directions = c(1, -1)
  )
  oriented <- orient_pleiotropy_matrix(sp$x_matrices[["9002"]], "9002")$x_matrix
  anta <- oriented[
    sp$ground_truth$program_drivers[["9002:shared2"]],
    sp$ground_truth$program_snps[["9002:shared2"]]
  ]
  conc <- oriented[
    sp$ground_truth$program_drivers[["9002:shared1"]],
    sp$ground_truth$program_snps[["9002:shared1"]]
  ]
  expect_gt(mean(anta < 0, na.rm = TRUE), 0.9)
  expect_gt(mean(conc > 0, na.rm = TRUE), 0.9)
})

test_that("shared loci carry one locus effect, so target effects correlate across traits", {
  # This is what module_rg() reads: the correlation of the two targets' own
  # effects across a shared program's loci, signed by the planted direction.
  sp <- small_pair(
    K_shared = 2L, locus_overlap = 1, profile_overlap = 1, directions = c(1, -1)
  )
  for (k in 1:2) {
    snps <- sp$ground_truth$program_snps[[sprintf("9001:shared%d", k)]]
    r <- stats::cor(target_beta(sp, "9001", snps), target_beta(sp, "9002", snps))
    if (k == 1) {
      expect_gt(r, 0.7)
    } else {
      expect_lt(r, -0.7)
    }
  }
})

test_that("simulate_trait_pair rejects an over-subscribed shared locus pool", {
  expect_error(
    simulate_trait_pair(
      n_loci_per_trait = 60L, n_shared_loci = 10L, K_shared = 3L,
      module_size = 12L, locus_overlap = 1, seed = 1
    ),
    "shared locus pool too small"
  )
})

test_that("map_programs_to_planted maps a dense feature table by loading weight", {
  # The feature table carries a loading for every study on every program, so a
  # set overlap would score a perfect match at 8/200 and reject it. The weighted
  # share must recover it.
  sp <- small_pair(locus_overlap = 0.5, profile_overlap = 0.5)
  drivers <- sp$ground_truth$program_drivers[["9001:shared1"]]
  all_features <- as.character(c(drivers, 30001:30192))
  profiles <- data.frame(
    program_id = "9001:1", trait_id = "9001", trait_name = "t",
    program = 1L,
    feature_trait_id = all_features,
    feature_trait_name = all_features,
    # All the loading mass sits on the planted drivers.
    loading = c(rep(1, length(drivers)), rep(0.01, length(all_features) - length(drivers))),
    lfsr = 0.5, stringsAsFactors = FALSE
  )
  mapped <- map_programs_to_planted(profiles, sp$ground_truth, by = "drivers")
  expect_equal(mapped$planted_label, "shared1")
  expect_gt(mapped$fold, 3)
})

test_that("evaluate_multitrait_simulation credits recall once per planted program", {
  sp <- small_pair(locus_overlap = 0.5, profile_overlap = 0.8)
  mapping <- data.frame(
    program_id = c("9001:1", "9001:2", "9002:1"),
    trait_id = c("9001", "9001", "9002"),
    planted_label = c("shared1", "shared1", "shared1"),
    score = 1, fold = 10, stringsAsFactors = FALSE
  )
  # One planted program recovered through two links: recall must count it once,
  # and precision must still see two links.
  pairs <- data.frame(
    program_id_a = c("9001:1", "9001:2"),
    program_id_b = c("9002:1", "9002:1"),
    link_tier = "primary",
    direction = "concordant",
    stringsAsFactors = FALSE
  )
  ev <- evaluate_multitrait_simulation(pairs, mapping, sp$ground_truth)
  expect_equal(ev$n_links, 2L)
  expect_equal(ev$n_true_positive, 2L)
  expect_equal(ev$n_recovered_shared, 1L)
  expect_equal(ev$recall, 1 / 3)
  expect_equal(ev$precision, 1)
})

test_that("evaluate_multitrait_simulation counts unmatched links as false positives", {
  sp <- small_pair(locus_overlap = 0.5, profile_overlap = 0.8)
  mapping <- data.frame(
    program_id = c("9001:1", "9002:1"),
    trait_id = c("9001", "9002"),
    # Different planted programs: the link is real-looking but wrong.
    planted_label = c("shared1", "shared2"),
    score = 1, fold = 10, stringsAsFactors = FALSE
  )
  pairs <- data.frame(
    program_id_a = "9001:1", program_id_b = "9002:1",
    link_tier = "primary", direction = "concordant", stringsAsFactors = FALSE
  )
  ev <- evaluate_multitrait_simulation(pairs, mapping, sp$ground_truth)
  expect_equal(ev$n_true_positive, 0L)
  expect_equal(ev$n_false_positive, 1L)
  expect_equal(ev$precision, 0)
  expect_true(is.na(ev$sign_accuracy))
})

test_that("sign accuracy is scored against the planted direction", {
  sp <- small_pair(locus_overlap = 0.5, profile_overlap = 0.8)
  planted <- sp$ground_truth$correspondence
  label <- planted$program_label[planted$direction == "antagonistic"][1]
  mapping <- data.frame(
    program_id = c("9001:1", "9002:1"), trait_id = c("9001", "9002"),
    planted_label = label, score = 1, fold = 10, stringsAsFactors = FALSE
  )
  right <- data.frame(
    program_id_a = "9001:1", program_id_b = "9002:1", link_tier = "primary",
    direction = "antagonistic", stringsAsFactors = FALSE
  )
  wrong <- right
  wrong$direction <- "concordant"
  expect_equal(
    evaluate_multitrait_simulation(right, mapping, sp$ground_truth)$sign_accuracy, 1
  )
  expect_equal(
    evaluate_multitrait_simulation(wrong, mapping, sp$ground_truth)$sign_accuracy, 0
  )
})

test_that("evaluate_multitrait_simulation handles an empty link table", {
  sp <- small_pair(locus_overlap = 0.5, profile_overlap = 0.8)
  ev <- evaluate_multitrait_simulation(
    NULL, data.frame(), sp$ground_truth
  )
  expect_equal(ev$n_links, 0L)
  expect_equal(ev$recall, 0)
  expect_equal(ev$n_planted_shared, 3L)
})
