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

test_that("unshared_driver_rows = 'background' puts every driver in both matrices", {
  sp <- small_pair(
    locus_overlap = 0.5, profile_overlap = 0.5, unshared_driver_rows = "background"
  )
  all_drivers <- unique(unlist(sp$ground_truth$program_drivers))
  for (tid in names(sp$x_matrices)) {
    expect_true(all(all_drivers %in% rownames(sp$x_matrices[[tid]])))
  }
  expect_equal(sp$ground_truth$parameters$unshared_driver_rows, "background")
})

test_that("a driver private to one trait carries no planted signal in the other", {
  sp <- small_pair(
    locus_overlap = 0.5, profile_overlap = 0.5, unshared_driver_rows = "background"
  )
  gt <- sp$ground_truth
  private_a <- setdiff(
    gt$program_drivers[["9001:shared1"]], gt$program_drivers[["9002:shared1"]]
  )
  expect_gt(length(private_a), 0)
  snps_b <- gt$program_snps[["9002:shared1"]]
  cells <- sp$x_matrices[["9002"]][private_a, snps_b]
  # An unstructured background row: sparse, not a dense block proportional to
  # the planted effect as it is in its own trait.
  expect_lt(mean(!is.na(cells)), 0.2)
  own <- sp$x_matrices[["9001"]][private_a, gt$program_snps[["9001:shared1"]]]
  expect_true(all(!is.na(own)))
})

test_that("the default unshared_driver_rows reproduces the historical simulation", {
  default <- small_pair(locus_overlap = 0.5, profile_overlap = 0.5, seed = 7)
  absent <- small_pair(
    locus_overlap = 0.5, profile_overlap = 0.5, seed = 7,
    unshared_driver_rows = "absent"
  )
  expect_identical(default$x_matrices, absent$x_matrices)
  expect_identical(default$traits, absent$traits)
  # A trait's matrix holds only its own drivers.
  drivers_b_only <- setdiff(
    default$ground_truth$program_drivers[["9002:shared1"]],
    default$ground_truth$program_drivers[["9001:shared1"]]
  )
  expect_false(any(drivers_b_only %in% rownames(default$x_matrices[["9001"]])))
})

test_that("two shared programs fit at zero locus overlap with enough private loci", {
  expect_no_error(simulate_trait_pair(
    n_loci_per_trait = 240L, n_shared_loci = 80L, K_shared = 2L,
    K_specific = 2L, module_size = 36L, locus_overlap = 0,
    profile_overlap = 0.5, directions = c(1, -1),
    unshared_driver_rows = "background", seed = 1
  ))
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

test_that("map_programs_to_planted maps programs on their high-confidence SNPs", {
  sp <- small_pair(locus_overlap = 0.5, profile_overlap = 0.5)
  truth <- sp$ground_truth$program_snps[["9001:shared1"]]
  decoys <- paste0("decoy", seq_along(truth))
  loadings <- data.frame(
    program_id = rep(c("9001:1", "9001:2"), each = length(truth) + 2),
    trait_id = "9001",
    snp_id = c(truth, "x1", "x2", decoys, "x3", "x4"),
    # The two extra SNPs of each program are low-confidence and ignored.
    high_confidence = rep(c(rep(TRUE, length(truth)), FALSE, FALSE), 2),
    stringsAsFactors = FALSE
  )
  mapped <- map_programs_to_planted(loadings, sp$ground_truth)
  expect_setequal(names(mapped), c("program_id", "trait_id", "planted_label", "score"))
  expect_equal(mapped$planted_label[mapped$program_id == "9001:1"], "shared1")
  expect_equal(mapped$score[mapped$program_id == "9001:1"], 1)
  expect_true(is.na(mapped$planted_label[mapped$program_id == "9001:2"]))
})
