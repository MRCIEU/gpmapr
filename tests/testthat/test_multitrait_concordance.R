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

make_coloc_groups <- function(trait_id, loci, beta, se = 0.05) {
  return(data.frame(
    coloc_group_id = loci,
    variant_id = paste0("v", loci),
    trait_id = trait_id,
    beta = beta,
    se = se,
    stringsAsFactors = FALSE
  ))
}


test_that("a planted concordant module is significant and target-aligned", {
  set.seed(1)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 10), rep(0, 10))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))
  b <- make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))
  res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)

  expect_equal(nrow(res$pairs), 1L)
  expect_equal(res$pairs$n_loci_axis, 20L)
  expect_true(res$pairs$locus_concordance > 0)
  expect_true(res$pairs$concordance_z > 0)
  expect_equal(res$pairs$direction, "concordant")
  expect_true(res$pairs$shared)
  expect_true(all(res$matches$significant))
})


test_that("a planted antagonistic module is significant and opposed", {
  set.seed(2)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 10), rep(0, 10))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))
  b <- make_program("t2", 1, loci, -shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))
  res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)

  expect_true(res$pairs$locus_concordance < 0)
  expect_equal(res$pairs$direction, "antagonistic")
  expect_true(res$pairs$shared)
})


test_that("breadth of shared loading is not penalised", {
  # Regression test for the failure in the removed profile-correlation
  # statistic, where a 3-locus pair beat a 19-locus pair because correlating
  # background profiles rewards tightness rather than amount of shared signal.
  # Evidence must accumulate with the amount of shared loading instead.
  set.seed(3)
  loci <- paste0("g", 1:40)
  feats <- paste0("f", 1:10)

  broad <- c(rep(1, 20), rep(0, 20))
  a_broad <- make_program("t1", 1, loci, broad + stats::rnorm(40, 0, 0.3), rep(TRUE, 40),
                          feats, stats::rnorm(10))
  b_broad <- make_program("t2", 1, loci, broad + stats::rnorm(40, 0, 0.3), rep(TRUE, 40),
                          feats, stats::rnorm(10))

  # Same per-locus magnitude, far fewer loci.
  spike <- c(rep(0, 37), rep(1, 3))
  a_spiky <- make_program("t1", 2, loci, spike, rep(TRUE, 40), feats, stats::rnorm(10))
  b_spiky <- make_program("t2", 2, loci, spike, rep(TRUE, 40), feats, stats::rnorm(10))

  program_data <- list(
    list(loadings = rbind(a_broad$loadings, a_spiky$loadings),
         profiles = rbind(a_broad$profiles, a_spiky$profiles)),
    list(loadings = rbind(b_broad$loadings, b_spiky$loadings),
         profiles = rbind(b_broad$profiles, b_spiky$profiles))
  )
  res <- compare_program_pairs_loadings(program_data, n_perm = 200, seed = 1)

  broad_row <- res$pairs$program_id_a == "t1:1" & res$pairs$program_id_b == "t2:1"
  spiky_row <- res$pairs$program_id_a == "t1:2" & res$pairs$program_id_b == "t2:2"

  # At comparable per-locus magnitude the raw score tracks how much shared
  # loading there is, so the broad pair scores far higher.
  expect_true(
    res$pairs$locus_concordance[broad_row] > res$pairs$locus_concordance[spiky_row]
  )
  # The broad pair is detected on its own merits, not squeezed out by the decoy.
  expect_true(res$pairs$p_concordance[broad_row] < 0.05)
  # n_eff_loci reports how much independent support each pair rests on.
  expect_true(res$pairs$n_eff_loci[broad_row] > res$pairs$n_eff_loci[spiky_row])
  # The scale-free cosine cannot tell the two apart -- it prefers the tighter,
  # narrower pair, which is exactly why it is a diagnostic and not the decision
  # statistic.
  expect_true(res$pairs$cosine_loci[spiky_row] > res$pairs$cosine_loci[broad_row])
})


test_that("independent loadings yield no significant correspondence", {
  set.seed(4)
  loci <- paste0("g", 1:30)
  feats <- paste0("f", 1:10)
  a <- make_program("t1", 1, loci, stats::rnorm(30), rep(TRUE, 30), feats, stats::rnorm(10))
  b <- make_program("t2", 1, loci, stats::rnorm(30), rep(TRUE, 30), feats, stats::rnorm(10))
  res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)

  expect_false(any(res$pairs$shared))
  expect_false(any(res$matches$significant))
})


test_that("links are reported in primary and candidate tiers", {
  set.seed(8)
  loci <- paste0("g", 1:30)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 15), rep(0, 15))
  mk <- function(tid, prog, ll) {
    make_program(tid, prog, loci, ll, rep(TRUE, 30), feats, stats::rnorm(10))
  }
  # Two programs in t2 both track the same t1 program, so only one can be the
  # reciprocal best match; the other should surface as a candidate.
  a <- mk("t1", 1, shared + stats::rnorm(30, 0, 0.05))
  b1 <- mk("t2", 1, shared + stats::rnorm(30, 0, 0.05))
  b2 <- mk("t2", 2, shared + stats::rnorm(30, 0, 0.20))
  program_data <- list(
    a,
    list(loadings = rbind(b1$loadings, b2$loadings),
         profiles = rbind(b1$profiles, b2$profiles))
  )
  res <- compare_program_pairs_loadings(program_data, n_perm = 200, seed = 1)

  expect_true("link_tier" %in% names(res$pairs))
  expect_true(all(res$pairs$link_tier[res$pairs$shared] == "primary"))
  expect_true(any(res$pairs$link_tier == "candidate", na.rm = TRUE))
  # Every primary link is also a reciprocal best match, and vice versa.
  expect_equal(sum(res$pairs$link_tier == "primary", na.rm = TRUE), sum(res$pairs$shared))
  # A candidate never carries the shared flag.
  cand <- !is.na(res$pairs$link_tier) & res$pairs$link_tier == "candidate"
  expect_false(any(res$pairs$shared[cand]))
})


test_that("module_rg recovers the planted direction with a CI", {
  set.seed(5)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 12), rep(0, 8))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))
  b <- make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))

  effect <- stats::rnorm(20, 0, 0.2)
  cg <- list(
    t1 = make_coloc_groups("t1", loci, effect),
    t2 = make_coloc_groups("t2", loci, effect + stats::rnorm(20, 0, 0.02))
  )
  res <- compare_program_pairs_loadings(
    list(a, b), n_perm = 200, seed = 1, coloc_groups = cg
  )

  expect_equal(nrow(res$module_rg), 1L)
  expect_true(res$module_rg$rg > 0.5)
  expect_equal(res$module_rg$rg_direction, "concordant")
  expect_true(res$module_rg$rg_ci_lower > 0)
  expect_equal(res$module_rg$n_allele_mismatch, 0L)
  expect_true(res$module_rg$n_loci_rg > 10L)

  # Opposed effects in the second trait must flip the sign.
  cg_opposed <- cg
  cg_opposed$t2$beta <- -cg_opposed$t2$beta
  rg_opposed <- module_rg(res, program_data = list(a, b), coloc_groups = cg_opposed)
  expect_true(rg_opposed$rg < 0)
  expect_equal(rg_opposed$rg_direction, "antagonistic")
})


test_that("build_program_families accepts the concordance result unchanged", {
  set.seed(6)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 10), rep(0, 10))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))
  b <- make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))
  res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)

  fam <- build_program_families(res)
  expect_true(all(c("nodes", "families") %in% names(fam)))
  expect_equal(nrow(fam$nodes), 2L)
  expect_equal(dplyr::n_distinct(fam$nodes$family), 1L)
  # locus_concordance now populates the family metric that was silently NA.
  expect_true(is.finite(fam$families$mean_locus_concordance[1]))
})
