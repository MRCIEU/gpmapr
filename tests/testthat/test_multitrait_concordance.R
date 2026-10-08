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


# A trait profile for programs meant to describe the same biology: their trait
# loadings agree, so their profile congruence is ~1. Unrelated programs get
# random profiles.
prof_shared <- c(2, -1.5, 1, 0.5, -0.8, 1.2, 0.2, 0.3, -0.4, 0.9)


test_that("a planted concordant module is significant and target-aligned", {
  set.seed(1)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 10), rep(0, 10))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, prof_shared + stats::rnorm(10, 0, 0.05))
  b <- make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, prof_shared + stats::rnorm(10, 0, 0.05))
  res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)

  expect_equal(nrow(res$pairs), 1L)
  expect_equal(res$pairs$n_loci_axis, 20L)
  expect_true(res$pairs$locus_concordance > 0)
  expect_true(res$pairs$concordance_z > 0)
  expect_true(res$pairs$eligible)
  expect_equal(res$pairs$locus_status, "agree")
  expect_gt(res$pairs$phi_traits, 0.95)
  expect_equal(res$pairs$profile_strength, "strong")
  expect_equal(res$pairs$tier, "high")
  expect_equal(res$pairs$direction_profile, "concordant")
  expect_true(res$pairs$linked)
  expect_false(any(c("direction", "link_tier", "shared", "r2_loci", "significant",
                     "loci_ok", "evidence") %in% names(res$pairs)))
  expect_equal(res$programs$n_links, c(1L, 1L))
})


test_that("shared loci with unrelated profiles are low tier, not a link", {
  # The locus test is two-sided: opposite-signed SNP loadings at the same loci
  # are as much shared architecture as same-signed ones. But the two programs'
  # trait profiles are unrelated, so it is shared loci, different biology.
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
  expect_true(res$pairs$concordance_z < 0)
  expect_equal(res$pairs$locus_status, "agree")
  expect_lt(abs(res$pairs$phi_traits), 0.7)
  expect_equal(res$pairs$profile_strength, "weak")
  expect_equal(res$pairs$tier, "low")
  expect_false(res$pairs$linked)
})


test_that("correlated near-zero tails on unclaimed loci do not link", {
  # Each program claims only its own loci, which the other trait does not
  # carry, but both carry tiny, correlated loadings across a shared block --
  # as every program does through the target row. Over the whole shared axis a
  # shuffle null would call that significant; scored on claimed loci only,
  # there is nothing to score.
  set.seed(5)
  shared_loci <- paste0("s", 1:30)
  tail_pattern <- stats::rnorm(30)
  feats <- paste0("f", 1:10)
  one <- function(trait_id, own) {
    loci <- c(own, shared_loci)
    loading <- c(rep(1.5, 10), 1e-3 * (tail_pattern + stats::rnorm(30, 0, 0.1)))
    return(make_program(
      trait_id, 1, loci, loading, c(rep(TRUE, 10), rep(FALSE, 30)),
      feats, stats::rnorm(10), lfsr = c(rep(0.001, 10), rep(0.6, 30))
    ))
  }
  a <- one("t1", paste0("a", 1:10))
  b <- one("t2", paste0("b", 1:10))
  res <- compare_program_pairs_loadings(
    list(a, b), n_perm = 200, seed = 1, min_shared_loci = 0L
  )

  expect_equal(res$pairs$n_loci_axis, 0L)
  expect_true(is.na(res$pairs$concordance_z))
  expect_equal(res$pairs$locus_status, "too_few_loci")
})


test_that("two programs claiming exactly the same loci still link", {
  # Uniform loadings on an identical claimed set: nothing to distinguish within
  # the claimed loci, so the null must draw from the whole shared universe.
  set.seed(1)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 10), rep(0, 10))
  claimed <- c(rep(TRUE, 10), rep(FALSE, 10))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), claimed,
                    feats, prof_shared + stats::rnorm(10, 0, 0.05))
  b <- make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), claimed,
                    feats, prof_shared + stats::rnorm(10, 0, 0.05))
  res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)

  expect_equal(res$pairs$n_loci_axis, 10L)
  expect_gt(res$pairs$concordance_z, 3)
  expect_equal(res$pairs$locus_status, "agree")
  expect_true(res$pairs$linked)
})


test_that("a pair below min_shared_loci is scored but cannot link", {
  set.seed(1)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 10), rep(0, 10))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))
  # Only two of the co-loading loci are high-confidence in the second program.
  b <- make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05),
                    c(TRUE, TRUE, rep(FALSE, 18)), feats, stats::rnorm(10))

  gated <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)
  expect_equal(gated$pairs$n_loci_shared, 2L)
  expect_true(gated$pairs$concordance_z > 0)
  expect_false(gated$pairs$eligible)
  expect_equal(gated$pairs$locus_status, "too_few_loci")
  # Still part of the multiple-testing correction: eligibility is not
  # independent of the statistic, so it must not choose the tests corrected.
  expect_true(is.finite(gated$pairs$q_concordance))

  open <- compare_program_pairs_loadings(
    list(a, b), n_perm = 200, seed = 1, min_shared_loci = 0L
  )
  expect_equal(open$pairs$locus_status, "agree")
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
  # Both pairs share all of their loading, so both are fully explained.
  expect_gt(res$pairs$alignment[broad_row], 0.6)
  expect_gt(res$pairs$alignment[spiky_row], 0.8)
})


test_that("alignment compares the shared loci; comparability says how much is shared", {
  set.seed(31)
  n <- 60
  loci <- paste0("g", seq_len(n))
  feats <- paste0("f", 1:10)
  # Loadings 1.5 on the claimed loci, near zero elsewhere.
  mk <- function(tid, claimed, extra = character(0)) {
    ll <- ifelse(seq_len(n) %in% claimed, 1.5, 0) + stats::rnorm(n, 0, 0.02)
    p <- make_program(tid, 1, loci, ll, seq_len(n) %in% claimed, feats, prof_shared)
    if (length(extra) > 0) {
      # Loci the other trait does not carry at all.
      p2 <- make_program(tid, 1, extra, rep(1.5, length(extra)), rep(TRUE, length(extra)),
                         feats, prof_shared)
      p$loadings <- rbind(p$loadings, p2$loadings)
    }
    return(p)
  }
  pairs_of <- function(a, b, ...) {
    res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1, ...)
    return(res$pairs)
  }
  same <- pairs_of(mk("t1", 1:20), mk("t2", 1:20))
  expect_gt(same$alignment, 0.95)
  expect_gt(same$comparability, 0.95)
  expect_equal(same$tier, "high")

  # 5 loci nested in a 30-locus program: enough shared loci to test, but the
  # loadings disagree beyond those 5 (alignment about 5/30). A strong profile
  # whose loci were testable and disagree is demoted to low, not linked.
  nested <- pairs_of(mk("t1", 1:30), mk("t2", 1:5))
  expect_equal(nested$alignment, 5 / 30, tolerance = 0.1)
  expect_equal(nested$profile_strength, "strong")
  expect_equal(nested$locus_status, "disagree")
  expect_equal(nested$tier, "low")
  expect_false(nested$linked)

  # Same 10 comparable loci, but trait 1's program has 40 more loci trait 2
  # does not carry: the loci that can be compared agree fully, and
  # comparability records that only a fifth of trait 1's program is visible.
  private <- pairs_of(mk("t1", 1:10, extra = paste0("own", 1:40)), mk("t2", 1:10))
  expect_gt(private$alignment, 0.95)
  expect_equal(private$comparability, 10 / 50, tolerance = 0.1)
})


test_that("the lfsr weight is (1 - 2 lfsr)_+ by default, 1 - lfsr as an option", {
  with_weight <- function(method, code) {
    old <- options(gpmapr.lfsr_weight = method)
    on.exit(options(old))
    return(code)
  }
  expect_equal(gpmapr:::.confidence_factor(c(0, 0.25, 0.5, 0.8, NA)), c(1, 0.5, 0, 0, 1))
  expect_equal(
    with_weight("one_sided", gpmapr:::.confidence_factor(c(0, 0.25, 0.5, NA))),
    c(1, 0.75, 0.5, 1)
  )
  expect_error(with_weight("bogus", gpmapr:::.confidence_factor(0.1)))

  # Unrelated, uncertain tails (lfsr 0.4) dilute alignment less when a
  # coin-flip-ish sign gets little weight.
  set.seed(32)
  n <- 200
  loci <- paste0("g", seq_len(n))
  feats <- paste0("f", 1:10)
  mk <- function(tid) {
    ll <- ifelse(seq_len(n) <= 20, 1.5, stats::rnorm(n, 0, 0.25))
    p <- make_program(tid, 1, loci, ll, seq_len(n) <= 20, feats, prof_shared)
    p$loadings$lfsr <- ifelse(seq_len(n) <= 20, 0.001, 0.4)
    return(p)
  }
  pd <- list(mk("t1"), mk("t2"))
  two <- compare_program_pairs_loadings(pd, n_perm = 50, seed = 1)$pairs$alignment
  one <- with_weight(
    "one_sided", compare_program_pairs_loadings(pd, n_perm = 50, seed = 1)$pairs$alignment
  )
  expect_gt(two, one)
})


test_that("profile congruence is signed, and reports how many traits carry it", {
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  ll <- c(rep(1.5, 10), rep(0, 10))
  a <- make_program("t1", 1, loci, ll, ll > 0, feats, prof_shared)
  b_same <- make_program("t2", 1, loci, ll, ll > 0, feats, prof_shared)
  b_mirror <- make_program("t2", 1, loci, ll, ll > 0, feats, -prof_shared)
  same <- compare_program_pairs_loadings(list(a, b_same), n_perm = 50, seed = 1)$pairs
  mirror <- compare_program_pairs_loadings(list(a, b_mirror), n_perm = 50, seed = 1)$pairs
  expect_equal(same$phi_traits, 1, tolerance = 1e-8)
  expect_equal(same$direction_profile, "concordant")
  expect_equal(mirror$phi_traits, -1, tolerance = 1e-8)
  expect_equal(mirror$direction_profile, "antagonistic")
  expect_equal(mirror$tier, "high")
  expect_equal(same$n_traits_shared, 10L)

  # One trait dominates: congruence still high, but n_eff_traits near 1.
  dom <- c(10, rep(0.01, 9))
  a_dom <- make_program("t1", 1, loci, ll, ll > 0, feats, dom)
  b_dom <- make_program("t2", 1, loci, ll, ll > 0, feats, dom + c(0, stats::rnorm(9, 0, 0.01)))
  d <- compare_program_pairs_loadings(list(a_dom, b_dom), n_perm = 50, seed = 1)$pairs
  expect_gt(d$phi_traits, 0.99)
  expect_lt(d$n_eff_traits, 1.2)
})


test_that("anchor traits and excluded proxies are left out of the profile comparison", {
  loci <- paste0("g", 1:20)
  ll <- c(rep(1.5, 10), rep(0, 10))
  # The two programs agree only on the anchors' own rows and a proxy ("px");
  # on the remaining traits they are unrelated (orthogonal).
  feats <- c("t1", "t2", "px", paste0("f", 1:4))
  la <- c(5, 5, 5, 1, 0, 1, 0)
  lb <- c(5, 5, 5, 0, 1, 0, 1)
  a <- make_program("t1", 1, loci, ll, ll > 0, feats, la)
  b <- make_program("t2", 1, loci, ll, ll > 0, feats, lb)
  naive <- sum(la * lb) / sqrt(sum(la^2) * sum(lb^2))
  expect_gt(naive, 0.9)
  res <- compare_program_pairs_loadings(
    list(a, b), n_perm = 50, seed = 1, exclude_features = list(t1 = "px")
  )$pairs
  expect_equal(res$n_traits_shared, 4L)
  expect_equal(res$phi_traits, 0, tolerance = 1e-8)
})


test_that("independent loadings yield no significant correspondence", {
  set.seed(4)
  loci <- paste0("g", 1:30)
  feats <- paste0("f", 1:10)
  a <- make_program("t1", 1, loci, stats::rnorm(30), rep(TRUE, 30), feats, stats::rnorm(10))
  b <- make_program("t2", 1, loci, stats::rnorm(30), rep(TRUE, 30), feats, stats::rnorm(10))
  res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)

  expect_false(any(res$pairs$locus_status == "agree"))
})


test_that("one program can link to several programs in another trait", {
  set.seed(8)
  loci <- paste0("g", 1:30)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 15), rep(0, 15))
  mk <- function(tid, prog, ll) {
    make_program(tid, prog, loci, ll, rep(TRUE, 30), feats, prof_shared + stats::rnorm(10, 0, 0.05))
  }
  # Two programs in t2 both track the same t1 program, as when EBMF splits one
  # module into two factors: both are links, not one best match.
  a <- mk("t1", 1, shared + stats::rnorm(30, 0, 0.05))
  b1 <- mk("t2", 1, shared + stats::rnorm(30, 0, 0.05))
  b2 <- mk("t2", 2, shared + stats::rnorm(30, 0, 0.20))
  program_data <- list(
    a,
    list(loadings = rbind(b1$loadings, b2$loadings),
         profiles = rbind(b1$profiles, b2$profiles))
  )
  res <- compare_program_pairs_loadings(program_data, n_perm = 200, seed = 1)

  expect_true(all(res$pairs$linked))
  expect_equal(res$programs$n_links[res$programs$program_id == "t1:1"], 2L)
  expect_true(all(res$pairs$locus_status == "agree"))
})


test_that("module_rg recovers the planted direction with a CI", {
  set.seed(5)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 12), rep(0, 8))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, prof_shared + stats::rnorm(10, 0, 0.05))
  b <- make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, prof_shared + stats::rnorm(10, 0, 0.05))

  effect <- stats::rnorm(20, 0, 0.2)
  cg <- list(
    t1 = make_coloc_groups("t1", loci, effect),
    t2 = make_coloc_groups("t2", loci, effect + stats::rnorm(20, 0, 0.02))
  )
  res <- compare_program_pairs_loadings(
    list(a, b), n_perm = 200, seed = 1, coloc_groups = cg
  )

  expect_equal(nrow(res$module_rg), 1L)
  expect_false("link_tier" %in% names(res$module_rg))
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

  # Half the loci agree and half oppose: an estimate near zero whose CI spans
  # zero, so no direction is called.
  cg_null <- cg
  cg_null$t2$beta <- cg$t1$beta * rep(c(1, -1), 10)
  rg_null <- module_rg(res, program_data = list(a, b), coloc_groups = cg_null)
  expect_true(is.finite(rg_null$rg))
  expect_true(rg_null$rg_ci_lower < 0 && rg_null$rg_ci_upper > 0)
  expect_equal(rg_null$rg_direction, "undetermined")

  # linked_only = TRUE skips a pair that is not linked.
  unlinked <- res$pairs
  unlinked$linked <- FALSE
  expect_equal(nrow(module_rg(unlinked, list(a, b), cg)), 0L)
  expect_equal(nrow(module_rg(unlinked, list(a, b), cg, linked_only = FALSE)), 1L)
})


test_that("tiers follow profile strength x locus status", {
  pairs <- data.frame(
    phi_traits = c(0.9, 0.9, 0.9, -0.8, -0.8, -0.8, 0.3, 0.3),
    n_loci_shared = c(10, 2, 10, 10, 2, 10, 10, 10),
    n_eff_loci = c(8, 2, 8, 8, 2, 4, 8, 8),
    alignment = c(0.9, 0.9, 0.1, 0.9, 0.9, 0.9, 0.9, 0.9),
    q_concordance = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.5)
  )
  out <- gpmapr:::.assign_link_tiers(pairs)
  expect_equal(
    out$profile_strength,
    c("strong", "strong", "strong", "moderate", "moderate", "moderate", "weak", "weak")
  )
  # Row 6: enough shared loci but n_eff_loci below 5, so still too few to test.
  expect_equal(
    out$locus_status,
    c("agree", "too_few_loci", "disagree", "agree", "too_few_loci", "too_few_loci",
      "agree", "disagree")
  )
  expect_equal(
    out$tier,
    c("high", "medium", "low", "medium", "low", "low", "low", "none")
  )
  expect_equal(out$linked, out$tier %in% c("high", "medium"))

  # A moderate profile whose loci disagree is no link at all.
  md <- gpmapr:::.assign_link_tiers(data.frame(
    phi_traits = 0.75, n_loci_shared = 10, n_eff_loci = 8,
    alignment = 0.1, q_concordance = 0.01
  ))
  expect_equal(md$tier, "none")
})


test_that("direction_check compares rg's call with the profile's sign", {
  expect_equal(
    gpmapr:::.direction_check(
      c("concordant", "concordant", "antagonistic", "concordant", "concordant"),
      c("concordant", "antagonistic", "antagonistic", "undetermined", NA)
    ),
    c("agree", "conflict", "agree", "no call", "no call")
  )

  set.seed(5)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 12), rep(0, 8))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, prof_shared + stats::rnorm(10, 0, 0.05))
  b <- make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, prof_shared + stats::rnorm(10, 0, 0.05))
  effect <- stats::rnorm(20, 0, 0.2)
  cg <- list(
    t1 = make_coloc_groups("t1", loci, effect),
    t2 = make_coloc_groups("t2", loci, effect + stats::rnorm(20, 0, 0.02))
  )
  res <- compare_program_pairs_loadings(
    list(a, b), n_perm = 200, seed = 1, coloc_groups = cg
  )
  expect_equal(res$pairs$direction_profile, "concordant")
  expect_equal(res$module_rg$direction_check, "agree")

  # The traits' own effects oppose while the profile says concordant.
  cg$t2$beta <- -cg$t2$beta
  expect_equal(module_rg(res, list(a, b), cg)$direction_check, "conflict")

  # Pairs without direction_profile get no direction_check column.
  bare <- res$pairs[, c("program_id_a", "program_id_b", "trait_id_a", "trait_id_b", "linked")]
  expect_false("direction_check" %in% names(module_rg(bare, list(a, b), cg)))
})


test_that("profile_contributions names the traits carrying phi", {
  loci <- paste0("g", 1:20)
  ll <- c(rep(1.5, 10), rep(0, 10))
  feats <- c("t1", "t2", "px", "dom", "f2", "f3", "f4")
  # "dom" carries most of the agreement; "f4" works against it; the anchors
  # and the excluded proxy "px" agree strongly but must not be counted.
  la <- c(5, 5, 5, 3, 1, 0.5, 1)
  lb <- c(5, 5, 5, 3, 1, 0.5, -0.5)
  a <- make_program("t1", 1, loci, ll, ll > 0, feats, la)
  b <- make_program("t2", 1, loci, ll, ll > 0, feats, lb)
  res <- compare_program_pairs_loadings(
    list(a, b), n_perm = 50, seed = 1, exclude_features = list(t1 = "px")
  )
  expect_equal(res$settings$exclude_features, list(t1 = "px"))

  contrib <- profile_contributions(res, list(a, b), top_n = 10)
  traits <- contrib$traits
  expect_setequal(traits$feature_trait_id, c("dom", "f2", "f3", "f4"))
  expect_equal(traits$feature_trait_id[1], "dom")
  expect_equal(sum(abs(traits$share)), 1, tolerance = 1e-8)
  expect_lt(traits$share[traits$feature_trait_id == "f4"], 0)
  expect_equal(nrow(contrib$pairs), 1L)
  expect_match(contrib$pairs$top_traits, "^dom 84%")

  top1 <- profile_contributions(res, list(a, b), top_n = 1)
  expect_equal(nrow(top1$traits), 1L)
})


test_that("module_rg does not depend on which allele each locus is coded on", {
  set.seed(5)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 12), rep(0, 8))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))
  b <- make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, stats::rnorm(10))
  effect <- stats::rnorm(20, 0.1, 0.2)
  cg <- list(
    t1 = make_coloc_groups("t1", loci, effect),
    t2 = make_coloc_groups("t2", loci, effect + stats::rnorm(20, 0, 0.05))
  )
  pairs <- data.frame(
    program_id_a = "t1:1", program_id_b = "t2:1", trait_id_a = "t1",
    trait_id_b = "t2", linked = TRUE, stringsAsFactors = FALSE
  )
  base <- module_rg(pairs, list(a, b), cg, n_boot = 0)
  # Recode half the loci on the other allele: both traits' effects flip.
  flip <- seq(1, 20, by = 2)
  cg$t1$beta[flip] <- -cg$t1$beta[flip]
  cg$t2$beta[flip] <- -cg$t2$beta[flip]
  recoded <- module_rg(pairs, list(a, b), cg, n_boot = 0)
  expect_equal(recoded$rg, base$rg)
})


test_that("build_program_families accepts the concordance result unchanged", {
  set.seed(6)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 10), rep(0, 10))
  a <- make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, prof_shared + stats::rnorm(10, 0, 0.05))
  b <- make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    feats, prof_shared + stats::rnorm(10, 0, 0.05))
  res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)

  fam <- build_program_families(res)
  expect_true(all(c("nodes", "families", "family_loci") %in% names(fam)))
  expect_equal(nrow(fam$nodes), 2L)
  expect_equal(dplyr::n_distinct(fam$nodes$family), 1L)
  expect_true(is.finite(fam$families$mean_phi[1]))
})

test_that("module_rg uses only the module's own loci", {
  # Two programs claiming 4 of the 20 shared loci. The old behaviour used all
  # 20 for every pair, so n_loci_rg was identical across pairs regardless of
  # what the programs actually claimed.
  loci <- paste0("cg", 1:20)
  mk <- function(tid, claimed) {
    data.frame(
      program_id = paste0(tid, ":1"), trait_id = tid,
      trait_name = paste0("t", tid), program = 1L,
      snp_id = loci, coloc_group_id = loci, chr = 1L, bp = seq_along(loci),
      loading = 1, abs_loading = 1, lfsr = 0.01,
      high_confidence = loci %in% claimed, stringsAsFactors = FALSE
    )
  }
  ld <- dplyr::bind_rows(mk("1", loci[1:4]), mk("2", loci[3:6]))
  axis <- gpmapr:::.module_locus_axis(ld, "1:1", "2:1", "1", "2")
  expect_setequal(axis, loci[1:6])
  expect_lt(length(axis), length(gpmapr:::.shared_locus_axis(ld, "1", "2")))
})


# --- LD-based locus bridging for uploaded traits -------------------------
#
# An uploaded GWAS assigns its own coloc_group_id values independently of
# GPMap's shared numbering (confirmed on real data: an uploaded T2D trait
# shared 1 coloc_group_id with an existing triglycerides trait by direct
# equality, ~40 by variant_id/LD). These tests use a literal GUID-format
# trait id so is_guid() -- and therefore .needs_locus_bridge() -- is TRUE,
# with deliberately disjoint coloc_group_id strings on each side, mirroring
# that real-world mismatch.

test_that(".needs_locus_bridge is TRUE only when a trait is an uploaded GUID", {
  guid <- "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
  expect_true(gpmapr:::.needs_locus_bridge(guid, "1992"))
  expect_true(gpmapr:::.needs_locus_bridge("1992", guid))
  expect_false(gpmapr:::.needs_locus_bridge("1992", "1993"))
})


test_that(".ld_best_match finds exact and LD-tier matches, dropping the rest", {
  fake_ld <- function(variant_ids) {
    # v2 (queried) is in LD with two targets; keep the higher-r one (21).
    # v3 is only in LD with 30, below the 0.5 threshold -> excluded.
    # v4 has no proxy information returned at all -> no match.
    list(lds = data.frame(
      lead_variant_id = c(2, 2, 3),
      proxy_variant_id = c(20, 21, 30),
      ld_block_id = c(1, 1, 2),
      r = c(0.8, 0.95, 0.3),
      stringsAsFactors = FALSE
    ))
  }
  m <- gpmapr:::.ld_best_match(
    query_variant_ids = c("1", "2", "3", "4"),
    target_variant_ids = c("1", "20", "21", "30"),
    min_r = 0.5, ld_proxies_fn = fake_ld
  )
  exact <- m[m$query_variant_id == "1", ]
  expect_equal(exact$target_variant_id, "1")
  expect_equal(exact$r, 1)
  ld2 <- m[m$query_variant_id == "2", ]
  expect_equal(nrow(ld2), 1L)
  expect_equal(ld2$target_variant_id, "21")
  expect_equal(ld2$r, 0.95)
  expect_equal(nrow(m[m$query_variant_id == "3", ]), 0L)
  expect_equal(nrow(m[m$query_variant_id == "4", ]), 0L)
})


test_that(".ld_best_match degrades to exact-match-only when ld_proxies_fn errors", {
  failing_ld <- function(variant_ids) stop("network down")
  m <- gpmapr:::.ld_best_match(
    query_variant_ids = c("1", "2"), target_variant_ids = c("1", "99"),
    min_r = 0.5, ld_proxies_fn = failing_ld
  )
  expect_equal(nrow(m), 1L)
  expect_equal(m$query_variant_id, "1")
  expect_equal(m$r, 1)
})


test_that(".build_locus_bridge translates coloc_group_id via variant_id then LD", {
  cg_a <- data.frame(
    trait_id = "upload", coloc_group_id = c("u1", "u2", "u3"),
    variant_id = c("100", "200", "300"), stringsAsFactors = FALSE
  )
  cg_b <- data.frame(
    trait_id = "existing", coloc_group_id = c("e1", "e2"),
    variant_id = c("100", "250"), stringsAsFactors = FALSE
  )
  fake_ld <- function(variant_ids) {
    list(lds = data.frame(
      lead_variant_id = "200", proxy_variant_id = "250", ld_block_id = 1, r = 0.9,
      stringsAsFactors = FALSE
    ))
  }
  bridge <- gpmapr:::.build_locus_bridge(
    cg_a, cg_b, "upload", "existing", min_r = 0.5, ld_proxies_fn = fake_ld
  )
  expect_equal(nrow(bridge), 2L)
  m1 <- bridge[bridge$coloc_group_id_a == "u1", ]
  expect_equal(m1$coloc_group_id_b, "e1")
  expect_equal(m1$r, 1)
  m2 <- bridge[bridge$coloc_group_id_a == "u2", ]
  expect_equal(m2$coloc_group_id_b, "e2")
  expect_equal(m2$r, 0.9)
  expect_false("u3" %in% bridge$coloc_group_id_a)
})


test_that(".translate_locus_ids and its positional variant respect direction and no-ops", {
  guid <- "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
  bridge <- data.frame(
    trait_a = guid, trait_b = "1992",
    coloc_group_id_a = c("u1", "u2"), coloc_group_id_b = c("e1", "e2"),
    r = c(1, 0.9), stringsAsFactors = FALSE
  )
  bridges <- stats::setNames(list(bridge), gpmapr:::.locus_bridge_key(guid, "1992"))

  expect_equal(
    sort(gpmapr:::.translate_locus_ids(c("u1", "u2", "u3"), guid, "1992", bridges)),
    c("e1", "e2")
  )
  expect_equal(
    gpmapr:::.translate_locus_ids_positional(c("u1", "u3", "u2"), guid, "1992", bridges),
    c("e1", NA_character_, "e2")
  )
  expect_equal(
    gpmapr:::.translate_locus_ids(c("e1", "e2"), "1992", guid, bridges), c("u1", "u2")
  )
  # No-ops: no bridges at all, or no entry for this specific pair.
  expect_equal(gpmapr:::.translate_locus_ids("u1", guid, "1992", NULL), "u1")
  expect_equal(gpmapr:::.translate_locus_ids("x1", guid, "1993", bridges), "x1")
})


test_that(".shared_locus_axis and .locus_pair_stats bridge a GUID trait's own id space", {
  guid <- "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
  loci_up <- paste0("u", 1:5)
  loci_ex <- paste0("e", 1:5)
  ld <- dplyr::bind_rows(
    data.frame(
      program_id = paste0(guid, ":1"), trait_id = guid, trait_name = guid,
      program = 1L, snp_id = loci_up, coloc_group_id = loci_up,
      chr = NA_character_, bp = NA_real_, loading = 1, abs_loading = 1,
      lfsr = 0.01, high_confidence = TRUE, stringsAsFactors = FALSE
    ),
    data.frame(
      program_id = "1992:1", trait_id = "1992", trait_name = "1992",
      program = 1L, snp_id = loci_ex, coloc_group_id = loci_ex,
      chr = NA_character_, bp = NA_real_, loading = 1, abs_loading = 1,
      lfsr = 0.01, high_confidence = TRUE, stringsAsFactors = FALSE
    )
  )
  # Naive intersection finds nothing -- the premise of the whole fix.
  expect_equal(gpmapr:::.shared_locus_axis(ld, guid, "1992"), character(0))

  bridge <- data.frame(
    trait_a = guid, trait_b = "1992",
    coloc_group_id_a = loci_up[1:3], coloc_group_id_b = loci_ex[1:3],
    r = c(1, 1, 0.8), stringsAsFactors = FALSE
  )
  bridges <- stats::setNames(list(bridge), gpmapr:::.locus_bridge_key(guid, "1992"))
  expect_equal(
    sort(gpmapr:::.shared_locus_axis(ld, guid, "1992", bridges = bridges)),
    sort(loci_ex[1:3])
  )

  locus_sets <- gpmapr:::.program_locus_sets(ld, "high_confidence", 25L)
  universe <- gpmapr:::.locus_universe_by_trait(ld)
  key <- data.frame(
    program_id_a = paste0(guid, ":1"), program_id_b = "1992:1",
    trait_id_a = guid, trait_id_b = "1992", program_a = 1L, program_b = 1L,
    stringsAsFactors = FALSE
  )
  unbridged <- gpmapr:::.locus_pair_stats(key, locus_sets, universe, bridges = NULL)
  expect_equal(unbridged$n_loci_shared, 0L)
  bridged <- gpmapr:::.locus_pair_stats(key, locus_sets, universe, bridges = bridges)
  expect_equal(bridged$n_loci_shared, 3L)
  # Native sizes are unaffected by the translation.
  expect_equal(bridged$n_loci_a, 5L)
  expect_equal(bridged$n_loci_b, 5L)
})


# An upload's coloc table carries the upload's own rows and, in the same coloc
# groups, the rows of every existing trait that colocalised with it -- at the
# upload's lead variant and on the upload's allele. GPMap's own trait tables
# can report the same variant on the other allele, so module_rg() reads both
# effects of an upload pair from the upload's table.
upload_coloc_table <- function(guid, partner, loci_up, beta_up, beta_partner) {
  return(data.frame(
    coloc_group_id = c(loci_up, loci_up),
    variant_id = paste0("v", c(seq_along(loci_up), seq_along(loci_up))),
    trait_id = c(rep(guid, length(loci_up)), rep(partner, length(loci_up))),
    beta = c(beta_up, beta_partner), se = 0.05, stringsAsFactors = FALSE
  ))
}


test_that("module_rg reads an upload pair's effects from the upload's own table", {
  guid <- "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
  loci_up <- paste0("u", 1:12)
  loci_ex <- paste0("e", 1:12)
  set.seed(9)
  shared <- rep(1.5, 12)
  a <- make_program(guid, 1, loci_up, shared + stats::rnorm(12, 0, 0.05), rep(TRUE, 12),
                    paste0("f", 1:10), stats::rnorm(10))
  b <- make_program("1992", 1, loci_ex, shared + stats::rnorm(12, 0, 0.05), rep(TRUE, 12),
                    paste0("f", 1:10), stats::rnorm(10))

  effect <- stats::rnorm(12, 0, 0.2)
  partner <- effect + stats::rnorm(12, 0, 0.02)
  cg <- list()
  cg[[guid]] <- upload_coloc_table(guid, "1992", loci_up, effect, partner)
  # The partner's own table: same variants, but half of them reported on the
  # other allele. Reading the partner's effect from here would scramble rg.
  cg[["1992"]] <- data.frame(
    coloc_group_id = loci_ex, variant_id = paste0("v", 1:12), trait_id = "1992",
    beta = partner * rep(c(1, -1), 6), se = 0.05, stringsAsFactors = FALSE
  )

  bridge <- data.frame(
    trait_a = guid, trait_b = "1992",
    coloc_group_id_a = loci_up, coloc_group_id_b = loci_ex,
    r = 1, stringsAsFactors = FALSE
  )
  bridges <- stats::setNames(list(bridge), gpmapr:::.locus_bridge_key(guid, "1992"))
  pairs <- data.frame(
    program_id_a = paste0(guid, ":1"), program_id_b = "1992:1",
    trait_id_a = guid, trait_id_b = "1992", linked = TRUE,
    stringsAsFactors = FALSE
  )

  rg <- module_rg(pairs, program_data = list(a, b), coloc_groups = cg, bridges = bridges)
  expect_equal(nrow(rg), 1L)
  expect_true(rg$rg > 0.9)
  expect_equal(rg$n_allele_mismatch, 0L)
  expect_equal(rg$n_loci_rg, 12L)

  # Same answer with the upload on the B side.
  pairs_rev <- data.frame(
    program_id_a = "1992:1", program_id_b = paste0(guid, ":1"),
    trait_id_a = "1992", trait_id_b = guid, linked = TRUE,
    stringsAsFactors = FALSE
  )
  rg_rev <- module_rg(pairs_rev, program_data = list(b, a), coloc_groups = cg, bridges = bridges)
  expect_equal(rg_rev$rg, rg$rg, tolerance = 1e-8)

  # An empty (attempted-but-found-nothing) bridge: no loci match, empty row.
  # Passed explicitly (rather than bridges = NULL) so this stays network-free
  # -- module_rg() only auto-builds a bridge when none was supplied at all.
  empty_bridge <- bridge[0, ]
  bridges_empty <- stats::setNames(list(empty_bridge), gpmapr:::.locus_bridge_key(guid, "1992"))
  rg_unbridged <- module_rg(pairs, program_data = list(a, b), coloc_groups = cg, bridges = bridges_empty)
  expect_true(is.na(rg_unbridged$rg))
  expect_equal(rg_unbridged$n_loci_rg, 0L)
})


test_that("compare_program_pairs_loadings links a GUID trait via variant_id, no LD call needed", {
  guid <- "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
  set.seed(11)
  loci_up <- paste0("u", 1:20)
  loci_ex <- paste0("e", 1:20)
  shared <- c(rep(1.5, 12), rep(0, 8))
  a <- make_program(guid, 1, loci_up, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    paste0("f", 1:10), prof_shared + stats::rnorm(10, 0, 0.05))
  b <- make_program("1992", 1, loci_ex, shared + stats::rnorm(20, 0, 0.05), rep(TRUE, 20),
                    paste0("f", 1:10), prof_shared + stats::rnorm(10, 0, 0.05))

  # Disjoint coloc_group_id strings (as an upload's own numbering is), but
  # the same variant_id on both sides -- the free, network-free tier.
  beta <- stats::rnorm(20, 0, 0.2)
  cg <- list()
  cg[[guid]] <- upload_coloc_table(
    guid, "1992", loci_up, beta, beta + stats::rnorm(20, 0, 0.02)
  )
  cg[["1992"]] <- data.frame(
    coloc_group_id = loci_ex, variant_id = paste0("v", 1:20), trait_id = "1992",
    beta = beta + stats::rnorm(20, 0, 0.02), se = 0.05, stringsAsFactors = FALSE
  )

  failing_ld <- function(variant_ids) stop("should not be called: all matches are exact")
  res <- compare_program_pairs_loadings(
    list(a, b), n_perm = 200, seed = 1, coloc_groups = cg, ld_proxies_fn = failing_ld
  )

  expect_equal(nrow(res$pairs), 1L)
  expect_equal(res$pairs$n_loci_shared, 20L)
  expect_gt(res$pairs$n_loci_axis, 0)
  expect_true(res$pairs$linked)
  expect_true(gpmapr:::.locus_bridge_key(guid, "1992") %in% names(res$bridges))
  expect_equal(nrow(res$module_rg), 1L)
  expect_true(res$module_rg$rg > 0.5)
})


test_that(".ld_best_match looks proxies up from the target side and keeps negative r", {
  # Only targets have proxy lists (as for GPMap's own leads); query "1" is
  # reachable only in reverse, query "2" only via a negative r.
  table <- data.frame(
    lead_variant_id = c(10, 20, 20),
    proxy_variant_id = c(1, 2, 3),
    ld_block_id = 1,
    r = c(0.95, -0.97, 0.4),
    stringsAsFactors = FALSE
  )
  fake_ld <- function(variant_ids) {
    list(lds = table[table$lead_variant_id %in% variant_ids, , drop = FALSE])
  }
  m <- gpmapr:::.ld_best_match(
    c("1", "2", "3"), c("10", "20"), min_r = 0.5, ld_proxies_fn = fake_ld
  )
  expect_equal(m$target_variant_id[m$query_variant_id == "1"], "10")
  expect_equal(m$r[m$query_variant_id == "1"], 0.95)
  expect_equal(m$target_variant_id[m$query_variant_id == "2"], "20")
  expect_equal(m$r[m$query_variant_id == "2"], -0.97)
  expect_false("3" %in% m$query_variant_id)
})


test_that(".build_locus_bridge matches an upload to an existing trait by its own coloc record", {
  guid <- "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
  # Upload loci u1-u3; the upload's coloc run recorded trait 1992 in u1 and u2
  # (existing_study_extraction_id 501, 502). u1's leads are in negative LD,
  # u2's leads have no proxy row, u3 is only reachable by exact variant.
  cg_up <- data.frame(
    trait_id = c(guid, guid, guid, "1992", "1992"),
    coloc_group_id = c("u1", "u2", "u3", "u1", "u2"),
    variant_id = c("100", "200", "300", "100", "200"),
    existing_study_extraction_id = c(NA, NA, NA, 501, 502),
    study_extraction_id = c(9001, 9001, 9001, NA, NA),
    stringsAsFactors = FALSE
  )
  cg_ex <- data.frame(
    trait_id = "1992",
    coloc_group_id = c("e1", "e2", "e3"),
    variant_id = c("110", "210", "300"),
    study_extraction_id = c(501, 502, 503),
    stringsAsFactors = FALSE
  )
  fake_ld <- function(variant_ids) {
    t <- data.frame(lead_variant_id = 110, proxy_variant_id = 100, ld_block_id = 1, r = -0.93)
    list(lds = t[t$lead_variant_id %in% variant_ids, , drop = FALSE])
  }
  bridge <- gpmapr:::.build_locus_bridge(cg_up, cg_ex, guid, "1992", ld_proxies_fn = fake_ld)
  expect_equal(nrow(bridge), 3L)
  expect_equal(bridge$coloc_group_id_b[bridge$coloc_group_id_a == "u1"], "e1")
  expect_equal(bridge$r[bridge$coloc_group_id_a == "u1"], -0.93)
  expect_equal(bridge$coloc_group_id_b[bridge$coloc_group_id_a == "u2"], "e2")
  expect_true(is.na(bridge$r[bridge$coloc_group_id_a == "u2"]))
  expect_equal(bridge$coloc_group_id_b[bridge$coloc_group_id_a == "u3"], "e3")
  expect_equal(bridge$r[bridge$coloc_group_id_a == "u3"], 1)

  # Same matches with the upload on the B side.
  flipped <- gpmapr:::.build_locus_bridge(cg_ex, cg_up, "1992", guid, ld_proxies_fn = fake_ld)
  expect_equal(
    flipped$coloc_group_id_a[match(c("u1", "u2", "u3"), flipped$coloc_group_id_b)],
    c("e1", "e2", "e3")
  )
})


test_that("module_rg drops upload loci the partner has no row for in the upload's table", {
  guid <- "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
  loci_up <- paste0("u", 1:16)
  loci_ex <- paste0("e", 1:16)
  set.seed(5)
  shared <- rep(1.5, 16)
  a <- make_program(guid, 1, loci_up, shared + stats::rnorm(16, 0, 0.05), rep(TRUE, 16),
                    paste0("f", 1:10), stats::rnorm(10))
  b <- make_program("1992", 1, loci_ex, shared + stats::rnorm(16, 0, 0.05), rep(TRUE, 16),
                    paste0("f", 1:10), stats::rnorm(10))
  effect <- stats::rnorm(16, 0, 0.2)
  cg <- list()
  cg[[guid]] <- upload_coloc_table(guid, "1992", loci_up, effect, effect)
  # The partner has its own effect at every locus, but no row in the upload's
  # table at the last four, so those cannot be put on the upload's allele.
  missing <- cg[[guid]]$trait_id == "1992" & cg[[guid]]$coloc_group_id %in% loci_up[13:16]
  cg[[guid]] <- cg[[guid]][!missing, ]
  cg[["1992"]] <- data.frame(
    coloc_group_id = loci_ex, variant_id = paste0("b", 1:16), trait_id = "1992",
    beta = effect, se = 0.05, stringsAsFactors = FALSE
  )
  pairs <- data.frame(
    program_id_a = paste0(guid, ":1"), program_id_b = "1992:1",
    trait_id_a = guid, trait_id_b = "1992", linked = TRUE,
    stringsAsFactors = FALSE
  )
  bridges <- stats::setNames(
    list(data.frame(
      trait_a = guid, trait_b = "1992",
      coloc_group_id_a = loci_up, coloc_group_id_b = loci_ex, r = NA_real_,
      stringsAsFactors = FALSE
    )),
    gpmapr:::.locus_bridge_key(guid, "1992")
  )

  res <- module_rg(pairs, list(a, b), cg, bridges = bridges)
  expect_equal(res$n_loci_rg, 12L)
  expect_equal(res$n_allele_mismatch, 4L)
  expect_true(res$rg > 0.9)
})
