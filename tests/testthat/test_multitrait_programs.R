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

# Trait profiles for programs meant to describe the same biology (their trait
# loadings agree); unrelated programs get random ones.
prof_shared <- c(2, -1.5, 1, 0.5, -0.8, 1.2, 0.2, 0.3, -0.4, 0.9)
prof_other <- c(-0.3, 0.4, 2, -1.2, 0.1, -0.9, 1.5, -0.6, 0.8, 0.2)
prof15 <- c(prof_shared, 1, -1, 0.5, 0.7, -0.2)

test_that("build_program_families groups linked programs and drops unlinked ones", {
  set.seed(6)
  loci <- paste0("g", 1:20)
  feats <- paste0("f", 1:15)
  hc <- rep(TRUE, 20)
  # The signal is planted on the locus axis, which is what the concordance
  # statistic reads: t1:1 and t2:1 load on the same loci, the others do not.
  shared <- c(rep(1.5, 10), rep(0, 10))
  programs <- list(
    make_program("t1", 1, loci, shared + stats::rnorm(20, 0, 0.05), hc,
                 feats, prof15 + stats::rnorm(15, 0, 0.05)),
    make_program("t1", 2, loci, stats::rnorm(20), hc, feats, stats::rnorm(15)),
    make_program("t2", 1, loci, shared + stats::rnorm(20, 0, 0.05), hc,
                 feats, prof15 + stats::rnorm(15, 0, 0.05)),
    make_program("t2", 2, loci, stats::rnorm(20), hc, feats, stats::rnorm(15))
  )
  res <- compare_program_pairs_loadings(programs, n_perm = 200, seed = 1)
  fam <- build_program_families(res)
  fam_of <- function(pid) fam$nodes$family[fam$nodes$program_id == pid]
  expect_setequal(fam$nodes$program_id[fam$nodes$family == fam_of("t1:1")], c("t1:1", "t2:1"))
  expect_false("t1:2" %in% fam$nodes$program_id[fam$nodes$family == fam_of("t1:1")])
  expect_equal(fam$n_unlinked, 4L - nrow(fam$nodes))
  expect_true(all(igraph::degree(fam$graph) > 0))
})

test_that("build_program_families handles no links", {
  set.seed(7)
  loci <- paste0("g", 1:10)
  feats <- paste0("f", 1:10)
  a <- make_program("t1", 1, loci, stats::rnorm(10), rep(TRUE, 10), feats, stats::rnorm(10))
  b <- make_program("t2", 1, loci, stats::rnorm(10), rep(TRUE, 10), feats, stats::rnorm(10))
  res <- compare_program_pairs_loadings(list(a, b), n_perm = 200, seed = 1)
  fam <- build_program_families(res)
  expect_null(fam$graph)
  expect_equal(nrow(fam$families), 0L)
  expect_equal(nrow(fam$nodes), 0L)
  expect_equal(fam$n_unlinked, 2L)
})

# Three traits, one shared module: every pair links, so the three programs are
# one family rather than three separate pairs.
triangle_fixture <- function(seed = 12) {
  set.seed(seed)
  loci <- paste0("g", 1:30)
  feats <- paste0("f", 1:10)
  shared <- c(rep(1.5, 12), rep(0, 18))
  claimed <- c(rep(TRUE, 12), rep(FALSE, 18))
  mk <- function(tid, prog, ll) {
    make_program(tid, prog, loci, ll, claimed, feats, prof_shared + stats::rnorm(10, 0, 0.05))
  }
  mk_other <- function(tid) {
    make_program(tid, 2, loci, other + stats::rnorm(30, 0, 0.05), other_claimed,
                 feats, prof_other + stats::rnorm(10, 0, 0.05))
  }
  bind <- function(p1, p2) {
    return(list(
      loadings = rbind(p1$loadings, p2$loadings),
      profiles = rbind(p1$profiles, p2$profiles)
    ))
  }
  other <- c(rep(0, 15), rep(1.5, 12), rep(0, 3))
  other_claimed <- c(rep(FALSE, 15), rep(TRUE, 12), rep(FALSE, 3))
  list(
    bind(mk("t1", 1, shared + stats::rnorm(30, 0, 0.05)), mk_other("t1")),
    mk("t2", 1, shared + stats::rnorm(30, 0, 0.05)),
    bind(mk("t3", 1, shared + stats::rnorm(30, 0, 0.05)), mk_other("t3"))
  )
}

test_that("a module shared by three traits forms one family of three", {
  pd <- triangle_fixture()
  res <- compare_program_pairs_loadings(pd, n_perm = 200, seed = 1)
  fam <- build_program_families(res, program_data = pd)
  fam_of <- function(pid) fam$nodes$family[fam$nodes$program_id == pid]
  expect_equal(fam_of("t1:1"), fam_of("t2:1"))
  expect_equal(fam_of("t1:1"), fam_of("t3:1"))
  row <- fam$families[fam$families$family == fam_of("t1:1"), ]
  expect_equal(row$n_programs, 3L)
  expect_equal(row$n_traits, 3L)
  # All three claim g1-g12.
  expect_equal(row$n_core_loci, 12L)
  expect_equal(row$n_union_loci, 12L)
  core <- fam$family_loci[fam$family_loci$family == fam_of("t1:1"), ]
  expect_setequal(core$locus, paste0("g", 1:12))
  expect_true(all(core$n_traits == 3L))
})

test_that("family locus sets separate core and union loci", {
  nodes <- data.frame(
    program_id = c("t1:1", "t2:1", "t3:1"), trait_id = c("t1", "t2", "t3"),
    trait_name = c("t1", "t2", "t3"), family = "family_1", stringsAsFactors = FALSE
  )
  mk <- function(pid, tid, loci) {
    data.frame(
      program_id = pid, trait_id = tid, trait_name = tid, program = 1L,
      snp_id = paste0("rs_", loci), coloc_group_id = loci, loading = 1,
      abs_loading = 1, lfsr = 0.01, high_confidence = TRUE, stringsAsFactors = FALSE
    )
  }
  loadings <- rbind(
    mk("t1:1", "t1", c("a", "b", "c")),
    mk("t2:1", "t2", c("a", "b", "d")),
    mk("t3:1", "t3", c("a", "e"))
  )
  fl <- gpmapr:::.family_loci(nodes, loadings)
  m <- gpmapr:::.family_metrics(nodes, loadings[0, ], fl)
  expect_equal(m$n_core_loci, 1L)   # a
  expect_equal(m$n_union_loci, 5L)  # a-e
})

test_that("mutual_best keeps only each program's closest partner per trait", {
  # t1:1 links to both t2:1 (same profile) and t2:2 (a similar profile, same
  # loci with uneven loadings): both are links, but only t2:1 is t1:1's
  # closest partner in t2.
  set.seed(21)
  loci <- paste0("g", 1:100)
  feats <- paste0("f", 1:10)
  core <- c(rep(1.5, 16), rep(0, 84))
  half <- core * rep(c(1.3, 0.7), 50)
  mk <- function(tid, prog, ll, prof) {
    make_program(tid, prog, loci, ll + stats::rnorm(100, 0, 0.05), abs(ll) > 0,
                 feats, prof)
  }
  t2 <- list(
    mk("t2", 1, core, prof_shared),
    mk("t2", 2, half, prof_shared + c(0.4, -0.3, 0.3, 0, 0.2, -0.2, 0, 0, 0.1, 0))
  )
  pd <- list(
    mk("t1", 1, core, prof_shared),
    list(loadings = rbind(t2[[1]]$loadings, t2[[2]]$loadings),
         profiles = rbind(t2[[1]]$profiles, t2[[2]]$profiles))
  )
  res <- compare_program_pairs_loadings(pd, n_perm = 200, seed = 1)
  expect_true(all(res$pairs$linked))
  phi <- stats::setNames(res$pairs$phi_traits, res$pairs$program_id_b)
  expect_gt(phi[["t2:1"]], phi[["t2:2"]])

  all_links <- build_program_families(res)
  expect_equal(nrow(all_links$edges), 2L)
  best <- build_program_families(res, mutual_best = TRUE)
  expect_equal(nrow(best$edges), 1L)
  expect_equal(best$edges$program_id_b, "t2:1")
  expect_false("t2:2" %in% best$nodes$program_id)
})

test_that("an upload's loci are put on the shared id space through its bridge", {
  guid <- "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
  bridge <- data.frame(
    trait_a = guid, trait_b = "1992",
    coloc_group_id_a = c("u1", "u2"), coloc_group_id_b = c("e1", "e2"),
    r = 1, stringsAsFactors = FALSE
  )
  bridges <- stats::setNames(list(bridge), gpmapr:::.locus_bridge_key(guid, "1992"))
  expect_equal(
    gpmapr:::.canonical_locus_ids(c("u1", "u2", "u9"), guid, bridges),
    c("e1", "e2", paste0(guid, ":u9"))
  )
  expect_equal(gpmapr:::.canonical_locus_ids(c("e1", "e5"), "1992", bridges), c("e1", "e5"))

  nodes <- data.frame(
    program_id = c(paste0(guid, ":1"), "1992:1"), trait_id = c(guid, "1992"),
    trait_name = c("up", "ex"), family = "family_1", stringsAsFactors = FALSE
  )
  mk <- function(pid, tid, loci) {
    data.frame(
      program_id = pid, trait_id = tid, trait_name = tid, program = 1L,
      snp_id = loci, coloc_group_id = loci, loading = 1, abs_loading = 1,
      lfsr = 0.01, high_confidence = TRUE, stringsAsFactors = FALSE
    )
  }
  loadings <- rbind(mk(paste0(guid, ":1"), guid, c("u1", "u2", "u9")),
                    mk("1992:1", "1992", c("e1", "e2", "e5")))
  fl <- gpmapr:::.family_loci(nodes, loadings, bridges)
  m <- gpmapr:::.family_metrics(nodes, loadings[0, ], fl)
  expect_equal(m$n_core_loci, 2L)   # e1, e2 (u1, u2 translated)
  expect_equal(m$n_union_loci, 4L)  # e1, e2, e5, upload-only u9
})

test_that("family_module_rg pools a family's loci and reads the direction", {
  pd <- triangle_fixture()
  res <- compare_program_pairs_loadings(pd, n_perm = 200, seed = 1)
  fam <- build_program_families(res, program_data = pd)
  set.seed(3)
  loci <- paste0("g", 1:30)
  effect <- stats::rnorm(30, 0, 0.2)
  cg <- list(
    t1 = data.frame(coloc_group_id = loci, variant_id = loci, trait_id = "t1",
                    beta = effect, se = 0.05, stringsAsFactors = FALSE),
    t2 = data.frame(coloc_group_id = loci, variant_id = loci, trait_id = "t2",
                    beta = effect + stats::rnorm(30, 0, 0.02), se = 0.05,
                    stringsAsFactors = FALSE),
    t3 = data.frame(coloc_group_id = loci, variant_id = loci, trait_id = "t3",
                    beta = -effect + stats::rnorm(30, 0, 0.02), se = 0.05,
                    stringsAsFactors = FALSE)
  )
  frg <- family_module_rg(fam, pd, cg)
  fam_id <- fam$nodes$family[fam$nodes$program_id == "t1:1"]
  frg <- frg[frg$family == fam_id, ]
  expect_equal(nrow(frg), 3L)
  dir <- stats::setNames(frg$rg_direction, paste(frg$trait_id_a, frg$trait_id_b))
  expect_equal(unname(dir["t1 t2"]), "concordant")
  expect_equal(unname(dir["t1 t3"]), "antagonistic")
  expect_equal(unname(dir["t2 t3"]), "antagonistic")
  expect_true(all(frg$n_loci_rg == 12L))
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
  valid_ids <- summary$programs$program[
    as.character(summary$programs$confidence_tier) %in% c("high", "medium")
  ]
  expect_setequal(unique(ex$loadings$program), valid_ids)
  expect_setequal(unique(ex$profiles$program), valid_ids)
  expect_true(all(abs(ex$loadings$loading) == ex$loadings$abs_loading))
  expect_true(all(ex$programs$n_high_confidence >= 0))
  expect_true(all(ex$programs$n_features > 0))
})
