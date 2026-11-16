library(testthat)

test_that("anchor_proxies finds the duplicate group and the proxy group", {
  p <- anchor_proxies("100", duplicates = fake_duplicates())
  expect_setequal(p$trait_id, c("101", "102", "103", "104"))
  expect_equal(p$source[p$trait_id == "101"], "duplicate_of_anchor")
  expect_equal(p$rule[p$trait_id == "102"], "firth_spa")
  expect_equal(p$source[p$trait_id == "104"], "proxy_group")
  expect_equal(p$rule[p$trait_id == "104"], "wrapper;concept_map")
  expect_true(all(p$anchor_id == "100"))
  expect_equal(nrow(anchor_proxies("200", duplicates = fake_duplicates())), 1L)
  expect_equal(nrow(anchor_proxies("999", duplicates = fake_duplicates())), 0L)
})

test_that("an uploaded GWAS lists the trait it measures as well as its groups", {
  p <- .proxy_table("100", fake_duplicates(), include_self = TRUE)
  expect_equal(p$source[p$trait_id == "100"], "anchor_trait")
  expect_equal(nrow(p), 5L)
})

test_that("the older duplicates shape falls back to the duplicate group", {
  old <- data.frame(
    trait_id = c(101L, 102L, 201L), trait_name = c("A (UKB)", "A (SPA)", "B copy"),
    parent_trait_id = c(100L, 100L, 200L), parent_trait_name = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  expect_setequal(anchor_proxies("100", duplicates = old)$trait_id, c("101", "102"))
  expect_setequal(anchor_proxies("101", duplicates = old)$trait_id, c("100", "102"))
  from_child <- anchor_proxies("101", duplicates = old)
  expect_equal(from_child$trait_name[from_child$trait_id == "100"], "A")
})

test_that("proxy evidence matches the multi-trait vignette's duplicate_evidence()", {
  cg <- data.frame(
    coloc_group_id = rep(1:4, each = 3),
    trait_id = rep(c("100", "101", "104"), times = 4),
    beta = c(1, 2, -0.5, -1, -1.5, 0.5, 2, 1, 1, 0.5, 0.4, NA),
    se = rep(0.1, 12),
    min_p = rep(1e-8, 12),
    stringsAsFactors = FALSE
  )
  p <- anchor_proxies("100", duplicates = fake_duplicates(), coloc_groups = cg)
  # 101: same sign at every anchor locus, z oriented = 20, 15, 10, 4 vs |z| 10, 10, 20, 5
  r101 <- p[p$trait_id == "101", ]
  expect_equal(r101$n_loci, 4L)
  expect_equal(r101$coverage, 1)
  expect_equal(r101$pct_same_sign, 100)
  z <- c(20, 15, 10, 4)
  za <- c(10, 10, 20, 5)
  expect_equal(r101$r_z, sum(z * za) / sqrt(sum(z^2) * sum(za^2)))
  r104 <- p[p$trait_id == "104", ]
  expect_equal(r104$n_loci, 3L)
  expect_equal(r104$coverage, 0.75)
  expect_equal(r104$pct_same_sign, 100 * 1 / 3)
  expect_equal(p$n_loci[p$trait_id == "102"], 0L)
})
