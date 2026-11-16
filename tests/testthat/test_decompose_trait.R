library(testthat)

test_that("decompose_trait returns a classed, annotated decomposition", {
  d <- simulated_decompositions()[[1]]
  expect_s3_class(d, "gpmap_decomposition")
  expect_true(all(c(
    "trait", "programs", "memberships", "loadings", "profiles", "extraction", "annotation",
    "excluded", "exclude_trait_ids", "fit", "summary", "coloc_groups", "settings"
  ) %in% names(d)))
  expect_equal(d$trait$id, "9001")
  expect_s3_class(d$settings, "gpmap_settings")
  expect_true(all(c("program_id", "trait_id", "trait_name", "confidence_tier") %in% names(d$programs)))
  expect_true(all(startsWith(d$programs$program_id, "9001:")))
  expect_true(all(as.character(d$memberships$program) %in% as.character(d$programs$program)))
  expect_setequal(names(d$annotation), c(
    "specific_traits", "categories", "tissues", "pathways", "pathways_status",
    "distinctness", "genes"
  ))
  expect_equal(d$annotation$pathways_status, "ok")
  expect_identical(summary(d), d$programs)
})

test_that("programs(), annotation() and program_cards() return tidy tables", {
  d <- simulated_decompositions()[[1]]
  kept <- programs(d)
  expect_true(all(kept$confidence_tier %in% c("high", "medium")))
  expect_equal(nrow(programs(d, tiers = NULL)), nrow(d$programs))
  cats <- annotation(d, "categories")
  expect_true(all(c("program", "trait_category", "fdr", "r_squared", "significant") %in% names(cats)))
  expect_true(all(cats$program %in% as.character(d$programs$program)))
  expect_true(all(c("program", "trait_name", "specificity") %in% names(annotation(d, "specific_traits"))))
  expect_true(all(c("program_a", "program_b", "phi") %in% names(annotation(d, "distinctness"))))
  cards <- program_cards(d)
  expect_equal(nrow(cards), nrow(kept))
  expect_true(all(c(
    "program_id", "program_label", "confidence_tier", "specific_traits", "categories",
    "tissues", "reactome", "kegg", "genes"
  ) %in% names(cards)))
  expect_true(all(cards$categories != "not annotated"))
})

test_that("decompose_trait removes proxies, hand additions and background duplicates", {
  sim <- simulated_pair()
  obj <- sim$traits[[1]]
  counts <- table(as.character(obj$coloc_groups$trait_id))
  background <- setdiff(names(sort(counts, decreasing = TRUE)), "9001")[1:6]
  dups <- data.frame(
    trait_id = c("9001", background[1], background[2], background[3], background[4]),
    trait_name = c("anchor", "dup of anchor", "proxy of anchor", "bg", "bg copy"),
    duplicate_group = c(1L, 1L, 2L, 3L, 3L),
    is_representative = c(TRUE, FALSE, TRUE, TRUE, FALSE),
    match_rules = c(NA, "synonym", NA, NA, "ukb_field"),
    proxy_group_id = c(5L, 5L, 5L, NA, NA),
    proxy_group_name = c("anchor", "anchor", "anchor", NA, NA),
    proxy_rules = c(NA, NA, "icd10", NA, NA),
    stringsAsFactors = FALSE
  )
  d <- decompose_trait(
    obj, settings = fast_settings(), duplicates = dups, exclude_trait_ids = background[5],
    annotate = FALSE, cores = 1, verbose = FALSE
  )
  ex <- excluded(d)
  expect_equal(ex$source[ex$trait_id == background[1]], "duplicate_of_anchor")
  expect_equal(ex$source[ex$trait_id == background[2]], "proxy_group")
  expect_equal(ex$source[ex$trait_id == background[5]], "user")
  expect_equal(ex$source[ex$trait_id == background[4]], "background_duplicate")
  expect_false(background[3] %in% ex$trait_id)
  expect_setequal(d$exclude_trait_ids, background[c(1, 2, 5)])
  expect_true(all(ex$n_loci[ex$source != "background_duplicate"] > 0))
  fitted <- as.character(d$fit$trait_info$trait_id)
  expect_false(any(background[c(1, 2, 4, 5)] %in% fitted))
  expect_null(d$annotation)
  expect_equal(program_cards(d)$categories[1], "not annotated")
})

test_that("an uploaded GWAS without proxy_trait_id warns, and with it lists its trait", {
  guid <- "c6d2e20f-6ece-d093-78ee-37534547c466"
  cg <- data.frame(
    coloc_group_id = 1:2, trait_id = c(guid, "101"), trait_name = c("upload", "Anchor (UKB)"),
    beta = c(1, 1), se = c(0.1, 0.1), min_p = c(1e-8, 1e-8), stringsAsFactors = FALSE
  )
  expect_warning(
    ex <- .anchor_exclusions(guid, cg, fake_duplicates()),
    "proxy_trait_id"
  )
  expect_equal(nrow(ex), 0L)
  ex <- .anchor_exclusions(guid, cg, fake_duplicates(), proxy_trait_id = "100", exclude_trait_ids = "300")
  expect_setequal(ex$trait_id, c("100", "101", "102", "103", "104", "300"))
  expect_equal(ex$source[ex$trait_id == "100"], "anchor_trait")
  expect_true(all(ex$anchor_id == guid))
})

test_that("print and plots of a decomposition", {
  d <- simulated_decompositions()[[1]]
  expect_snapshot(print(d))
  for (type in c("programs", "similarity", "enrichment")) {
    expect_s3_class(plot(d, type = type), "ggplot")
  }
  expect_s3_class(plot_program_profile(d, programs(d)$program[1]), "ggplot")
  expect_s3_class(plot_program_enrichment(d, "pathways"), "ggplot")
  expect_error(plot_program_profile(d, 999), "no extracted profile")
  expect_error(plot_families(d), "trait_correspondence")
})

test_that("results do not depend on cores or the caller's random state", {
  skip_on_cran()
  skip_on_os("windows")
  obj <- simulated_pair()$traits[[1]]
  run <- function(cores, seed) {
    set.seed(seed)
    return(decompose_trait(
      obj, settings = fast_settings(n_rep = 2), proxies = FALSE, remove_duplicates = FALSE,
      annotate = FALSE, cores = cores, verbose = FALSE
    ))
  }
  a <- run(1, 1)
  b <- run(2, 99)
  expect_equal(a$programs, b$programs)
  expect_equal(a$loadings, b$loadings)
})
