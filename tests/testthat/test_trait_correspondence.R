library(testthat)

test_that("trait_correspondence returns classed tables at every level", {
  x <- simulated_correspondence()
  expect_s3_class(x, "gpmap_correspondence")
  expect_equal(x$traits$trait_id, c("9001", "9002"))
  expect_equal(nrow(x$trait_pairs), 1L)
  expect_true(x$trait_pairs$n_loci_shared > 0)
  expect_equal(
    sum(x$trait_pairs[, c("n_high", "n_medium", "n_low", "n_none")]), nrow(x$pairs)
  )
  expect_true(all(c("rg", "direction_check", "top_traits", "trait_name_a", "trait_name_b") %in% names(x$pairs)))
  expect_true(nrow(families(x)) >= 1)
  expect_true(all(c("family", "locus") %in% names(family_loci(x))))
  expect_equal(nrow(family_loci(x, family = "family_1")), sum(x$family_loci$family == "family_1"))
  expect_true(all(x$loci$trait_id %in% c("9001", "9002")))
  expect_identical(summary(x), x$trait_pairs)
  p <- programs(x)
  expect_true("family" %in% names(p))
  expect_setequal(unique(p$trait_id), c("9001", "9002"))
  expect_equal(nrow(excluded(x)), 0L)
})

test_that("links() filters and labels the scored pairs", {
  x <- simulated_correspondence()
  l <- links(x)
  expect_true(nrow(l) >= 1)
  expect_true(all(l$tier %in% c("high", "medium")))
  expect_true(all(l$linked))
  expect_equal(names(l)[1:2], c("label_a", "label_b"))
  expect_true(all(grepl("^Simulated trait 900[12] P[0-9]+$", c(l$label_a, l$label_b))))
  expect_equal(nrow(links(x, tiers = NULL)), nrow(x$pairs))
})

test_that("the planted shared programs are linked", {
  x <- simulated_correspondence()
  truth <- simulated_pair()$ground_truth$correspondence
  expect_true(nrow(links(x)) >= 1)
  expect_true(nrow(links(x)) <= nrow(truth) + 1)
})

test_that("family enrichment is reported per family and trait", {
  x <- simulated_correspondence()
  s <- annotation(x, "summary")
  expect_equal(nrow(s), nrow(x$nodes[!duplicated(x$nodes[, c("family", "trait_id")]), ]))
  expect_true(all(c("family", "trait_name", "programs", "n_snps", "categories", "tissues",
                    "reactome", "kegg") %in% names(s)))
  expect_true(all(c("family", "trait_id", "trait_category", "fdr") %in% names(annotation(x, "categories"))))
})

test_that("relink re-tiers without rescoring and rejects settings that need it", {
  x <- simulated_correspondence()
  same <- relink(x)
  expect_equal(same$pairs, x$pairs)
  expect_equal(same$families, x$families)
  loose <- relink(x, list(phi_high = 0.3, phi_low = 0.2))
  expect_true(sum(loose$pairs$linked) >= sum(x$pairs$linked))
  expect_equal(loose$settings$phi_high, 0.3)
  expect_equal(loose$pairs$p_concordance, x$pairs$p_concordance)
  expect_error(relink(x, list(n_perm = 10)), "trait_correspondence\\(x\\$decompositions")
  expect_error(relink(x, list(ebmf_prior = "point_normal")), "decompose_trait")
  mutual <- relink(x, gpmap_settings(unclass(x$settings), mutual_best = TRUE))
  expect_true(mutual$settings$mutual_best)
})

test_that("trait_correspondence checks its inputs", {
  ds <- simulated_decompositions()
  expect_error(trait_correspondence(ds[1], verbose = FALSE), "two traits")
  expect_error(trait_correspondence(list(ds[[1]], ds[[1]]), verbose = FALSE), "only once")
  other <- ds[[2]]
  other$settings$ebmf_prior <- "point_normal"
  expect_error(
    trait_correspondence(list(ds[[1]], other), settings = fast_settings(), annotate = FALSE, verbose = FALSE),
    "different EBMF settings"
  )
})

test_that("a permutation floor close to the FDR threshold warns", {
  ds <- simulated_decompositions()
  expect_warning(
    x <- trait_correspondence(ds, settings = fast_settings(n_perm = 10), annotate = FALSE, verbose = FALSE),
    "increase n_perm"
  )
  expect_true(is.finite(x$q_floor))
  expect_output(print(x), "permutation floor")
})

test_that("print and plots of a correspondence", {
  x <- simulated_correspondence()
  expect_snapshot(print(x))
  for (type in c("families", "phi")) {
    expect_s3_class(plot(x, type = type), "ggplot")
  }
  l <- links(x)
  expect_s3_class(plot_link(x, l$program_id_a[1], l$program_id_b[1]), "ggplot")
  expect_s3_class(plot_family_loci(x, "family_1"), "ggplot")
  expect_s3_class(plot_program_profile(x, l$program_id_a[1]), "ggplot")
  expect_s3_class(program_cards(x), "data.frame")
  expect_error(plot_programs(x), "decompose_trait")
})
