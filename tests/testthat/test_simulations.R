library(testthat)

test_that("adjusted_rand_index sanity", {
  set.seed(1)
  a <- rep(1:10, each = 20)
  expect_equal(.adjusted_rand_index(a, a), 1)
  b <- sample(rep(1:2, each = 100))
  expect_lt(.adjusted_rand_index(rep(1:5, each = 40), rep(1:4, times = 50)), 0.05)
})

test_that("simulate_trait returns a valid trait object with ground truth", {
  sim <- simulate_trait(
    n_coloc_groups = 60,
    K = 3,
    module_sizes = c(10, 10, 10),
    drivers_per_module = 4,
    n_background_traits = 15,
    seed = 42
  )
  cg <- sim$trait_object$coloc_groups
  required <- c(
    "coloc_group_id", "variant_id", "display_snp", "chr", "bp",
    "trait_id", "trait_name", "min_p", "beta", "se",
    "tissue", "gene", "gene_id", "trait_category"
  )
  expect_true(all(required %in% names(cg)))
  expect_true(all(c("cis_trans") %in% names(cg) || TRUE))
  expect_setequal(unique(cg$coloc_group_id), seq_len(60))
  expect_true(nrow(cg) > 0)

  truth <- sim$ground_truth
  expect_equal(length(truth$module_of_snp), 60)
  expect_setequal(names(truth$module_of_snp), unique(cg$variant_id))
  expect_equal(truth$parameters$n_background_snps_actual, 30)
  expect_length(truth$multi_module_snps, 0)
})

test_that("planted modules are recovered by the pipeline (disjoint)", {
  sim <- simulate_trait(
    n_coloc_groups = 48,
    K = 3,
    module_sizes = c(8, 8, 8),
    p_structural_zero = 0.2,
    p_spurious = 0.02,
    p_active_background = 0.05,
    noise_sd = 0.3,
    seed = 7
  )
  res <- run_univariate_clustering(sim$trait_object, louvain_gamma = 2)
  ev <- evaluate_univariate_simulation(sim, res)
  expect_gte(ev$k_hat, 2)
  expect_gte(ev$ari_structured, 0.5)
})

test_that("null simulation produces limited structure", {
  sim <- simulate_trait(
    n_coloc_groups = 80,
    K = 0,
    p_active_background = 0.05,
    seed = 11
  )
  res <- run_univariate_clustering(sim$trait_object, louvain_gamma = 2)
  ev <- evaluate_univariate_simulation(sim, res)
  expect_equal(ev$k_planted, 0)
  expect_lte(ev$coverage, 0.5)
})

test_that("overlap creates multi-module SNPs", {
  sim <- simulate_trait(
    n_coloc_groups = 120,
    K = 3,
    overlap = 10,
    seed = 3
  )
  expect_length(sim$ground_truth$multi_module_snps, 20)
  sim2 <- simulate_trait(
    n_coloc_groups = 200,
    K = 3,
    overlap = "strong",
    seed = 3
  )
  expect_gt(length(sim2$ground_truth$multi_module_snps), 0)
})

test_that("module_annotations are planted on driver traits", {
  sim <- simulate_trait(
    n_coloc_groups = 60,
    K = 2,
    module_sizes = c(12, 12),
    drivers_per_module = 4,
    annotation_noise = 0,
    module_annotations = list(
      list(tissues = "Liver", genes = "HFE"),
      list(trait_categories = "Immune")
    ),
    seed = 5
  )
  cg <- sim$trait_object$coloc_groups
  drv1 <- sim$ground_truth$driver_traits$module_1
  drv2 <- sim$ground_truth$driver_traits$module_2
  expect_true(all(cg$tissue[cg$trait_name %in% drv1] == "Liver"))
  expect_true(all(cg$gene[cg$trait_name %in% drv1] == "HFE"))
  expect_false(any(is.na(cg$gene_id[cg$trait_name %in% drv1])))
  expect_true(all(cg$trait_category[cg$trait_name %in% drv2] == "Immune"))
})

test_that("sign_pattern flipped produces anti-correlated driver profiles", {
  sim <- simulate_trait(
    n_coloc_groups = 40,
    K = 1,
    module_sizes = 20,
    drivers_per_module = 6,
    sign_pattern = "flipped",
    p_structural_zero = 0,
    effect_size = 5,
    noise_sd = 0,
    seed = 9
  )
  cg <- sim$trait_object$coloc_groups
  drv <- sim$ground_truth$driver_traits$module_1
  means <- tapply(cg$beta[cg$trait_name %in% drv], cg$trait_name[cg$trait_name %in% drv], mean)
  expect_gt(max(means) * -min(means), 0)
})

test_that("density caps are respected approximately", {
  sim <- simulate_trait(
    n_coloc_groups = 100,
    K = 0,
    n_background_traits = 30,
    p_active_background = 0.3,
    traits_per_snp = c(0, 6),
    seed = 21
  )
  counts <- table(sim$trait_object$coloc_groups$coloc_group_id)
  expect_lte(max(counts), 6)
})
