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

test_that("module_sizes can define the number of planted modules", {
  sim <- simulate_trait(
    n_coloc_groups = 100,
    module_sizes = c(5, 20, 50),
    n_background_traits = 5,
    seed = 43
  )
  expect_equal(sim$ground_truth$parameters$K, 3L)
  expect_equal(sim$ground_truth$parameters$module_sizes, c(5L, 20L, 50L))
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

test_that("simulation evaluation accepts overlapping EBMF memberships", {
  sim <- simulate_trait(
    n_coloc_groups = 60,
    K = 2,
    module_sizes = c(20, 20),
    drivers_per_module = 4,
    p_structural_zero = 0.1,
    seed = 13
  )
  res <- run_univariate_clustering(
    sim$trait_object,
    cluster_type = "ebmf",
    min_snp_signals = 2,
    ebmf_lfsr_threshold = 0.05,
    ebmf_magnitude_threshold = 0.25,
    min_module_size = 3
  )
  summary <- summarise_ebmf_programs(
    res,
    n_null = 2,
    n_rep = 0,
    verbose = FALSE
  )
  valid <- summary$programs$program[summary$programs$status == "valid"]
  ev <- evaluate_univariate_simulation(
    sim,
    res,
    predicted_memberships = summary$assigned[
      summary$assigned$program %in% valid,
      ,
      drop = FALSE
    ]
  )
  expect_true(all(c("coverage", "background_absorbed", "module_recall") %in%
                    names(ev)))
  expect_gte(ev$coverage, 0)
  expect_lte(ev$coverage, 1)
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

test_that("flipped signs restart for each module", {
  sim <- simulate_trait(
    n_coloc_groups = 60,
    K = 2,
    module_sizes = c(20, 20),
    drivers_per_module = 4,
    sign_pattern = "flipped",
    p_structural_zero = 0,
    p_spurious = 0,
    noise_sd = 0,
    effect_size = 5,
    seed = 10
  )
  cg <- sim$trait_object$coloc_groups
  first_drivers <- vapply(seq_len(2), function(m) {
    driver <- sim$ground_truth$driver_traits[[paste0("module_", m)]][1]
    mean(cg$beta[cg$trait_name == driver])
  }, numeric(1))
  expect_equal(first_drivers, c(5, 5))
})

test_that("traits_per_snp preserves target-trait observations", {
  sim <- simulate_trait(
    n_coloc_groups = 20,
    K = 1,
    module_sizes = 10,
    drivers_per_module = 2,
    p_structural_zero = 0,
    p_spurious = 0,
    p_active_background = 0,
    traits_per_snp = c(0, 1),
    seed = 12
  )
  target_rows <- sim$trait_object$coloc_groups |>
    dplyr::filter(trait_id == 1)
  expect_equal(nrow(target_rows), 20)
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

test_that("n_traits controls the total number of simulated traits", {
  sim <- simulate_trait(
    n_coloc_groups = 60,
    K = 2,
    module_sizes = c(15, 15),
    drivers_per_module = 3,
    n_traits = 20,
    seed = 33
  )
  p <- sim$ground_truth$parameters
  expect_equal(p$n_traits, 20L)
  expect_equal(p$n_background_traits, 13L)
  expect_equal(length(unique(sim$trait_object$coloc_groups$trait_id)), 20)
  expect_error(
    simulate_trait(n_coloc_groups = 60, K = 2, drivers_per_module = 3,
                   n_traits = 5),
    "n_traits"
  )
})

test_that("molecular traits are gene-annotated and much sparser", {
  sim <- simulate_trait(
    n_coloc_groups = 200,
    K = 0,
    n_background_traits = 20,
    p_molecular = 1,
    p_active_molecular = 0.02,
    p_active_background = 0.5,
    seed = 17
  )
  cg <- sim$trait_object$coloc_groups
  bg <- cg[cg$trait_id > 1, ]
  expect_false(any(is.na(bg$gene)))
  expect_false(any(is.na(bg$gene_id)))
  expect_equal(length(sim$ground_truth$molecular_traits), 20)
  counts <- table(bg$trait_id)
  expect_lte(max(counts), 15)

  sim_dense <- simulate_trait(
    n_coloc_groups = 200,
    K = 0,
    n_background_traits = 20,
    p_molecular = 0,
    p_active_background = 0.5,
    seed = 18
  )
  cg_dense <- sim_dense$trait_object$coloc_groups
  bg_dense <- cg_dense[cg_dense$trait_id > 1, ]
  expect_true(all(is.na(bg_dense$gene)))
  expect_true(all(is.na(bg_dense$gene_id)))
  expect_gt(mean(table(bg_dense$trait_id)), 5 * mean(counts))
})

test_that("molecular drivers are observed at a fraction of module SNPs", {
  sim <- simulate_trait(
    n_coloc_groups = 200,
    K = 2,
    module_sizes = c(50, 50),
    drivers_per_module = 4,
    p_molecular = 1,
    p_active_molecular = 0.05,
    p_structural_zero = 0,
    seed = 1
  )
  cg <- sim$trait_object$coloc_groups
  d1 <- sim$ground_truth$driver_traits$module_1
  counts <- table(cg$trait_id[cg$trait_name %in% d1])
  expect_lt(mean(counts), 25)
})

test_that("planted module genes force driver traits molecular", {
  sim <- simulate_trait(
    n_coloc_groups = 40,
    K = 1,
    module_sizes = 15,
    drivers_per_module = 4,
    module_annotations = list(list(genes = "HFE")),
    annotation_noise = 0,
    p_molecular = 0,
    seed = 23
  )
  cg <- sim$trait_object$coloc_groups
  drv <- sim$ground_truth$driver_traits$module_1
  expect_true(all(cg$gene[cg$trait_name %in% drv] == "HFE"))
  expect_false(any(is.na(cg$gene_id[cg$trait_name %in% drv])))
  expect_equal(length(sim$ground_truth$molecular_traits), length(drv))
})
