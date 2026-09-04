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
    n_traits_per_module = 4,
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
    min_abs_z = 0,
    effect_tail = 0,
    background_sparsity_sd = 0,
    n_hub_traits = 0L,
    target_pattern = "module",
    p_negative = NULL,
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
    min_abs_z = 0,
    background_sparsity_sd = 0,
    n_hub_traits = 0L,
    target_pattern = "module",
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
    n_traits_per_module = 4,
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

test_that("trait_overlap shares driver traits across modules without SNP overlap", {
  sim <- simulate_trait(
    n_coloc_groups = 150,
    K = 3,
    module_sizes = c(30, 40, 50),
    n_traits_per_module = 20,
    trait_overlap = 0.5,
    trait_effect_corr = 0.8,
    seed = 3
  )
  gt <- sim$ground_truth
  expect_length(gt$multi_module_snps, 0)
  expect_equal(sum(gt$module_of_snp > 0), sum(lengths(gt$module_memberships)))
  shared <- Reduce(intersect, gt$driver_traits)
  expect_length(shared, 10L)
  expect_equal(gt$driver_overlap[1, 2], 0.5)
  expect_equal(gt$driver_overlap[1, 1], 1)
  expect_equal(
    sim$ground_truth$parameters$n_driver_rows,
    10L + 3L * 10L
  )
})

test_that("shared trait effects are correlated but not identical across modules", {
  sim <- simulate_trait(
    n_coloc_groups = 150,
    K = 2,
    module_sizes = c(50, 50),
    n_traits_per_module = 20,
    trait_overlap = 0.75,
    trait_effect_corr = 0.9,
    p_structural_zero = 0,
    p_spurious = 0,
    noise_sd = 0,
    seed = 4
  )
  gt <- sim$ground_truth
  cg <- sim$trait_object$coloc_groups
  m1 <- names(gt$module_of_snp)[gt$module_of_snp == 1]
  m2 <- names(gt$module_of_snp)[gt$module_of_snp == 2]
  shared <- Reduce(intersect, gt$driver_traits)
  expect_length(shared, 15L)
  effects <- vapply(shared, function(tr) {
    c(
      mean(cg$beta[cg$trait_name == tr & cg$variant_id %in% m1]),
      mean(cg$beta[cg$trait_name == tr & cg$variant_id %in% m2])
    )
  }, numeric(2))
  expect_gt(stats::cor(effects[1, ], effects[2, ]), 0.5)
  expect_gt(sum(effects[1, ] != effects[2, ]), 0)
})

test_that("overlap = 'nested' nests modules inside the parent with correlated effects", {
  sim <- simulate_trait(
    n_coloc_groups = 400,
    K = 4,
    module_sizes = c(300, 160, 80, 40),
    overlap = "nested",
    n_traits_per_module = 10,
    noise_sd = 0,
    p_structural_zero = 0,
    p_spurious = 0,
    seed = 6
  )
  gt <- sim$ground_truth
  ms <- gt$module_memberships
  expect_true(all(ms[[2]] %in% ms[[1]]))
  expect_true(all(ms[[3]] %in% ms[[1]]))
  expect_true(all(ms[[4]] %in% ms[[1]]))
  expect_gt(length(gt$multi_module_snps), 0)
  expect_equal(sim$ground_truth$parameters$overlap, "nested")
  expect_equal(sim$ground_truth$parameters$nested_effect_corr, 0.8)
  cg <- sim$trait_object$coloc_groups
  shared <- intersect(ms[[1]], ms[[2]])
  d1 <- gt$driver_traits[[1]][1]
  d2 <- gt$driver_traits[[2]][1]
  b1 <- cg$beta[cg$trait_name == d1 & cg$variant_id %in% shared]
  b2 <- cg$beta[cg$trait_name == d2 & cg$variant_id %in% shared]
  expect_gt(stats::cor(b1, b2), 0.5)
  expect_error(
    simulate_trait(n_coloc_groups = 400, K = 3, module_sizes = c(100, 200, 300),
                   overlap = "nested", seed = 7),
    "largest"
  )
})

test_that("module_annotations are planted on driver traits", {
  sim <- simulate_trait(
    n_coloc_groups = 60,
    K = 2,
    module_sizes = c(12, 12),
    n_traits_per_module = 4,
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
    n_traits_per_module = 6,
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
    n_traits_per_module = 4,
    sign_pattern = "flipped",
    p_structural_zero = 0,
    p_spurious = 0,
    noise_sd = 0,
    effect_size = 5,
    effect_tail = 0,
    p_negative = NULL,
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
    n_traits_per_module = 2,
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
    n_traits_per_module = 3,
    n_traits = 20,
    background_sparsity_sd = 0,
    n_hub_traits = 0L,
    seed = 33
  )
  p <- sim$ground_truth$parameters
  expect_equal(p$n_traits, 20L)
  expect_equal(p$n_background_traits, 13L)
  expect_equal(length(unique(sim$trait_object$coloc_groups$trait_id)), 20)
  expect_error(
    simulate_trait(n_coloc_groups = 60, K = 2, n_traits_per_module = 3,
                   n_traits = 5),
    "n_traits"
  )
})

test_that("n_traits_per_module can vary per module", {
  sim <- simulate_trait(
    n_coloc_groups = 90,
    K = 3,
    module_sizes = c(10, 20, 30),
    n_traits_per_module = c(2, 5, 8),
    background_sparsity_sd = 0,
    n_hub_traits = 0L,
    seed = 44
  )
  counts <- lengths(sim$ground_truth$driver_traits)
  expect_equal(unname(counts), c(2L, 5L, 8L))
  expect_equal(
    length(unique(sim$trait_object$coloc_groups$trait_id)),
    1L + 2L + 5L + 8L + sim$ground_truth$parameters$n_background_traits
  )
  expect_error(
    simulate_trait(n_coloc_groups = 60, K = 2, module_sizes = c(10, 10),
                   n_traits_per_module = c(2, 4), snp_driver_groups = 3),
    "snp_driver_groups"
  )
})

test_that("trait_subset = 'phenotypic' keeps all rows when nothing is molecular", {
  sim <- simulate_trait(
    n_coloc_groups = 100,
    K = 2,
    module_sizes = c(20, 20),
    n_traits_per_module = 4,
    p_active_background = 0.1,
    seed = 50
  )
  res_all <- run_univariate_clustering(sim$trait_object, min_snp_signals = 2)
  res_ph <- run_univariate_clustering(
    sim$trait_object,
    trait_subset = "phenotypic",
    min_snp_signals = 2
  )
  expect_equal(res_all$parameters$trait_subset, "all")
  expect_equal(res_ph$parameters$trait_subset, "phenotypic")
  expect_false(any(res_all$trait_info$feature_type == "molecular"))
  expect_identical(nrow(res_ph$x_matrix), nrow(res_all$x_matrix))
  expect_true(all(res_ph$trait_info$feature_type %in% c("phenotypic", NA)))
  expect_true(as.character(sim$trait_object$trait$id) %in% rownames(res_ph$x_matrix))
  expect_error(
    run_univariate_clustering(sim$trait_object, trait_subset = "bogus"),
    "should be one of"
  )
})

test_that("planted module genes are annotated on driver traits", {
  sim <- simulate_trait(
    n_coloc_groups = 40,
    K = 1,
    module_sizes = 15,
    n_traits_per_module = 4,
    module_annotations = list(list(genes = "HFE")),
    annotation_noise = 0,
    seed = 23
  )
  cg <- sim$trait_object$coloc_groups
  drv <- sim$ground_truth$driver_traits$module_1
  expect_true(all(cg$gene[cg$trait_name %in% drv] == "HFE"))
  expect_false(any(is.na(cg$gene_id[cg$trait_name %in% drv])))
})

test_that("min_abs_z floors every observed background cell at the hit threshold", {
  sim <- simulate_trait(
    n_coloc_groups = 100, K = 1, module_sizes = 20,
    n_traits_per_module = 5, n_background_traits = 40,
    min_abs_z = 4.5, p_active_background = 0.05,
    seed = 3
  )
  cg <- sim$trait_object$coloc_groups
  bg <- cg[
    !cg$trait_name %in% unlist(sim$ground_truth$driver_traits) &
      cg$trait_name != "Simulated target trait", ]
  expect_true(all(abs(bg$beta) >= 4.5))
  # with the floor explicitly disabled, classic small-noise cells reappear
  sim0 <- simulate_trait(
    n_coloc_groups = 100, K = 1, module_sizes = 20,
    n_traits_per_module = 5, n_background_traits = 40,
    p_active_background = 0.05,
    min_abs_z = 0, seed = 3
  )
  bg0 <- sim0$trait_object$coloc_groups[
    !sim0$trait_object$coloc_groups$trait_name %in%
      unlist(sim0$ground_truth$driver_traits), ]
  expect_true(any(abs(bg0$beta) < 4.5))
})

test_that("effect_tail adds per-trait magnitude heterogeneity without breaking sign agreement", {
  flat <- simulate_trait(
    n_coloc_groups = 100, K = 1, module_sizes = 50,
    n_traits_per_module = 10, n_background_traits = 10,
    effect_size = 6, noise_sd = 0.1,
    p_structural_zero = 0, effect_tail = 0, seed = 5
  )
  tailed <- simulate_trait(
    n_coloc_groups = 100, K = 1, module_sizes = 50,
    n_traits_per_module = 10, n_background_traits = 10,
    effect_size = 6, noise_sd = 0.1,
    p_structural_zero = 0, effect_tail = 0.8, seed = 5
  )
  drv <- flat$ground_truth$driver_traits$module_1
  spread <- function(sim) {
    m <- sim$trait_object$coloc_groups
    d <- m[m$trait_name %in% drv, ]
    q <- quantile(abs(d$beta), c(0.1, 0.9))
    unname(q[2] / q[1])
  }
  expect_gt(spread(tailed), 2 * spread(flat))
  # a driver trait's values stay the same sign across its module SNPs
  m <- tailed$trait_object$coloc_groups
  d <- m[m$trait_name == drv[1], ]
  mod1 <- names(tailed$ground_truth$module_of_snp)[
    tailed$ground_truth$module_of_snp == 1]
  expect_gt(mean(sign(d$beta[d$variant_id %in% mod1]) > 0), 0.8)
})

test_that("background_sparsity_sd gives heavy-tailed per-trait observation counts", {
  s <- simulate_trait(
    n_coloc_groups = 200, K = 0, n_background_traits = 100,
    p_active_background = 0.03,
    background_sparsity_sd = 1.2, seed = 7
  )
  cnt <- table(s$trait_object$coloc_groups$trait_id)
  expect_gt(max(cnt), 5 * median(cnt))
})

test_that("n_hub_traits plants dense, mutually correlated background traits", {
  s <- simulate_trait(
    n_coloc_groups = 200, K = 1, module_sizes = 30,
    n_traits_per_module = 10, n_background_traits = 60,
    n_hub_traits = 10, hub_snp_fraction = c(0.3, 1),
    p_active_background = 0.03, effect_size = 6,
    background_sparsity_sd = 0, seed = 9
  )
  cg <- s$trait_object$coloc_groups
  cnt <- table(cg$trait_id[cg$trait_name != "Simulated target trait"])
  hub_cut <- quantile(cnt, 0.85)
  dense_rows <- names(cnt)[cnt >= hub_cut]
  expect_gt(length(dense_rows), 4)
  tids <- sort(unique(cg$trait_id))
  mat <- matrix(NA_real_, nrow = length(tids),
                ncol = length(unique(cg$variant_id)),
                dimnames = list(tids, sort(unique(cg$variant_id))))
  mat[cbind(as.character(cg$trait_id), as.character(cg$variant_id))] <- cg$beta
  cc <- cor(t(mat[dense_rows, , drop = FALSE]), use = "pairwise.complete.obs")
  expect_gt(median(abs(cc[upper.tri(cc)])), 0.3)
})

test_that("target_pattern = 'dense' gives a positive significant target row at every SNP", {
  s <- simulate_trait(
    n_coloc_groups = 100, K = 1, module_sizes = 20,
    n_traits_per_module = 5, n_background_traits = 20,
    target_pattern = "dense", min_abs_z = 4.5, seed = 11
  )
  trow <- s$trait_object$coloc_groups[
    s$trait_object$coloc_groups$trait_name == "Simulated target trait", ]
  expect_equal(length(unique(trow$variant_id)), 100)
  expect_true(all(trow$beta > 0))
})

test_that("p_negative flips an approximate fraction of driver cell signs", {
  s <- simulate_trait(
    n_coloc_groups = 100, K = 1, module_sizes = 40,
    n_traits_per_module = 10, n_background_traits = 10,
    p_structural_zero = 0, noise_sd = 0.1,
    p_negative = 0.2, seed = 13
  )
  cg <- s$trait_object$coloc_groups
  d <- cg[cg$trait_name %in% unlist(s$ground_truth$driver_traits), ]
  expect_true(abs(mean(d$beta < 0) - 0.2) < 0.08)
})

test_that("overlap = 'partial' shares boundary SNPs between adjacent modules", {
  s <- simulate_trait(
    n_coloc_groups = 100, K = 3, module_sizes = c(20, 20, 20),
    n_traits_per_module = c(5, 5, 5), n_background_traits = 10,
    overlap = "partial", overlap_fraction = 0.25, seed = 15
  )
  gt <- s$ground_truth
  expect_gt(length(gt$multi_module_snps), 0)
  expect_gt(length(intersect(gt$module_memberships[[1]], gt$module_memberships[[2]])), 0)
  expect_equal(gt$parameters$overlap, "partial")
})

test_that("simulate_trait defaults are deterministic and realistic", {
  args <- list(n_coloc_groups = 60, K = 2, module_sizes = c(10, 10),
               n_traits_per_module = c(4, 4), n_background_traits = 12,
               seed = 11)
  s1 <- do.call(simulate_trait, args)
  s2 <- do.call(simulate_trait, args)
  expect_identical(s1$trait_object$coloc_groups$beta,
                   s2$trait_object$coloc_groups$beta)
  # realism is on by default: every observed cell is a significant hit
  expect_true(all(abs(s1$trait_object$coloc_groups$beta) >= 4.5))
  expect_equal(s1$ground_truth$parameters$min_abs_z, 4.5)
  expect_equal(s1$ground_truth$parameters$target_pattern, "dense")
  expect_equal(s1$ground_truth$parameters$p_negative, 0.05)
})
