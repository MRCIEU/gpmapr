library(testthat)

make_fake_ebmf_result <- function(F_pm, lfsr = 0.001) {
  F_lfsr <- matrix(
    lfsr, nrow = nrow(F_pm), ncol = ncol(F_pm), dimnames = dimnames(F_pm)
  )
  list(
    parameters = list(method = "ebmf"),
    cluster_details = list(flash_fit = list(
      F_pm = F_pm, F_lfsr = F_lfsr, n_factors = ncol(F_pm)
    ))
  )
}

make_planted_fixture <- function(seed = 42) {
  set.seed(seed)
  snp_ids <- paste0("snp", 1:30)
  loading_p1 <- c(
    stats::rnorm(10, 0.8, 0.05), stats::rnorm(10, 0.02, 0.05),
    stats::rnorm(10, 0.02, 0.05)
  )
  loading_p2 <- c(
    stats::rnorm(10, 0.02, 0.05), stats::rnorm(10, 0.8, 0.05),
    stats::rnorm(10, 0.02, 0.05)
  )
  F_pm <- cbind(loading_p1, loading_p2)
  rownames(F_pm) <- snp_ids

  cg <- data.frame(
    variant_id = snp_ids,
    coloc_group_id = 1:30,
    gene_id = c(rep(101L, 10), rep(102L, 10), rep(103L, 10)),
    gene = c(rep("GENE1", 10), rep("GENE2", 10), rep("GENE3", 10)),
    tissue = c(rep("TissueA", 10), rep("TissueB", 10), rep("TissueC", 10)),
    trait_category = c(rep("A", 10), rep("B", 10), rep("C", 10)),
    stringsAsFactors = FALSE
  )

  mappings <- list(mappings = data.frame(
    gene_id = c(101L, 102L, 103L),
    term_id = c("path1", "path2", "path3"),
    source = c("KEGG", "KEGG", "KEGG"),
    description = c("Pathway One", "Pathway Two", "Pathway Three"),
    stringsAsFactors = FALSE
  ))

  list(
    clustering_result = make_fake_ebmf_result(F_pm),
    coloc_groups = cg,
    mappings = mappings
  )
}


test_that("enrich_program_loadings_trait_categories recovers a planted signal", {
  fx <- make_planted_fixture()
  res <- enrich_program_loadings_trait_categories(fx$clustering_result, fx$coloc_groups)

  expect_setequal(names(res), c("by_program", "summary"))
  expect_equal(nrow(res$summary), 2L)
  expect_equal(res$summary$top_category, c("A", "B"))

  comp1 <- res$by_program[[1]]$comparison
  expect_equal(comp1$trait_category[1], "A")
  expect_true(comp1$enrichment[1] > 0)
  expect_true(comp1$p[1] < 1e-10)
  expect_true(all(diff(comp1$fdr) >= -1e-12))
})


test_that("enrich_program_loadings_tissues recovers a planted signal", {
  fx <- make_planted_fixture()
  res <- enrich_program_loadings_tissues(fx$clustering_result, fx$coloc_groups)

  expect_equal(res$summary$top_category, c("TissueA", "TissueB"))
  expect_true(res$by_program[[2]]$comparison$enrichment[1] > 0)
})


test_that("enrich_program_loadings_pathways recovers a planted signal via genes", {
  fx <- make_planted_fixture()
  res <- enrich_program_loadings_pathways(
    fx$clustering_result, fx$coloc_groups,
    mappings = fx$mappings
  )

  expect_setequal(names(res), c("by_program", "summary", "mappings"))
  expect_equal(res$summary$top_pathway, c(
    "KEGG: Pathway One", "KEGG: Pathway Two"
  ))
  comp1 <- res$by_program[[1]]$comparison
  expect_equal(comp1$term_id[1], "path1")
  expect_setequal(
    names(comp1),
    c("term_id", "source", "description", "enrichment", "se", "z", "p", "fdr",
      "n_snps", "n_category_snps")
  )
})


test_that("min_category_size filters out small categories", {
  fx <- make_planted_fixture()
  res <- enrich_program_loadings_trait_categories(
    fx$clustering_result, fx$coloc_groups,
    min_category_size = 15L
  )
  expect_true(all(vapply(res$by_program, function(x) nrow(x$comparison) == 0L, logical(1))))
})


test_that("continuous trait-category enrichment errors without a trait_category column", {
  fx <- make_planted_fixture()
  cg <- fx$coloc_groups
  cg$trait_category <- NULL
  expect_error(
    enrich_program_loadings_trait_categories(fx$clustering_result, cg),
    "trait_category"
  )
})


test_that("min_loading_magnitude gates diluted members and strengthens the contrast", {
  set.seed(1)
  n_strong <- 5L
  n_weak <- 50L
  n_bg <- 20L
  snp_ids <- paste0("snp", seq_len(n_strong + n_weak + n_bg))
  loading <- c(
    stats::rnorm(n_strong, 0.8, 0.02),
    stats::rnorm(n_weak, 0, 0.005),
    stats::rnorm(n_bg, 0, 0.005)
  )
  F_pm <- matrix(loading, ncol = 1, dimnames = list(snp_ids, NULL))
  clustering_result <- make_fake_ebmf_result(F_pm)
  cg <- data.frame(
    variant_id = snp_ids,
    coloc_group_id = seq_along(snp_ids),
    trait_category = c(rep("A", n_strong + n_weak), rep("B", n_bg)),
    stringsAsFactors = FALSE
  )

  gated <- enrich_program_loadings_trait_categories(
    clustering_result, cg,
    min_loading_magnitude = 0.02
  )$by_program[[1]]$comparison
  ungated <- enrich_program_loadings_trait_categories(
    clustering_result, cg,
    min_loading_magnitude = 0
  )$by_program[[1]]$comparison

  gated_a <- gated[gated$trait_category == "A", ]
  ungated_a <- ungated[ungated$trait_category == "A", ]

  # Ungated: the 50 no-signal "A"-linked SNPs dilute the group mean enough
  # that the true 5-SNP signal isn't detected.
  expect_equal(ungated_a$n_category_snps, n_strong + n_weak)
  expect_true(ungated_a$p > 0.05)

  # Gated: those 50 SNPs fall below the magnitude threshold and are
  # relabelled x = 0, leaving only the 5 true signal SNPs as members.
  expect_equal(gated_a$n_category_snps, n_strong)
  expect_true(gated_a$p < 1e-10)
  expect_true(gated_a$enrichment > ungated_a$enrichment)

  # The background pool (n_snps) is unaffected by gating -- SNPs are
  # relabelled, not dropped.
  expect_equal(gated_a$n_snps, ungated_a$n_snps)
})


test_that("enrich_program_loadings_pathways pools KEGG and Reactome into one FDR family", {
  set.seed(7)
  snp_ids <- paste0("snp", 1:20)
  loading <- c(stats::rnorm(5, 0.8, 0.02), stats::rnorm(15, 0, 0.01))
  F_pm <- matrix(loading, ncol = 1, dimnames = list(snp_ids, NULL))
  clustering_result <- make_fake_ebmf_result(F_pm)
  cg <- data.frame(
    variant_id = snp_ids,
    coloc_group_id = 1:20,
    gene_id = c(rep(1L, 5), 2:16),
    gene = c(rep("G1", 5), paste0("G", 2:16)),
    stringsAsFactors = FALSE
  )
  kegg_map <- data.frame(
    gene_id = 1L, term_id = "kegg_true", source = "KEGG",
    description = "KEGG true", stringsAsFactors = FALSE
  )
  set.seed(11)
  reactome_map <- do.call(rbind, lapply(1:20, function(i) {
    data.frame(
      gene_id = sample(2:16, 5), term_id = paste0("react", i), source = "Reactome",
      description = paste0("React ", i), stringsAsFactors = FALSE
    )
  }))
  mappings <- list(mappings = rbind(kegg_map, reactome_map))

  res <- enrich_program_loadings_pathways(
    clustering_result, cg,
    mappings = mappings, sources = c("KEGG", "Reactome"),
    min_category_size = 5L, min_loading_magnitude = 0
  )
  comp <- res$by_program[[1]]$comparison

  kegg_rows <- comp[comp$source == "KEGG", ]
  reactome_rows <- comp[comp$source == "Reactome", ]
  expect_equal(nrow(kegg_rows), 1L)
  expect_true(nrow(reactome_rows) > 1L)

  # KEGG and Reactome share one BH family: reconstructing fdr from the pooled
  # p-values across both sources reproduces every reported fdr exactly. (The
  # planted signal is strong enough that kegg_rows$p underflows to exactly 0,
  # so fdr == p there is expected and not itself evidence of per-source
  # correction -- the pooled reconstruction below is the real check.)
  pooled_fdr <- stats::p.adjust(comp$p, method = "BH")
  expect_equal(sort(comp$fdr), sort(pooled_fdr))
})


test_that("continuous enrichment pools FDR across programs, not within one program", {
  # Two programs, each with its own distinct planted category -- program 1's
  # SNPs only carry signal for category A, program 2's only for category B.
  set.seed(21)
  snp_ids <- paste0("snp", 1:60)
  loading_p1 <- c(stats::rnorm(10, 0.8, 0.02), stats::rnorm(50, 0, 0.01))
  loading_p2 <- c(stats::rnorm(10, 0, 0.01), stats::rnorm(10, 0.8, 0.02), stats::rnorm(40, 0, 0.01))
  F_pm <- cbind(loading_p1, loading_p2)
  rownames(F_pm) <- snp_ids
  clustering_result <- make_fake_ebmf_result(F_pm)
  cg <- data.frame(
    variant_id = snp_ids,
    coloc_group_id = 1:60,
    trait_category = c(rep("A", 10), rep("B", 10), rep("C", 40)),
    stringsAsFactors = FALSE
  )

  res <- enrich_program_loadings_trait_categories(
    clustering_result, cg,
    min_category_size = 5L, min_loading_magnitude = 0
  )
  pooled <- dplyr::bind_rows(lapply(res$by_program, function(x) x$comparison))
  expect_equal(sort(pooled$fdr), sort(stats::p.adjust(pooled$p, method = "BH")))

  # The pooled family spans both programs' tests, so it is strictly larger
  # than either program's own comparison table.
  expect_true(nrow(pooled) > nrow(res$by_program[[1]]$comparison))
})


test_that("continuous enrichment functions return an empty result for an empty posterior", {
  empty_fit <- make_fake_ebmf_result(matrix(numeric(0), nrow = 0, ncol = 0))
  fx <- make_planted_fixture()

  res_cat <- enrich_program_loadings_trait_categories(empty_fit, fx$coloc_groups)
  expect_equal(res_cat$by_program, list())
  expect_equal(nrow(res_cat$summary), 0L)

  res_path <- enrich_program_loadings_pathways(
    empty_fit, fx$coloc_groups,
    mappings = fx$mappings
  )
  expect_equal(res_path$by_program, list())
  expect_equal(nrow(res_path$summary), 0L)
})
