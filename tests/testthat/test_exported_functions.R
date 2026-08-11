library(testthat)

test_that("health() returns expected output", {
  result <- health_api()
  expect_type(result, "list")
  expect_true(result$status %in% c("healthy", "unhealthy"))
})

test_that("search_gpmap() returns expected output", {
  result <- search_gpmap("haemoglobin")
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
  expected_names <- c(
    "call", "name", "type", "type_id", "num_coloc_groups",
    "num_coloc_studies", "num_rare_results", "num_study_extractions"
  )
  expect_true(all(expected_names %in% names(result)))
})

test_that("trait() returns expected output", {
  trait_id <- 5020
  result <- trait(trait_id)
  expect_type(result, "list")
  expect_true(result$trait$id == trait_id)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("trait() returns full_associations when requested", {
  trait_id <- 5020
  result <- trait(trait_id, include_full_associations = TRUE)
  expect_type(result, "list")
  expect_true(is.data.frame(result$full_associations))
  expect_true(nrow(result$full_associations) > 0)
  expected_names <- c("variant_id", "study_id", "beta", "se", "p", "eaf", "imputed")
  expect_true(all(expected_names %in% names(result$full_associations)))
})

test_that("build_ebmf_matrix() builds oriented z-score features with gene collapse", {
  skip_if_not(requireNamespace("flashier", quietly = TRUE))
  trait_id <- 5020
  ebmf_data <- build_ebmf_matrix(
    trait_id = trait_id,
    snp_key = "variant_id",
    min_snp_signals = 5L
  )
  expect_true(is.matrix(ebmf_data$beta_matrix))
  expect_true(is.matrix(ebmf_data$se_matrix))
  expect_identical(dim(ebmf_data$beta_matrix), dim(ebmf_data$se_matrix))
  expect_gt(nrow(ebmf_data$beta_matrix), 1)
  expect_gt(ncol(ebmf_data$beta_matrix), 1)
  expect_true(all(ebmf_data$se_matrix > 0, na.rm = TRUE))
  expect_true(all(c("trait", "gene") %in% ebmf_data$trait_info$feature_type) ||
    "trait" %in% ebmf_data$trait_info$feature_type)
  expect_false("trait_annotations" %in% names(ebmf_data))
  expect_false("snp_annotations" %in% names(ebmf_data))
  expect_true(is.matrix(ebmf_data$verification_trait_matrix))
  expect_equal(
    rownames(ebmf_data$verification_trait_matrix),
    colnames(ebmf_data$beta_matrix)
  )
  expect_true(!is.null(ebmf_data$prep_summary))
  expect_equal(ebmf_data$prep_summary$compress_method, "asinh")
  expect_lt(
    ebmf_data$prep_summary$max_abs_used,
    ebmf_data$prep_summary$max_abs_raw
  )
})

test_that("compress_effect_matrix() soft-compresses extremes without hard clip", {
  x <- matrix(
    c(1, 5, 1000, -1000),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("f1", "f2"))
  )
  y <- compress_effect_matrix(x, method = "asinh")
  expect_equal(dim(y), dim(x))
  expect_equal(sign(y), sign(x))
  expect_lt(max(abs(y)), max(abs(x)))
  expect_gt(abs(y[1, 2]), abs(y[1, 1]))
  expect_identical(compress_effect_matrix(x, method = "none"), x)
})

test_that("variant() returns expected output", {
  variant_id <- 5553693
  result <- variant(variant_id)
  expect_type(result, "list")
  expect_true(result$variant$id == variant_id)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("gene() returns expected output", {
  gene_id <- "WNT7B"
  result <- gene(gene_id)
  expect_type(result, "list")
  expect_true(result$gene$gene == gene_id)
  expect_true(nrow(result$coloc_groups) > 0)
})

test_that("associations() returns expected output", {
  variant_ids <- c(80750)
  study_ids <- c(5020, 4870)
  result <- associations(variant_ids, study_ids)
  expect_type(result, "list")
  expect_true(all(result$variant_id %in% variant_ids))
  expect_true(all(result$study_id %in% study_ids))
})

test_that("get_all_gene_pleiotropies() returns expected output", {
  result <- get_all_gene_pleiotropies()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("get_all_variant_pleiotropies() returns expected output", {
  result <- get_all_variant_pleiotropies()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("all_genes() returns expected output", {
  result <- all_genes()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("all_traits() returns expected output", {
  result <- all_traits()
  expect_type(result, "list")
  expect_true(nrow(result) > 0)
})

test_that("traits(trait_ids) returns expected output", {
  trait_ids <- c(4405, 4872)
  result <- traits(trait_ids = trait_ids, include_associations = TRUE)
  expect_type(result, "list")
  expect_true(all(result$traits$id %in% trait_ids))
  expect_true(nrow(result$coloc_groups) > 0)
  expect_true(!is.null(result$study_extractions))
  expect_true(nrow(result$study_extractions) > 0)
})

test_that("genes(gene_ids) returns expected output", {
  gene_ids <- c("WNT7B", "WNT7A")
  result <- genes(gene_ids = gene_ids, include_associations = TRUE)
  expect_type(result, "list")
  expect_true(all(result$genes$gene %in% gene_ids))
  expect_true(nrow(result$coloc_groups) > 0)
  expect_true(!is.null(result$study_extractions))
  expect_true(nrow(result$study_extractions) > 0)
})

test_that("variants() returns expected output", {
  variant_ids <- c(5553693, 5553694)
  result <- variants(variants = variant_ids, include_associations = TRUE, expand = TRUE)
  expect_type(result, "list")
  expect_true(all(result$variants$variant_id %in% variant_ids))
  expect_true(nrow(result$coloc_groups) > 0)
  expect_true(!is.null(result$study_extractions))
  expect_true(nrow(result$study_extractions) > 0)
})

test_that("pathway_enrichment() returns expected output with gene IDs", {
  genes <- all_genes()
  gene_ids <- genes[genes$gene %in% c("TREM2", "APOE"), "id"]

  result <- pathway_enrichment(gene_ids, minimum_count_in_network = 1)
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_true("input_gene_count" %in% names(result))
  expect_true("matched_gene_count" %in% names(result))
  expect_true("p_value_threshold" %in% names(result))
  expect_true("total_terms_tested" %in% names(result))
  expect_equal(result$input_gene_count, length(gene_ids))
  expect_true(is.data.frame(result$results) || is.null(result$results))
  if (is.data.frame(result$results) && nrow(result$results) > 0) {
    expected_cols <- c(
      "term_id", "source", "description", "pathway_size",
      "background_size", "overlap", "p_value", "fdr", "gene_ids",
      "pathway_gene_ids"
    )
    expect_true(all(expected_cols %in% names(result$results)))
    expect_true(all(result$results$fdr <= 0.05))
    expect_true(is.list(result$results$gene_ids))
    expect_true(is.list(result$results$pathway_gene_ids))
    expect_true(length(result$results$pathway_gene_ids[[1]]) >=
                  length(result$results$gene_ids[[1]]))
  }
})

test_that("pathway_enrichment() accepts gene names", {
  result <- pathway_enrichment(
    c("TREM2", "APOE"),
    minimum_count_in_network = 1
  )
  expect_type(result, "list")
  expect_equal(result$input_gene_count, 2L)
})

test_that("pathway_enrichment() validates inputs", {
  expect_error(pathway_enrichment(NULL), "genes is required")
  expect_error(pathway_enrichment(c(1, NA)), "genes must not contain NA")
  expect_error(pathway_enrichment(c(TRUE)), "genes must be numeric gene IDs or character gene names")
  expect_error(pathway_enrichment(1, source = "Invalid"), "source must be one of")
  expect_error(pathway_enrichment(1, p_value_threshold = 2), "p_value_threshold must be")
  expect_error(pathway_enrichment(1, minimum_count_in_network = 0),
               "minimum_count_in_network must be a positive integer")
})

test_that("genes_at_snps() maps molecular QTL genes to SNPs", {
  coloc_groups <- data.frame(
    variant_id = c("v1", "v1", "v2", "v3"),
    coloc_group_id = c(10L, 10L, 11L, 12L),
    gene_id = c(1L, NA, 2L, 3L),
    gene = c("GENEA", NA, "GENEB", "GENEC"),
    stringsAsFactors = FALSE
  )
  out <- genes_at_snps(c("v1", "v2"), coloc_groups)
  expect_equal(nrow(out), 2L)
  expect_equal(sort(out$gene), c("GENEA", "GENEB"))
  expect_true(all(c("snp_id", "coloc_group_id", "gene_id", "gene") %in% names(out)))
})

test_that("compare_group_pathways() flags recovered and split pathways", {
  trait_enrichment <- list(
    pathways = data.frame(
      source = c("KEGG", "Reactome"),
      term_id = c("k1", "r1"),
      description = c("Path A", "Path B"),
      overlap = c(3L, 4L),
      p_value = c(0.01, 0.02),
      fdr = c(0.03, 0.04),
      input_genes = c("1, 2", "3, 4"),
      stringsAsFactors = FALSE
    )
  )
  group_enrichment <- list(
    by_group = list(
      list(
        group = "1",
        pathways = data.frame(
          source = c("KEGG", "KEGG"),
          term_id = c("k1", "k2"),
          description = c("Path A", "Path C"),
          overlap = c(2L, 2L),
          p_value = c(0.01, 0.02),
          fdr = c(0.03, 0.04),
          input_genes = c("1", "5"),
          stringsAsFactors = FALSE
        )
      ),
      list(
        group = "2",
        pathways = data.frame(
          source = "KEGG",
          term_id = "k1",
          description = "Path A",
          overlap = 2L,
          p_value = 0.01,
          fdr = 0.03,
          input_genes = "2",
          stringsAsFactors = FALSE
        )
      )
    )
  )
  cmp <- compare_group_pathways(trait_enrichment, group_enrichment)
  expect_true("split_across_groups" %in% cmp$pathway_status$status)
  expect_true("baseline_only" %in% cmp$pathway_status$status)
  expect_true("group_specific" %in% cmp$pathway_status$status)
  expect_equal(cmp$group_overlap$n_shared_with_baseline[1], 1L)
  expect_equal(cmp$group_overlap$n_novel[1], 1L)
})

test_that("extract_ebmf_clusters defaults to stricter magnitude_threshold", {
  expect_equal(formals(extract_ebmf_clusters)$magnitude_threshold, 0.25)
  expect_equal(formals(run_ebmf_comparison)$magnitude_threshold, 0.25)
  expect_null(formals(run_ebmf_comparison)$label_schemes)
  expect_equal(formals(build_ebmf_matrix)$min_snp_signals, 5L)
  expect_equal(eval(formals(build_ebmf_matrix)$compress_method), c("asinh", "none"))
  expect_equal(formals(select_ebmf_comparison_run)$min_total_pve, 0.01)
  expect_equal(eval(formals(run_ebmf)$se_mode), c("unit", "matrix"))
})

test_that("zscore_feature_columns() standardises columns", {
  x <- matrix(
    c(1, 2, 3, 10, 20, 30),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), c("f1", "f2"))
  )
  z <- zscore_feature_columns(x)
  expect_equal(unname(round(colMeans(z), 10)), c(0, 0))
  expect_equal(unname(round(apply(z, 2, sd), 10)), c(1, 1))
})

test_that("scale_feature_columns() can skip centering", {
  x <- matrix(
    c(1, 2, 3, 10, 20, 30),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), c("f1", "f2"))
  )
  scaled <- scale_feature_columns(x, center = FALSE)
  expect_false(isTRUE(all.equal(unname(colMeans(scaled)), c(0, 0))))
  expect_equal(unname(round(apply(scaled, 2, sd), 10)), c(1, 1))
})

test_that("filter_perturbation_features() drops sparse columns", {
  trait_matrix <- matrix(
    c(1, 2, 3, 3, NA, 4),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), c("t1", "t2"))
  )
  gene_matrix <- matrix(
    c(1, 2, 3, NA, NA, 1),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), c("g1", "g2"))
  )
  out <- filter_perturbation_features(trait_matrix, gene_matrix, min_snps = 3)
  expect_equal(out$n_traits_kept, 1L)
  expect_equal(out$n_genes_kept, 1L)
  expect_equal(colnames(out$trait_matrix), "t1")
  expect_equal(colnames(out$gene_matrix), "g1")
})

test_that("orient_perturbation_matrices() flips SNP rows by target sign", {
  trait_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("t1", "t2"))
  )
  gene_matrix <- matrix(
    c(5, 6, 7, 8),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("g1", "g2"))
  )
  z_target <- c(s1 = 1.5, s2 = -2)
  out <- orient_perturbation_matrices(trait_matrix, gene_matrix, z_target)
  expect_equal(out$trait_matrix["s2", ], -trait_matrix["s2", ])
  expect_equal(out$gene_matrix["s1", ], gene_matrix["s1", ])
})

test_that("weighted_snp_similarity() combines trait and gene cosines", {
  trait_matrix <- matrix(
    c(1, 0, 0, 1, 1, 0),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), c("t1", "t2"))
  )
  gene_matrix <- matrix(
    c(1, 1, 0, 0, 1, 1),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), c("g1", "g2"))
  )
  sim <- weighted_snp_similarity(trait_matrix, gene_matrix, w_gene = 0.5, w_trait = 0.5)
  expect_equal(dim(sim$s_matrix), c(3L, 3L))
  expect_equal(unname(diag(sim$s_matrix)), rep(1, 3), tolerance = 1e-8)
  expect_equal(sim$s_matrix, 0.5 * sim$s_gene + 0.5 * sim$s_trait)
})

test_that("summarise_baseline_traits() ranks by abs mean z", {
  trait_matrix <- matrix(
    c(1, 2, NA, 0, 3, 4),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), c("t1", "t2"))
  )
  out <- summarise_baseline_traits(
    trait_matrix,
    top_n = 2,
    min_snps_with_signal = 1L
  )
  expect_true(is.list(out))
  expect_equal(out$traits$trait_id[1], "t2")
  expect_true(out$traits$abs_mean_z[1] >= out$traits$abs_mean_z[2])
})

test_that("summarise_snp_group_traits() reports top trait categories", {
  trait_matrix <- matrix(
    c(2, 1, 0, 3, 1, 0, 0, 0, 4),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), c("10", "20", "30"))
  )
  coloc_groups <- data.frame(
    trait_id = c(10L, 20L, 30L),
    trait_category = c("Anthropometric Measures", "Anthropometric Measures", "Metabolic Disease"),
    stringsAsFactors = FALSE
  )
  groups <- c(s1 = 1, s2 = 1, s3 = 1)
  out <- summarise_snp_group_traits(
    trait_matrix = trait_matrix,
    groups = groups,
    coloc_groups = coloc_groups,
    min_group_size = 1,
    top_n = 3,
    min_snps_with_signal = 1L,
    n_categories = 2L
  )
  expect_true("top_trait_categories" %in% names(out$summary))
  expect_true(grepl("Anthropometric", out$summary$top_trait_categories[1]))
  expect_equal(nrow(out$by_group[[1]]$categories) <= 2, TRUE)
})

test_that("summarise_module_specific_traits() ranks by module-vs-rest", {
  trait_matrix <- matrix(
    c(
      5, 0.2, 0.1,
      5, 0.2, 0.1,
      0.2, 4, 0.1,
      0.2, 4, 0.1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("s1", "s2", "s3", "s4"), c("10", "20", "30"))
  )
  coloc_groups <- data.frame(
    trait_id = c(10L, 20L, 30L),
    trait_category = c("Lipids", "Glycemic", "Physiological Measures"),
    stringsAsFactors = FALSE
  )
  groups <- c(s1 = 1, s2 = 1, s3 = 2, s4 = 2)
  out <- summarise_module_specific_traits(
    trait_matrix = trait_matrix,
    groups = groups,
    coloc_groups = coloc_groups,
    min_group_size = 1,
    top_n = 3,
    min_snps_with_signal = 1L,
    n_categories = 1L,
    min_specificity = 1.1
  )
  expect_equal(length(out$by_group), 2L)
  expect_true("specificity" %in% names(out$by_group[[1]]$traits))
  expect_true("passes_min_specificity" %in% names(out$by_group[[1]]$traits))
  expect_equal(out$by_group[[1]]$traits$trait_id[1], "10")
  expect_equal(out$by_group[[2]]$traits$trait_id[1], "20")
  expect_true(grepl("Lipids", out$summary$top_trait_categories[1]))
  expect_true(all(out$summary$n_traits_ranked > 0))
})

test_that("enrich_snp_group_tissues() flags enriched tissues", {
  coloc_groups <- data.frame(
    variant_id = c("s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"),
    gene_id = 1:8,
    tissue = c(
      "Adipose", "Adipose", "Adipose", "Blood",
      "Blood", "Blood", "Blood", "Brain"
    ),
    stringsAsFactors = FALSE
  )
  groups <- c(s1 = 1, s2 = 1, s3 = 2, s4 = 2)
  out <- enrich_snp_group_tissues(
    groups = groups,
    coloc_groups = coloc_groups,
    min_group_size = 1,
    background = "rest"
  )
  expect_true(nrow(out$summary) >= 1)
  expect_true("by_group" %in% names(out))
  expect_true("background" %in% names(out))
  expect_equal(out$background_mode, "rest")
  expect_true("fold_enrichment" %in% names(out$by_group[[1]]$enriched) ||
    nrow(out$by_group[[1]]$enriched) == 0)
})


test_that("build_pleiotropy_matrix filters to target trait SNPs from shared coloc_groups", {
  coloc_groups <- data.frame(
    trait_id = c(1L, 1L, 1L, 2L, 2L, 99L, 99L),
    trait_name = c("A", "A", "A", "B", "B", "bg", "bg"),
    coloc_group_id = c(10L, 11L, 12L, 20L, 21L, 10L, 20L),
    variant_id = c("v1", "v2", "v3", "v4", "v5", "v1", "v4"),
    display_snp = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs1", "rs4"),
    chr = 1L,
    bp = seq_len(7),
    beta = c(1, 1, 1, 1, 1, 0.5, 0.5),
    se = rep(0.1, 7),
    min_p = 1e-10,
    stringsAsFactors = FALSE
  )

  p1 <- build_pleiotropy_matrix(1L, coloc_groups = coloc_groups)
  p2 <- build_pleiotropy_matrix(2L, coloc_groups = coloc_groups)
  expect_equal(ncol(p1$x_matrix), 3L)
  expect_equal(ncol(p2$x_matrix), 2L)
  expect_true(setequal(colnames(p1$x_matrix), c("v1", "v2", "v3")))
  expect_true(setequal(colnames(p2$x_matrix), c("v4", "v5")))

  bivariate <- build_bivariate_pleiotropy_matrices(
    trait_id_1 = 1L,
    trait_id_2 = 2L,
    coloc_groups = coloc_groups
  )
  expect_equal(ncol(bivariate$x1_matrix), 3L)
  expect_equal(ncol(bivariate$x2_matrix), 2L)
  expect_true(setequal(colnames(bivariate$x1_matrix), c("v1", "v2", "v3")))
  expect_true(setequal(colnames(bivariate$x2_matrix), c("v4", "v5")))
})

test_that("cluster_cross_trait_snps() hierarchical method is backward compatible", {
  set.seed(1)
  s_star <- matrix(stats::runif(30, min = -1, max = 1), nrow = 5, ncol = 6)
  rownames(s_star) <- paste0("r", seq_len(5))
  colnames(s_star) <- paste0("c", seq_len(6))

  result <- cluster_cross_trait_snps(
    s_star,
    method = "hierarchical",
    linkage = "average",
    k1 = 2,
    k2 = 2
  )

  expect_equal(result$method, "hierarchical")
  expect_length(result$clusters_rows, 5)
  expect_length(result$clusters_cols, 6)
  expect_null(result$profile_similarity_rows)
})

test_that("cluster_cross_trait_snps() louvain clusters cross-trait profiles", {
  set.seed(1)
  s_star <- matrix(stats::runif(30, min = -1, max = 1), nrow = 5, ncol = 6)
  rownames(s_star) <- paste0("r", seq_len(5))
  colnames(s_star) <- paste0("c", seq_len(6))

  result <- cluster_cross_trait_snps(
    s_star,
    method = "louvain",
    similarity_threshold = 0.2,
    gamma = 1
  )

  expect_equal(result$method, "louvain")
  expect_length(result$clusters_rows, 5)
  expect_length(result$clusters_cols, 6)
  expect_equal(nrow(result$profile_similarity_rows), 5)
  expect_equal(nrow(result$profile_similarity_cols), 6)
  expect_true(all(abs(diag(result$profile_similarity_rows) - 1) < 1e-10))
})