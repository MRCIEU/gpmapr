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

  y_gentle <- compress_effect_matrix(x, method = "asinh", asinh_scale = 2)
  expect_equal(y_gentle, 2 * asinh(x / 2))
  expect_gt(abs(y_gentle[1, 2]), abs(y[1, 2]))
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

test_that("GWAS upload own coloc rows are stamped with the lookup GUID", {
  guid <- "7a289615-c1b4-91f3-3d97-887f60de9155"
  # As returned by the API: gwas_upload_id is set on every row, background
  # studies included. The last row is a background study whose GPMap trait id
  # happens to equal the upload id and must not be stamped.
  cg <- data.frame(
    study_id = c(10, NA, 10, NA),
    existing_study_id = c(NA, 99, NA, 98),
    gwas_upload_id = c(55, 55, 55, 55),
    trait_id = c(NA, 2527, 55, 55),
    trait_name = c(NA, "Allergic rhinitis", NA, "Other trait"),
    variant_id = c(1, 1, 2, 2),
    stringsAsFactors = FALSE
  )
  associations <- data.frame(
    study_id = c(10, 10),
    existing_study_id = c(NA, NA),
    variant_id = c(1, 2),
    beta = c(0.2, -0.1),
    se = c(0.05, 0.04),
    stringsAsFactors = FALSE
  )
  upload <- list(
    trait = list(id = 55, name = "Atopic"),
    coloc_groups = cg,
    associations = associations
  )
  decorated <- .decorate_gwas_upload(upload, guid, include_associations = TRUE)
  expect_true(all(c("beta", "se") %in% names(decorated$coloc_groups)))
  own <- !is.na(decorated$coloc_groups$study_id) &
    decorated$coloc_groups$study_id == 10
  expect_equal(decorated$coloc_groups$trait_id[own], c(guid, guid))
  expect_equal(
    decorated$coloc_groups$trait_id[is.na(decorated$coloc_groups$study_id)],
    c("2527", "55")
  )
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

test_that("extract_ebmf_clusters defaults to stricter magnitude_threshold", {
  expect_equal(formals(extract_ebmf_clusters)$magnitude_threshold, 0.25)
  expect_equal(formals(run_ebmf)$missing_se, 1)
  expect_false(any(c("se_matrix", "se_mode", "var_type") %in% names(formals(run_ebmf))))
})

test_that("EBMF inputs encode unobserved cells as zeros", {
  beta <- matrix(
    c(2, NA, 3, 1, NA, 4),
    nrow = 2,
    dimnames = list(c("t1", "t2"), c("s1", "s2", "s3"))
  )

  hard <- gpmapr:::.prepare_ebmf_flash_inputs(beta)
  expect_equal(hard$data["t2", "s1"], 0)
  expect_false(anyNA(hard$data))
  expect_equal(hard$S, 1)

  soft <- gpmapr:::.prepare_ebmf_flash_inputs(beta, missing_se = 2.5)
  expect_equal(soft$data["t2", "s1"], 0)
  expect_equal(soft$S["t2", "s1"], 2.5)
  expect_equal(soft$S["t1", "s1"], 1)

  # A row or column with no observed non-zero value is removed, with a warning.
  empty_col <- cbind(beta, s4 = c(NA, 0))
  expect_warning(
    trimmed <- gpmapr:::.prepare_ebmf_flash_inputs(empty_col),
    "1 all-zero/all-missing column"
  )
  expect_equal(colnames(trimmed$data), c("s1", "s2", "s3"))

  expect_error(run_ebmf(beta, missing_se = 0), "missing_se")
})

test_that("EBMF refits default to SE 1 for results saved before ebmf_missing_se", {
  expect_equal(gpmapr:::.ebmf_missing_se(list()), 1)
  expect_equal(gpmapr:::.ebmf_missing_se(list(ebmf_missing_se = 2.5)), 2.5)
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
    min_specificity = 1.1
  )
  expect_equal(length(out$by_group), 2L)
  expect_true("specificity" %in% names(out$by_group[[1]]$traits))
  expect_true("passes_min_specificity" %in% names(out$by_group[[1]]$traits))
  expect_equal(out$by_group[[1]]$traits$trait_id[1], "10")
  expect_equal(out$by_group[[2]]$traits$trait_id[1], "20")
  expect_equal(out$by_group[[1]]$traits$trait_category[1], "Lipids")
  expect_true(all(vapply(out$by_group, function(x) nrow(x$traits), integer(1)) > 0))
})

test_that("summarise_module_specific_traits() returns no groups below min_group_size", {
  trait_matrix <- matrix(
    c(1, 2, 3),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), "10")
  )
  groups <- c(s1 = 1, s2 = 2, s3 = 3)
  out <- summarise_module_specific_traits(
    trait_matrix = trait_matrix,
    groups = groups,
    min_group_size = 2L,
    top_n = 3,
    min_snps_with_signal = 1L
  )
  expect_equal(length(out$by_group), 0L)
  expect_false("summary" %in% names(out))

  keep_three <- summarise_module_specific_traits(
    trait_matrix = trait_matrix,
    groups = c(s1 = 1, s2 = 1, s3 = 1),
    min_group_size = 3L,
    top_n = 3,
    min_snps_with_signal = 1L
  )
  expect_equal(length(keep_three$by_group), 1L)
  expect_equal(keep_three$by_group[[1]]$group, "1")
})


test_that("enrich_snp_set() tests tissues against the rest of the trait's SNPs", {
  coloc_groups <- data.frame(
    variant_id = c("s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"),
    gene_id = c(1:7, NA),
    tissue = c(
      "Adipose", "Adipose", "Adipose", "Blood",
      "Blood", "Blood", "Blood", "Brain"
    ),
    stringsAsFactors = FALSE
  )
  out <- enrich_snp_set(
    c("s1", "s2"), c("s1", "s2", "s3", "s4"), coloc_groups,
    annotation = "tissue", min_set_size = 1
  )
  expect_true(all(c(
    "tissue", "n_set", "n_background", "n_set_total", "n_background_total",
    "frac_set", "frac_background", "fold_enrichment", "p_value", "fdr"
  ) %in% names(out)))
  # Distinct SNP-tissue links: the set has Adipose x2 and Blood x1; the
  # background (s3, s4) has Blood x2, and s4's Brain row has no gene_id so it
  # is not a tissue link.
  adipose <- out[out$tissue == "Adipose", ]
  expect_equal(adipose$n_set, 2L)
  expect_equal(adipose$n_background, 0L)
  expect_equal(adipose$n_set_total, 3L)
  expect_equal(adipose$n_background_total, 2L)
  expect_false("Brain" %in% out$tissue)
  expect_equal(out$tissue[1], "Adipose")
})


test_that("enrich_snp_set() removes the set from the background", {
  coloc_groups <- data.frame(
    variant_id = c("s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"),
    trait_category = c(
      "Metabolic", "Metabolic", "Metabolic", "Blood",
      "Blood", "Blood", "Blood", "Neurological"
    ),
    stringsAsFactors = FALSE
  )
  # Passing the set inside the background gives the same result as leaving it
  # out: the set is always compared with the trait's other SNPs.
  with_set <- enrich_snp_set(c("s1", "s2"), c("s1", "s2", "s3", "s4"),
                             coloc_groups, min_set_size = 1)
  without_set <- enrich_snp_set(c("s1", "s2"), c("s3", "s4"),
                                coloc_groups, min_set_size = 1)
  expect_equal(with_set, without_set)
  expect_true("trait_category" %in% names(with_set))
  expect_equal(nrow(enrich_snp_set("s1", c("s1", "s2"), coloc_groups)), 0L)
})


# The traits x SNPs matrix run_univariate_clustering() builds before filtering.
pleiotropy_matrix_for <- function(trait_id, coloc_groups) {
  locus_data <- gpmapr:::.prepare_pleiotropy_locus_data(trait_id, coloc_groups)
  return(gpmapr:::.finalize_pleiotropy_from_locus(locus_data, as.character(trait_id)))
}

test_that("the pleiotropy matrix keeps the target trait's SNPs from shared coloc_groups", {
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

  p1 <- pleiotropy_matrix_for(1L, coloc_groups)
  p2 <- pleiotropy_matrix_for(2L, coloc_groups)
  expect_equal(ncol(p1$x_matrix), 3L)
  expect_equal(ncol(p2$x_matrix), 2L)
  expect_true(setequal(colnames(p1$x_matrix), c("v1", "v2", "v3")))
  expect_true(setequal(colnames(p2$x_matrix), c("v4", "v5")))
  expect_true("feature_type" %in% names(p1$trait_info))
})

test_that("univariate input can exclude trans markers", {
  coloc_groups <- data.frame(
    trait_id = c(1L, 1L, 2L, 2L),
    trait_name = c("target", "target", "other", "other"),
    variant_id = c(101L, 102L, 101L, 102L),
    min_p = c(1e-10, 1e-10, 1e-10, 1e-10),
    beta = c(1, 1, 1, 1),
    se = c(1, 1, 1, 1),
    cis_trans = c("cis", "trans", "cis", "trans"),
    coloc_group_id = c(11L, 12L, 11L, 12L),
    display_snp = c("101", "102", "101", "102"),
    chr = c(1L, 1L, 1L, 1L),
    bp = c(101L, 102L, 101L, 102L),
    stringsAsFactors = FALSE
  )

  cis_only <- coloc_groups |>
    dplyr::filter(is.na(cis_trans) | tolower(cis_trans) != "trans")
  out <- pleiotropy_matrix_for(1L, cis_only)

  expect_identical(colnames(out$x_matrix), "101")
})

test_that("snp_similarity_matrix treats missing cells as zero effects", {
  x_matrix <- matrix(
    c(
      1, 1, NA,
      1, NA, 1
    ),
    nrow = 3,
    dimnames = list(
      paste0("trait", 1:3),
      c("snp1", "snp2")
    )
  )
  out <- snp_similarity_matrix(x_matrix)
  expect_equal(out$s_matrix["snp1", "snp2"], 0.5)
  expect_equal(out$overlap_matrix["snp1", "snp2"], 1)
  expect_true(out$eligible_matrix["snp1", "snp2"])
  expect_setequal(names(out), c("s_matrix", "col_norms", "overlap_matrix", "eligible_matrix"))
})

test_that("enrich_snp_set() is a one-sided Fisher test with BH FDR", {
  coloc_groups <- data.frame(
    variant_id = c("s1", "s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4", "s5", "s6"),
    gene_id = 1:11,
    tissue = c(
      "Adipose", "Adipose", "Brain", "Adipose", "Blood",
      "Blood", "Blood", "Blood", "Brain", "Brain", "Brain"
    ),
    trait_category = c(
      "Metabolic", "Metabolic", "Neuro", "Metabolic", "Blood",
      "Blood", "Blood", "Blood", "Neuro", "Neuro", "Neuro"
    ),
    stringsAsFactors = FALSE
  )
  set <- c("s1", "s2", "s3")
  out <- enrich_snp_set(set, paste0("s", 1:6), coloc_groups, min_set_size = 1)
  for (i in seq_len(nrow(out))) {
    tab <- matrix(c(
      out$n_set[i], out$n_set_total[i] - out$n_set[i],
      out$n_background[i], out$n_background_total[i] - out$n_background[i]
    ), nrow = 2)
    expect_equal(
      out$p_value[i],
      stats::fisher.test(tab, alternative = "greater")$p.value
    )
  }
  expect_equal(out$fdr, stats::p.adjust(out$p_value, method = "BH"))
  expect_false(any(c("fe_ci_lower", "fe_ci_upper") %in% names(out)))

  # A set with more links than its background still gets a p-value.
  big <- enrich_snp_set(paste0("s", 1:5), paste0("s", 1:6), coloc_groups,
                        min_set_size = 1)
  expect_true(all(is.finite(big$p_value)))
})

test_that("genes_at_snps can include situated genes", {
  coloc_groups <- data.frame(
    variant_id = "s1",
    coloc_group_id = 1L,
    gene_id = NA_integer_,
    gene = NA_character_,
    situated_gene_id = 42L,
    situated_gene = "GENE42",
    stringsAsFactors = FALSE
  )

  ordinary <- genes_at_snps("s1", coloc_groups)
  situated <- genes_at_snps("s1", coloc_groups, include_situated_gene = TRUE)

  expect_equal(nrow(ordinary), 0L)
  expect_equal(situated$gene_id, 42L)
  expect_equal(situated$gene, "GENE42")
})

test_that("enrich_snp_set_pathways() reports the top term per source", {
  genes <- all_genes()
  gene_ids <- genes[genes$gene %in% c("TREM2", "APOE", "FTO", "MC4R", "LPL"), c("id", "gene")]
  coloc_groups <- data.frame(
    variant_id = c("s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4", "s5", "s5"),
    coloc_group_id = 1:10,
    gene_id = rep(gene_ids$id, 2),
    gene = rep(gene_ids$gene, 2),
    stringsAsFactors = FALSE
  )

  out <- enrich_snp_set_pathways(paste0("s", 1:5), coloc_groups)
  expect_setequal(names(out), c("n_snps", "genes", "pathways", "top_term"))
  expect_equal(out$n_snps, 5L)
  expect_equal(names(out$top_term), c("KEGG", "Reactome"))
  expect_true(all(out$pathways$source %in% c("KEGG", "Reactome")))
  expect_true(is.na(out$top_term[["KEGG"]]) || grepl("^KEGG:", out$top_term[["KEGG"]]))
  expect_true(is.na(out$top_term[["Reactome"]]) || grepl("^Reactome:", out$top_term[["Reactome"]]))

  too_small <- enrich_snp_set_pathways(c("s1", "s2"), coloc_groups)
  expect_equal(nrow(too_small$pathways), 0L)
})

test_that("search_gpmap() orders results by coloc groups and rare results", {
  result <- search_gpmap("haemoglobin")
  importance <- result$num_coloc_groups + result$num_rare_results
  expect_false(is.unsorted(rev(importance)))
})

test_that("trait_duplicates() returns expected output", {
  result <- trait_duplicates()
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  # Older API versions name each duplicate's parent; newer ones list duplicate
  # groups with a representative. remove_duplicate_traits() reads both.
  parent_shape <- c("trait_id", "trait_name", "parent_trait_id", "parent_trait_name")
  group_shape <- c("trait_id", "trait_name", "duplicate_group", "is_representative")
  expect_true(all(parent_shape %in% names(result)) || all(group_shape %in% names(result)))
})

test_that("remove_duplicate_traits() drops duplicates and keeps parents", {
  duplicates <- data.frame(
    trait_id = c(11L, 12L, 21L),
    trait_name = c("A copy 1", "A copy 2", "B copy"),
    parent_trait_id = c(1L, 1L, 2L),
    parent_trait_name = c("A", "A", "B")
  )
  data <- data.frame(trait_id = c(1, 11, 12, 2, 21, 3, 1), value = 1:7)

  result <- remove_duplicate_traits(data, duplicates)
  expect_equal(result$trait_id, c(1, 2, 3, 1))

  # a kept duplicate stays, and its parent goes instead
  result <- remove_duplicate_traits(data, duplicates, keep_trait_ids = "21")
  expect_equal(result$trait_id, c(1, 21, 3, 1))

  expect_error(remove_duplicate_traits(data.frame(id = 1), duplicates), "trait_id")
})

test_that("remove_duplicate_traits() reads duplicate groups with a representative", {
  # The same duplicates as above, in the duplicate-group shape: each group's
  # representative is the parent, and trait 3 is in a group of its own.
  duplicates <- data.frame(
    trait_id = c(1L, 11L, 12L, 2L, 21L, 3L),
    trait_name = c("A", "A copy 1", "A copy 2", "B", "B copy", "C"),
    duplicate_group = c(1L, 1L, 1L, 2L, 2L, 3L),
    is_representative = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE)
  )
  data <- data.frame(trait_id = c(1, 11, 12, 2, 21, 3, 1), value = 1:7)

  result <- remove_duplicate_traits(data, duplicates)
  expect_equal(result$trait_id, c(1, 2, 3, 1))
  result <- remove_duplicate_traits(data, duplicates, keep_trait_ids = "21")
  expect_equal(result$trait_id, c(1, 21, 3, 1))

  expect_error(
    remove_duplicate_traits(data, data.frame(trait_id = 1L, group = 1L)),
    "duplicate_group"
  )
})

test_that("delete_gwas() validates its inputs", {
  guid <- "00000000-0000-0000-0000-000000000000"
  expect_error(delete_gwas("not-a-guid", "user@example.com"), "GUID")
  expect_error(delete_gwas(guid), "email is required")
  expect_error(delete_gwas(guid, ""), "email is required")
})

test_that("delete_gwas() surfaces API errors", {
  expect_error(
    delete_gwas("00000000-0000-0000-0000-000000000000", "user@example.com"),
    "not found"
  )
})

test_that("upload_gwas() validates column arguments", {
  file <- tempfile(fileext = ".tsv")
  writeLines(c("chr\tpos\tea\toa\tp\tbeta\tse", "1\t100\tA\tG\t1e-9\t0.1\t0.01"), file)
  on.exit(unlink(file))
  upload <- function(...) {
    upload_gwas(file, name = "test", email = "user@example.com", sample_size = 1000, ...)
  }

  expect_error(
    upload(chr_col = "chr", bp_col = "pos", ea_col = "ea", oa_col = "oa", beta_col = "beta", se_col = "se"),
    "Missing required column arguments: p_col"
  )
  expect_error(
    upload(chr_col = "chr", bp_col = "pos", ea_col = "ea", oa_col = "oa", p_col = "p", beta_col = "beta"),
    "Either beta_col and se_col"
  )
  expect_error(
    upload(chr_col = "chr", bp_col = "pos", ea_col = "ea", oa_col = "oa", p_col = "pval",
      beta_col = "beta", se_col = "se"),
    "Columns not found in file: pval"
  )
  expect_error(
    upload(column_names = list(CHR = "chr"), chr_col = "chr"),
    "not both"
  )
  expect_error(
    suppressWarnings(upload(column_names = list(CHR = "chr", POS = "pos"))),
    "Unknown column_names: POS"
  )
  expect_warning(
    expect_error(upload(column_names = list(CHR = "chr", BP = "pos")), "Missing required"),
    "deprecated"
  )
})

test_that("legacy column_names are mapped to API column names", {
  result <- suppressWarnings(build_gwas_column_names(
    list(chr = "chr", BP = "pos", EA = "ea", OA = "oa", P = "p", OR = "or", LB = "lb", UB = "ub"),
    list()
  ))
  expect_equal(names(result), c("CHR", "BP", "EA", "OA", "P", "OR", "OR_LB", "OR_UB"))
  expect_equal(result$OR_LB, "lb")
})
