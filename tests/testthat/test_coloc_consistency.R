library(testthat)

# Tests of the coloc-consistency research functions (EBMF, enrichment, pathways, pleiotropy matrices,
# GWAS-upload stamping). The pathway tests call the GPMap API and skip without it (helper-live-api.R).

pleiotropy_matrix_for <- function(trait_id, coloc_groups) {
  locus_data <- gpmapr:::.prepare_pleiotropy_locus_data(trait_id, coloc_groups)
  return(gpmapr:::.finalize_pleiotropy_from_locus(locus_data, as.character(trait_id)))
}

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

test_that("GWAS upload own coloc rows are stamped with the lookup GUID", {
  guid <- "7a289615-c1b4-91f3-3d97-887f60de9155"
  # As returned by API 1.1.0: the upload's own rows have source "upload" and its gwas_upload_id (55); a
  # compared upload's rows have its own id (77); existing studies have their trait_id. The last existing
  # row's GPMap trait id happens to equal the upload id and must not be stamped.
  cg <- data.frame(
    source = c("upload", "existing", "upload", "existing", "upload"),
    gwas_upload_id = c(55, NA, 55, NA, 77),
    study_id = c(NA, 99, NA, 98, NA),
    study_extraction_id = c(1, 5, 2, 6, 3),
    trait_id = c(NA, 2527, NA, 55, NA),
    trait_name = c(NA, "Allergic rhinitis", NA, "Other trait", "Compared"),
    variant_id = c(1, 1, 2, 2, 1),
    stringsAsFactors = FALSE
  )
  associations <- data.frame(
    source = c("upload", "upload"),
    study_id = c(NA, NA),
    gwas_upload_id = c(55, 55),
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
  decorated <- .decorate_gwas_upload(upload, guid)
  groups <- decorated$coloc_groups
  expect_null(decorated$associations)
  expect_true(all(c("beta", "se") %in% names(groups)))
  own <- groups$source == "upload" & groups$gwas_upload_id %in% 55
  expect_equal(groups$trait_id[own], c(guid, guid))
  expect_equal(groups$trait_name[own], c("Atopic", "Atopic"))
  expect_equal(groups$beta[own], c(0.2, -0.1))
  expect_equal(groups$trait_id[groups$source == "existing"], c("2527", "55"))
  expect_true(is.na(groups$trait_id[groups$gwas_upload_id %in% 77]))
  expect_equal(decorated$trait$guid, guid)
})

test_that("pathway_enrichment() returns expected output with gene IDs", {
  skip_if_no_api()
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
  skip_if_no_api()
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
  skip_if_no_api()
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
