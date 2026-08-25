#!/usr/bin/env Rscript
# Recompute the parameter-grid summary that backs the "Appendix: Choosing
# clustering parameters" sections of the investigation vignettes.
#
# Mirrors the compute of vignettes/investigation-univariate-analysis.Rmd (steps
# 1-4.5 plus module enrichments) for each grid combo and prints one TSV row per
# run with:
#   run gamma mss thr n_modules n_reliable reliable_snps mean_internal
#   mean_separation mean_connectedness n_category_modules n_tissue_modules
#   n_pathway_modules n_phenotype_modules
# followed by a per-module detail table (group, n_snps, top category / tissue /
# pathway / phenotype labels) used to write the "Notable modules" narrative.
#
# The pleiotropy matrix is built once and reused across combos (it only depends
# on associations, not on min_snp_signals / gamma / similarity_threshold).

pkg_root <- if (file.exists("../DESCRIPTION")) ".." else "."
suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_root, quiet = TRUE)
  } else {
    library(gpmapr)
  }
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
combo_spec <- if (length(args) > 0) args[1] else "all"

select_api("production")

target_trait_id <- 1992
min_module_size <- 3
max_snp_fraction <- 0.8
min_mean_internal <- 0.3
min_connectedness <- 0.5
compress_method <- "asinh"
associations <- "coloc"

combo_list <- list(
  c(gamma = 2, mss = 3, thr = 0),
  c(gamma = 2, mss = 3, thr = 0.2),
  c(gamma = 2, mss = 3, thr = 0.5),
  c(gamma = 2, mss = 5, thr = 0),
  c(gamma = 2, mss = 5, thr = 0.5),
  c(gamma = 3, mss = 3, thr = 0),
  c(gamma = 3, mss = 3, thr = 0.5),
  c(gamma = 3, mss = 5, thr = 0),
  c(gamma = 3, mss = 5, thr = 0.5)
)

if (combo_spec != "all") {
  parts <- strsplit(combo_spec, "_")[[1]]
  combo_list <- combo_list[vapply(
    combo_list,
    function(cb) {
      cb["gamma"] == as.numeric(sub("^g", "", parts[1])) &&
        cb["mss"] == as.numeric(sub("^mss", "", parts[2])) &&
        cb["thr"] == as.numeric(sub("^thr", "", parts[3]))
    },
    logical(1)
  )]
}

cat("step: fetching coloc groups\n", file = stderr())
cg <- trait(target_trait_id, include_associations = TRUE)$coloc_groups
cat("step: building pleiotropy matrix\n", file = stderr())
pleiotropy <- build_pleiotropy_matrix(
  trait_id = target_trait_id,
  coloc_groups = cg,
  snp_key = "variant_id",
  associations = associations
)
X_full <- pleiotropy$x_matrix

detail_rows <- list()

for (cb in combo_list) {
  gamma <- as.numeric(cb["gamma"])
  mss <- as.integer(cb["mss"])
  thr <- as.numeric(cb["thr"])
  run <- paste0("g", gamma, " mss", mss, " thr", thr)
  cat(">>> run:", run, "\n", file = stderr())

  X <- X_full
  trait_snp_counts <- rowSums(!is.na(X))
  trait_snp_frac <- trait_snp_counts / ncol(X)
  sparse_trait_ids <- names(trait_snp_counts)[trait_snp_counts < mss]
  ubiquitous_trait_ids <- names(trait_snp_frac)[trait_snp_frac > max_snp_fraction]
  drop_trait_ids <- setdiff(
    union(sparse_trait_ids, ubiquitous_trait_ids),
    as.character(target_trait_id)
  )
  X <- X[!rownames(X) %in% drop_trait_ids, , drop = FALSE]

  X_star <- orient_pleiotropy_matrix(X, target_trait_id = target_trait_id)$x_matrix
  X_star <- compress_effect_matrix(X_star, method = compress_method)

  similarity <- snp_similarity_matrix(X_star)
  S <- similarity$s_matrix

  clusters_louvain <- cluster_snp_profiles_louvain(
    S,
    similarity_threshold = thr,
    gamma = gamma
  )

  module_quality <- summarise_snp_module_quality(
    S,
    clusters_louvain$cluster,
    edge_threshold = thr,
    min_module_size = min_module_size,
    min_mean_internal = min_mean_internal,
    min_connectedness = min_connectedness
  )
  q <- module_quality |> dplyr::filter(n_snps >= min_module_size)
  reliable_q <- q |> dplyr::filter(reliable)

  reliable_ids <- module_quality$cluster[module_quality$reliable]
  clusters_analysis <- clusters_louvain$cluster[
    clusters_louvain$cluster %in% reliable_ids
  ]

  n_pw <- 0L
  n_ph <- 0L
  if (length(clusters_analysis) > 0) {
    baseline_snp_ids <- names(clusters_louvain$cluster)

    group_categories <- enrich_snp_group_trait_categories(
      groups = clusters_analysis,
      coloc_groups = cg,
      snp_key = "variant_id",
      min_group_size = min_module_size,
      baseline_snp_ids = baseline_snp_ids
    )
    group_tissues <- enrich_snp_group_tissues(
      groups = clusters_analysis,
      coloc_groups = cg,
      snp_key = "variant_id",
      min_group_size = min_module_size,
      background = "trait",
      baseline_snp_ids = baseline_snp_ids
    )
    module_enrichment <- enrich_snp_group_pathways(
      groups = clusters_analysis,
      coloc_groups = cg,
      min_group_size = min_module_size,
      sources = c("KEGG", "Reactome", "HP")
    )

    n_cat <- sum(!is.na(group_categories$summary$top_enriched_category_any))
    n_tis <- sum(!is.na(group_tissues$summary$top_enriched_tissue_any))
    n_pw <- sum(!is.na(module_enrichment$summary$top_enriched_pathway))
    n_ph <- sum(!is.na(module_enrichment$summary$top_enriched_phenotype))

    detail <- reliable_q |>
      dplyr::transmute(run = run, group = as.character(cluster), n_snps = n_snps) |>
      dplyr::full_join(
        group_categories$summary |>
          dplyr::mutate(group = as.character(group)) |>
          dplyr::select(group, top_enriched_category_any),
        by = "group"
      ) |>
      dplyr::full_join(
        group_tissues$summary |>
          dplyr::mutate(group = as.character(group)) |>
          dplyr::select(group, top_enriched_tissue_any),
        by = "group"
      ) |>
      dplyr::full_join(
        module_enrichment$summary |>
          dplyr::mutate(group = as.character(group)) |>
          dplyr::select(group, top_enriched_pathway, top_enriched_phenotype),
        by = "group"
      ) |>
      dplyr::mutate(n_snps = as.integer(n_snps))
    detail_rows[[length(detail_rows) + 1L]] <- detail
  } else {
    n_cat <- 0L
    n_tis <- 0L
  }

  cat(
    run,
    gamma, mss, thr,
    nrow(q), sum(q$reliable, na.rm = TRUE),
    sum(q$n_snps[q$reliable], na.rm = TRUE),
    round(mean(reliable_q$mean_internal_similarity, na.rm = TRUE), 3),
    round(mean(reliable_q$separation, na.rm = TRUE), 3),
    round(mean(reliable_q$connectedness, na.rm = TRUE), 3),
    n_cat, n_tis, n_pw, n_ph,
    sep = "\t"
  )
  cat("\n")
}

cat("\n### DETAIL ###\n")
detail_all <- dplyr::bind_rows(detail_rows)
cat("run\tgroup\tn_snps\ttop_category\ttop_tissue\ttop_pathway\ttop_phenotype\n")
for (i in seq_len(nrow(detail_all))) {
  d <- detail_all[i, ]
  fmt <- function(x) ifelse(is.na(x), "", x)
  cat(
    d$run, d$group, d$n_snps,
    fmt(d$top_enriched_category_any),
    fmt(d$top_enriched_tissue_any),
    fmt(d$top_enriched_pathway),
    fmt(d$top_enriched_phenotype),
    sep = "\t"
  )
  cat("\n")
}
