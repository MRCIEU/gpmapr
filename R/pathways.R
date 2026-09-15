#' @title Pathway Enrichment
#' @description Test a set of genes for pathway enrichment using the GPMap pathway database.
#' @param genes A vector of numeric gene IDs (from `all_genes()$id` or coloc group `gene_id`
#'   values) or gene names (e.g. `"APOE"`).
#' @param source Optional pathway source to filter by: `"Reactome"`, `"KEGG"`, or `"HP"`.
#' @param p_value_threshold FDR-adjusted p-value threshold for filtering results. Defaults to 0.05.
#' @param minimum_count_in_network Optional minimum number of input genes that must overlap a
#'   pathway for it to be tested. If `NULL`, the API default is used.
#' @return A list with:
#'   \itemize{
#'     \item results: a dataframe of enriched pathways with columns term_id, source,
#'       description, pathway_size, background_size, overlap, p_value, fdr, gene_ids
#'       (input genes overlapping the pathway), and pathway_gene_ids (all genes in the pathway)
#'     \item input_gene_count: number of genes submitted
#'     \item matched_gene_count: number of genes matched in the pathway database
#'     \item source: pathway source filter applied (if any)
#'     \item p_value_threshold: FDR threshold used
#'     \item minimum_count_in_network: minimum overlap threshold used (if returned by API)
#'     \item total_terms_tested: total pathway terms tested
#'   }
#' @export
pathway_enrichment <- function(genes,
                               source = NULL,
                               p_value_threshold = 0.05,
                               minimum_count_in_network = NULL) {
  if (is.null(genes) || length(genes) == 0) {
    stop("genes is required")
  }
  if (any(is.na(genes))) {
    stop("genes must not contain NA values")
  }
  if (!is.numeric(genes) && !is.character(genes)) {
    stop("genes must be numeric gene IDs or character gene names")
  }
  if (is.character(genes) && any(genes == "")) {
    stop("genes must not contain empty strings")
  }

  genes <- unique(genes)

  if (!is.null(source)) {
    valid_sources <- c("Reactome", "KEGG", "HP")
    if (!source %in% valid_sources) {
      stop("source must be one of: ", paste(valid_sources, collapse = ", "))
    }
  }

  if (!is.numeric(p_value_threshold) || p_value_threshold <= 0 || p_value_threshold > 1) {
    stop("p_value_threshold must be a number between 0 and 1")
  }

  if (!is.null(minimum_count_in_network)) {
    if (!is.numeric(minimum_count_in_network) ||
        minimum_count_in_network < 1 ||
        minimum_count_in_network != as.integer(minimum_count_in_network)) {
      stop("minimum_count_in_network must be a positive integer")
    }
  }

  return(pathway_enrichment_api(
    genes = genes,
    source = source,
    p_value_threshold = p_value_threshold,
    minimum_count_in_network = minimum_count_in_network
  ))
}


#' @title Pathway Gene Mappings
#' @description Fetch the full gene x pathway membership universe (pathway
#' sizes and every gene belonging to each pathway), for local continuous
#' enrichment tests that need the complete pathway gene sets rather than only
#' the genes returned by `pathway_enrichment()`'s overrepresentation test.
#' @param source Optional pathway source to filter by: `"Reactome"`, `"KEGG"`, or `"HP"`.
#' @return A list with:
#'   \itemize{
#'     \item sizes: dataframe of term_id, source, description, pathway_size, background_size
#'     \item mappings: dataframe of gene_id, term_id, source, description
#'       (one row per gene-pathway membership)
#'   }
#' @export
pathway_mappings <- function(source = NULL) {
  if (!is.null(source)) {
    valid_sources <- c("Reactome", "KEGG", "HP")
    if (!source %in% valid_sources) {
      stop("source must be one of: ", paste(valid_sources, collapse = ", "))
    }
  }

  response <- pathway_mappings_api(source = source)
  sizes <- response$sizes
  mappings <- response$mappings

  if (!is.data.frame(sizes) || nrow(sizes) == 0) {
    sizes <- data.frame(
      term_id = character(0), source = character(0), description = character(0),
      pathway_size = integer(0), background_size = integer(0),
      stringsAsFactors = FALSE
    )
  }
  if (!is.data.frame(mappings) || nrow(mappings) == 0) {
    mappings <- data.frame(
      gene_id = integer(0), term_id = character(0), source = character(0),
      description = character(0), stringsAsFactors = FALSE
    )
  } else {
    mappings$term_id <- as.character(mappings$term_id)
  }

  return(list(sizes = sizes, mappings = mappings))
}


#' @title Genes Linked To SNPs Via Coloc Groups
#' @description Return molecular QTL genes mapped to a set of SNPs through
#' `coloc_groups` rows with non-missing `gene_id`.
#' @param snp_ids Character or numeric SNP identifiers matching `snp_key`.
#' @param coloc_groups Coloc-group dataframe from `trait()` /
#'   `traits(..., include_associations = TRUE)`.
#' @param snp_key Column used to match `snp_ids`: `"variant_id"`, `"display_snp"`,
#'   or `"coloc_group_id"`. Defaults to `"variant_id"`.
#' @param include_situated_gene Include `situated_gene_id`/`situated_gene` links
#'   in addition to `gene_id`/`gene`. Defaults to `FALSE`.
#' @return A dataframe with `snp_id`, `coloc_group_id`, `gene_id`, and `gene`
#'   (one row per SNP-gene link).
#' @export
genes_at_snps <- function(snp_ids,
                          coloc_groups,
                          snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                          include_situated_gene = FALSE) {
  snp_key <- match.arg(snp_key)
  if (is.null(snp_ids) || length(snp_ids) == 0) {
    return(data.frame(
      snp_id = character(0),
      coloc_group_id = integer(0),
      gene_id = integer(0),
      gene = character(0),
      stringsAsFactors = FALSE
    ))
  }
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("coloc_groups is required")
  }
  if (!snp_key %in% names(coloc_groups)) {
    stop("coloc_groups must include column: ", snp_key)
  }
  if (!all(c("gene_id", "gene", "coloc_group_id") %in% names(coloc_groups))) {
    stop("coloc_groups must include gene_id, gene, and coloc_group_id")
  }
  if (!is.logical(include_situated_gene) || length(include_situated_gene) != 1L ||
      is.na(include_situated_gene)) {
    stop("include_situated_gene must be TRUE or FALSE")
  }

  snp_ids <- unique(as.character(snp_ids))
  snp_col <- as.character(coloc_groups[[snp_key]])

  out <- coloc_groups |>
    dplyr::mutate(snp_id = snp_col) |>
    dplyr::filter(snp_id %in% snp_ids, !is.na(gene_id)) |>
    dplyr::distinct(snp_id, coloc_group_id, gene_id, gene) |>
    dplyr::arrange(snp_id, gene)

  if (include_situated_gene && all(c("situated_gene_id", "situated_gene") %in% names(coloc_groups))) {
    situated <- coloc_groups |>
      dplyr::mutate(snp_id = snp_col) |>
      dplyr::filter(snp_id %in% snp_ids, !is.na(situated_gene_id)) |>
      dplyr::transmute(
        snp_id, coloc_group_id,
        gene_id = situated_gene_id,
        gene = situated_gene
      ) |>
      dplyr::distinct()
    out <- dplyr::bind_rows(out, situated) |>
      dplyr::distinct(snp_id, coloc_group_id, gene_id, gene) |>
      dplyr::arrange(snp_id, gene)
  }

  return(out)
}


#' @title Pathway Enrichment For SNP Groups
#' @description Verification helper: for each SNP grouping (EBMF program, etc.)
#' with more than `min_group_size` SNPs, collect molecular
#' QTL genes at those SNPs and test pathway enrichment. Compare the per-group
#' results to a trait-level baseline to see whether baseline pathways reappear
#' and whether they split cleanly across groups.
#' @param groups Either a named vector (`names` = SNP ids, values = group ids)
#'   or a dataframe with `snp_id` plus a group column (`group`, `cluster`, or
#'   `program`).
#' @param coloc_groups Coloc-group dataframe used to map SNPs to genes.
#' @param min_group_size Only enrich groups with at least this many SNPs.
#'   Defaults to 5.
#' @param snp_key Column used to match SNP ids in `coloc_groups`.
#' @param include_situated_gene Include situated-gene links in addition to
#'   ordinary gene links. Defaults to `FALSE`.
#' @param sources Pathway sources to query, e.g. `c("KEGG", "Reactome")`.
#'   Reactome and
#'   KEGG hits are summarised separately (and together as `top_enriched_pathway`
#'   for backwards compatibility); HP (Human Phenotype Ontology) hits are
#'   summarised separately as phenotypes.
#' @param p_value_threshold FDR threshold passed to `pathway_enrichment()`.
#' @param minimum_count_in_network Minimum overlap passed to `pathway_enrichment()`.
#' @return A list with:
#'   \itemize{
#'     \item by_group: list of per-group results (`group`, `n_snps`, `genes`,
#'       `snp_genes`, `pathways`)
#'     \item summary: one row per enriched group. `n_enriched_pathways` /
#'       `top_enriched_pathway` combine the non-HP sources (KEGG + Reactome);
#'       `n_enriched_reactome` / `top_enriched_reactome` and
#'       `n_enriched_kegg` / `top_enriched_kegg` report each source separately;
#'       `n_enriched_phenotypes` / `top_enriched_phenotype` refer to HP terms
#'       only (NA / 0 when `sources` has no HP).
#'   }
#' @export
enrich_snp_group_pathways <- function(groups,
                                      coloc_groups,
                                      min_group_size = 5L,
                                      snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                      sources = c("KEGG", "Reactome"),
                                      p_value_threshold = 0.05,
                                      minimum_count_in_network = 2L,
                                      include_situated_gene = FALSE) {
  snp_key <- match.arg(snp_key)
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("coloc_groups is required")
  }
  if (!is.numeric(min_group_size) || min_group_size < 0) {
    stop("min_group_size must be a non-negative number")
  }

  group_df <- .normalize_snp_groups(groups)
  group_sizes <- table(group_df$group)
  large_groups <- names(group_sizes)[group_sizes >= min_group_size]

  empty_summary <- data.frame(
    group = character(0),
    n_snps = integer(0),
    n_genes = integer(0),
    n_enriched_pathways = integer(0),
    top_enriched_pathway = character(0),
    n_enriched_reactome = integer(0),
    top_enriched_reactome = character(0),
    n_enriched_kegg = integer(0),
    top_enriched_kegg = character(0),
    n_enriched_phenotypes = integer(0),
    top_enriched_phenotype = character(0),
    stringsAsFactors = FALSE
  )

  if (length(large_groups) == 0) {
    return(list(
      by_group = list(),
      summary = empty_summary,
      min_group_size = as.integer(min_group_size)
    ))
  }

  by_group <- lapply(large_groups, function(grp) {
    snp_ids <- group_df$snp_id[group_df$group == grp]
    snp_genes <- genes_at_snps(
      snp_ids = snp_ids,
      coloc_groups = coloc_groups,
      snp_key = snp_key,
      include_situated_gene = include_situated_gene
    )
    genes <- snp_genes |>
      dplyr::distinct(gene_id, gene) |>
      dplyr::arrange(gene)
    pathways <- .enrich_pathway_sources(
      gene_ids = genes$gene_id,
      sources = sources,
      p_value_threshold = p_value_threshold,
      minimum_count_in_network = minimum_count_in_network
    )
    list(
      group = grp,
      n_snps = length(unique(snp_ids)),
      genes = genes,
      snp_genes = snp_genes,
      pathways = pathways
    )
  })

  summary_df <- dplyr::bind_rows(lapply(by_group, function(x) {
    subsets <- .pathway_source_subsets(x$pathways)
    data.frame(
      group = x$group,
      n_snps = x$n_snps,
      n_genes = nrow(x$genes),
      n_enriched_pathways = nrow(subsets$reactome) + nrow(subsets$kegg),
      top_enriched_pathway = .top_pathway_label(
        rbind(subsets$reactome, subsets$kegg)
      ),
      n_enriched_reactome = nrow(subsets$reactome),
      top_enriched_reactome = .top_pathway_label(subsets$reactome),
      n_enriched_kegg = nrow(subsets$kegg),
      top_enriched_kegg = .top_pathway_label(subsets$kegg),
      n_enriched_phenotypes = nrow(subsets$phenotype),
      top_enriched_phenotype = .top_pathway_label(subsets$phenotype),
      stringsAsFactors = FALSE
    )
  }))

  return(list(
    by_group = by_group,
    summary = summary_df,
    min_group_size = as.integer(min_group_size)
  ))
}


#' @title Continuous Pathway Enrichment From EBMF Loadings
#' @description Parallel to `enrich_snp_group_pathways()`, but without a hard
#' lFSR/magnitude cutoff: every SNP with a finite EBMF loading contributes to
#' every program's test, weighted by its squared loading, and is regressed
#' against each pathway's local 0/1 gene-membership indicator
#' (`lm(loading^2 ~ pathway_membership)`), following the approach in
#' `scripts/enrichment_via_loadings_example.r`. Pathway gene sets come from
#' `pathway_mappings()` (the full local gene x pathway universe) rather than
#' the overrepresentation-only `pathway_enrichment()` API. Every pathway term
#' tested for every program and every kept `source` (KEGG and Reactome pooled
#' together) is one Benjamini-Hochberg family, so `fdr` reflects the true
#' number of regressions run for this trait -- not one program's or one
#' source's slice of them. Trait-category, tissue, and pathway enrichment
#' (`enrich_program_loadings_trait_categories()`,
#' `enrich_program_loadings_tissues()`, this function) are each their own
#' separate family; none of the three is pooled with another.
#'
#' IMPORTANT: as in that script, SNPs are treated as independent
#' observations; for formal inference, LD should be accounted for (e.g. SNP
#' -> gene aggregation, or cluster-robust SEs / permutation using LD blocks).
#' The enrichment effect estimates are nevertheless useful for comparison
#' against the hard-cutoff hypergeometric pathway test.
#' @param clustering_result Result of `run_univariate_clustering()`.
#' @param coloc_groups Coloc-group dataframe used to map SNPs to genes.
#' @param mappings Optional pre-fetched `pathway_mappings()` result, reused
#'   across programs/traits instead of re-querying the API. Defaults to
#'   calling `pathway_mappings(source = NULL)`.
#' @param snp_key Column used to match SNP ids in `coloc_groups`.
#' @param sources Pathway sources to keep, e.g. `c("KEGG", "Reactome")`. `NULL`
#'   keeps every source present in `mappings`. Sources are pooled into one FDR
#'   family, not corrected independently.
#' @param min_category_size Only test a pathway if at least this many SNPs
#'   map (via `genes_at_snps()`) to one of its genes. Defaults to 5.
#' @param min_loading_magnitude Soft membership gate: a SNP only counts
#'   towards a pathway's `x = 1` group if its absolute loading on that
#'   program exceeds this value; otherwise it falls back to `x = 0` for that
#'   program's test (it stays in the analysis, just not counted as a
#'   member). This keeps near-zero, shrinkage-only loadings from diluting
#'   the mean squared-loading contrast without imposing a hard
#'   lFSR/magnitude gate. Defaults to 0.02.
#' @param include_situated_gene Include situated-gene links in addition to
#'   ordinary gene links when mapping SNPs to genes. Defaults to `FALSE`.
#' @return A list with:
#'   \itemize{
#'     \item by_program: list of per-program results (`program`, `n_snps`,
#'       `comparison` with columns term_id, source, description, enrichment,
#'       se, z, p, fdr, n_snps, n_category_snps; `fdr` is corrected across
#'       every program and source tested -- `source` is a label only)
#'     \item summary: one row per program (`n_pathways_tested`, `n_enriched`,
#'       `top_pathway`)
#'     \item mappings: the `pathway_mappings()` result used
#'   }
#' @export
enrich_program_loadings_pathways <- function(clustering_result,
                                             coloc_groups,
                                             mappings = NULL,
                                             snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                             sources = c("KEGG", "Reactome"),
                                             min_category_size = 5L,
                                             min_loading_magnitude = 0,
                                             include_situated_gene = FALSE) {
  snp_key <- match.arg(snp_key)
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("coloc_groups is required")
  }
  if (is.null(mappings)) {
    mappings <- pathway_mappings(source = NULL)
  }
  mapping_rows <- mappings$mappings
  if (!is.null(sources)) {
    mapping_rows <- mapping_rows[mapping_rows$source %in% sources, , drop = FALSE]
  }

  posterior <- ebmf_posterior_table(clustering_result)
  posterior <- posterior[is.finite(posterior$loading), , drop = FALSE]

  empty_summary <- data.frame(
    program = integer(0), n_snps = integer(0), n_pathways_tested = integer(0),
    n_enriched = integer(0), top_pathway = character(0),
    stringsAsFactors = FALSE
  )
  if (nrow(posterior) == 0 || nrow(mapping_rows) == 0) {
    return(list(by_program = list(), summary = empty_summary, mappings = mappings))
  }

  snp_ids <- unique(posterior$snp_id)
  snp_pathway <- .snp_pathway_matrix(
    snp_ids = snp_ids,
    coloc_groups = coloc_groups,
    snp_key = snp_key,
    mappings = mapping_rows,
    include_situated_gene = include_situated_gene
  )
  pathway_labels <- mapping_rows |>
    dplyr::distinct(term_id, source, description)
  split_sources <- if (!is.null(sources)) sources else sort(unique(mapping_rows$source))

  empty_comparison <- data.frame(
    term_id = character(0), source = character(0), description = character(0),
    enrichment = numeric(0), se = numeric(0), z = numeric(0), p = numeric(0),
    fdr = numeric(0), n_snps = integer(0), n_category_snps = integer(0),
    stringsAsFactors = FALSE
  )

  # Every program's every source is tested first, without correcting for
  # multiplicity, so the pooled BH family below spans every (program x source
  # x term) test actually run -- KEGG and Reactome share one family, and that
  # family spans every program, not just one program's or one source's slice.
  programs <- sort(unique(posterior$program))
  raw_by_program <- lapply(programs, function(pg) {
    prog_loadings <- posterior[posterior$program == pg, , drop = FALSE]
    y <- prog_loadings$loading[match(snp_ids, prog_loadings$snp_id)]

    by_source <- lapply(split_sources, function(src) {
      src_cols <- intersect(
        colnames(snp_pathway),
        pathway_labels$term_id[pathway_labels$source == src]
      )
      if (length(src_cols) == 0) {
        return(NULL)
      }
      res <- .continuous_link_enrichment(
        y, snp_pathway[, src_cols, drop = FALSE],
        min_category_size = min_category_size,
        min_loading_magnitude = min_loading_magnitude
      )
      if (nrow(res) == 0) {
        return(NULL)
      }
      names(res)[names(res) == "value"] <- "term_id"
      res$source <- src
      return(res)
    })
    comparison <- dplyr::bind_rows(by_source)
    if (nrow(comparison) > 0) {
      comparison <- comparison |>
        dplyr::left_join(pathway_labels, by = c("term_id", "source")) |>
        dplyr::select(
          term_id, source, description, enrichment, se, z, p,
          n_snps, n_category_snps
        )
    } else {
      comparison <- empty_comparison[, setdiff(names(empty_comparison), "fdr")]
    }
    return(list(program = pg, n_snps = sum(is.finite(y)), comparison = comparison))
  })

  pooled_p <- unlist(lapply(raw_by_program, function(x) x$comparison$p))
  pooled_fdr <- stats::p.adjust(pooled_p, method = "BH")
  offset <- 0L
  by_program <- lapply(raw_by_program, function(x) {
    n <- nrow(x$comparison)
    comparison <- x$comparison
    if (n > 0) {
      comparison$fdr <- pooled_fdr[(offset + 1L):(offset + n)]
      comparison <- comparison |>
        dplyr::select(
          term_id, source, description, enrichment, se, z, p, fdr,
          n_snps, n_category_snps
        ) |>
        dplyr::arrange(fdr, p)
    } else {
      comparison$fdr <- numeric(0)
    }
    offset <<- offset + n
    return(list(program = x$program, n_snps = x$n_snps, comparison = comparison))
  })

  summary_df <- dplyr::bind_rows(lapply(by_program, function(x) {
    top_pathway <- if (nrow(x$comparison) > 0) {
      paste0(x$comparison$source[1], ": ", x$comparison$description[1])
    } else {
      NA_character_
    }
    data.frame(
      program = x$program,
      n_snps = x$n_snps,
      n_pathways_tested = nrow(x$comparison),
      n_enriched = sum(x$comparison$fdr <= 0.05, na.rm = TRUE),
      top_pathway = top_pathway,
      stringsAsFactors = FALSE
    )
  }))

  return(list(by_program = by_program, summary = summary_df, mappings = mappings))
}


#' @title Compare Trait Vs Group Pathway Enrichments
#' @description Summarise whether baseline trait pathways reappear in SNP-group
#' enrichments, whether they split across multiple groups, and which pathways
#' are group-specific.
#' @param trait_enrichment Trait-level enrichment list with a `pathways`
#'   dataframe (columns `source`, `term_id`, `description`).
#' @param group_enrichment Output of `enrich_snp_group_pathways()`.
#' @return A list with:
#'   \itemize{
#'     \item pathway_status: one row per pathway key seen in baseline and/or
#'       groups (`in_baseline`, `n_groups`, `groups`, `status`)
#'     \item group_overlap: per-group counts of shared vs novel pathways
#'   }
#' @export
compare_group_pathways <- function(trait_enrichment, group_enrichment) {
  if (is.null(trait_enrichment$pathways) || is.null(group_enrichment$by_group)) {
    stop("trait_enrichment and group_enrichment must be enrichment helper outputs")
  }

  baseline_keys <- .pathway_keys(trait_enrichment$pathways)
  group_rows <- dplyr::bind_rows(lapply(group_enrichment$by_group, function(x) {
    if (nrow(x$pathways) == 0) {
      return(data.frame(
        group = character(0),
        pathway_key = character(0),
        description = character(0),
        stringsAsFactors = FALSE
      ))
    }
    x$pathways |>
      dplyr::mutate(
        group = as.character(x$group),
        pathway_key = paste(source, term_id, sep = ":")
      ) |>
      dplyr::select(group, pathway_key, description)
  }))

  if (nrow(group_rows) == 0) {
    group_rows <- data.frame(
      group = character(0),
      pathway_key = character(0),
      description = character(0),
      stringsAsFactors = FALSE
    )
  }

  all_keys <- union(baseline_keys, unique(group_rows$pathway_key))
  if (length(all_keys) == 0) {
    return(list(
      pathway_status = data.frame(
        pathway_key = character(0),
        description = character(0),
        in_baseline = logical(0),
        n_groups = integer(0),
        groups = character(0),
        status = character(0),
        stringsAsFactors = FALSE
      ),
      group_overlap = dplyr::bind_rows(lapply(group_enrichment$by_group, function(x) {
        data.frame(
          group = x$group,
          n_enriched_pathways = 0L,
          n_shared_with_baseline = 0L,
          n_novel = 0L,
          stringsAsFactors = FALSE
        )
      }))
    ))
  }

  pathway_status <- dplyr::bind_rows(lapply(all_keys, function(key) {
    groups_hit <- sort(unique(group_rows$group[group_rows$pathway_key == key]))
    in_baseline <- key %in% baseline_keys
    n_groups <- length(groups_hit)
    status <- if (in_baseline && n_groups == 0) {
      "baseline_only"
    } else if (in_baseline && n_groups == 1) {
      "recovered_one_group"
    } else if (in_baseline && n_groups > 1) {
      "split_across_groups"
    } else if (!in_baseline && n_groups == 1) {
      "group_specific"
    } else {
      "group_specific_multi"
    }
    desc <- group_rows$description[group_rows$pathway_key == key][1]
    if (is.na(desc) || identical(desc, character(0))) {
      bp <- trait_enrichment$pathways
      if (nrow(bp) > 0) {
        bp_key <- paste(bp$source, bp$term_id, sep = ":")
        desc <- bp$description[bp_key == key][1]
      }
    }
    data.frame(
      pathway_key = key,
      description = desc,
      in_baseline = in_baseline,
      n_groups = n_groups,
      groups = paste(groups_hit, collapse = ", "),
      status = status,
      stringsAsFactors = FALSE
    )
  }))

  if (nrow(pathway_status) > 0) {
    pathway_status <- pathway_status |>
      dplyr::arrange(status, dplyr::desc(n_groups), pathway_key)
  }

  group_overlap <- dplyr::bind_rows(lapply(group_enrichment$by_group, function(x) {
    keys <- .pathway_keys(x$pathways)
    data.frame(
      group = x$group,
      n_enriched_pathways = length(keys),
      n_shared_with_baseline = sum(keys %in% baseline_keys),
      n_novel = sum(!keys %in% baseline_keys),
      stringsAsFactors = FALSE
    )
  }))

  return(list(
    pathway_status = pathway_status,
    group_overlap = group_overlap
  ))
}


.normalize_snp_groups <- function(groups) {
  if (is.data.frame(groups)) {
    snp_col <- intersect(c("snp_id", "variant_id"), names(groups))
    group_col <- intersect(c("group", "cluster", "program"), names(groups))
    if (length(snp_col) == 0 || length(group_col) == 0) {
      stop(
        "groups dataframe must include snp_id (or variant_id) and ",
        "group, cluster, or program"
      )
    }
    return(data.frame(
      snp_id = as.character(groups[[snp_col[1]]]),
      group = as.character(groups[[group_col[1]]]),
      stringsAsFactors = FALSE
    ))
  }

  if (is.null(names(groups))) {
    stop("groups vector must be named by SNP id")
  }

  return(data.frame(
    snp_id = as.character(names(groups)),
    group = as.character(groups),
    stringsAsFactors = FALSE
  ))
}


.enrich_pathway_sources <- function(gene_ids,
                                    sources,
                                    p_value_threshold,
                                    minimum_count_in_network) {
  if (length(gene_ids) == 0) {
    return(.empty_pathway_results())
  }

  if (is.null(sources)) {
    enrichment <- pathway_enrichment(
      genes = gene_ids,
      source = NULL,
      p_value_threshold = p_value_threshold,
      minimum_count_in_network = minimum_count_in_network
    )
    return(.format_pathway_results(enrichment$results))
  }

  rows <- lapply(sources, function(src) {
    enrichment <- pathway_enrichment(
      genes = gene_ids,
      source = src,
      p_value_threshold = p_value_threshold,
      minimum_count_in_network = minimum_count_in_network
    )
    .format_pathway_results(enrichment$results)
  })

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    return(out)
  }
  return(dplyr::arrange(out, fdr, source, term_id))
}


.format_pathway_results <- function(results) {
  if (!is.data.frame(results) || nrow(results) == 0) {
    return(.empty_pathway_results())
  }

  out <- results |>
    dplyr::mutate(
      input_genes = vapply(
        gene_ids,
        function(g) paste(g, collapse = ", "),
        character(1)
      )
    ) |>
    dplyr::select(
      source, term_id, description, overlap, p_value, fdr, input_genes
    )

  return(out)
}


.empty_pathway_results <- function() {
  return(data.frame(
    source = character(0),
    term_id = character(0),
    description = character(0),
    overlap = integer(0),
    p_value = numeric(0),
    fdr = numeric(0),
    input_genes = character(0),
    stringsAsFactors = FALSE
  ))
}


.pathway_source_subsets <- function(pathways) {
  if (is.null(pathways) || nrow(pathways) == 0) {
    empty <- .empty_pathway_results()
    return(list(reactome = empty, kegg = empty, phenotype = empty))
  }
  list(
    reactome = pathways[pathways$source == "Reactome", , drop = FALSE],
    kegg = pathways[pathways$source == "KEGG", , drop = FALSE],
    phenotype = pathways[pathways$source == "HP", , drop = FALSE]
  )
}


.top_pathway_label <- function(pathways) {
  if (!is.data.frame(pathways) || nrow(pathways) == 0) {
    return(NA_character_)
  }
  return(paste0(pathways$source[1], ": ", pathways$description[1]))
}


.pathway_keys <- function(pathways) {
  if (!is.data.frame(pathways) || nrow(pathways) == 0) {
    return(character(0))
  }
  return(paste(pathways$source, pathways$term_id, sep = ":"))
}


# Gene x pathway 0/1 membership matrix restricted to `gene_ids`, built from a
# `pathway_mappings()$mappings`-shaped dataframe (gene_id, term_id, ...).
.gene_pathway_matrix <- function(gene_ids, mappings) {
  gene_ids <- unique(gene_ids)
  pathway_ids <- unique(mappings$term_id)
  mat <- matrix(
    0L,
    nrow = length(gene_ids), ncol = length(pathway_ids),
    dimnames = list(as.character(gene_ids), pathway_ids)
  )
  if (length(gene_ids) == 0 || length(pathway_ids) == 0) {
    return(mat)
  }
  hits <- mappings[mappings$gene_id %in% gene_ids, , drop = FALSE]
  if (nrow(hits) > 0) {
    mat[cbind(as.character(hits$gene_id), hits$term_id)] <- 1L
  }
  return(mat)
}


# SNP x pathway 0/1 membership matrix: a SNP is "in" a pathway if any gene
# mapped to it (via genes_at_snps()) belongs to that pathway.
.snp_pathway_matrix <- function(snp_ids,
                                coloc_groups,
                                snp_key,
                                mappings,
                                include_situated_gene = FALSE) {
  snp_ids <- unique(as.character(snp_ids))
  pathway_ids <- unique(mappings$term_id)
  empty <- matrix(
    0L,
    nrow = length(snp_ids), ncol = length(pathway_ids),
    dimnames = list(snp_ids, pathway_ids)
  )
  if (length(snp_ids) == 0 || length(pathway_ids) == 0) {
    return(empty)
  }

  snp_genes <- genes_at_snps(
    snp_ids = snp_ids,
    coloc_groups = coloc_groups,
    snp_key = snp_key,
    include_situated_gene = include_situated_gene
  )
  if (nrow(snp_genes) == 0) {
    return(empty)
  }

  gene_ids <- unique(snp_genes$gene_id)
  gene_mat <- .gene_pathway_matrix(gene_ids, mappings)

  snp_gene <- snp_genes |>
    dplyr::distinct(snp_id, gene_id) |>
    dplyr::mutate(gene_id = as.character(gene_id))
  snp_gene_mat <- matrix(
    0L,
    nrow = length(snp_ids), ncol = length(gene_ids),
    dimnames = list(snp_ids, as.character(gene_ids))
  )
  snp_gene_mat[cbind(snp_gene$snp_id, snp_gene$gene_id)] <- 1L

  mat <- (snp_gene_mat %*% gene_mat > 0) * 1L
  dimnames(mat) <- list(snp_ids, pathway_ids)
  return(mat)
}
