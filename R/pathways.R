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


#' @title Genes Linked To SNPs Via Coloc Groups
#' @description Return molecular QTL genes mapped to a set of SNPs through
#' `coloc_groups` rows with non-missing `gene_id`.
#' @param snp_ids Character or numeric SNP identifiers matching `snp_key`.
#' @param coloc_groups Coloc-group dataframe from `trait()` /
#'   `traits(..., include_associations = TRUE)`.
#' @param snp_key Column used to match `snp_ids`: `"variant_id"`, `"display_snp"`,
#'   or `"coloc_group_id"`. Defaults to `"variant_id"`.
#' @return A dataframe with `snp_id`, `coloc_group_id`, `gene_id`, and `gene`
#'   (one row per SNP-gene link).
#' @export
genes_at_snps <- function(snp_ids,
                          coloc_groups,
                          snp_key = c("variant_id", "display_snp", "coloc_group_id")) {
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

  snp_ids <- unique(as.character(snp_ids))
  snp_col <- as.character(coloc_groups[[snp_key]])

  out <- coloc_groups |>
    dplyr::mutate(snp_id = snp_col) |>
    dplyr::filter(snp_id %in% snp_ids, !is.na(gene_id)) |>
    dplyr::distinct(snp_id, coloc_group_id, gene_id, gene) |>
    dplyr::arrange(snp_id, gene)

  return(out)
}


#' @title Baseline Pathway Enrichment For A Trait
#' @description Verification helper: collect molecular QTL genes at genome-wide
#' significant loci for a trait and test pathway enrichment (default KEGG +
#' Reactome). Use before clustering / EBMF to see whether the trait gene set
#' has recoverable pathway signal, and which SNPs those genes map to.
#' @param trait_id Numeric trait ID.
#' @param coloc_groups Optional coloc-group dataframe. If `NULL`, fetched via
#'   `trait(trait_id, include_associations = TRUE)$coloc_groups`.
#' @param p_threshold P-value threshold for target-trait SNPs. Defaults to 5e-8.
#' @param snp_key Column used to identify SNPs. Defaults to `"variant_id"`.
#' @param sources Character vector of pathway sources to query separately and
#'   bind. Defaults to `c("KEGG", "Reactome")`. Use `NULL` for a single
#'   unfiltered `pathway_enrichment()` call (all sources).
#' @param p_value_threshold FDR threshold passed to `pathway_enrichment()`.
#' @param minimum_count_in_network Minimum overlap passed to `pathway_enrichment()`.
#' @return A list with:
#'   \itemize{
#'     \item trait_id, n_snps, n_genes
#'     \item snp_genes: SNP-gene links from `genes_at_snps()`
#'     \item genes: distinct genes tested
#'     \item pathways: enriched pathway dataframe (possibly empty)
#'     \item summary: one-row overview
#'   }
#' @export
enrich_trait_pathways <- function(trait_id,
                                  coloc_groups = NULL,
                                  p_threshold = 5e-8,
                                  snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                  sources = c("KEGG", "Reactome"),
                                  p_value_threshold = 0.05,
                                  minimum_count_in_network = 2L) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }
  snp_key <- match.arg(snp_key)

  if (is.null(coloc_groups)) {
    coloc_groups <- trait(trait_id, include_associations = TRUE)$coloc_groups
  }
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("No coloc_groups data available")
  }

  target_id <- trait_id
  snp_col <- as.character(coloc_groups[[snp_key]])
  target_snps <- coloc_groups |>
    dplyr::mutate(snp_id = snp_col) |>
    dplyr::filter(
      trait_id == target_id,
      if (!is.null(p_threshold)) min_p <= p_threshold else TRUE
    ) |>
    dplyr::distinct(snp_id)

  snp_genes <- genes_at_snps(
    snp_ids = target_snps$snp_id,
    coloc_groups = coloc_groups,
    snp_key = snp_key
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

  summary_df <- data.frame(
    trait_id = trait_id,
    n_snps = nrow(target_snps),
    n_genes = nrow(genes),
    n_enriched_pathways = nrow(pathways),
    top_enriched_pathway = .top_pathway_label(pathways),
    stringsAsFactors = FALSE
  )

  return(list(
    trait_id = trait_id,
    n_snps = nrow(target_snps),
    n_genes = nrow(genes),
    snp_genes = snp_genes,
    genes = genes,
    pathways = pathways,
    summary = summary_df
  ))
}


#' @title Pathway Enrichment For SNP Groups
#' @description Verification helper: for each SNP grouping (Louvain module,
#' EBMF program, etc.) with more than `min_group_size` SNPs, collect molecular
#' QTL genes at those SNPs and test pathway enrichment. Compare against
#' `enrich_trait_pathways()` to see whether baseline pathways reappear and
#' whether they split cleanly across groups.
#' @param groups Either a named vector (`names` = SNP ids, values = group ids)
#'   or a dataframe with `snp_id` plus a group column (`group`, `cluster`, or
#'   `program`).
#' @param coloc_groups Coloc-group dataframe used to map SNPs to genes.
#' @param min_group_size Only enrich groups with more than this many SNPs.
#'   Defaults to 5.
#' @param snp_key Column used to match SNP ids in `coloc_groups`.
#' @param sources Pathway sources; see `enrich_trait_pathways()`.
#' @param p_value_threshold FDR threshold passed to `pathway_enrichment()`.
#' @param minimum_count_in_network Minimum overlap passed to `pathway_enrichment()`.
#' @return A list with:
#'   \itemize{
#'     \item by_group: list of per-group results (`group`, `n_snps`, `genes`,
#'       `snp_genes`, `pathways`)
#'     \item summary: one row per enriched group
#'   }
#' @export
enrich_snp_group_pathways <- function(groups,
                                      coloc_groups,
                                      min_group_size = 5L,
                                      snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                      sources = c("KEGG", "Reactome"),
                                      p_value_threshold = 0.05,
                                      minimum_count_in_network = 2L) {
  snp_key <- match.arg(snp_key)
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("coloc_groups is required")
  }
  if (!is.numeric(min_group_size) || min_group_size < 0) {
    stop("min_group_size must be a non-negative number")
  }

  group_df <- .normalize_snp_groups(groups)
  group_sizes <- table(group_df$group)
  large_groups <- names(group_sizes)[group_sizes > min_group_size]

  by_group <- lapply(large_groups, function(grp) {
    snp_ids <- group_df$snp_id[group_df$group == grp]
    snp_genes <- genes_at_snps(
      snp_ids = snp_ids,
      coloc_groups = coloc_groups,
      snp_key = snp_key
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
    data.frame(
      group = x$group,
      n_snps = x$n_snps,
      n_genes = nrow(x$genes),
      n_enriched_pathways = nrow(x$pathways),
      top_enriched_pathway = .top_pathway_label(x$pathways),
      stringsAsFactors = FALSE
    )
  }))

  return(list(
    by_group = by_group,
    summary = summary_df,
    min_group_size = as.integer(min_group_size)
  ))
}


#' @title Compare Trait Vs Group Pathway Enrichments
#' @description Summarise whether baseline trait pathways reappear in SNP-group
#' enrichments, whether they split across multiple groups, and which pathways
#' are group-specific.
#' @param trait_enrichment Output of `enrich_trait_pathways()`.
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
