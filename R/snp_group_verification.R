#' @title Baseline Top Associated Traits
#' @description Rank complex-trait columns by mean absolute effect across a set
#' of SNPs. Use as a verification baseline when pathway enrichment is empty or
#' weak: do SNP groups recover distinct trait drivers?
#' @param trait_matrix SNPs x traits numeric matrix (e.g. oriented trait matrix
#'   from EBMF / verification helpers).
#' @param snp_ids Optional SNP ids to include. Defaults to all rows.
#' @param trait_info Optional dataframe with `feature_id`/`trait_id` and
#'   `feature_name`/`trait_name` columns for labels.
#' @param coloc_groups Optional coloc-group dataframe used to attach
#'   `trait_category` labels (and for category summaries).
#' @param exclude_trait_ids Optional trait ids to drop (e.g. the target trait).
#' @param top_n Maximum traits to return. Defaults to 20.
#' @param min_snps_with_signal Drop traits observed in fewer than this many SNPs
#'   in the summarised set. Defaults to 5.
#' @param n_categories Number of leading trait categories to report (by summed
#'   `abs_mean_z`). Defaults to 2.
#' @return A list with:
#'   \itemize{
#'     \item traits: ranked trait dataframe (includes `trait_category` when available)
#'     \item categories: top trait categories for the SNP set
#'   }
#' @export
summarise_baseline_traits <- function(trait_matrix,
                                      snp_ids = NULL,
                                      trait_info = NULL,
                                      coloc_groups = NULL,
                                      exclude_trait_ids = NULL,
                                      top_n = 20L,
                                      min_snps_with_signal = 5L,
                                      n_categories = 2L) {
  traits <- .summarise_traits_for_snps(
    trait_matrix = trait_matrix,
    snp_ids = snp_ids,
    group = "baseline",
    trait_info = trait_info,
    coloc_groups = coloc_groups,
    exclude_trait_ids = exclude_trait_ids,
    top_n = NULL,
    min_snps_with_signal = min_snps_with_signal
  )
  categories <- .top_trait_categories(traits, n = n_categories)
  traits <- utils::head(traits, top_n)
  return(list(
    traits = traits,
    categories = categories,
    top_categories = .category_label(categories)
  ))
}


#' @title Top Associated Traits Per SNP Group
#' @description For each SNP group of size at least `min_group_size`, rank trait
#' columns by mean absolute effect among member SNPs and report the leading
#' trait categories. Prefer \code{summarise_module_specific_traits()} when
#' common trait categories dominate magnitude rankings.
#' @inheritParams summarise_baseline_traits
#' @param groups Named vector (`names` = SNP ids) or dataframe with `snp_id`
#'   plus `group` / `cluster` / `program`.
#' @param min_group_size Only summarise groups with at least this many SNPs.
#' @param top_n Maximum traits per group. Defaults to 10.
#' @param min_snps_with_signal Drop traits observed in fewer than this many SNPs
#'   within the group. Defaults to 3.
#' @param n_categories Number of leading trait categories per group. Defaults to 2.
#' @return A list with `by_group` (per-group trait tables and categories) and
#'   `summary` (top trait and top categories per group).
#' @export
summarise_snp_group_traits <- function(trait_matrix,
                                       groups,
                                       trait_info = NULL,
                                       coloc_groups = NULL,
                                       exclude_trait_ids = NULL,
                                       min_group_size = 5L,
                                       top_n = 10L,
                                       min_snps_with_signal = 3L,
                                       n_categories = 2L) {
  group_df <- .normalize_snp_groups(groups)
  group_sizes <- table(group_df$group)
  large_groups <- names(group_sizes)[group_sizes >= min_group_size]

  empty_summary <- data.frame(
    group = character(0),
    n_snps = integer(0),
    n_traits_ranked = integer(0),
    top_trait = character(0),
    top_abs_mean_z = numeric(0),
    top_trait_categories = character(0),
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
    traits_all <- .summarise_traits_for_snps(
      trait_matrix = trait_matrix,
      snp_ids = snp_ids,
      group = grp,
      trait_info = trait_info,
      coloc_groups = coloc_groups,
      exclude_trait_ids = exclude_trait_ids,
      top_n = NULL,
      min_snps_with_signal = min_snps_with_signal
    )
    categories <- .top_trait_categories(traits_all, n = n_categories)
    list(
      group = grp,
      n_snps = length(unique(snp_ids)),
      traits = utils::head(traits_all, top_n),
      categories = categories,
      top_categories = .category_label(categories)
    )
  })

  summary_df <- dplyr::bind_rows(lapply(by_group, function(x) {
    top_trait <- if (nrow(x$traits) > 0) {
      x$traits$trait_name[1]
    } else {
      NA_character_
    }
    data.frame(
      group = x$group,
      n_snps = x$n_snps,
      n_traits_ranked = nrow(x$traits),
      top_trait = top_trait,
      top_abs_mean_z = if (nrow(x$traits) > 0) x$traits$abs_mean_z[1] else NA_real_,
      top_trait_categories = x$top_categories,
      stringsAsFactors = FALSE
    )
  }))

  return(list(
    by_group = by_group,
    summary = summary_df,
    min_group_size = as.integer(min_group_size)
  ))
}


#' @title Module-Specific Trait Drivers
#' @description For each SNP group of size at least `min_group_size`, rank traits by
#' how much stronger their mean absolute effect is inside the module than in the
#' remaining SNPs:
#' \deqn{\mathrm{specificity} = \bar{|z|}_{\mathrm{module}} /
#' (\bar{|z|}_{\mathrm{rest}} + \varepsilon)}.
#' Trait categories are ranked by the same module-vs-rest contrast so common
#' categories (e.g. Physiological Measures) do not dominate unless they are
#' specifically elevated in the module.
#'
#' Every large module always returns up to `top_n` traits (and category labels
#' derived from them). `min_specificity` is only a soft flag
#' (`passes_min_specificity`); it does not drop rows.
#' @inheritParams summarise_snp_group_traits
#' @param specificity_eps Pseudocount in the specificity denominator. Defaults
#'   to 0.1.
#' @param min_specificity Soft threshold used to set
#'   `passes_min_specificity` on ranked traits (default 1.25). Does not filter
#'   the returned table.
#' @return A list with `by_group` and `summary` (same shape as
#'   `summarise_snp_group_traits()`, plus specificity columns).
#' @export
summarise_module_specific_traits <- function(trait_matrix,
                                             groups,
                                             trait_info = NULL,
                                             coloc_groups = NULL,
                                             exclude_trait_ids = NULL,
                                             min_group_size = 5L,
                                             top_n = 10L,
                                             min_snps_with_signal = 3L,
                                             n_categories = 2L,
                                             specificity_eps = 0.1,
                                             min_specificity = 1.25) {
  if (!is.matrix(trait_matrix)) {
    stop("trait_matrix must be a matrix")
  }
  if (!is.numeric(specificity_eps) || specificity_eps < 0) {
    stop("specificity_eps must be a non-negative number")
  }
  if (!is.numeric(min_specificity) || min_specificity < 0) {
    stop("min_specificity must be a non-negative number")
  }

  group_df <- .normalize_snp_groups(groups)
  group_sizes <- table(group_df$group)
  large_groups <- names(group_sizes)[group_sizes >= min_group_size]
  all_snps <- rownames(trait_matrix)

  empty_summary <- data.frame(
    group = character(0),
    n_snps = integer(0),
    n_traits_ranked = integer(0),
    n_traits_passing_min_specificity = integer(0),
    top_trait = character(0),
    top_specificity = numeric(0),
    top_abs_mean_z = numeric(0),
    top_trait_categories = character(0),
    stringsAsFactors = FALSE
  )

  if (length(large_groups) == 0) {
    return(list(
      by_group = list(),
      summary = empty_summary,
      min_group_size = as.integer(min_group_size),
      min_specificity = min_specificity,
      specificity_eps = specificity_eps
    ))
  }

  by_group <- lapply(large_groups, function(grp) {
    snp_ids <- unique(group_df$snp_id[group_df$group == grp])
    rest_ids <- setdiff(all_snps, snp_ids)
    traits_mod <- .summarise_traits_for_snps(
      trait_matrix = trait_matrix,
      snp_ids = snp_ids,
      group = grp,
      trait_info = trait_info,
      coloc_groups = coloc_groups,
      exclude_trait_ids = exclude_trait_ids,
      top_n = NULL,
      min_snps_with_signal = min_snps_with_signal
    )
    traits_rest <- .summarise_traits_for_snps(
      trait_matrix = trait_matrix,
      snp_ids = rest_ids,
      group = "rest",
      trait_info = trait_info,
      coloc_groups = coloc_groups,
      exclude_trait_ids = exclude_trait_ids,
      top_n = NULL,
      min_snps_with_signal = 1L
    )

    rest_lookup <- stats::setNames(traits_rest$abs_mean_z, traits_rest$trait_id)
    traits_all <- traits_mod |>
      dplyr::mutate(
        abs_mean_z_rest = as.numeric(rest_lookup[as.character(trait_id)]),
        abs_mean_z_rest = dplyr::if_else(
          is.na(abs_mean_z_rest),
          0,
          abs_mean_z_rest
        ),
        specificity = abs_mean_z / (abs_mean_z_rest + specificity_eps),
        abs_mean_z_delta = abs_mean_z - abs_mean_z_rest,
        passes_min_specificity = is.finite(specificity) &
          specificity >= min_specificity
      ) |>
      dplyr::filter(is.finite(specificity)) |>
      dplyr::arrange(
        dplyr::desc(specificity),
        dplyr::desc(abs_mean_z_delta),
        dplyr::desc(abs_mean_z)
      )

    traits_top <- utils::head(traits_all, top_n)
    categories <- .top_specific_trait_categories(traits_top, n = n_categories)
    list(
      group = grp,
      n_snps = length(snp_ids),
      n_rest_snps = length(rest_ids),
      traits = traits_top,
      categories = categories,
      top_categories = .category_label(categories),
      n_traits_passing_min_specificity = sum(traits_all$passes_min_specificity)
    )
  })

  summary_df <- dplyr::bind_rows(lapply(by_group, function(x) {
    top_trait <- if (nrow(x$traits) > 0) {
      x$traits$trait_name[1]
    } else {
      NA_character_
    }
    data.frame(
      group = x$group,
      n_snps = x$n_snps,
      n_traits_ranked = nrow(x$traits),
      n_traits_passing_min_specificity = x$n_traits_passing_min_specificity,
      top_trait = top_trait,
      top_specificity = if (nrow(x$traits) > 0) x$traits$specificity[1] else NA_real_,
      top_abs_mean_z = if (nrow(x$traits) > 0) x$traits$abs_mean_z[1] else NA_real_,
      top_trait_categories = x$top_categories,
      stringsAsFactors = FALSE
    )
  }))

  return(list(
    by_group = by_group,
    summary = summary_df,
    min_group_size = as.integer(min_group_size),
    min_specificity = min_specificity,
    specificity_eps = specificity_eps
  ))
}


#' @title Compare Baseline Vs Group Trait Drivers
#' @description Mark whether each group's top traits overlap the baseline top
#' trait set, and whether baseline drivers split across multiple groups.
#' @param baseline_traits Output of `summarise_baseline_traits()`.
#' @param group_traits Output of `summarise_snp_group_traits()` or
#'   `summarise_module_specific_traits()`.
#' @param top_n Number of top traits per group/baseline to compare.
#' @return A list with `trait_status` and `group_overlap`.
#' @export
compare_group_traits <- function(baseline_traits, group_traits, top_n = 10L) {
  baseline_df <- .as_trait_table(baseline_traits)
  if (is.null(baseline_df) || is.null(group_traits$by_group)) {
    stop("baseline_traits and group_traits must be verification helper outputs")
  }

  baseline_ids <- utils::head(as.character(baseline_df$trait_id), top_n)
  group_rows <- dplyr::bind_rows(lapply(group_traits$by_group, function(x) {
    traits <- utils::head(x$traits, top_n)
    if (nrow(traits) == 0) {
      return(data.frame(
        group = character(0),
        trait_id = character(0),
        trait_name = character(0),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      group = x$group,
      trait_id = as.character(traits$trait_id),
      trait_name = as.character(traits$trait_name),
      stringsAsFactors = FALSE
    )
  }))

  all_ids <- union(baseline_ids, unique(group_rows$trait_id))
  if (length(all_ids) == 0) {
    return(list(
      trait_status = data.frame(
        trait_id = character(0),
        trait_name = character(0),
        in_baseline = logical(0),
        n_groups = integer(0),
        groups = character(0),
        status = character(0),
        stringsAsFactors = FALSE
      ),
      group_overlap = data.frame(
        group = character(0),
        n_top_traits = integer(0),
        n_shared_with_baseline = integer(0),
        n_novel = integer(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  name_lookup <- c(
    stats::setNames(as.character(baseline_df$trait_name), baseline_df$trait_id),
    stats::setNames(group_rows$trait_name, group_rows$trait_id)
  )

  trait_status <- dplyr::bind_rows(lapply(all_ids, function(id) {
    groups_hit <- sort(unique(group_rows$group[group_rows$trait_id == id]))
    in_baseline <- id %in% baseline_ids
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
    data.frame(
      trait_id = id,
      trait_name = unname(name_lookup[id][1]),
      in_baseline = in_baseline,
      n_groups = n_groups,
      groups = paste(groups_hit, collapse = ", "),
      status = status,
      stringsAsFactors = FALSE
    )
  })) |>
    dplyr::arrange(status, dplyr::desc(n_groups), trait_name)

  group_overlap <- dplyr::bind_rows(lapply(group_traits$by_group, function(x) {
    ids <- utils::head(as.character(x$traits$trait_id), top_n)
    data.frame(
      group = x$group,
      n_top_traits = length(ids),
      n_shared_with_baseline = sum(ids %in% baseline_ids),
      n_novel = sum(!ids %in% baseline_ids),
      stringsAsFactors = FALSE
    )
  }))

  return(list(
    trait_status = trait_status,
    group_overlap = group_overlap
  ))
}


#' @title Baseline Tissue Composition
#' @description Count tissues linked to a SNP set via molecular QTL coloc rows.
#' @param snp_ids SNP identifiers matching `snp_key`.
#' @param coloc_groups Coloc-group dataframe with a `tissue` column.
#' @param snp_key Column used to match SNPs. Defaults to `"variant_id"`.
#' @return A dataframe with `tissue`, `n_links`, `n_snps`, and `frac_links`.
#' @export
summarise_baseline_tissues <- function(snp_ids,
                                       coloc_groups,
                                       snp_key = c(
                                         "variant_id",
                                         "display_snp",
                                         "coloc_group_id"
                                       )) {
  snp_key <- match.arg(snp_key)
  links <- .snp_tissue_links(snp_ids, coloc_groups, snp_key)
  return(.summarise_tissue_counts(links))
}


#' @title Tissue Enrichment Per SNP Group
#' @description For each SNP group of size at least `min_group_size`, count
#' tissues at member SNPs (from `coloc_groups`) and compare observed vs expected
#' composition against a trait-level baseline by default. Baseline links come
#' only from SNPs attributed to the analysed trait (`baseline_snp_ids`), not the
#' full GPMap universe. For `background = "trait"`, each module is compared to
#' the rest of the trait SNPs (leave-one-out within the trait).
#' @param groups Named SNP group vector or group dataframe.
#' @param coloc_groups Coloc-group dataframe with `tissue`.
#' @param snp_key Column used to match SNPs.
#' @param min_group_size Only test groups with at least this many SNPs.
#' @param p_value_threshold FDR threshold used to populate `enriched`. The full
#'   observed-vs-expected table is always returned in `comparison`.
#' @param background `"trait"` (default; leave-one-out vs other trait SNPs),
#'   `"rest"` (leave-one-out vs other large groups), or `"pooled"` (all SNPs in
#'   large groups, including the tested group).
#' @param baseline_snp_ids SNP ids defining the trait baseline when
#'   `background = "trait"`. Defaults to all SNPs in `groups` when `NULL`.
#' @param bootstrap Compute bootstrap 95% confidence intervals for the fold
#'   enrichment (`fe_ci_lower`/`fe_ci_upper`). Defaults to `TRUE`.
#' @param permutations Number of link-label permutations for the empirical
#'   p-values (`p_perm`). Set to `0` to skip. Defaults to `1000L`.
#' @param seed RNG seed for the bootstrap FE confidence intervals and the
#'   permutation p-values. Defaults to `1L`.
#' @return A list with `by_group` (each with `tissues`, `comparison`,
#'   `enriched`, `enriched_any`), `summary`, and `background` (baseline tissue
#'   counts). `comparison` adds three alternative significance summaries
#'   alongside the hypergeometric `fdr`: `fe_ci_lower`/`fe_ci_upper` (bootstrap
#'   95% CI for `fold_enrichment`), `p_perm`/`fdr_perm` (permutation p and BH
#'   FDR), and `p_modules`/`fdr_modules` (one-sided Fisher against the pooled
#'   other modules). `sig_any` marks rows significant by any method; `enriched`
#'   is the FDR-only subset and `enriched_any` the any-method subset.
#' @export
enrich_snp_group_tissues <- function(groups,
                                     coloc_groups,
                                     snp_key = c(
                                       "variant_id",
                                       "display_snp",
                                       "coloc_group_id"
                                     ),
                                     min_group_size = 5L,
                                     p_value_threshold = 0.05,
                                     background = c("trait", "rest", "pooled"),
                                     baseline_snp_ids = NULL,
                                     bootstrap = TRUE,
                                     permutations = 1000L,
                                     seed = 1L) {
  snp_key <- match.arg(snp_key)
  background <- match.arg(background)
  group_df <- .normalize_snp_groups(groups)
  group_sizes <- table(group_df$group)
  large_groups <- names(group_sizes)[group_sizes >= min_group_size]
  all_module_snps <- unique(group_df$snp_id[group_df$group %in% large_groups])

  empty_summary <- data.frame(
    group = character(0),
    n_snps = integer(0),
    n_tissue_links = integer(0),
    n_enriched_tissues = integer(0),
    top_tissue = character(0),
    top_enriched_tissue = character(0),
    n_enriched_tissues_any = integer(0),
    top_enriched_tissue_any = character(0),
    n_enriched_tissues_perm = integer(0),
    n_enriched_tissues_modules = integer(0),
    stringsAsFactors = FALSE
  )

  if (length(large_groups) == 0) {
    return(list(
      by_group = list(),
      summary = empty_summary,
      background = .summarise_tissue_counts(
        .snp_tissue_links(character(0), coloc_groups, snp_key)
      ),
      n_background_links = 0L,
      background_mode = background,
      min_group_size = as.integer(min_group_size)
    ))
  }

  pooled_snps <- unique(group_df$snp_id[group_df$group %in% large_groups])
  if (identical(background, "trait")) {
    if (is.null(baseline_snp_ids)) {
      baseline_snp_ids <- unique(group_df$snp_id)
    }
    baseline_snp_ids <- unique(as.character(baseline_snp_ids))
    universe_links <- .snp_tissue_links(
      baseline_snp_ids,
      coloc_groups,
      snp_key
    )
  } else {
    universe_links <- .snp_tissue_links(pooled_snps, coloc_groups, snp_key)
  }
  universe_counts <- .summarise_tissue_counts(universe_links)

  by_group <- lapply(seq_along(large_groups), function(ii) {
    grp <- large_groups[ii]
    snp_ids <- unique(group_df$snp_id[group_df$group == grp])
    links <- .snp_tissue_links(snp_ids, coloc_groups, snp_key)
    counts <- .summarise_tissue_counts(links)

    if (identical(background, "trait")) {
      bg_snps <- setdiff(baseline_snp_ids, snp_ids)
      bg_links <- .snp_tissue_links(bg_snps, coloc_groups, snp_key)
      n_bg_snps <- length(bg_snps)
    } else if (identical(background, "rest")) {
      bg_snps <- setdiff(pooled_snps, snp_ids)
      bg_links <- .snp_tissue_links(bg_snps, coloc_groups, snp_key)
      n_bg_snps <- length(bg_snps)
    } else {
      bg_links <- universe_links
      n_bg_snps <- length(pooled_snps)
    }

    comparison <- .hypergeometric_link_enrichment(
      group_values = links$tissue,
      universe_values = bg_links$tissue,
      p_value_threshold = NULL,
      bootstrap = bootstrap,
      permutations = permutations,
      seed = seed + ii
    )
    names(comparison)[names(comparison) == "value"] <- "tissue"

    other_module_snps <- setdiff(all_module_snps, snp_ids)
    other_links <- .snp_tissue_links(other_module_snps, coloc_groups, snp_key)
    p_modules <- .fisher_vs_pool(
      group_values = links$tissue,
      other_values = other_links$tissue,
      values = comparison$tissue
    )
    comparison$p_modules <- p_modules[comparison$tissue]
    comparison$fdr_modules <- stats::p.adjust(comparison$p_modules, method = "BH")

    comparison$sig_any <-
      (!is.na(comparison$fdr) & comparison$fdr <= p_value_threshold) |
      (!is.na(comparison$fdr_perm) & comparison$fdr_perm <= p_value_threshold) |
      (!is.na(comparison$fdr_modules) & comparison$fdr_modules <= p_value_threshold)

    enriched <- comparison
    if (!is.null(p_value_threshold) && nrow(enriched) > 0) {
      enriched <- enriched[!is.na(enriched$fdr) & enriched$fdr <= p_value_threshold, , drop = FALSE]
    }
    enriched_any <- comparison
    if (nrow(enriched_any) > 0) {
      enriched_any <- enriched_any[
        enriched_any$sig_any & enriched_any$fold_enrichment >= 1,
        ,
        drop = FALSE
      ]
    }

    list(
      group = grp,
      n_snps = length(snp_ids),
      n_baseline_snps = n_bg_snps,
      n_tissue_links = nrow(links),
      tissues = counts,
      comparison = comparison,
      enriched = enriched,
      enriched_any = enriched_any
    )
  })

  summary_df <- dplyr::bind_rows(lapply(by_group, function(x) {
    top_tissue <- if (nrow(x$tissues) > 0) x$tissues$tissue[1] else NA_character_
    top_enriched <- if (nrow(x$enriched) > 0) {
      paste0(
        x$enriched$tissue[1],
        " (FDR=", signif(x$enriched$fdr[1], 3),
        ", FE=", signif(x$enriched$fold_enrichment[1], 3), ")"
      )
    } else {
      NA_character_
    }
    any_x <- if (nrow(x$enriched_any) > 0) {
      best <- x$enriched_any[which.min(x$enriched_any$fdr), , drop = FALSE]
      tag <- .method_tag(
        best$fdr, best$fdr_perm, best$fdr_modules, p_value_threshold
      )
      paste0(
        best$tissue[1],
        " (FE=", signif(best$fold_enrichment[1], 2),
        if (nzchar(tag)) paste0(", ", tag) else "",
        ")"
      )
    } else {
      NA_character_
    }
    data.frame(
      group = x$group,
      n_snps = x$n_snps,
      n_tissue_links = x$n_tissue_links,
      n_enriched_tissues = nrow(x$enriched),
      top_tissue = top_tissue,
      top_enriched_tissue = top_enriched,
      n_enriched_tissues_any = nrow(x$enriched_any),
      top_enriched_tissue_any = any_x,
      n_enriched_tissues_perm = sum(
        !is.na(x$comparison$fdr_perm) & x$comparison$fdr_perm <= p_value_threshold,
        na.rm = TRUE
      ),
      n_enriched_tissues_modules = sum(
        !is.na(x$comparison$fdr_modules) &
          x$comparison$fdr_modules <= p_value_threshold,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))

  return(list(
    by_group = by_group,
    summary = summary_df,
    background = universe_counts,
    n_background_links = nrow(universe_links),
    background_mode = background,
    min_group_size = as.integer(min_group_size)
  ))
}


#' @title Trait-Category Enrichment Per SNP Group
#' @description For each SNP group, compare trait-category composition in
#' `coloc_groups` to the rest of a trait-level SNP baseline (expected vs
#' observed). Uses distinct SNP-category links from coloc results attributed to
#' the analysed trait SNPs.
#' @inheritParams enrich_snp_group_tissues
#' @param baseline_snp_ids SNP ids defining the trait baseline. Defaults to all
#'   SNPs in `groups` when `NULL`.
#' @return A list with `by_group` (each with `categories`, `comparison`,
#'   `enriched`, `enriched_any`), `summary`, and `background` category counts.
#'   `comparison` and `summary` add the alternative significance summaries
#'   described in `enrich_snp_group_tissues()`.
#' @export
enrich_snp_group_trait_categories <- function(groups,
                                              coloc_groups,
                                              snp_key = c(
                                                "variant_id",
                                                "display_snp",
                                                "coloc_group_id"
                                              ),
                                              min_group_size = 5L,
                                              p_value_threshold = 0.05,
                                              baseline_snp_ids = NULL,
                                              bootstrap = TRUE,
                                              permutations = 1000L,
                                              seed = 1L) {
  snp_key <- match.arg(snp_key)
  group_df <- .normalize_snp_groups(groups)
  group_sizes <- table(group_df$group)
  large_groups <- names(group_sizes)[group_sizes >= min_group_size]
  all_module_snps <- unique(group_df$snp_id[group_df$group %in% large_groups])

  empty_summary <- data.frame(
    group = character(0),
    n_snps = integer(0),
    n_category_links = integer(0),
    n_enriched_categories = integer(0),
    top_category = character(0),
    top_enriched_category = character(0),
    n_enriched_categories_any = integer(0),
    top_enriched_category_any = character(0),
    n_enriched_categories_perm = integer(0),
    n_enriched_categories_modules = integer(0),
    stringsAsFactors = FALSE
  )

  if (is.null(baseline_snp_ids)) {
    baseline_snp_ids <- unique(group_df$snp_id)
  }
  baseline_snp_ids <- unique(as.character(baseline_snp_ids))

  if (length(large_groups) == 0) {
    return(list(
      by_group = list(),
      summary = empty_summary,
      background = .summarise_category_counts(
        .snp_category_links(character(0), coloc_groups, snp_key)
      ),
      n_background_links = 0L,
      background_mode = "trait",
      min_group_size = as.integer(min_group_size)
    ))
  }

  universe_links <- .snp_category_links(
    baseline_snp_ids,
    coloc_groups,
    snp_key
  )
  universe_counts <- .summarise_category_counts(universe_links)

  by_group <- lapply(seq_along(large_groups), function(ii) {
    grp <- large_groups[ii]
    snp_ids <- unique(group_df$snp_id[group_df$group == grp])
    links <- .snp_category_links(snp_ids, coloc_groups, snp_key)
    counts <- .summarise_category_counts(links)
    bg_snps <- setdiff(baseline_snp_ids, snp_ids)
    bg_links <- .snp_category_links(bg_snps, coloc_groups, snp_key)
    comparison <- .hypergeometric_link_enrichment(
      group_values = links$trait_category,
      universe_values = bg_links$trait_category,
      p_value_threshold = NULL,
      bootstrap = bootstrap,
      permutations = permutations,
      seed = seed + ii
    )
    names(comparison)[names(comparison) == "value"] <- "trait_category"

    other_module_snps <- setdiff(all_module_snps, snp_ids)
    other_links <- .snp_category_links(other_module_snps, coloc_groups, snp_key)
    p_modules <- .fisher_vs_pool(
      group_values = links$trait_category,
      other_values = other_links$trait_category,
      values = comparison$trait_category
    )
    comparison$p_modules <- p_modules[comparison$trait_category]
    comparison$fdr_modules <- stats::p.adjust(comparison$p_modules, method = "BH")

    comparison$sig_any <-
      (!is.na(comparison$fdr) & comparison$fdr <= p_value_threshold) |
      (!is.na(comparison$fdr_perm) & comparison$fdr_perm <= p_value_threshold) |
      (!is.na(comparison$fdr_modules) & comparison$fdr_modules <= p_value_threshold)

    enriched <- comparison
    if (!is.null(p_value_threshold) && nrow(enriched) > 0) {
      enriched <- enriched[!is.na(enriched$fdr) & enriched$fdr <= p_value_threshold, , drop = FALSE]
    }
    enriched_any <- comparison
    if (nrow(enriched_any) > 0) {
      enriched_any <- enriched_any[
        enriched_any$sig_any & enriched_any$fold_enrichment >= 1,
        ,
        drop = FALSE
      ]
    }

    list(
      group = grp,
      n_snps = length(snp_ids),
      n_baseline_snps = length(baseline_snp_ids),
      n_category_links = nrow(links),
      categories = counts,
      comparison = comparison,
      enriched = enriched,
      enriched_any = enriched_any
    )
  })

  summary_df <- dplyr::bind_rows(lapply(by_group, function(x) {
    top_category <- if (nrow(x$categories) > 0) {
      x$categories$trait_category[1]
    } else {
      NA_character_
    }
    top_enriched <- if (nrow(x$enriched) > 0) {
      paste0(
        x$enriched$trait_category[1],
        " (FDR=", signif(x$enriched$fdr[1], 3),
        ", FE=", signif(x$enriched$fold_enrichment[1], 3), ")"
      )
    } else {
      NA_character_
    }
    any_x <- if (nrow(x$enriched_any) > 0) {
      best <- x$enriched_any[which.min(x$enriched_any$fdr), , drop = FALSE]
      tag <- .method_tag(
        best$fdr, best$fdr_perm, best$fdr_modules, p_value_threshold
      )
      paste0(
        best$trait_category[1],
        " (FE=", signif(best$fold_enrichment[1], 2),
        if (nzchar(tag)) paste0(", ", tag) else "",
        ")"
      )
    } else {
      NA_character_
    }
    data.frame(
      group = x$group,
      n_snps = x$n_snps,
      n_category_links = x$n_category_links,
      n_enriched_categories = nrow(x$enriched),
      top_category = top_category,
      top_enriched_category = top_enriched,
      n_enriched_categories_any = nrow(x$enriched_any),
      top_enriched_category_any = any_x,
      n_enriched_categories_perm = sum(
        !is.na(x$comparison$fdr_perm) & x$comparison$fdr_perm <= p_value_threshold,
        na.rm = TRUE
      ),
      n_enriched_categories_modules = sum(
        !is.na(x$comparison$fdr_modules) &
          x$comparison$fdr_modules <= p_value_threshold,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))

  return(list(
    by_group = by_group,
    summary = summary_df,
    background = universe_counts,
    n_background_links = nrow(universe_links),
    background_mode = "trait",
    min_group_size = as.integer(min_group_size)
  ))
}


.summarise_traits_for_snps <- function(trait_matrix,
                                       snp_ids,
                                       group,
                                       trait_info,
                                       coloc_groups = NULL,
                                       exclude_trait_ids,
                                       top_n,
                                       min_snps_with_signal = 1L) {
  if (!is.matrix(trait_matrix)) {
    stop("trait_matrix must be a matrix")
  }
  if (is.null(snp_ids)) {
    snp_ids <- rownames(trait_matrix)
  }
  snp_ids <- as.character(snp_ids)
  row_idx <- match(snp_ids, rownames(trait_matrix))
  row_idx <- row_idx[!is.na(row_idx)]
  if (length(row_idx) == 0 || ncol(trait_matrix) == 0) {
    return(data.frame(
      group = character(0),
      trait_id = character(0),
      trait_name = character(0),
      trait_category = character(0),
      mean_z = numeric(0),
      abs_mean_z = numeric(0),
      n_snps_with_signal = integer(0),
      n_snps = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  mat <- trait_matrix[row_idx, , drop = FALSE]
  mean_z <- apply(mat, 2, function(z) mean(z, na.rm = TRUE))
  abs_mean_z <- apply(mat, 2, function(z) mean(abs(z), na.rm = TRUE))
  n_signal <- colSums(is.finite(mat) & mat != 0)

  out <- data.frame(
    group = group,
    trait_id = colnames(trait_matrix),
    trait_name = colnames(trait_matrix),
    trait_category = NA_character_,
    mean_z = as.numeric(mean_z),
    abs_mean_z = as.numeric(abs_mean_z),
    n_snps_with_signal = as.integer(n_signal),
    n_snps = length(row_idx),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$abs_mean_z), , drop = FALSE]
  out <- out[out$n_snps_with_signal >= min_snps_with_signal, , drop = FALSE]

  if (!is.null(exclude_trait_ids)) {
    out <- out[!out$trait_id %in% as.character(exclude_trait_ids), , drop = FALSE]
  }

  out <- .join_trait_labels(out, trait_info)
  out <- .join_trait_categories(out, coloc_groups)
  out <- out[order(-out$abs_mean_z, -out$n_snps_with_signal), , drop = FALSE]
  if (!is.null(top_n)) {
    out <- utils::head(out, top_n)
  }
  rownames(out) <- NULL
  return(out)
}


.as_trait_table <- function(baseline_traits) {
  if (is.data.frame(baseline_traits)) {
    return(baseline_traits)
  }
  if (is.list(baseline_traits) && is.data.frame(baseline_traits$traits)) {
    return(baseline_traits$traits)
  }
  return(NULL)
}


.top_trait_categories <- function(traits, n = 2L) {
  empty <- data.frame(
    trait_category = character(0),
    score = numeric(0),
    n_traits = integer(0),
    stringsAsFactors = FALSE
  )
  if (!is.data.frame(traits) || nrow(traits) == 0) {
    return(empty)
  }
  if (!"trait_category" %in% names(traits)) {
    return(empty)
  }

  scored <- traits |>
    dplyr::filter(
      !is.na(trait_category),
      trait_category != "",
      is.finite(abs_mean_z)
    )
  if (nrow(scored) == 0) {
    return(empty)
  }

  out <- scored |>
    dplyr::group_by(trait_category) |>
    dplyr::summarise(
      score = sum(abs_mean_z),
      n_traits = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(score), dplyr::desc(n_traits), trait_category)

  return(utils::head(out, n))
}


.top_specific_trait_categories <- function(traits, n = 2L) {
  empty <- data.frame(
    trait_category = character(0),
    score = numeric(0),
    mean_specificity = numeric(0),
    n_traits = integer(0),
    stringsAsFactors = FALSE
  )
  if (!is.data.frame(traits) || nrow(traits) == 0) {
    return(empty)
  }
  needed <- c("trait_category", "specificity", "abs_mean_z_delta")
  if (!all(needed %in% names(traits))) {
    return(empty)
  }

  scored <- traits |>
    dplyr::filter(
      !is.na(trait_category),
      trait_category != "",
      is.finite(specificity),
      is.finite(abs_mean_z_delta)
    )
  if (nrow(scored) == 0) {
    return(empty)
  }

  out <- scored |>
    dplyr::group_by(trait_category) |>
    dplyr::summarise(
      score = sum(pmax(abs_mean_z_delta, 0) * specificity),
      mean_specificity = mean(specificity),
      n_traits = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(score),
      dplyr::desc(mean_specificity),
      dplyr::desc(n_traits),
      trait_category
    )

  return(utils::head(out, n))
}


.category_label <- function(categories) {
  if (!is.data.frame(categories) || nrow(categories) == 0) {
    return(NA_character_)
  }
  return(paste(categories$trait_category, collapse = " / "))
}


.join_trait_categories <- function(out, coloc_groups) {
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    return(out)
  }
  if (!all(c("trait_id", "trait_category") %in% names(coloc_groups))) {
    return(out)
  }

  cats <- coloc_groups |>
    dplyr::mutate(trait_id = as.character(trait_id)) |>
    dplyr::filter(!is.na(trait_category), trait_category != "") |>
    dplyr::group_by(trait_id) |>
    dplyr::summarise(
      trait_category = dplyr::first(trait_category),
      .groups = "drop"
    )

  out <- out |>
    dplyr::select(-trait_category) |>
    dplyr::left_join(cats, by = "trait_id")

  return(out)
}


.join_trait_labels <- function(out, trait_info) {
  if (is.null(trait_info) || nrow(trait_info) == 0) {
    return(out)
  }

  info <- trait_info
  if ("feature_id" %in% names(info) && !"trait_id" %in% names(info)) {
    info <- dplyr::rename(info, trait_id = feature_id)
  }
  if ("feature_name" %in% names(info) && !"trait_name" %in% names(info)) {
    info <- dplyr::rename(info, trait_name = feature_name)
  }
  if (!all(c("trait_id", "trait_name") %in% names(info))) {
    return(out)
  }

  info <- info |>
    dplyr::mutate(trait_id = as.character(trait_id)) |>
    dplyr::distinct(trait_id, .keep_all = TRUE) |>
    dplyr::select(trait_id, trait_name)

  out <- out |>
    dplyr::select(-trait_name) |>
    dplyr::left_join(info, by = "trait_id") |>
    dplyr::mutate(
      trait_name = dplyr::if_else(is.na(trait_name), trait_id, trait_name)
    )

  return(out)
}


.snp_tissue_links <- function(snp_ids, coloc_groups, snp_key) {
  if (is.null(snp_ids) || length(snp_ids) == 0) {
    return(data.frame(
      snp_id = character(0),
      tissue = character(0),
      stringsAsFactors = FALSE
    ))
  }
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("coloc_groups is required")
  }
  if (!"tissue" %in% names(coloc_groups)) {
    stop("coloc_groups must include a tissue column")
  }
  if (!snp_key %in% names(coloc_groups)) {
    stop("coloc_groups must include column: ", snp_key)
  }

  snp_ids <- unique(as.character(snp_ids))
  snp_col <- as.character(coloc_groups[[snp_key]])
  out <- coloc_groups |>
    dplyr::mutate(snp_id = snp_col) |>
    dplyr::filter(
      snp_id %in% snp_ids,
      !is.na(tissue),
      tissue != "",
      !is.na(gene_id)
    ) |>
    dplyr::distinct(snp_id, tissue)

  return(out)
}


.snp_category_links <- function(snp_ids, coloc_groups, snp_key) {
  if (is.null(snp_ids) || length(snp_ids) == 0) {
    return(data.frame(
      snp_id = character(0),
      trait_category = character(0),
      stringsAsFactors = FALSE
    ))
  }
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("coloc_groups is required")
  }
  if (!"trait_category" %in% names(coloc_groups)) {
    stop("coloc_groups must include a trait_category column")
  }
  if (!snp_key %in% names(coloc_groups)) {
    stop("coloc_groups must include column: ", snp_key)
  }

  snp_ids <- unique(as.character(snp_ids))
  snp_col <- as.character(coloc_groups[[snp_key]])
  out <- coloc_groups |>
    dplyr::mutate(snp_id = snp_col) |>
    dplyr::filter(
      snp_id %in% snp_ids,
      !is.na(trait_category),
      trait_category != ""
    ) |>
    dplyr::distinct(snp_id, trait_category)

  return(out)
}


.summarise_tissue_counts <- function(links) {
  if (nrow(links) == 0) {
    return(data.frame(
      tissue = character(0),
      n_links = integer(0),
      n_snps = integer(0),
      frac_links = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  total <- nrow(links)
  out <- links |>
    dplyr::group_by(tissue) |>
    dplyr::summarise(
      n_links = dplyr::n(),
      n_snps = dplyr::n_distinct(snp_id),
      .groups = "drop"
    ) |>
    dplyr::mutate(frac_links = n_links / total) |>
    dplyr::arrange(dplyr::desc(n_links), tissue)

  return(out)
}


.summarise_category_counts <- function(links) {
  if (nrow(links) == 0) {
    return(data.frame(
      trait_category = character(0),
      n_links = integer(0),
      n_snps = integer(0),
      frac_links = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  total <- nrow(links)
  out <- links |>
    dplyr::group_by(trait_category) |>
    dplyr::summarise(
      n_links = dplyr::n(),
      n_snps = dplyr::n_distinct(snp_id),
      .groups = "drop"
    ) |>
    dplyr::mutate(frac_links = n_links / total) |>
    dplyr::arrange(dplyr::desc(n_links), trait_category)

  return(out)
}


.hypergeometric_link_enrichment <- function(group_values,
                                            universe_values,
                                            p_value_threshold = NULL,
                                            bootstrap = TRUE,
                                            permutations = 1000L,
                                            seed = 1L) {
  empty <- data.frame(
    value = character(0),
    n_group = integer(0),
    n_universe = integer(0),
    n_group_total = integer(0),
    n_universe_total = integer(0),
    frac_group = numeric(0),
    frac_universe = numeric(0),
    fold_enrichment = numeric(0),
    p_value = numeric(0),
    fdr = numeric(0),
    fe_ci_lower = numeric(0),
    fe_ci_upper = numeric(0),
    p_perm = numeric(0),
    fdr_perm = numeric(0),
    stringsAsFactors = FALSE
  )
  group_values <- as.character(group_values)
  universe_values <- as.character(universe_values)
  group_values <- group_values[!is.na(group_values) & group_values != ""]
  universe_values <- universe_values[
    !is.na(universe_values) & universe_values != ""
  ]

  if (length(group_values) == 0 || length(universe_values) == 0) {
    return(empty)
  }

  n_group_total <- length(group_values)
  n_universe_total <- length(universe_values)
  group_counts <- table(group_values)
  universe_counts <- table(universe_values)
  values <- names(group_counts)

  rows <- lapply(values, function(val) {
    k <- as.integer(group_counts[[val]])
    k_u <- if (val %in% names(universe_counts)) {
      as.integer(universe_counts[[val]])
    } else {
      0L
    }
    frac_group <- k / n_group_total
    frac_u <- k_u / n_universe_total
    fold_enrichment <- if (frac_u > 0) frac_group / frac_u else Inf
    p_value <- if (n_group_total > n_universe_total) {
      NA_real_
    } else {
      stats::phyper(
        q = k - 1L,
        m = k_u,
        n = n_universe_total - k_u,
        k = n_group_total,
        lower.tail = FALSE
      )
    }
    data.frame(
      value = val,
      n_group = k,
      n_universe = k_u,
      n_group_total = n_group_total,
      n_universe_total = n_universe_total,
      frac_group = frac_group,
      frac_universe = frac_u,
      fold_enrichment = fold_enrichment,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    return(empty)
  }
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")

  if (bootstrap || permutations > 0) {
    set.seed(seed)
  }
  if (bootstrap) {
    ci <- .fe_bootstrap_ci(
      n_group = out$n_group,
      n_universe = out$n_universe,
      n_group_total = out$n_group_total,
      n_universe_total = out$n_universe_total,
      B = 1000L
    )
    out$fe_ci_lower <- ci$lower
    out$fe_ci_upper <- ci$upper
  } else {
    out$fe_ci_lower <- NA_real_
    out$fe_ci_upper <- NA_real_
  }
  if (permutations > 0) {
    p_perm <- .permutation_link_pvalues(
      group_values = group_values,
      universe_values = universe_values,
      values = values,
      B = permutations
    )
    out$p_perm <- p_perm[out$value]
    out$fdr_perm <- stats::p.adjust(out$p_perm, method = "BH")
  } else {
    out$p_perm <- NA_real_
    out$fdr_perm <- NA_real_
  }

  if (!is.null(p_value_threshold)) {
    out <- out[out$fdr <= p_value_threshold, , drop = FALSE]
  }
  out <- out[order(
    out$fdr,
    -out$fold_enrichment,
    out$p_value,
    out$value
  ), , drop = FALSE]
  rownames(out) <- NULL
  return(out)
}


# Percentile-bootstrap CI for fold enrichment. The group and universe link
# counts are drawn as binomials (parametric bootstrap) from the observed
# proportions, so the CI reflects the small-count uncertainty of the ratio.
.fe_bootstrap_ci <- function(n_group,
                             n_universe,
                             n_group_total,
                             n_universe_total,
                             B = 1000L) {
  n <- length(n_group)
  lower <- upper <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    k <- n_group[i]
    k_u <- n_universe[i]
    ng <- n_group_total[i]
    nu <- n_universe_total[i]
    if (is.na(k) || is.na(k_u) || ng <= 0 || nu <= 0) next
    if (k_u <= 0) next
    if (k == 0) {
      lower[i] <- 0
      upper[i] <- 0
      next
    }
    kg <- stats::rbinom(B, size = ng, prob = k / ng)
    ku_draw <- stats::rbinom(B, size = nu, prob = k_u / nu)
    fe <- (kg / ng) / (ku_draw / nu)
    fe <- fe[is.finite(fe)]
    if (length(fe) == 0) next
    lower[i] <- unname(stats::quantile(fe, 0.025, na.rm = TRUE))
    upper[i] <- unname(stats::quantile(fe, 0.975, na.rm = TRUE))
  }
  return(list(lower = lower, upper = upper))
}


# Permutation p-values for module-vs-rest enrichment. Link labels are shuffled
# (draw `n_group_total` links from the pooled baseline links, without
# replacement), so the null distribution matches the hypergeometric test but is
# empirical rather than asymptotic.
.permutation_link_pvalues <- function(group_values,
                                      universe_values,
                                      values,
                                      B = 1000L) {
  p <- stats::setNames(rep(NA_real_, length(values)), values)
  n_group_total <- length(group_values)
  pool <- c(group_values, universe_values)
  n_pool <- length(pool)
  group_counts <- table(group_values)
  if (n_group_total == 0 || n_pool == 0 || n_group_total >= n_pool) {
    return(p)
  }
  for (v in values) {
    k_obs <- if (v %in% names(group_counts)) {
      as.integer(group_counts[[v]])
    } else {
      0L
    }
    if (k_obs == 0) {
      p[[v]] <- 1
      next
    }
    cnt <- 0L
    for (b in seq_len(B)) {
      idx <- sample.int(n_pool, n_group_total, replace = FALSE)
      if (sum(pool[idx] == v) >= k_obs) {
        cnt <- cnt + 1L
      }
    }
    p[[v]] <- (cnt + 1) / (B + 1)
  }
  return(p)
}


# One-sided Fisher test of the group against a pooled comparison set (used for
# module-vs-other-modules): is the group enriched for `value` relative to the
# pooled links of the other modules?
.fisher_vs_pool <- function(group_values, other_values, values) {
  n_g <- length(group_values)
  n_o <- length(other_values)
  p <- stats::setNames(rep(NA_real_, length(values)), values)
  if (n_g == 0 || n_o == 0) {
    return(p)
  }
  gc <- table(group_values)
  oc <- table(other_values)
  for (v in values) {
    k <- if (v %in% names(gc)) as.integer(gc[[v]]) else 0L
    k_o <- if (v %in% names(oc)) as.integer(oc[[v]]) else 0L
    if (k == 0) {
      p[[v]] <- 1
      next
    }
    tab <- matrix(c(k, n_g - k, k_o, n_o - k_o), nrow = 2, byrow = TRUE)
    p[[v]] <- stats::fisher.test(tab, alternative = "greater")$p.value
  }
  return(p)
}


# Compact labels of which significance methods pass the threshold for a row.
.method_tag <- function(fdr, fdr_perm, fdr_modules, threshold = 0.05) {
  tags <- character(0)
  if (!is.na(fdr) && fdr <= threshold) {
    tags <- c(tags, paste0("FDR=", signif(fdr, 3)))
  }
  if (!is.na(fdr_perm) && fdr_perm <= threshold) {
    tags <- c(tags, paste0("perm=", signif(fdr_perm, 3)))
  }
  if (!is.na(fdr_modules) && fdr_modules <= threshold) {
    tags <- c(tags, paste0("vsmod=", signif(fdr_modules, 3)))
  }
  return(paste(tags, collapse = ", "))
}


# Backward-compatible wrapper used by older call sites / tests if needed.
.hypergeometric_tissue_enrichment <- function(group_links,
                                              background_links,
                                              p_value_threshold) {
  out <- .hypergeometric_link_enrichment(
    group_values = group_links$tissue,
    universe_values = background_links$tissue,
    p_value_threshold = p_value_threshold
  )
  names(out)[names(out) == "value"] <- "tissue"
  names(out)[names(out) == "n_universe"] <- "n_background"
  names(out)[names(out) == "n_universe_total"] <- "n_background_total"
  return(out)
}


#' @title Pairwise Module-vs-Module Enrichment
#' @description For each trait-category (or tissue) value, compare every pair of
#' SNP modules with a two-sided Fisher exact test on link counts, then apply a
#' global Benjamini-Hochberg correction across all module-pair x value tests.
#' This complements the module-vs-rest hypergeometric enrichment: instead of
#' asking whether a module differs from the trait background, it asks which
#' modules differ from *each other*.
#' @param groups Named vector (`names` = SNP ids) or dataframe of module
#'   assignments, as used by `enrich_snp_group_trait_categories()`.
#' @param coloc_groups A coloc_groups dataframe (from
#'   `trait(..., include_associations = TRUE)`).
#' @param value_col Which annotation to compare: `"trait_category"` or
#'   `"tissue"`.
#' @param snp_key SNP id column in `coloc_groups`: `"variant_id"`,
#'   `"display_snp"`, or `"coloc_group_id"`.
#' @param min_group_size Only compare modules with at least this many SNPs.
#' @param p_value_threshold FDR threshold for flagging a pairwise difference.
#'   Defaults to 0.05.
#' @return A list with:
#'   \itemize{
#'     \item summary: one row per module with the number of significant pairwise
#'       differences, the top differing value, its FDR, and the top differing
#'       partner module
#'     \item pairwise: long dataframe of all module-pair x value tests with
#'       link counts, rates, `p_value`, and the global `fdr`
#'     \item by_value: per-value counts of significant module pairs
#'   }
#' @export
compare_snp_group_modules <- function(groups,
                                      coloc_groups,
                                      value_col = c("trait_category", "tissue"),
                                      snp_key = c(
                                        "variant_id",
                                        "display_snp",
                                        "coloc_group_id"
                                      ),
                                      min_group_size = 5L,
                                      p_value_threshold = 0.05) {
  value_col <- match.arg(value_col)
  snp_key <- match.arg(snp_key)
  group_df <- .normalize_snp_groups(groups)
  group_sizes <- table(group_df$group)
  large_groups <- names(group_sizes)[group_sizes >= min_group_size]

  empty <- list(
    summary = data.frame(
      module = character(0),
      n_pairwise_sig = integer(0),
      top_pairwise_value = character(0),
      top_pairwise_fdr = numeric(0),
      top_pairwise_partner = character(0),
      stringsAsFactors = FALSE
    ),
    pairwise = data.frame(
      value = character(0),
      module_a = character(0),
      module_b = character(0),
      n_a = integer(0),
      k_a = integer(0),
      n_b = integer(0),
      k_b = integer(0),
      rate_a = numeric(0),
      rate_b = numeric(0),
      p_value = numeric(0),
      fdr = numeric(0),
      stringsAsFactors = FALSE
    ),
    by_value = data.frame(
      value = character(0),
      n_sig_pairs = integer(0),
      stringsAsFactors = FALSE
    )
  )
  if (length(large_groups) < 2) {
    return(empty)
  }

  links_fun <- if (identical(value_col, "tissue")) {
    .snp_tissue_links
  } else {
    .snp_category_links
  }
  link_col <- if (identical(value_col, "tissue")) "tissue" else "trait_category"

  links_by_module <- lapply(large_groups, function(grp) {
    snp_ids <- unique(group_df$snp_id[group_df$group == grp])
    return(links_fun(snp_ids, coloc_groups, snp_key)[[link_col]])
  })
  names(links_by_module) <- large_groups

  all_values <- sort(unique(unlist(links_by_module, use.names = FALSE)))
  if (length(all_values) == 0) {
    return(empty)
  }

  pair_idx <- combn(seq_along(large_groups), 2L, simplify = FALSE)
  pair_rows <- lapply(pair_idx, function(ij) {
    i <- ij[1]
    j <- ij[2]
    vi <- links_by_module[[i]]
    vj <- links_by_module[[j]]
    ni <- length(vi)
    nj <- length(vj)
    if (ni == 0 || nj == 0) {
      return(NULL)
    }
    rows <- lapply(all_values, function(v) {
      ki <- sum(vi == v)
      kj <- sum(vj == v)
      tab <- matrix(c(ki, ni - ki, kj, nj - kj), nrow = 2, byrow = TRUE)
      p_value <- tryCatch(
        stats::fisher.test(tab, alternative = "two.sided")$p.value,
        error = function(e) NA_real_
      )
      return(data.frame(
        value = v,
        module_a = large_groups[i],
        module_b = large_groups[j],
        n_a = ni,
        k_a = ki,
        n_b = nj,
        k_b = kj,
        rate_a = ki / ni,
        rate_b = kj / nj,
        p_value = p_value,
        stringsAsFactors = FALSE
      ))
    })
    return(dplyr::bind_rows(rows))
  })
  pairwise <- dplyr::bind_rows(pair_rows)
  if (nrow(pairwise) == 0) {
    return(empty)
  }
  pairwise$fdr <- stats::p.adjust(pairwise$p_value, method = "BH")
  sig <- pairwise[
    !is.na(pairwise$fdr) & pairwise$fdr <= p_value_threshold,
    ,
    drop = FALSE
  ]

  by_value <- if (nrow(sig) > 0) {
    sig |>
      dplyr::count(value, name = "n_sig_pairs") |>
      dplyr::arrange(dplyr::desc(n_sig_pairs), value)
  } else {
    empty$by_value
  }

  summary_df <- if (nrow(sig) > 0) {
    lapply(large_groups, function(m) {
      rows <- sig[sig$module_a == m | sig$module_b == m, , drop = FALSE]
      if (nrow(rows) == 0) {
        return(NULL)
      }
      best <- rows[which.min(rows$fdr), , drop = FALSE]
      partner <- if (best$module_a == m) best$module_b else best$module_a
      return(data.frame(
        module = m,
        n_pairwise_sig = nrow(rows),
        top_pairwise_value = best$value,
        top_pairwise_fdr = best$fdr,
        top_pairwise_partner = partner,
        stringsAsFactors = FALSE
      ))
    }) |>
      dplyr::bind_rows() |>
      dplyr::arrange(dplyr::desc(n_pairwise_sig), top_pairwise_fdr)
  } else {
    empty$summary
  }

  return(list(
    summary = summary_df,
    pairwise = pairwise,
    by_value = by_value,
    value_col = value_col,
    min_group_size = as.integer(min_group_size),
    p_value_threshold = p_value_threshold
  ))
}
