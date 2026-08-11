#' @title Baseline Top Associated Traits
#' @description Rank complex-trait columns by mean absolute effect across a set
#' of SNPs. Use as a verification baseline when pathway enrichment is empty or
#' weak: do SNP groups recover distinct trait drivers?
#' @param trait_matrix SNPs x traits numeric matrix (e.g. from
#'   `build_perturbation_matrices()$trait_matrix`).
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
#' @description For each SNP group larger than `min_group_size`, rank trait
#' columns by mean absolute effect among member SNPs and report the leading
#' trait categories. Prefer \code{summarise_module_specific_traits()} when
#' common trait categories dominate magnitude rankings.
#' @inheritParams summarise_baseline_traits
#' @param groups Named vector (`names` = SNP ids) or dataframe with `snp_id`
#'   plus `group` / `cluster` / `program`.
#' @param min_group_size Only summarise groups with more than this many SNPs.
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
  large_groups <- names(group_sizes)[group_sizes > min_group_size]

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
#' @description For each SNP group larger than `min_group_size`, rank traits by
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
  large_groups <- names(group_sizes)[group_sizes > min_group_size]
  all_snps <- rownames(trait_matrix)

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
#' @description For each SNP group larger than `min_group_size`, count tissues at
#' member SNPs and test enrichment (hypergeometric on unique SNP-tissue links).
#' By default the background is leave-one-out (`background = "rest"`): SNPs in
#' other large groups, so reported tissues are module-specific rather than
#' merely abundant in the pooled set.
#' @param groups Named SNP group vector or group dataframe.
#' @param coloc_groups Coloc-group dataframe with `tissue`.
#' @param snp_key Column used to match SNPs.
#' @param min_group_size Only test groups with more than this many SNPs.
#' @param p_value_threshold FDR threshold for reporting enriched tissues.
#' @param background `"rest"` (default; leave-one-out) or `"pooled"` (all SNPs
#'   in large groups, including the tested group).
#' @return A list with `by_group`, `summary`, and `background` (pooled counts
#'   across large groups, for reference).
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
                                     background = c("rest", "pooled")) {
  snp_key <- match.arg(snp_key)
  background <- match.arg(background)
  group_df <- .normalize_snp_groups(groups)
  group_sizes <- table(group_df$group)
  large_groups <- names(group_sizes)[group_sizes > min_group_size]

  pooled_snps <- unique(group_df$snp_id[group_df$group %in% large_groups])
  pooled_links <- .snp_tissue_links(pooled_snps, coloc_groups, snp_key)
  pooled_counts <- .summarise_tissue_counts(pooled_links)

  by_group <- lapply(large_groups, function(grp) {
    snp_ids <- unique(group_df$snp_id[group_df$group == grp])
    if (identical(background, "rest")) {
      bg_snps <- setdiff(pooled_snps, snp_ids)
    } else {
      bg_snps <- pooled_snps
    }
    links <- .snp_tissue_links(snp_ids, coloc_groups, snp_key)
    bg_links <- .snp_tissue_links(bg_snps, coloc_groups, snp_key)
    counts <- .summarise_tissue_counts(links)
    enriched <- .hypergeometric_tissue_enrichment(
      group_links = links,
      background_links = bg_links,
      p_value_threshold = p_value_threshold
    )
    list(
      group = grp,
      n_snps = length(snp_ids),
      n_rest_snps = length(bg_snps),
      n_tissue_links = nrow(links),
      tissues = counts,
      enriched = enriched
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
    data.frame(
      group = x$group,
      n_snps = x$n_snps,
      n_tissue_links = x$n_tissue_links,
      n_enriched_tissues = nrow(x$enriched),
      top_tissue = top_tissue,
      top_enriched_tissue = top_enriched,
      stringsAsFactors = FALSE
    )
  }))

  return(list(
    by_group = by_group,
    summary = summary_df,
    background = pooled_counts,
    n_background_links = nrow(pooled_links),
    background_mode = background,
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


.hypergeometric_tissue_enrichment <- function(group_links,
                                              background_links,
                                              p_value_threshold) {
  if (nrow(group_links) == 0 || nrow(background_links) == 0) {
    return(data.frame(
      tissue = character(0),
      n_group = integer(0),
      n_background = integer(0),
      n_group_total = integer(0),
      n_background_total = integer(0),
      fold_enrichment = numeric(0),
      p_value = numeric(0),
      fdr = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  n_group_total <- nrow(group_links)
  n_background_total <- nrow(background_links)
  background_counts <- table(background_links$tissue)
  group_counts <- table(group_links$tissue)
  # Include tissues present only in the module (absent from rest/background).
  tissues <- names(group_counts)

  rows <- lapply(tissues, function(tissue) {
    k <- as.integer(group_counts[[tissue]])
    k_bg <- if (tissue %in% names(background_counts)) {
      as.integer(background_counts[[tissue]])
    } else {
      0L
    }
    frac_group <- k / n_group_total
    frac_bg <- k_bg / n_background_total
    fold_enrichment <- if (frac_bg > 0) frac_group / frac_bg else Inf
    p_value <- stats::phyper(
      q = k - 1L,
      m = k_bg,
      n = n_background_total - k_bg,
      k = n_group_total,
      lower.tail = FALSE
    )
    data.frame(
      tissue = tissue,
      n_group = k,
      n_background = k_bg,
      n_group_total = n_group_total,
      n_background_total = n_background_total,
      fold_enrichment = fold_enrichment,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    return(out)
  }
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[out$fdr <= p_value_threshold, , drop = FALSE]
  out <- out[order(
    out$fdr,
    -out$fold_enrichment,
    out$p_value,
    out$tissue
  ), , drop = FALSE]
  rownames(out) <- NULL
  return(out)
}
