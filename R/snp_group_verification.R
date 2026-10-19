#' @title Module-Specific Trait Drivers
#' @description For each SNP group of size at least `min_group_size`, rank traits by
#' how much stronger their mean absolute effect is inside the module than in the
#' remaining SNPs:
#' \deqn{\mathrm{specificity} = \bar{|z|}_{\mathrm{module}} /
#' (\bar{|z|}_{\mathrm{rest}} + \varepsilon)}.
#'
#' Every large module always returns up to `top_n` traits. `min_specificity` is
#' only a soft flag (`passes_min_specificity`); it does not drop rows.
#' @param trait_matrix SNPs x traits numeric matrix (e.g.
#'   `run_univariate_clustering()$trait_matrix`).
#' @param groups Named vector (`names` = SNP ids) or dataframe with `snp_id`
#'   plus `group` / `cluster` / `program`.
#' @param trait_info Optional dataframe with `feature_id`/`trait_id` and
#'   `feature_name`/`trait_name` columns for labels.
#' @param coloc_groups Optional coloc-group dataframe used to attach
#'   `trait_category` labels.
#' @param exclude_trait_ids Optional trait ids to drop (e.g. the target trait).
#' @param min_group_size Only summarise groups with at least this many SNPs.
#' @param top_n Maximum traits per group. Defaults to 10.
#' @param min_snps_with_signal Drop traits observed in fewer than this many SNPs
#'   within the group. Defaults to 3.
#' @param specificity_eps Pseudocount in the specificity denominator. Defaults
#'   to 0.1.
#' @param min_rest_snps_with_signal Minimum number of SNPs outside the module
#'   at which a trait must be observed before it gets a specificity score. A
#'   trait seen only inside the module has an essentially empty denominator, so
#'   the ratio measures sparsity rather than specificity; such traits get `NA`.
#'   Defaults to `2`.
#' @param min_specificity Soft threshold used to set
#'   `passes_min_specificity` on ranked traits (default 1.25). Does not filter
#'   the returned table.
#' @return A list with `by_group`, one entry per large group with `group`,
#'   `n_snps`, `n_rest_snps`, `traits` (the top `top_n` traits with
#'   `abs_mean_z`, `abs_mean_z_rest`, `specificity`, `abs_mean_z_delta` and
#'   `passes_min_specificity`) and `n_traits_passing_min_specificity`, plus the
#'   settings used.
#' @export
summarise_module_specific_traits <- function(trait_matrix,
                                             groups,
                                             trait_info = NULL,
                                             coloc_groups = NULL,
                                             exclude_trait_ids = NULL,
                                             min_group_size = 5L,
                                             top_n = 10L,
                                             min_snps_with_signal = 3L,
                                             specificity_eps = 0.1,
                                             min_specificity = 1.25,
                                             min_rest_snps_with_signal = 2L) {
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

  if (length(large_groups) == 0) {
    return(list(
      by_group = list(),
      min_group_size = as.integer(min_group_size),
      min_specificity = min_specificity,
      specificity_eps = specificity_eps,
      min_rest_snps_with_signal = as.integer(min_rest_snps_with_signal)
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
    rest_n_lookup <- stats::setNames(
      traits_rest$n_snps_with_signal, traits_rest$trait_id
    )
    traits_all <- traits_mod |>
      dplyr::mutate(
        abs_mean_z_rest = as.numeric(rest_lookup[as.character(trait_id)]),
        abs_mean_z_rest = dplyr::if_else(
          is.na(abs_mean_z_rest),
          0,
          abs_mean_z_rest
        ),
        n_snps_with_signal_rest = as.integer(
          rest_n_lookup[as.character(trait_id)]
        ),
        n_snps_with_signal_rest = dplyr::if_else(
          is.na(n_snps_with_signal_rest),
          0L,
          n_snps_with_signal_rest
        ),
        # A trait observed at few or no SNPs outside the module has an
        # essentially empty denominator, so abs_mean_z / specificity_eps
        # returns a huge ratio that reflects sparsity, not specificity. Such
        # traits get NA rather than a score that would top every ranking.
        specificity = dplyr::if_else(
          n_snps_with_signal_rest >= min_rest_snps_with_signal,
          abs_mean_z / (abs_mean_z_rest + specificity_eps),
          NA_real_
        ),
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

    list(
      group = grp,
      n_snps = length(snp_ids),
      n_rest_snps = length(rest_ids),
      traits = utils::head(traits_all, top_n),
      n_traits_passing_min_specificity = sum(traits_all$passes_min_specificity)
    )
  })

  return(list(
    by_group = by_group,
    min_group_size = as.integer(min_group_size),
    min_specificity = min_specificity,
    specificity_eps = specificity_eps,
    min_rest_snps_with_signal = as.integer(min_rest_snps_with_signal)
  ))
}


#' @title Annotation Enrichment of a SNP Set
#' @description Test whether a set of SNPs is enriched for a trait category or
#' a tissue compared with the rest of a trait's SNPs. Each SNP is linked to the
#' distinct categories (or tissues) of the traits it colocalises with in
#' `coloc_groups`; tissue links come from rows with a `gene_id` (molecular QTL
#' studies). For every annotation value seen in the set, a one-sided Fisher
#' exact test compares the set's links (with vs without the value) against the
#' background's links, where the background is `background_snp_ids` minus the
#' set. P-values are BH-adjusted across the values tested.
#'
#' Links rather than SNPs are counted, so a SNP colocalising with many traits
#' contributes several links; links from one SNP, and SNPs in LD, are treated
#' as independent.
#' @param snp_ids Variant ids of the SNP set.
#' @param background_snp_ids Variant ids of the trait's SNPs (the set may be
#'   included; it is removed before testing).
#' @param coloc_groups Coloc-group dataframe with `variant_id` and the
#'   annotation column (`trait_category`, or `tissue` and `gene_id`).
#' @param annotation `"trait_category"` (default) or `"tissue"`.
#' @param min_set_size Sets with fewer distinct SNPs return an empty table.
#'   Defaults to `5`.
#' @return A dataframe with one row per annotation value seen in the set,
#'   ordered by `fdr`: the value (in a column named after `annotation`),
#'   `n_set` and `n_background` (links with the value), `n_set_total` and
#'   `n_background_total` (all links), `frac_set`, `frac_background`,
#'   `fold_enrichment` (`frac_set / frac_background`), `p_value` and `fdr`.
#' @export
enrich_snp_set <- function(snp_ids,
                           background_snp_ids,
                           coloc_groups,
                           annotation = c("trait_category", "tissue"),
                           min_set_size = 5L) {
  annotation <- match.arg(annotation)
  snp_ids <- unique(as.character(snp_ids))
  empty <- .link_enrichment_test(character(0), character(0))
  names(empty)[names(empty) == "value"] <- annotation
  if (length(snp_ids) < min_set_size) {
    return(empty)
  }
  background <- setdiff(unique(as.character(background_snp_ids)), snp_ids)
  set_links <- .snp_annotation_links(snp_ids, coloc_groups, annotation)
  background_links <- .snp_annotation_links(background, coloc_groups, annotation)
  out <- .link_enrichment_test(set_links[[annotation]], background_links[[annotation]])
  names(out)[names(out) == "value"] <- annotation
  return(out)
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


# Distinct SNP-annotation links: one row per SNP and trait category (or
# tissue) among the coloc rows at that SNP. Tissue links come only from rows
# with a gene_id (molecular QTL studies).
.snp_annotation_links <- function(snp_ids,
                                  coloc_groups,
                                  annotation = c("trait_category", "tissue"),
                                  snp_key = "variant_id") {
  annotation <- match.arg(annotation)
  empty <- data.frame(snp_id = character(0), value = character(0), stringsAsFactors = FALSE)
  names(empty)[2] <- annotation
  if (is.null(snp_ids) || length(snp_ids) == 0) {
    return(empty)
  }
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("coloc_groups is required")
  }
  needed <- c(snp_key, annotation, if (annotation == "tissue") "gene_id")
  missing_cols <- setdiff(needed, names(coloc_groups))
  if (length(missing_cols) > 0) {
    stop("coloc_groups must include column(s): ", paste(missing_cols, collapse = ", "))
  }

  snp_ids <- unique(as.character(snp_ids))
  values <- as.character(coloc_groups[[annotation]])
  keep <- as.character(coloc_groups[[snp_key]]) %in% snp_ids &
    !is.na(values) & values != ""
  if (annotation == "tissue") {
    keep <- keep & !is.na(coloc_groups$gene_id)
  }
  out <- data.frame(
    snp_id = as.character(coloc_groups[[snp_key]][keep]),
    value = values[keep],
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out), , drop = FALSE]
  names(out)[2] <- annotation
  rownames(out) <- NULL
  return(out)
}


# One-sided Fisher exact test of each annotation value seen in a set's links
# against a disjoint background's links: set links with / without the value vs
# background links with / without it. This is the hypergeometric tail with the
# set and background pooled as the universe. BH across the values tested.
.link_enrichment_test <- function(set_values, background_values) {
  empty <- data.frame(
    value = character(0),
    n_set = integer(0),
    n_background = integer(0),
    n_set_total = integer(0),
    n_background_total = integer(0),
    frac_set = numeric(0),
    frac_background = numeric(0),
    fold_enrichment = numeric(0),
    p_value = numeric(0),
    fdr = numeric(0),
    stringsAsFactors = FALSE
  )
  set_values <- as.character(set_values)
  background_values <- as.character(background_values)
  set_values <- set_values[!is.na(set_values) & set_values != ""]
  background_values <- background_values[!is.na(background_values) & background_values != ""]
  if (length(set_values) == 0 || length(background_values) == 0) {
    return(empty)
  }

  n_set_total <- length(set_values)
  n_background_total <- length(background_values)
  set_counts <- table(set_values)
  background_counts <- table(background_values)
  values <- names(set_counts)
  k <- as.integer(set_counts[values])
  k_b <- as.integer(background_counts[values])
  k_b[is.na(k_b)] <- 0L

  frac_set <- k / n_set_total
  frac_background <- k_b / n_background_total
  out <- data.frame(
    value = values,
    n_set = k,
    n_background = k_b,
    n_set_total = n_set_total,
    n_background_total = n_background_total,
    frac_set = frac_set,
    frac_background = frac_background,
    fold_enrichment = ifelse(frac_background > 0, frac_set / frac_background, Inf),
    p_value = stats::phyper(
      q = k - 1L,
      m = k + k_b,
      n = n_set_total + n_background_total - k - k_b,
      k = n_set_total,
      lower.tail = FALSE
    ),
    stringsAsFactors = FALSE
  )
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[order(out$fdr, -out$fold_enrichment, out$p_value, out$value), , drop = FALSE]
  rownames(out) <- NULL
  return(out)
}
