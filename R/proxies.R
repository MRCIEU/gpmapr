#' @title Proxy Traits of an Anchor Trait
#' @description List the GPMap traits that measure the same thing as an anchor
#' (target) trait, from GPMap's own duplicate and proxy groups
#' ([trait_duplicates()]):
#' \describe{
#'   \item{`duplicate_of_anchor`}{other traits in the anchor's duplicate group:
#'   other GWAS of the same measurement (`rule` is the API's `match_rules`).}
#'   \item{`proxy_group`}{other traits in the anchor's proxy group: a wider
#'   grouping, for example T2D GWAS together with diabetes diagnosis,
#'   medication, family history and complications (`rule` is the API's
#'   `proxy_rules`).}
#' }
#' After a trait's SNPs are oriented to its own risk allele, a proxy of the
#' anchor is positive at almost every SNP, so it adds the same component to
#' every SNP's profile and hides the programs that differ between SNPs.
#' [decompose_trait()] removes them before fitting.
#'
#' With `coloc_groups` (the anchor's `trait()$coloc_groups`), each proxy is
#' also described by how closely it tracks the anchor at the anchor's own
#' loci: `n_loci` (loci where it is observed), `coverage` (share of the
#' anchor's loci), `pct_same_sign` (share with an effect in the same direction
#' after orienting to the anchor's risk allele) and `r_z` (uncentred
#' correlation of its oriented z-scores with the anchor's \eqn{|z|}; 1 means
#' the same effects up to scale). A true proxy is near 100% and near 1.
#' @param trait_id The anchor trait id.
#' @param duplicates The output of [trait_duplicates()]. Fetched from the API
#'   if not supplied. With the older shape that has no proxy-group columns,
#'   only the duplicate group is used.
#' @param coloc_groups Optional coloc-group table holding the anchor's rows,
#'   for the evidence columns.
#' @return A data.frame with one row per proxy trait: `anchor_id`, `trait_id`,
#'   `trait_name`, `source`, `rule`, and the evidence columns when
#'   `coloc_groups` is supplied.
#' @examplesIf gpmap_available()
#' anchor_proxies("1631")
#' @concept trait decomposition
#' @export
anchor_proxies <- function(trait_id, duplicates = trait_duplicates(), coloc_groups = NULL) {
  if (length(trait_id) != 1 || is.na(trait_id)) {
    stop("trait_id must be a single trait id", call. = FALSE)
  }
  trait_id <- as.character(trait_id)
  out <- .proxy_table(trait_id, duplicates)
  if (!is.null(coloc_groups)) {
    out <- .add_proxy_evidence(out, coloc_groups, trait_id)
  }
  return(out)
}


# The anchor's duplicate group and proxy group from trait_duplicates(), in
# either of its shapes. With include_self, the anchor itself is listed too
# (source "anchor_trait"), for an uploaded GWAS that measures it.
.proxy_table <- function(trait_id, duplicates, include_self = FALSE) {
  trait_id <- as.character(trait_id)
  out <- .empty_proxy_table()
  if (is.null(duplicates) || !is.data.frame(duplicates) || nrow(duplicates) == 0) {
    return(out)
  }
  duplicates$trait_id <- as.character(duplicates$trait_id)
  name_of <- stats::setNames(as.character(duplicates$trait_name), duplicates$trait_id)
  rows <- list()

  if ("duplicate_group" %in% names(duplicates)) {
    group <- duplicates$duplicate_group[duplicates$trait_id == trait_id]
    group <- group[!is.na(group)]
    if (length(group) > 0) {
      members <- duplicates[
        duplicates$duplicate_group %in% group & duplicates$trait_id != trait_id, ,
        drop = FALSE
      ]
      rule <- if ("match_rules" %in% names(members)) members$match_rules else NA_character_
      rows$dup <- data.frame(
        trait_id = members$trait_id, trait_name = as.character(members$trait_name),
        source = "duplicate_of_anchor", rule = as.character(rule),
        stringsAsFactors = FALSE
      )
    }
  } else {
    pairs <- .duplicate_parent_pairs(duplicates)
    parent <- c(pairs$parent_trait_id[pairs$trait_id == trait_id], trait_id)
    ids <- unique(c(
      pairs$trait_id[pairs$parent_trait_id %in% parent], parent
    ))
    ids <- setdiff(ids, trait_id)
    if (length(ids) > 0) {
      parent_names <- if ("parent_trait_name" %in% names(duplicates)) {
        stats::setNames(
          as.character(duplicates$parent_trait_name), as.character(duplicates$parent_trait_id)
        )
      } else {
        character(0)
      }
      names_all <- c(name_of, parent_names)
      rows$dup <- data.frame(
        trait_id = ids, trait_name = unname(names_all[ids]),
        source = "duplicate_of_anchor", rule = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  if ("proxy_group_id" %in% names(duplicates)) {
    group <- duplicates$proxy_group_id[duplicates$trait_id == trait_id]
    group <- group[!is.na(group)]
    if (length(group) > 0) {
      members <- duplicates[
        duplicates$proxy_group_id %in% group & duplicates$trait_id != trait_id, ,
        drop = FALSE
      ]
      members <- members[!members$trait_id %in% rows$dup$trait_id, , drop = FALSE]
      rule <- if ("proxy_rules" %in% names(members)) members$proxy_rules else NA_character_
      rows$proxy <- data.frame(
        trait_id = members$trait_id, trait_name = as.character(members$trait_name),
        source = "proxy_group", rule = as.character(rule),
        stringsAsFactors = FALSE
      )
    }
  }

  if (include_self) {
    self_name <- unname(name_of[trait_id])
    rows <- c(list(self = data.frame(
      trait_id = trait_id,
      trait_name = if (length(self_name) == 0) NA_character_ else self_name,
      source = "anchor_trait", rule = NA_character_,
      stringsAsFactors = FALSE
    )), rows)
  }
  if (length(rows) == 0) {
    return(out)
  }
  out <- dplyr::bind_rows(rows)
  out <- out[!duplicated(out$trait_id), , drop = FALSE]
  out <- data.frame(anchor_id = trait_id, out, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  return(out)
}


.empty_proxy_table <- function() {
  return(data.frame(
    anchor_id = character(0), trait_id = character(0), trait_name = character(0),
    source = character(0), rule = character(0), stringsAsFactors = FALSE
  ))
}


# How closely each listed trait tracks the anchor at the anchor's own loci.
.add_proxy_evidence <- function(table, coloc_groups, anchor_id) {
  ev <- .proxy_evidence(coloc_groups, anchor_id, table$trait_id)
  if (nrow(table) == 0) {
    return(cbind(table, ev[0, c("n_loci", "coverage", "pct_same_sign", "r_z")]))
  }
  out <- dplyr::left_join(table, ev[, c("trait_id", "n_loci", "coverage", "pct_same_sign", "r_z")],
    by = "trait_id"
  )
  return(out)
}


# n_loci, coverage, pct_same_sign and r_z of each of `ids` against the anchor,
# at the anchor's own loci (its strongest row per coloc group).
.proxy_evidence <- function(coloc_groups, anchor_id, ids) {
  ids <- unique(as.character(ids))
  empty <- data.frame(
    trait_id = ids, n_loci = rep(0L, length(ids)), coverage = rep(0, length(ids)),
    pct_same_sign = rep(NA_real_, length(ids)), r_z = rep(NA_real_, length(ids)),
    stringsAsFactors = FALSE
  )
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0 || length(ids) == 0 ||
        !all(c("trait_id", "coloc_group_id", "beta", "se", "min_p") %in% names(coloc_groups))) {
    return(empty)
  }
  cg <- coloc_groups[!is.na(coloc_groups$beta) & !is.na(coloc_groups$se) & coloc_groups$se > 0, ,
    drop = FALSE
  ]
  cg$trait_id <- as.character(cg$trait_id)
  anchor <- cg[cg$trait_id == as.character(anchor_id), , drop = FALSE] |>
    dplyr::group_by(coloc_group_id) |>
    dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(coloc_group_id, z_anchor = beta / se)
  hits <- cg[cg$trait_id %in% ids, , drop = FALSE] |>
    dplyr::group_by(trait_id, coloc_group_id) |>
    dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::inner_join(anchor, by = "coloc_group_id") |>
    dplyr::mutate(z = beta / se * sign(z_anchor), za = abs(z_anchor)) |>
    dplyr::group_by(trait_id) |>
    dplyr::summarise(
      n_loci = dplyr::n(),
      pct_same_sign = 100 * mean(z > 0),
      r_z = sum(z * za) / sqrt(sum(z^2) * sum(za^2)),
      .groups = "drop"
    )
  out <- data.frame(trait_id = ids, stringsAsFactors = FALSE) |>
    dplyr::left_join(hits, by = "trait_id") |>
    dplyr::mutate(
      n_loci = dplyr::coalesce(as.integer(n_loci), 0L),
      coverage = n_loci / max(nrow(anchor), 1)
    ) |>
    dplyr::select(trait_id, n_loci, coverage, pct_same_sign, r_z)
  return(as.data.frame(out, stringsAsFactors = FALSE))
}
