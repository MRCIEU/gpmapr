# Display tables for the create_report() templates (inst/reports/): rounding,
# ordering and column choice only, no analysis.

.report_settings_table <- function(settings) {
  groups <- .gpmap_settings_groups()
  return(dplyr::bind_rows(lapply(names(groups), function(g) {
    keys <- groups[[g]]
    return(data.frame(
      group = g, setting = keys,
      value = vapply(keys, function(k) paste(format(settings[[k]]), collapse = ", "), character(1)),
      stringsAsFactors = FALSE
    ))
  })))
}


.report_excluded_table <- function(excluded, include_background = FALSE) {
  ex <- excluded
  if (is.null(ex) || nrow(ex) == 0) {
    return(data.frame())
  }
  if (!include_background) {
    ex <- ex[ex$source != "background_duplicate", , drop = FALSE]
  }
  source_order <- c("anchor_trait", "duplicate_of_anchor", "proxy_group", "user", "background_duplicate")
  ex <- ex[order(ex$anchor_id, match(ex$source, source_order), ex$trait_name), , drop = FALSE]
  out <- data.frame(
    anchor_id = ex$anchor_id,
    source = ex$source,
    trait_id = ex$trait_id,
    trait_name = ex$trait_name,
    rule = ifelse(is.na(ex$rule), "", ex$rule),
    n_loci = ex$n_loci,
    coverage = round(ex$coverage, 2),
    pct_same_sign = round(ex$pct_same_sign),
    r_z = round(ex$r_z, 2),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  return(out)
}


.report_background_counts <- function(excluded) {
  ex <- excluded[excluded$source == "background_duplicate", , drop = FALSE]
  if (nrow(ex) == 0) {
    return(0L)
  }
  return(nrow(ex))
}


.report_fit_table <- function(d) {
  X <- d$fit$x_matrix
  p <- d$programs
  tiers <- table(factor(as.character(p$confidence_tier), levels = c("high", "medium", "single_trait", "low")))
  return(data.frame(
    n_snps = ncol(X),
    n_background_traits = nrow(X),
    fraction_observed = round(mean(!is.na(X)), 3),
    n_ebmf_factors = d$fit$cluster_details$n_programs,
    n_programs = nrow(p),
    n_dropped_small = if (is.null(d$summary$dropped_programs)) 0L else nrow(d$summary$dropped_programs),
    n_high = as.integer(tiers[["high"]]),
    n_medium = as.integer(tiers[["medium"]]),
    n_single_trait = as.integer(tiers[["single_trait"]]),
    n_low = as.integer(tiers[["low"]]),
    background_mean_S = round(d$summary$background$mean_similarity, 3),
    stringsAsFactors = FALSE
  ))
}


.report_program_table <- function(d) {
  p <- programs(d, tiers = NULL)
  if (nrow(p) == 0) {
    return(p)
  }
  col <- function(name) {
    return(if (name %in% names(p)) p[[name]] else rep(NA, nrow(p)))
  }
  return(data.frame(
    program = col("program_label"),
    confidence_tier = as.character(p$confidence_tier),
    n_snps_filtered = col("n_snps_filtered"),
    effective_n_traits = round(col("effective_n_traits"), 2),
    mean_internal_similarity = round(col("mean_internal_similarity"), 3),
    connectedness = round(col("connectedness"), 3),
    similarity_q = signif(col("similarity_q"), 2),
    similarity_z = round(col("similarity_z"), 1),
    factor_pve_pct = round(100 * col("factor_pve"), 2),
    replication = round(col("replication"), 2),
    parent_program = col("parent_program"),
    stringsAsFactors = FALSE
  ))
}


.report_cards_table <- function(d) {
  cards <- program_cards(d)
  if (nrow(cards) == 0) {
    return(cards)
  }
  text_cols <- c("specific_traits", "categories", "tissues", "reactome", "kegg", "genes")
  for (col in text_cols) {
    cards[[col]] <- gsub("; ", "<br>", .html_escape(cards[[col]]), fixed = TRUE)
  }
  return(data.frame(
    program = cards$program_label,
    confidence_tier = cards$confidence_tier,
    n_snps = cards$n_snps_filtered,
    effective_n_traits = round(cards$effective_n_traits, 2),
    cards[, text_cols],
    stringsAsFactors = FALSE
  ))
}


.report_links_table <- function(x, tiers = c("high", "medium", "low")) {
  l <- links(x, tiers = tiers)
  if (nrow(l) == 0) {
    return(l)
  }
  num <- function(name, digits, fn = round) {
    return(if (name %in% names(l)) fn(l[[name]], digits) else rep(NA, nrow(l)))
  }
  return(data.frame(
    program_a = l$label_a,
    program_b = l$label_b,
    tier = l$tier,
    profile_strength = l$profile_strength,
    phi_traits = num("phi_traits", 3),
    n_eff_traits = num("n_eff_traits", 1),
    top_traits = .html_escape(ifelse(is.na(l$top_traits), "", l$top_traits)),
    locus_status = l$locus_status,
    n_loci_shared = l$n_loci_shared,
    n_eff_loci = num("n_eff_loci", 1),
    alignment = num("alignment", 3),
    comparability = num("comparability", 3),
    q_concordance = num("q_concordance", 3, signif),
    direction_profile = l$direction_profile,
    rg = num("rg", 2),
    rg_ci = ifelse(
      is.finite(l$rg_ci_lower),
      sprintf("[%.2f, %.2f]", l$rg_ci_lower, l$rg_ci_upper), ""
    ),
    direction_check = ifelse(is.na(l$direction_check), "", l$direction_check),
    stringsAsFactors = FALSE
  ))
}


.report_tier_table <- function(x) {
  tp <- x$trait_pairs
  return(data.frame(
    trait_pair = paste(tp$trait_name_a, "~", tp$trait_name_b),
    n_loci_shared = tp$n_loci_shared,
    high = tp$n_high, medium = tp$n_medium, low = tp$n_low, none = tp$n_none,
    stringsAsFactors = FALSE
  ))
}


.report_direction_table <- function(x) {
  l <- links(x)
  if (nrow(l) == 0) {
    return(data.frame())
  }
  l$trait_pair <- paste(l$trait_name_a, "~", l$trait_name_b)
  l$rg_check <- ifelse(is.na(l$direction_check), "no call", l$direction_check)
  out <- l |>
    dplyr::count(trait_pair, direction_profile, rg_check) |>
    tidyr::pivot_wider(names_from = rg_check, values_from = n, values_fill = 0)
  return(as.data.frame(out, stringsAsFactors = FALSE))
}


.report_family_table <- function(x) {
  f <- families(x)
  if (nrow(f) == 0) {
    return(f)
  }
  members <- vapply(f$family, function(fam) {
    m <- x$nodes[x$nodes$family == fam, , drop = FALSE]
    return(paste(.program_display_label(m$program_id, m$trait_name), collapse = "<br>"))
  }, character(1))
  f$mean_phi <- round(f$mean_phi, 2)
  f$members <- members
  f$traits <- NULL
  front <- c("family", "members")
  return(f[, c(front, setdiff(names(f), front)), drop = FALSE])
}


.html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  return(gsub(">", "&gt;", x, fixed = TRUE))
}
