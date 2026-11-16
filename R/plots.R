#' @title Plot a Decomposition or Correspondence
#' @description Draw one of the standard figures of a [decompose_trait()] or
#' [trait_correspondence()] result. Each `type` is also its own function,
#' which takes more options.
#' @param x A [decompose_trait()] or [trait_correspondence()] result.
#' @param type Which figure: for a decomposition `"programs"`
#'   ([plot_programs()]), `"similarity"` ([plot_similarity_heatmap()]),
#'   `"enrichment"` ([plot_program_enrichment()]) or `"profile"`
#'   ([plot_program_profile()]); for a correspondence `"families"`
#'   ([plot_families()]), `"phi"` ([plot_phi()]), `"link"` ([plot_link()]) or
#'   `"family_loci"` ([plot_family_loci()]).
#' @param ... Passed to the plotting function (e.g. `program` for
#'   `"profile"`, `program_a` and `program_b` for `"link"`).
#' @param y Unused.
#' @return A ggplot object.
#' @name plot.gpmap
NULL


#' @rdname plot.gpmap
#' @export
plot.gpmap_decomposition <- function(x, y, type = c("programs", "similarity", "enrichment", "profile"),
                                     ...) {
  type <- match.arg(type)
  fn <- switch(type,
    programs = plot_programs,
    similarity = plot_similarity_heatmap,
    enrichment = plot_program_enrichment,
    profile = plot_program_profile
  )
  return(fn(x, ...))
}


#' @rdname plot.gpmap
#' @export
plot.gpmap_correspondence <- function(x, y, type = c("families", "phi", "link", "family_loci"), ...) {
  type <- match.arg(type)
  fn <- switch(type,
    families = plot_families,
    phi = plot_phi,
    link = plot_link,
    family_loci = plot_family_loci
  )
  return(fn(x, ...))
}


.tier_colours <- function() {
  return(c(high = "#2e7d32", medium = "#ef8f00", single_trait = "#8a8a8a", low = "#c62828", none = "#cfcfcf"))
}


.direction_colours <- function() {
  return(c(concordant = "#2a78d6", antagonistic = "#eb6834"))
}


.assert_class <- function(x, cls, fn) {
  if (!inherits(x, cls)) {
    needed <- switch(cls,
      gpmap_decomposition = "decompose_trait()",
      gpmap_correspondence = "trait_correspondence()"
    )
    stop(fn, "() needs a ", needed, " result", call. = FALSE)
  }
  return(invisible(TRUE))
}


.empty_plot <- function(message, title = NULL) {
  return(
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = message, colour = "grey40") +
      ggplot2::theme_void() +
      ggplot2::labs(title = title)
  )
}


#' @title Program Sizes by Confidence Tier
#' @description One bar per program of a [decompose_trait()] result: how many
#' SNPs pass its membership gate (`n_snps_filtered`), coloured by confidence
#' tier, in confidence order.
#' @param x A [decompose_trait()] result.
#' @param tiers Tiers to show; `NULL` (default) shows every program.
#' @return A ggplot object.
#' @concept trait decomposition
#' @export
plot_programs <- function(x, tiers = NULL) {
  .assert_class(x, "gpmap_decomposition", "plot_programs")
  p <- programs(x, tiers = tiers)
  if (nrow(p) == 0) {
    return(.empty_plot("No programs", x$trait$name))
  }
  p$label <- paste0("P", p$program)
  p$label <- factor(p$label, levels = rev(unique(p$label)))
  p$confidence_tier <- factor(as.character(p$confidence_tier), levels = names(.tier_colours()))
  return(
    ggplot2::ggplot(p, ggplot2::aes(x = .data$n_snps_filtered, y = .data$label, fill = .data$confidence_tier)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::scale_fill_manual(values = .tier_colours(), drop = TRUE, name = "confidence tier") +
      ggplot2::labs(
        title = paste0(x$trait$name, ": programs"),
        x = "member SNPs (lFSR / magnitude gate)", y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  )
}


#' @title Background-trait Profile of a Program
#' @description The background traits with the largest loadings on one
#' program (its trait loadings, `L_pm`, sign-oriented to the anchor), with
#' transparency showing confidence (`1 - lfsr`). Only high- and
#' medium-confidence programs have extracted profiles.
#' @param x A [decompose_trait()] result, or a [trait_correspondence()]
#'   result (then `program` is a program id such as `"1631:3"`).
#' @param program The program number (or program id).
#' @param top_n Number of traits shown.
#' @return A ggplot object.
#' @concept trait decomposition
#' @export
plot_program_profile <- function(x, program, top_n = 20L) {
  if (inherits(x, "gpmap_correspondence")) {
    tid <- sub(":[^:]*$", "", as.character(program))
    x <- x$decompositions[[tid]]
    if (is.null(x)) {
      stop("No trait ", tid, " in this correspondence", call. = FALSE)
    }
  }
  .assert_class(x, "gpmap_decomposition", "plot_program_profile")
  k <- sub("^.*:", "", as.character(program))
  prof <- x$profiles[as.character(x$profiles$program) == k, , drop = FALSE]
  if (nrow(prof) == 0) {
    stop("Program ", k, " has no extracted profile (only high and medium programs do)", call. = FALSE)
  }
  prof <- prof[as.character(prof$feature_trait_id) != x$trait$id, , drop = FALSE]
  prof <- utils::head(prof[order(-abs(prof$loading)), , drop = FALSE], top_n)
  prof$confidence <- if ("lfsr" %in% names(prof)) pmax(0, 1 - prof$lfsr) else 1
  prof$feature_trait_name <- .truncate(prof$feature_trait_name, 60)
  prof$feature_trait_name <- factor(prof$feature_trait_name, levels = rev(unique(prof$feature_trait_name)))
  prof$sign <- ifelse(prof$loading >= 0, "positive", "negative")
  return(
    ggplot2::ggplot(prof, ggplot2::aes(
      x = .data$loading, y = .data$feature_trait_name, fill = .data$sign, alpha = .data$confidence
    )) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
      ggplot2::scale_fill_manual(values = c(positive = "#2a78d6", negative = "#eb6834"), guide = "none") +
      ggplot2::scale_alpha_continuous(range = c(0.25, 1), limits = c(0, 1), name = "1 - lfsr") +
      ggplot2::labs(
        title = paste0(x$trait$name, " P", k, ": background traits"),
        x = "trait loading (oriented to the anchor's risk allele)", y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 11)
  )
}


#' @title SNP Similarity Heatmap
#' @description The cosine similarity `S` between the SNPs of a
#' [decompose_trait()] result's programs, with SNPs ordered by the first
#' program each belongs to (in confidence order) and lines between programs.
#' Blocks along the diagonal are programs whose SNPs share a profile.
#' @param x A [decompose_trait()] result.
#' @param tiers Tiers whose programs are shown. Defaults to `c("high",
#'   "medium")`.
#' @return A ggplot object.
#' @concept trait decomposition
#' @export
plot_similarity_heatmap <- function(x, tiers = c("high", "medium")) {
  .assert_class(x, "gpmap_decomposition", "plot_similarity_heatmap")
  S <- x$fit$s_matrix
  p <- programs(x, tiers = tiers)
  m <- x$memberships
  m <- m[as.character(m$program) %in% as.character(p$program), , drop = FALSE]
  m$rank <- match(as.character(m$program), as.character(p$program))
  m <- m[order(m$rank, -m$abs_loading), , drop = FALSE]
  m <- m[!duplicated(m$snp_id) & m$snp_id %in% rownames(S), , drop = FALSE]
  if (is.null(S) || nrow(m) < 2) {
    return(.empty_plot("Fewer than 2 SNPs in the shown programs", x$trait$name))
  }
  snps <- as.character(m$snp_id)
  sub_s <- S[snps, snps, drop = FALSE]
  df <- data.frame(
    row = rep(seq_along(snps), times = length(snps)),
    col = rep(seq_along(snps), each = length(snps)),
    S = as.vector(sub_s)
  )
  bounds <- which(diff(m$rank) != 0) + 0.5
  mids <- tapply(seq_along(snps), m$rank, mean)
  return(
    ggplot2::ggplot(df, ggplot2::aes(x = .data$col, y = .data$row, fill = .data$S)) +
      ggplot2::geom_raster() +
      ggplot2::geom_hline(yintercept = bounds, colour = "grey30", linewidth = 0.2) +
      ggplot2::geom_vline(xintercept = bounds, colour = "grey30", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(
        low = "#eb6834", mid = "white", high = "#2a78d6", midpoint = 0, limits = c(-1, 1), name = "S",
        na.value = "white"
      ) +
      ggplot2::scale_x_continuous(
        breaks = as.numeric(mids), labels = paste0("P", p$program[as.integer(names(mids))]),
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_reverse(
        breaks = as.numeric(mids), labels = paste0("P", p$program[as.integer(names(mids))]),
        expand = c(0, 0)
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = paste0(x$trait$name, ": SNP similarity by program"), x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)
      )
  )
}


#' @title Loading-based Enrichment of Programs
#' @description The significant loading-based enrichment hits of a
#' [decompose_trait()] result's programs (see [annotation()]): one point per
#' program and term, sized by \eqn{R^2} and coloured by \eqn{-\log_{10}} FDR.
#' @param x A [decompose_trait()] result.
#' @param type `"categories"`, `"tissues"` or `"pathways"`.
#' @param tiers Tiers whose programs are shown.
#' @param top_n Maximum terms shown (most significant first).
#' @param fdr,min_r_squared Hit cutoffs.
#' @return A ggplot object.
#' @concept trait decomposition
#' @export
plot_program_enrichment <- function(x, type = c("categories", "tissues", "pathways"),
                                    tiers = c("high", "medium"), top_n = 25L,
                                    fdr = 0.05, min_r_squared = 0.05) {
  .assert_class(x, "gpmap_decomposition", "plot_program_enrichment")
  type <- match.arg(type)
  title <- paste0(x$trait$name, ": ", type, " enrichment")
  if (is.null(x$annotation)) {
    return(.empty_plot("Not annotated", title))
  }
  hits <- annotation(x, type, fdr = fdr, min_r_squared = min_r_squared)
  keep <- as.character(programs(x, tiers = tiers)$program)
  if (nrow(hits) > 0) {
    hits <- hits[hits$significant & hits$program %in% keep, , drop = FALSE]
  }
  if (nrow(hits) == 0) {
    return(.empty_plot("No significant hits", title))
  }
  term_col <- switch(type, categories = "trait_category", tissues = "tissue", pathways = "description")
  hits$term <- .truncate(hits[[term_col]], 60)
  top_terms <- utils::head(unique(hits$term[order(hits$fdr)]), top_n)
  hits <- hits[hits$term %in% top_terms, , drop = FALSE]
  hits$term <- factor(hits$term, levels = rev(top_terms))
  hits$program_label <- factor(paste0("P", hits$program), levels = paste0("P", keep))
  hits$neg_log10_fdr <- -log10(pmax(hits$fdr, 1e-300))
  return(
    ggplot2::ggplot(hits, ggplot2::aes(
      x = .data$program_label, y = .data$term, size = .data$r_squared, colour = .data$neg_log10_fdr
    )) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_gradient(low = "#9ecae1", high = "#08306b", name = "-log10 FDR") +
      ggplot2::scale_size_continuous(range = c(1.5, 6), name = "R2") +
      ggplot2::scale_x_discrete(drop = TRUE) +
      ggplot2::labs(title = title, x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 11)
  )
}


#' @title Multi-trait Program Families
#' @description The link graph of a [trait_correspondence()] result. Each node
#' is a linked program, coloured by trait and labelled with its program
#' number; shaded hulls are families. Edge width is \eqn{|\phi|}; edge colour is
#' the link's direction, from [module_rg()] where it calls one (solid) and
#' otherwise from the sign of the profile congruence (dashed). Programs with no
#' link are not drawn.
#' @param x A [trait_correspondence()] result.
#' @return A ggplot object.
#' @concept trait decomposition
#' @export
plot_families <- function(x) {
  .assert_class(x, "gpmap_correspondence", "plot_families")
  graph <- x$graph
  title <- "Multi-trait program families"
  if (is.null(graph) || igraph::vcount(graph) == 0) {
    return(.empty_plot("No links, so no families", title))
  }
  set.seed(x$settings$seed)
  layout <- igraph::layout_with_fr(graph)
  trait_ids <- x$traits$trait_id
  trait_names <- stats::setNames(x$traits$trait_name, trait_ids)
  nodes <- data.frame(
    name = igraph::V(graph)$name,
    x = layout[, 1], y = layout[, 2],
    trait = unname(trait_names[igraph::V(graph)$trait_id]),
    family = igraph::V(graph)$family,
    label = sub("^.*:", "P", igraph::V(graph)$name),
    stringsAsFactors = FALSE
  )
  ends <- igraph::ends(graph, igraph::E(graph), names = FALSE)
  rg_dir <- igraph::E(graph)$rg_direction
  rg_called <- !is.na(rg_dir) & rg_dir %in% c("concordant", "antagonistic")
  phi <- abs(igraph::E(graph)$phi_traits)
  edges <- data.frame(
    x = layout[ends[, 1], 1], y = layout[ends[, 1], 2],
    xend = layout[ends[, 2], 1], yend = layout[ends[, 2], 2],
    direction = ifelse(rg_called, rg_dir, igraph::E(graph)$direction_profile),
    source = ifelse(rg_called, "from rg", "from profile sign"),
    width = 0.4 + 2.5 * (phi - min(phi)) / max(1e-9, max(phi) - min(phi)),
    stringsAsFactors = FALSE
  )
  radius <- 0.06 * max(diff(range(layout[, 1])), diff(range(layout[, 2])), 1)
  angles <- seq(0, 2 * pi, length.out = 13)[-13]
  hulls <- dplyr::bind_rows(lapply(split(nodes, nodes$family), function(f) {
    if (nrow(f) < 2) {
      return(NULL)
    }
    pts <- data.frame(
      x = rep(f$x, each = length(angles)) + radius * cos(angles),
      y = rep(f$y, each = length(angles)) + radius * sin(angles)
    )
    pts <- pts[grDevices::chull(pts$x, pts$y), , drop = FALSE]
    pts$family <- f$family[1]
    return(pts)
  }))
  trait_pal <- stats::setNames(grDevices::hcl.colors(length(trait_ids), "Dark 3"), unname(trait_names))
  plot <- ggplot2::ggplot()
  if (nrow(hulls) > 0) {
    plot <- plot + ggplot2::geom_polygon(
      data = hulls, ggplot2::aes(x = .data$x, y = .data$y, group = .data$family),
      fill = "grey85", colour = "grey65", alpha = 0.35
    )
  }
  plot <- plot +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(
        x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
        colour = .data$direction, linetype = .data$source, linewidth = .data$width
      )
    ) +
    ggplot2::geom_point(
      data = nodes, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$trait),
      shape = 21, size = 7, colour = "white"
    ) +
    ggplot2::geom_text(
      data = nodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      colour = "white", size = 2.6
    ) +
    ggplot2::scale_colour_manual(values = .direction_colours(), na.value = "#9e9e9e", name = "link direction") +
    ggplot2::scale_linetype_manual(values = c("from rg" = "solid", "from profile sign" = "dashed"), name = NULL) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_fill_manual(values = trait_pal, name = "trait") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 11)
  return(plot)
}


#' @title Distribution of Profile Congruence
#' @description Histogram of \eqn{|\phi|} (`phi_traits`) over every scored
#' program pair of a [trait_correspondence()] result, filled by tier, with the
#' `phi_low` and `phi_high` cutoffs marked: how many pairs each cutoff admits.
#' @param x A [trait_correspondence()] result.
#' @param binwidth Histogram bin width.
#' @return A ggplot object.
#' @concept trait decomposition
#' @export
plot_phi <- function(x, binwidth = 0.025) {
  .assert_class(x, "gpmap_correspondence", "plot_phi")
  p <- x$pairs[is.finite(x$pairs$phi_traits), , drop = FALSE]
  if (nrow(p) == 0) {
    return(.empty_plot("No scored pairs", "Profile congruence"))
  }
  p$abs_phi <- abs(p$phi_traits)
  p$tier <- factor(p$tier, levels = c("none", "low", "medium", "high"))
  cutoffs <- data.frame(
    value = c(x$settings$phi_low, x$settings$phi_high),
    cutoff = c("phi_low", "phi_high")
  )
  return(
    ggplot2::ggplot(p, ggplot2::aes(x = .data$abs_phi, fill = .data$tier)) +
      ggplot2::geom_histogram(binwidth = binwidth, boundary = 0) +
      ggplot2::geom_vline(
        data = cutoffs, ggplot2::aes(xintercept = .data$value),
        linetype = "dashed", colour = "grey30"
      ) +
      ggplot2::annotate(
        "text", x = cutoffs$value, y = Inf, label = cutoffs$cutoff,
        vjust = 1.5, hjust = 1.05, size = 3, colour = "grey30"
      ) +
      ggplot2::scale_fill_manual(values = .tier_colours(), name = "tier", drop = FALSE) +
      ggplot2::labs(
        title = "Profile congruence of every program pair",
        x = "|phi| (Tucker congruence of trait loadings)", y = "program pairs"
      ) +
      ggplot2::theme_minimal(base_size = 11)
  )
}


#' @title Trait Profiles of Two Linked Programs
#' @description The two programs' trait loadings over the background traits
#' their profile congruence is computed on (both fits carry them; anchors and
#' proxies left out), one point per trait. Points along the diagonal are the
#' agreement \eqn{\phi} measures (along the anti-diagonal for an antagonistic
#' pair); the traits contributing most are labelled.
#' @param x A [trait_correspondence()] result.
#' @param program_a,program_b Program ids (e.g. `"1631:3"`), as in
#'   [links()].
#' @param n_labels Number of traits labelled.
#' @return A ggplot object.
#' @concept trait decomposition
#' @export
plot_link <- function(x, program_a, program_b, n_labels = 8L) {
  .assert_class(x, "gpmap_correspondence", "plot_link")
  profiles <- dplyr::bind_rows(lapply(x$decompositions, function(d) return(d$profiles)))
  pa <- profiles[profiles$program_id == program_a, , drop = FALSE]
  pb <- profiles[profiles$program_id == program_b, , drop = FALSE]
  if (nrow(pa) == 0 || nrow(pb) == 0) {
    stop("Both programs must be high or medium programs of this correspondence", call. = FALSE)
  }
  ta <- as.character(pa$trait_id[1])
  tb <- as.character(pb$trait_id[1])
  drop <- c(ta, tb, x$decompositions[[ta]]$exclude_trait_ids, x$decompositions[[tb]]$exclude_trait_ids)
  df <- dplyr::inner_join(
    pa[, c("feature_trait_id", "feature_trait_name", "loading")],
    pb[, c("feature_trait_id", "loading")],
    by = "feature_trait_id", suffix = c("_a", "_b")
  )
  df <- df[!as.character(df$feature_trait_id) %in% drop, , drop = FALSE]
  df$contribution <- abs(df$loading_a * df$loading_b)
  top <- utils::head(order(-df$contribution), n_labels)
  df$label <- ""
  df$label[top] <- .truncate(df$feature_trait_name[top], 40)
  in_pair <- (x$pairs$program_id_a == program_a & x$pairs$program_id_b == program_b) |
    (x$pairs$program_id_a == program_b & x$pairs$program_id_b == program_a)
  pair <- x$pairs[in_pair, , drop = FALSE]
  subtitle <- if (nrow(pair) == 1) {
    sprintf("phi = %.2f, tier %s, loci %s", pair$phi_traits, pair$tier, pair$locus_status)
  } else {
    NULL
  }
  label_a <- .program_display_label(program_a, x$decompositions[[ta]]$trait$name)
  label_b <- .program_display_label(program_b, x$decompositions[[tb]]$trait$name)
  text_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(ggplot2::aes(label = .data$label), size = 3, max.overlaps = Inf)
  } else {
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3, vjust = -0.6)
  }
  return(
    ggplot2::ggplot(df, ggplot2::aes(x = .data$loading_a, y = .data$loading_b)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
      ggplot2::geom_point(ggplot2::aes(size = .data$contribution), colour = "#2a78d6", alpha = 0.6) +
      text_layer +
      ggplot2::scale_size_continuous(range = c(0.8, 4), guide = "none") +
      ggplot2::labs(
        title = paste(label_a, "vs", label_b), subtitle = subtitle,
        x = paste(label_a, "trait loading"), y = paste(label_b, "trait loading")
      ) +
      ggplot2::theme_minimal(base_size = 11)
  )
}


#' @title Loci of a Program Family
#' @description Each member program's SNP loading at the loci its family
#' claims, on GPMap's shared locus numbering: one row per program, one column
#' per locus. Points mark loci a program claims as a high-confidence member.
#' Loci claimed by every member (the family's core) come first.
#' @param x A [trait_correspondence()] result.
#' @param family A family name, e.g. `"family_1"`.
#' @return A ggplot object.
#' @concept trait decomposition
#' @export
plot_family_loci <- function(x, family) {
  .assert_class(x, "gpmap_correspondence", "plot_family_loci")
  fl <- family_loci(x, family = family)
  members <- x$nodes[x$nodes$family == family, , drop = FALSE]
  if (nrow(fl) == 0 || nrow(members) == 0) {
    return(.empty_plot("No loci for this family", family))
  }
  loadings <- dplyr::bind_rows(lapply(x$decompositions, function(d) return(d$loadings)))
  loadings <- loadings[loadings$program_id %in% members$program_id & !is.na(loadings$coloc_group_id), , drop = FALSE]
  loadings$coloc_group_id <- as.character(loadings$coloc_group_id)
  loci_map <- dplyr::distinct(x$loci[, c("trait_id", "coloc_group_id", "locus")])
  df <- dplyr::inner_join(
    loadings[, c("program_id", "trait_id", "trait_name", "coloc_group_id", "loading", "high_confidence")],
    loci_map,
    by = c("trait_id", "coloc_group_id")
  )
  df <- df[df$locus %in% fl$locus, , drop = FALSE]
  if (nrow(df) == 0) {
    return(.empty_plot("No loadings at the family's loci", family))
  }
  locus_order <- fl$locus[order(-fl$n_programs, fl$locus)]
  df$locus <- factor(df$locus, levels = unique(locus_order))
  df$program <- .program_display_label(df$program_id, df$trait_name)
  return(
    ggplot2::ggplot(df, ggplot2::aes(x = .data$locus, y = .data$program, fill = .data$loading)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_point(data = df[df$high_confidence, , drop = FALSE], size = 0.8, colour = "grey15") +
      ggplot2::scale_fill_gradient2(low = "#eb6834", mid = "white", high = "#2a78d6", midpoint = 0, name = "loading") +
      ggplot2::labs(
        title = paste0(family, ": loadings at the family's loci"),
        subtitle = "points: loci the program claims as a high-confidence member",
        x = "locus (coloc group)", y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 6))
  )
}
