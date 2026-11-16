#' @title Programs of a Decomposition or Correspondence
#' @description One row per program, with its confidence tier and validation
#' metrics (see [summarise_ebmf_programs()]); for a [trait_correspondence()]
#' result, every trait's programs with the `family` each linked program
#' belongs to.
#' @param x A [decompose_trait()] or [trait_correspondence()] result.
#' @param tiers Confidence tiers to keep. Defaults to `c("high", "medium")`;
#'   `NULL` keeps every program.
#' @param ... Unused.
#' @return A data.frame.
#' @concept trait decomposition
#' @export
programs <- function(x, ...) {
  UseMethod("programs")
}


#' @rdname programs
#' @export
programs.gpmap_decomposition <- function(x, tiers = c("high", "medium"), ...) {
  p <- x$programs
  if (!is.null(tiers)) {
    p <- p[as.character(p$confidence_tier) %in% tiers, , drop = FALSE]
  }
  front <- c("program_id", "trait_id", "trait_name", "program", "program_label", "confidence_tier")
  p <- p[, c(intersect(front, names(p)), setdiff(names(p), front)), drop = FALSE]
  rownames(p) <- NULL
  return(p)
}


#' @rdname programs
#' @export
programs.gpmap_correspondence <- function(x, tiers = c("high", "medium"), ...) {
  p <- dplyr::bind_rows(lapply(x$decompositions, programs, tiers = tiers))
  family <- x$nodes$family[match(p$program_id, x$nodes$program_id)]
  p$family <- if (length(family) == 0) character(0) else family
  return(p)
}


#' @title Cross-trait Program Links
#' @description The scored program pairs of a [trait_correspondence()] result,
#' by tier and then by \eqn{|\phi|}, with display labels
#' (`"<trait name> P<k>"`). Links are pairs of tier high or medium.
#' @param x A [trait_correspondence()] result.
#' @param tiers Tiers to keep. Defaults to `c("high", "medium")`; `NULL` keeps
#'   every scored pair.
#' @param ... Unused.
#' @return A data.frame, one row per program pair: the scores and tiers of
#'   [compare_program_pairs_loadings()], the [module_rg()] estimate and
#'   `direction_check`, `top_traits` ([profile_contributions()]), and
#'   `label_a`, `label_b`.
#' @concept trait decomposition
#' @export
links <- function(x, ...) {
  UseMethod("links")
}


#' @rdname links
#' @export
links.gpmap_correspondence <- function(x, tiers = c("high", "medium"), ...) {
  p <- x$pairs
  if (!is.null(tiers)) {
    p <- p[p$tier %in% tiers, , drop = FALSE]
  }
  tier_levels <- c("high", "medium", "low", "none")
  p <- p[order(match(p$tier, tier_levels), -abs(p$phi_traits)), , drop = FALSE]
  p$label_a <- .program_display_label(p$program_id_a, p$trait_name_a)
  p$label_b <- .program_display_label(p$program_id_b, p$trait_name_b)
  front <- c("label_a", "label_b", "tier", "profile_strength", "locus_status", "phi_traits")
  p <- p[, c(front, setdiff(names(p), front)), drop = FALSE]
  rownames(p) <- NULL
  return(p)
}


#' @title Multi-trait Program Families
#' @description The families of a [trait_correspondence()] result (see
#' [build_program_families()]), and the loci their programs claim.
#' @param x A [trait_correspondence()] result.
#' @param family Optional family names (e.g. `"family_1"`) to keep.
#' @param ... Unused.
#' @return `families()`: one row per family. `family_loci()`: one row per
#'   family x locus.
#' @concept trait decomposition
#' @export
families <- function(x, ...) {
  UseMethod("families")
}


#' @rdname families
#' @export
families.gpmap_correspondence <- function(x, ...) {
  return(x$families)
}


#' @rdname families
#' @export
family_loci <- function(x, ...) {
  UseMethod("family_loci")
}


#' @rdname families
#' @export
family_loci.gpmap_correspondence <- function(x, family = NULL, ...) {
  out <- x$family_loci
  if (!is.null(family)) {
    out <- out[out$family %in% family, , drop = FALSE]
  }
  return(out)
}


#' @title Traits Removed Before Fitting
#' @description The proxy and duplicate traits removed from each anchor's
#' matrix (see [decompose_trait()]), with where each came from (`source`:
#' `duplicate_of_anchor`, `proxy_group`, `anchor_trait`, `user` or
#' `background_duplicate`) and, for the anchor's proxies, how closely it tracks
#' the anchor at the anchor's own loci (see [anchor_proxies()]).
#' @param x A [decompose_trait()] or [trait_correspondence()] result.
#' @param ... Unused.
#' @return A data.frame, one row per removed trait and anchor.
#' @concept trait decomposition
#' @export
excluded <- function(x, ...) {
  UseMethod("excluded")
}


#' @rdname excluded
#' @export
excluded.gpmap_decomposition <- function(x, ...) {
  return(x$excluded)
}


#' @rdname excluded
#' @export
excluded.gpmap_correspondence <- function(x, ...) {
  return(dplyr::bind_rows(lapply(x$decompositions, excluded)))
}


#' @title Biological Annotation of Programs or Families
#' @description Tidy tables of a result's annotation.
#'
#' For a [decompose_trait()] result, per program: `"specific_traits"`
#' (module-specific traits, [summarise_module_specific_traits()]),
#' `"categories"`, `"tissues"` and `"pathways"` (loading-based enrichment,
#' every tested term with `significant` marking FDR <= `fdr` and R^2 >=
#' `min_r_squared`), `"genes"` (colocalizing QTL genes at each program's SNPs)
#' and `"distinctness"` (pairwise program similarity, [program_distinctness()]).
#'
#' For a [trait_correspondence()] result, per family and trait: `"summary"`
#' (SNP counts and the top hits as text), `"categories"`, `"tissues"` and
#' `"pathways"` ([enrich_snp_set()], [enrich_snp_set_pathways()]).
#' @param x A [decompose_trait()] or [trait_correspondence()] result.
#' @param type Which table.
#' @param fdr,min_r_squared Cutoffs marking `significant` loading-based hits.
#' @param ... Unused.
#' @return A data.frame (empty when the result was not annotated).
#' @concept trait decomposition
#' @export
annotation <- function(x, ...) {
  UseMethod("annotation")
}


#' @rdname annotation
#' @export
annotation.gpmap_decomposition <- function(x,
                                           type = c(
                                             "specific_traits", "categories", "tissues",
                                             "pathways", "genes", "distinctness"
                                           ),
                                           fdr = 0.05,
                                           min_r_squared = 0.05,
                                           ...) {
  type <- match.arg(type)
  a <- x$annotation
  if (is.null(a)) {
    return(data.frame())
  }
  kept <- as.character(x$programs$program)
  out <- switch(type,
    specific_traits = {
      rows <- lapply(a$specific_traits$by_group, function(g) {
        if (is.null(g$traits) || nrow(g$traits) == 0) {
          return(NULL)
        }
        t <- g$traits
        t$group <- NULL
        return(data.frame(program = as.character(g$group), t, stringsAsFactors = FALSE))
      })
      dplyr::bind_rows(rows)
    },
    categories = .loading_hits(a$categories, kept, fdr, min_r_squared),
    tissues = .loading_hits(a$tissues, kept, fdr, min_r_squared),
    pathways = .loading_hits(a$pathways, kept, fdr, min_r_squared),
    genes = a$genes,
    distinctness = if (is.null(a$distinctness)) NULL else a$distinctness$pairs
  )
  if (is.null(out)) {
    return(data.frame())
  }
  rownames(out) <- NULL
  return(out)
}


#' @rdname annotation
#' @export
annotation.gpmap_correspondence <- function(x,
                                            type = c("summary", "categories", "tissues", "pathways"),
                                            fdr = 0.05,
                                            ...) {
  type <- match.arg(type)
  fe <- x$family_enrichment
  if (is.null(fe)) {
    return(data.frame())
  }
  if (type != "summary") {
    out <- fe[[type]]
    return(if (is.null(out)) data.frame() else out)
  }
  out <- fe$summary
  if (nrow(out) == 0) {
    return(out)
  }
  top_of <- function(table, value_col, n = 4L) {
    return(vapply(seq_len(nrow(out)), function(i) {
      if (out$status[i] == "too few high-confidence SNPs") {
        return("too few SNPs")
      }
      if (is.null(table) || nrow(table) == 0) {
        return("none")
      }
      keep <- table$family == out$family[i] & table$trait_id == out$trait_id[i] &
        !is.na(table$fdr) & table$fdr <= fdr
      hits <- table[keep, , drop = FALSE]
      if (nrow(hits) == 0) {
        return("none")
      }
      hits <- utils::head(hits[order(hits$fdr), , drop = FALSE], n)
      return(paste0(hits[[value_col]], " (FDR=", signif(hits$fdr, 2), ")", collapse = "; "))
    }, character(1)))
  }
  out$categories <- top_of(fe$categories, "trait_category")
  out$tissues <- top_of(fe$tissues, "tissue")
  pw <- fe$pathways
  for (src in c("Reactome", "KEGG")) {
    out[[tolower(src)]] <- top_of(if (is.null(pw)) NULL else pw[pw$source == src, , drop = FALSE],
      "description",
      n = 1L
    )
  }
  out[out$status == "pathways unavailable", c("reactome", "kegg")] <- "unavailable"
  return(out)
}


#' @title Program Summary Cards
#' @description One row per program, summarising what each one is about: its
#' confidence and size, the background traits most specific to it, its top
#' loading-based trait-category, tissue and pathway hits, and the colocalizing
#' QTL genes at its SNPs (strongest loading first). Every column is plain text
#' (items separated by `"; "`), so a vignette or report can add its own
#' columns by joining on `program_id`.
#' @param x A [decompose_trait()] or [trait_correspondence()] result.
#' @param tiers Confidence tiers to include. Defaults to `c("high",
#'   "medium")`; `NULL` includes every program.
#' @param n Number of items listed per annotation column.
#' @param n_genes Number of genes listed.
#' @param fdr,min_r_squared Cutoffs for a loading-based hit.
#' @return A data.frame with `program_id`, `trait_name`, `program`,
#'   `program_label`, `confidence_tier`, `n_snps_filtered`,
#'   `effective_n_traits`, `specific_traits`, `categories`, `tissues`,
#'   `reactome`, `kegg` and `genes`.
#' @concept trait decomposition
#' @export
program_cards <- function(x, tiers = c("high", "medium"), n = 3L, n_genes = 10L,
                          fdr = 0.05, min_r_squared = 0.05) {
  if (inherits(x, "gpmap_correspondence")) {
    return(dplyr::bind_rows(lapply(
      x$decompositions, program_cards,
      tiers = tiers, n = n, n_genes = n_genes, fdr = fdr, min_r_squared = min_r_squared
    )))
  }
  if (!inherits(x, "gpmap_decomposition")) {
    stop("x must come from decompose_trait() or trait_correspondence()", call. = FALSE)
  }
  p <- programs(x, tiers = tiers)
  cards <- data.frame(
    program_id = p$program_id,
    trait_name = p$trait_name,
    program = as.character(p$program),
    program_label = p$program_label,
    confidence_tier = as.character(p$confidence_tier),
    n_snps_filtered = p$n_snps_filtered,
    effective_n_traits = p$effective_n_traits,
    stringsAsFactors = FALSE
  )
  if (nrow(cards) == 0) {
    return(cards)
  }
  annotated <- !is.null(x$annotation)
  none <- if (annotated) "none" else "not annotated"
  collapse_by_program <- function(table, value_col, max_n, order_col = NULL) {
    if (!annotated || is.null(table) || nrow(table) == 0) {
      return(rep(none, nrow(cards)))
    }
    if (!is.null(order_col)) {
      table <- table[order(-table[[order_col]]), , drop = FALSE]
    }
    return(vapply(cards$program, function(k) {
      vals <- unique(as.character(table[[value_col]][as.character(table$program) == k]))
      vals <- utils::head(vals, max_n)
      return(if (length(vals) == 0) none else paste(vals, collapse = "; "))
    }, character(1), USE.NAMES = FALSE))
  }
  hits <- function(type) {
    h <- annotation(x, type, fdr = fdr, min_r_squared = min_r_squared)
    if (nrow(h) == 0) {
      return(h)
    }
    return(h[h$significant, , drop = FALSE])
  }
  st <- annotation(x, "specific_traits")
  if (nrow(st) > 0) {
    st <- st[st$passes_min_specificity, , drop = FALSE]
  }
  cards$specific_traits <- collapse_by_program(st, "trait_name", n)
  cards$categories <- collapse_by_program(hits("categories"), "trait_category", n, "enrichment")
  cards$tissues <- collapse_by_program(hits("tissues"), "tissue", n, "enrichment")
  pw <- hits("pathways")
  for (src in c("Reactome", "KEGG")) {
    sub_pw <- if (nrow(pw) == 0) pw else pw[pw$source == src, , drop = FALSE]
    cards[[tolower(src)]] <- if (annotated && identical(x$annotation$pathways_status, "unavailable")) {
      rep("unavailable", nrow(cards))
    } else {
      collapse_by_program(sub_pw, "description", n, "enrichment")
    }
  }
  genes <- annotation(x, "genes")
  cards$genes <- collapse_by_program(genes, "gene", n_genes)
  return(cards)
}


# Every tested term of a loading-based enrichment result, for the kept
# programs, with `significant` marking a hit.
.loading_hits <- function(res, kept, fdr, min_r_squared) {
  if (is.null(res) || length(res$by_program) == 0) {
    return(NULL)
  }
  rows <- lapply(res$by_program, function(b) {
    if (!as.character(b$program) %in% kept || is.null(b$comparison) || nrow(b$comparison) == 0) {
      return(NULL)
    }
    return(data.frame(program = as.character(b$program), b$comparison, stringsAsFactors = FALSE))
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    return(NULL)
  }
  out$significant <- is.finite(out$fdr) & out$fdr <= fdr &
    is.finite(out$r_squared) & out$r_squared >= min_r_squared
  return(out)
}


# The top module-specific traits of each program as one string, named by
# program.
.top_specific_traits <- function(x, n = 3L) {
  st <- annotation(x, "specific_traits")
  if (nrow(st) == 0) {
    return(stats::setNames(character(0), character(0)))
  }
  st <- st[st$passes_min_specificity, , drop = FALSE]
  progs <- unique(st$program)
  out <- vapply(progs, function(k) {
    return(paste(utils::head(st$trait_name[st$program == k], n), collapse = "; "))
  }, character(1))
  return(stats::setNames(out, progs))
}


# "<trait name> P<k>" for a program id "<trait_id>:<k>".
.program_display_label <- function(program_id, trait_name) {
  k <- sub("^.*:", "", as.character(program_id))
  return(paste0(trait_name, " P", k))
}
