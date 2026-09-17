#' @title Continuous Trait-Category Enrichment From EBMF Loadings
#' @description Parallel to `enrich_snp_group_trait_categories()`, but without
#' any hard lFSR/magnitude cutoff: every SNP with a finite EBMF loading
#' contributes to every program's test, weighted by its squared loading,
#' regressed against each trait category's local 0/1 membership indicator
#' (`lm(loading^2 ~ category_membership)`), following the approach in
#' `scripts/enrichment_via_loadings_example.r`. Every trait category tested for
#' every program is one Benjamini-Hochberg family, so `fdr` reflects the true
#' number of regressions run for this trait -- not one program's slice of
#' them. Tissue and pathway enrichment
#' (`enrich_program_loadings_tissues()`, `enrich_program_loadings_pathways()`)
#' are each their own separate family; neither is pooled with this one.
#'
#' IMPORTANT: as in that script, SNPs are treated as independent
#' observations; for formal inference, LD should be accounted for (e.g. SNP
#' -> gene aggregation, or cluster-robust SEs / permutation using LD blocks).
#' The enrichment effect estimates are nevertheless useful for comparison
#' against the hard-cutoff hypergeometric trait-category test.
#' @param clustering_result Result of `run_univariate_clustering()`.
#' @param coloc_groups Coloc-group dataframe with a `trait_category` column.
#' @param snp_key Column used to match SNP ids in `coloc_groups`.
#' @param min_category_size Only test a category if at least this many SNPs
#'   link to it. Defaults to 5.
#' @param min_loading_magnitude Soft membership gate: a SNP only counts
#'   towards a category's `x = 1` group if its absolute loading on that
#'   program exceeds this value; otherwise it falls back to `x = 0` for that
#'   program's test (it stays in the analysis, just not counted as a member).
#'   This keeps near-zero, shrinkage-only loadings from diluting the mean
#'   squared-loading contrast without imposing a hard lFSR/magnitude gate.
#'   Defaults to 0.02.
#' @return A list with:
#'   \itemize{
#'     \item by_program: list of per-program results (`program`, `n_snps`,
#'       `comparison` with columns trait_category, enrichment, se, z, p, fdr,
#'       r_squared, n_snps, n_category_snps; `fdr` is corrected across every
#'       program tested, not within this program alone. `r_squared` is the
#'       fraction of this program's squared-loading variance the category
#'       explains -- unlike `enrichment` it is bounded in [0, 1] and
#'       comparable across programs, which is what makes a fixed threshold on
#'       it (e.g. Cohen's 1988 benchmarks: 0.01 small, 0.09 medium, 0.25
#'       large) meaningful. It is not an independent effect-size-vs-power
#'       axis here: with `n_snps` roughly fixed per program, `r_squared` is
#'       close to a monotonic function of `p`.
#'     \item summary: one row per program (`n_categories_tested`,
#'       `n_enriched`, `top_category`)
#'   }
#' @export
enrich_program_loadings_trait_categories <- function(clustering_result,
                                                      coloc_groups,
                                                      snp_key = c(
                                                        "variant_id",
                                                        "display_snp",
                                                        "coloc_group_id"
                                                      ),
                                                      min_category_size = 5L,
                                                      min_loading_magnitude = 0) {
  snp_key <- match.arg(snp_key)
  return(.enrich_program_loadings_category(
    clustering_result = clustering_result,
    coloc_groups = coloc_groups,
    snp_key = snp_key,
    value_col = "trait_category",
    min_category_size = min_category_size,
    min_loading_magnitude = min_loading_magnitude
  ))
}


#' @title Continuous Tissue Enrichment From EBMF Loadings
#' @description Parallel to `enrich_snp_group_tissues()`, but without any hard
#' lFSR/magnitude cutoff: every SNP with a finite EBMF loading contributes to
#' every program's test, weighted by its squared loading, regressed against
#' each tissue's local 0/1 membership indicator
#' (`lm(loading^2 ~ tissue_membership)`), following the approach in
#' `scripts/enrichment_via_loadings_example.r`. Every tissue tested for every
#' program is one Benjamini-Hochberg family, so `fdr` reflects the true number
#' of regressions run for this trait -- not one program's slice of them.
#' Trait-category and pathway enrichment
#' (`enrich_program_loadings_trait_categories()`,
#' `enrich_program_loadings_pathways()`) are each their own separate family;
#' neither is pooled with this one.
#'
#' IMPORTANT: as in that script, SNPs are treated as independent
#' observations; for formal inference, LD should be accounted for (e.g. SNP
#' -> gene aggregation, or cluster-robust SEs / permutation using LD blocks).
#' The enrichment effect estimates are nevertheless useful for comparison
#' against the hard-cutoff hypergeometric tissue test.
#' @param clustering_result Result of `run_univariate_clustering()`.
#' @param coloc_groups Coloc-group dataframe with a `tissue` column.
#' @param snp_key Column used to match SNP ids in `coloc_groups`.
#' @param min_category_size Only test a tissue if at least this many SNPs
#'   link to it. Defaults to 5.
#' @param min_loading_magnitude Soft membership gate: a SNP only counts
#'   towards a tissue's `x = 1` group if its absolute loading on that program
#'   exceeds this value; otherwise it falls back to `x = 0` for that
#'   program's test (it stays in the analysis, just not counted as a member).
#'   This keeps near-zero, shrinkage-only loadings from diluting the mean
#'   squared-loading contrast without imposing a hard lFSR/magnitude gate.
#'   Defaults to 0.02.
#' @return A list with:
#'   \itemize{
#'     \item by_program: list of per-program results (`program`, `n_snps`,
#'       `comparison` with columns tissue, enrichment, se, z, p, fdr,
#'       r_squared, n_snps, n_category_snps; `fdr` is corrected across every
#'       program tested, not within this program alone. `r_squared` is the
#'       fraction of this program's squared-loading variance the tissue
#'       explains -- unlike `enrichment` it is bounded in [0, 1] and
#'       comparable across programs, which is what makes a fixed threshold on
#'       it (e.g. Cohen's 1988 benchmarks: 0.01 small, 0.09 medium, 0.25
#'       large) meaningful. It is not an independent effect-size-vs-power
#'       axis here: with `n_snps` roughly fixed per program, `r_squared` is
#'       close to a monotonic function of `p`.
#'     \item summary: one row per program (`n_categories_tested`,
#'       `n_enriched`, `top_category`)
#'   }
#' @export
enrich_program_loadings_tissues <- function(clustering_result,
                                            coloc_groups,
                                            snp_key = c(
                                              "variant_id",
                                              "display_snp",
                                              "coloc_group_id"
                                            ),
                                            min_category_size = 5L,
                                            min_loading_magnitude = 0) {
  snp_key <- match.arg(snp_key)
  return(.enrich_program_loadings_category(
    clustering_result = clustering_result,
    coloc_groups = coloc_groups,
    snp_key = snp_key,
    value_col = "tissue",
    min_category_size = min_category_size,
    min_loading_magnitude = min_loading_magnitude
  ))
}


.enrich_program_loadings_category <- function(clustering_result,
                                               coloc_groups,
                                               snp_key,
                                               value_col,
                                               min_category_size,
                                               min_loading_magnitude = 0) {
  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("coloc_groups is required")
  }
  if (!value_col %in% names(coloc_groups)) {
    stop("coloc_groups must include a ", value_col, " column")
  }

  posterior <- ebmf_posterior_table(clustering_result)
  posterior <- posterior[is.finite(posterior$loading), , drop = FALSE]

  empty_summary <- data.frame(
    program = integer(0), n_snps = integer(0), n_categories_tested = integer(0),
    n_enriched = integer(0), top_category = character(0),
    stringsAsFactors = FALSE
  )
  if (nrow(posterior) == 0) {
    return(list(by_program = list(), summary = empty_summary))
  }

  snp_ids <- unique(posterior$snp_id)
  membership <- .snp_category_matrix(
    snp_ids = snp_ids,
    coloc_groups = coloc_groups,
    snp_key = snp_key,
    value_col = value_col
  )

  # Every program's test is run first, without correcting for multiplicity, so
  # the BH family below spans every (program x category) test actually run for
  # this enrichment type -- not just one program's slice of it.
  programs <- sort(unique(posterior$program))
  raw_by_program <- lapply(programs, function(pg) {
    prog_loadings <- posterior[posterior$program == pg, , drop = FALSE]
    y <- prog_loadings$loading[match(snp_ids, prog_loadings$snp_id)]
    comparison <- .continuous_link_enrichment(
      y, membership,
      min_category_size = min_category_size,
      min_loading_magnitude = min_loading_magnitude
    )
    names(comparison)[names(comparison) == "value"] <- value_col
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
      comparison <- comparison[order(comparison$fdr, comparison$p), , drop = FALSE]
      rownames(comparison) <- NULL
    } else {
      comparison$fdr <- numeric(0)
    }
    offset <<- offset + n
    return(list(program = x$program, n_snps = x$n_snps, comparison = comparison))
  })

  summary_df <- dplyr::bind_rows(lapply(by_program, function(x) {
    top_category <- if (nrow(x$comparison) > 0) {
      x$comparison[[value_col]][1]
    } else {
      NA_character_
    }
    data.frame(
      program = x$program,
      n_snps = x$n_snps,
      n_categories_tested = nrow(x$comparison),
      n_enriched = sum(x$comparison$fdr <= 0.05, na.rm = TRUE),
      top_category = top_category,
      stringsAsFactors = FALSE
    )
  }))

  return(list(by_program = by_program, summary = summary_df))
}


# SNP x value 0/1 membership matrix (e.g. SNP x tissue, SNP x trait_category),
# pivoted from the long snp_id/value link tables already used by the
# hard-cutoff hypergeometric enrichment functions.
.snp_category_matrix <- function(snp_ids,
                                 coloc_groups,
                                 snp_key,
                                 value_col = c("tissue", "trait_category")) {
  value_col <- match.arg(value_col)
  snp_ids <- unique(as.character(snp_ids))
  links <- if (identical(value_col, "tissue")) {
    .snp_tissue_links(snp_ids, coloc_groups, snp_key)
  } else {
    .snp_category_links(snp_ids, coloc_groups, snp_key)
  }

  values <- sort(unique(links[[value_col]]))
  mat <- matrix(
    0L,
    nrow = length(snp_ids), ncol = length(values),
    dimnames = list(snp_ids, values)
  )
  if (nrow(links) > 0 && length(values) > 0) {
    mat[cbind(links$snp_id, links[[value_col]])] <- 1L
  }
  return(mat)
}


# Continuous analog of .hypergeometric_link_enrichment(): for each column of
# a SNP x value 0/1 membership matrix, test whether squared loading (`y`) is
# systematically greater for member SNPs via lm(y^2 ~ membership), as in
# scripts/enrichment_via_loadings_example.r. Returns one row per tested value
# with a "value" column (renamed by callers to trait_category/tissue/term_id)
# plus enrichment/se/z/p/r_squared and SNP counts -- NOT fdr. This tests one
# program (and, for pathways, one source); FDR correction is the caller's job
# once every program (and source) actually tested has been assembled, so the
# BH family reflects the true number of tests run rather than just this one
# call's slice of them. Columns not meeting `min_category_size`, or with no
# variation, are skipped.
#
# `r_squared` is the fraction of this program's squared-loading variance
# membership explains -- unlike `enrichment` (a raw mean-difference, in
# squared-loading units that vary by program and aren't comparable across
# programs) it is bounded in [0, 1] and comparable across programs, which is
# what makes a single fixed threshold on it meaningful. It is not an
# independent "effect size vs significance" axis here: with n_snps roughly
# fixed per program, r_squared is close to a monotonic function of p.
#
# `min_loading_magnitude` softly regates membership: a SNP linked to a value
# (x == 1) but with abs(loading) at or below this threshold is treated as
# x == 0 for that test instead, so shrinkage-only loadings near zero don't
# dilute the member-group mean. The SNP is not dropped from the analysis, so
# the background pool used for the residual-variance/SE estimate is
# unaffected.
.continuous_link_enrichment <- function(y, membership, min_category_size = 5L,
                                        min_loading_magnitude = 0) {
  empty <- data.frame(
    value = character(0), enrichment = numeric(0), se = numeric(0),
    z = numeric(0), p = numeric(0), r_squared = numeric(0),
    n_snps = integer(0), n_category_snps = integer(0),
    stringsAsFactors = FALSE
  )
  if (is.null(membership) || ncol(membership) == 0 || length(y) == 0) {
    return(empty)
  }

  y2 <- y^2
  keep <- is.finite(y2)
  y2 <- y2[keep]
  membership <- membership[keep, , drop = FALSE]
  if (length(y2) == 0) {
    return(empty)
  }

  values <- colnames(membership)
  rows <- lapply(seq_along(values), function(j) {
    x <- membership[, j]
    if (min_loading_magnitude > 0) {
      x[x == 1 & y2 <= min_loading_magnitude^2] <- 0L
    }
    n_category_snps <- sum(x == 1)
    if (n_category_snps < min_category_size || length(unique(x)) < 2) {
      return(NULL)
    }
    fit <- stats::lm(y2 ~ x)
    sm <- summary(fit)$coef
    if (!"x" %in% rownames(sm)) {
      return(NULL)
    }
    data.frame(
      value = values[j],
      enrichment = sm["x", "Estimate"],
      se = sm["x", "Std. Error"],
      z = sm["x", "t value"],
      p = sm["x", "Pr(>|t|)"],
      r_squared = summary(fit)$r.squared,
      n_snps = length(y2),
      n_category_snps = n_category_snps,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(empty)
  }

  out <- dplyr::bind_rows(rows)
  out <- out[order(out$p, out$value), , drop = FALSE]
  rownames(out) <- NULL
  return(out)
}
