#' @title Order EBMF Programs by Confidence, Children Under Parents
#' @description Order a program table (typically
#' `summarise_ebmf_programs()$programs`) for reporting. Top-level programs are
#' sorted by `confidence_tier` (high, medium, single_trait, low), then by
#' `similarity_z`
#' (descending) and `n_snps_filtered` (descending). A program with a
#' `parent_program` in the table is placed directly under its parent, and its
#' own children under it in turn, sorted the same way among siblings. A program
#' whose parent is not in the table is treated as top level.
#' @param programs Dataframe with at least `program` and `confidence_tier`;
#'   `parent_program`, `similarity_z` and `n_snps_filtered` are used when
#'   present.
#' @return `programs`, reordered, with two added columns: `depth` (0 for a
#'   top-level program, 1 for a child, and so on) and `program_label` (the
#'   program id, indented and prefixed with an arrow for children, for display
#'   in HTML tables).
#' @export
order_programs_by_confidence <- function(programs) {
  if (!is.data.frame(programs)) {
    stop("programs must be a dataframe")
  }
  if (!all(c("program", "confidence_tier") %in% names(programs))) {
    stop("programs must have program and confidence_tier columns")
  }
  n <- nrow(programs)
  if (n == 0) {
    programs$depth <- integer(0)
    programs$program_label <- character(0)
    return(programs)
  }

  tier_rank <- match(
    as.character(programs$confidence_tier),
    c("high", "medium", "single_trait", "low")
  )
  tier_rank[is.na(tier_rank)] <- 5L
  sim_z <- if ("similarity_z" %in% names(programs)) programs$similarity_z else rep(NA_real_, n)
  sim_z[!is.finite(sim_z)] <- -Inf
  size <- if ("n_snps_filtered" %in% names(programs)) programs$n_snps_filtered else rep(0, n)
  size[!is.finite(size)] <- 0
  base_order <- order(tier_rank, -sim_z, -size)

  ids <- as.character(programs$program)
  parent <- if ("parent_program" %in% names(programs)) {
    as.character(programs$parent_program)
  } else {
    rep(NA_character_, n)
  }
  parent[!parent %in% ids] <- NA_character_

  visited <- rep(FALSE, n)
  depth <- rep(0L, n)
  walk_order <- integer(0)
  visit <- function(i, d) {
    if (visited[i]) {
      return(invisible(NULL))
    }
    visited[i] <<- TRUE
    depth[i] <<- d
    walk_order <<- c(walk_order, i)
    children <- base_order[!is.na(parent[base_order]) & parent[base_order] == ids[i]]
    for (child in children) {
      visit(child, d + 1L)
    }
    return(invisible(NULL))
  }
  for (i in base_order[is.na(parent[base_order])]) {
    visit(i, 0L)
  }
  # A parent cycle cannot arise from summarise_ebmf_programs() (parents are
  # strictly larger), but keep anything left over rather than dropping rows.
  for (i in base_order[!visited[base_order]]) {
    visit(i, 0L)
  }

  out <- programs[walk_order, , drop = FALSE]
  out$depth <- depth[walk_order]
  indent <- vapply(out$depth, function(d) {
    return(strrep("  ", d))
  }, character(1))
  out$program_label <- paste0(
    indent,
    ifelse(out$depth > 0, "↳ ", ""),
    as.character(out$program)
  )
  rownames(out) <- NULL
  return(out)
}


#' @title HTML Table Coloured by Confidence Tier
#' @description Render a dataframe as an HTML `knitr::kable()` table with each
#' row's background coloured by its confidence tier: green for high, orange for
#' medium, red for low, grey for single_trait. Intended for vignettes and HTML
#' reports.
#' @param df Dataframe to render.
#' @param tier_col Name of the column holding the tier (`"high"`, `"medium"`,
#'   `"low"` or `"single_trait"`). Rows with another value are left uncoloured.
#' @param colours Named character vector of background colours for `high`,
#'   `medium`, `low` and `single_trait`.
#' @param drop_tier_col If `TRUE`, the tier column is used for the row colours
#'   but not shown as a column. Defaults to `FALSE`.
#' @param ... Passed to `knitr::kable()` (e.g. `caption`, `col.names`,
#'   `digits`).
#' @return A `knitr_kable` object in HTML format.
#' @export
kable_confidence <- function(df,
                             tier_col = "confidence_tier",
                             colours = c(
                               high = "#d9f2d9", medium = "#ffe5b4",
                               low = "#f8d0d0", single_trait = "#e6e6e6"
                             ),
                             drop_tier_col = FALSE,
                             ...) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required for kable_confidence()", call. = FALSE)
  }
  if (!tier_col %in% names(df)) {
    stop("tier_col not found in df: ", tier_col)
  }
  tiers <- as.character(df[[tier_col]])
  if (drop_tier_col) {
    df <- df[, names(df) != tier_col, drop = FALSE]
  }
  html <- knitr::kable(df, format = "html", escape = FALSE, row.names = FALSE, ...)
  lines <- strsplit(as.character(html), "\n", fixed = TRUE)[[1]]
  body_start <- grep("<tbody>", lines, fixed = TRUE)
  if (length(body_start) == 1 && length(tiers) > 0) {
    row_lines <- which(seq_along(lines) > body_start & grepl("^\\s*<tr>", lines))
    for (k in seq_along(row_lines)) {
      if (k > length(tiers)) {
        break
      }
      colour <- unname(colours[tiers[k]])
      if (length(colour) == 1 && !is.na(colour)) {
        lines[row_lines[k]] <- sub(
          "<tr>",
          paste0("<tr style=\"background-color:", colour, ";\">"),
          lines[row_lines[k]],
          fixed = TRUE
        )
      }
    }
  }
  return(structure(
    paste(lines, collapse = "\n"),
    format = "html",
    class = "knitr_kable"
  ))
}
