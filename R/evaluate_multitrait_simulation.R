#' @title Map Discovered Programs Onto Planted Programs
#' @description For each discovered program, find the planted program of the same
#' trait that it best overlaps on SNP membership. Used to translate EBMF's
#' arbitrary program numbering into the planted labels before cross-trait links
#' can be scored.
#' @param loadings A `extract_program_loadings()` loadings table.
#' @param ground_truth The `ground_truth` component of
#'   `simulate_trait_pair()`.
#' @param min_jaccard Minimum Jaccard overlap for a mapping to be accepted;
#'   below it the program is mapped to `NA`.
#' @return A dataframe with `program_id`, `trait_id`, `planted_label` and
#'   `score` (the best Jaccard overlap).
#' @export
map_programs_to_planted <- function(loadings,
                                    ground_truth,
                                    min_jaccard = 0.1) {
  empty <- data.frame(
    program_id = character(0), trait_id = character(0),
    planted_label = character(0), score = numeric(0),
    stringsAsFactors = FALSE
  )
  if (is.null(loadings) || nrow(loadings) == 0) {
    return(empty)
  }
  planted <- ground_truth$program_snps
  if (is.null(planted) || length(planted) == 0) {
    return(empty)
  }
  if (!"snp_id" %in% names(loadings)) {
    return(empty)
  }
  # Restrict to each program's high-confidence members when they are flagged.
  if ("high_confidence" %in% names(loadings)) {
    hc <- loadings$high_confidence
    hc[is.na(hc)] <- FALSE
    if (any(hc)) {
      loadings <- loadings[hc, , drop = FALSE]
    }
  }
  if (nrow(loadings) == 0) {
    return(empty)
  }

  by_program <- split(seq_len(nrow(loadings)), as.character(loadings$program_id))
  trait_of <- stats::setNames(
    as.character(loadings$trait_id), as.character(loadings$program_id)
  )

  rows <- lapply(names(by_program), function(pid) {
    idx <- by_program[[pid]]
    members <- unique(as.character(loadings$snp_id[idx]))
    tid <- trait_of[[pid]]
    candidates <- names(planted)[startsWith(names(planted), paste0(tid, ":"))]
    none <- data.frame(
      program_id = pid, trait_id = tid, planted_label = NA_character_,
      score = NA_real_, stringsAsFactors = FALSE
    )
    if (length(candidates) == 0 || length(members) == 0) {
      return(none)
    }

    scores <- vapply(candidates, function(cand) {
      truth <- planted[[cand]]
      uni <- length(union(members, truth))
      if (uni == 0) return(0)
      return(length(intersect(members, truth)) / uni)
    }, numeric(1))
    best <- which.max(scores)
    none$score <- scores[[best]]
    if (scores[[best]] >= min_jaccard) {
      none$planted_label <- sub("^[^:]*:", "", candidates[[best]])
    }
    return(none)
  })
  return(dplyr::bind_rows(rows))
}
