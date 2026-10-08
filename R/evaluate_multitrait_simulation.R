#' @title Map Discovered Programs Onto Planted Programs
#' @description For each discovered program, find the planted program of the same
#' trait that it best overlaps. Used to translate EBMF's arbitrary program
#' numbering into the planted labels before cross-trait links can be scored.
#'
#' Mapping is done on SNP membership by default. Set `by = "drivers"` to map on
#' the program's loading over background studies instead: a program whose loci
#' were planted privately still has the planted driver studies loading on it.
#' @param loadings A `extract_program_loadings()` loadings table (SNP axis) or
#'   profiles table (feature axis), depending on `by`.
#' @param ground_truth The `ground_truth` component of
#'   `simulate_trait_pair()`.
#' @param by `"snps"` (default) or `"drivers"`.
#' @param min_jaccard Minimum Jaccard overlap for a `by = "snps"` mapping to be
#'   accepted; below it the program is mapped to `NA`.
#' @param min_fold Minimum enrichment over chance for a `by = "drivers"`
#'   mapping. The feature table is dense -- every background study carries some
#'   loading on every program -- so a set overlap is the wrong measure there: the
#'   eight planted drivers of a program are only eight of several hundred
#'   features, and a Jaccard would reject a perfect match. Instead the score is
#'   the share of the program's total squared feature loading falling on the
#'   candidate driver set, divided by the share expected if loading were spread
#'   evenly. Defaults to `3`.
#' @return A dataframe with `program_id`, `trait_id`, `planted_label`, `score`
#'   (Jaccard for `"snps"`, weighted loading share for `"drivers"`) and `fold`
#'   (`NA` for `"snps"`).
#' @export
map_programs_to_planted <- function(loadings,
                                    ground_truth,
                                    by = c("snps", "drivers"),
                                    min_jaccard = 0.1,
                                    min_fold = 3) {
  by <- match.arg(by)
  empty <- data.frame(
    program_id = character(0), trait_id = character(0),
    planted_label = character(0), score = numeric(0), fold = numeric(0),
    stringsAsFactors = FALSE
  )
  if (is.null(loadings) || nrow(loadings) == 0) {
    return(empty)
  }
  planted <- if (by == "snps") {
    ground_truth$program_snps
  } else {
    ground_truth$program_drivers
  }
  if (is.null(planted) || length(planted) == 0) {
    return(empty)
  }

  id_col <- if (by == "snps") "snp_id" else "feature_trait_id"
  if (!id_col %in% names(loadings)) {
    return(empty)
  }
  # The SNP table can be restricted to a program's high-confidence members; the
  # feature table is dense by nature and is handled by weighting instead.
  if (by == "snps" && "high_confidence" %in% names(loadings)) {
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
    members <- as.character(loadings[[id_col]][idx])
    tid <- trait_of[[pid]]
    candidates <- names(planted)[startsWith(names(planted), paste0(tid, ":"))]
    none <- data.frame(
      program_id = pid, trait_id = tid, planted_label = NA_character_,
      score = NA_real_, fold = NA_real_, stringsAsFactors = FALSE
    )
    if (length(candidates) == 0 || length(members) == 0) {
      return(none)
    }

    if (by == "snps") {
      uniq <- unique(members)
      scores <- vapply(candidates, function(cand) {
        truth <- planted[[cand]]
        uni <- length(union(uniq, truth))
        if (uni == 0) return(0)
        return(length(intersect(uniq, truth)) / uni)
      }, numeric(1))
      best <- which.max(scores)
      none$score <- scores[[best]]
      if (scores[[best]] >= min_jaccard) {
        none$planted_label <- sub("^[^:]*:", "", candidates[[best]])
      }
      return(none)
    }

    # Feature axis: share of the program's squared feature loading landing on
    # each candidate driver set, against the evenly-spread expectation.
    w <- loadings$loading[idx]^2
    w[!is.finite(w) | w < 0] <- 0
    total <- sum(w)
    if (total <= 0) {
      return(none)
    }
    n_features <- length(unique(members))
    stats_by_cand <- vapply(candidates, function(cand) {
      truth <- planted[[cand]]
      share <- sum(w[members %in% truth]) / total
      expected <- length(intersect(unique(members), truth)) / n_features
      if (expected <= 0) return(c(share, 0))
      return(c(share, share / expected))
    }, numeric(2))
    best <- which.max(stats_by_cand[2, ])
    none$score <- stats_by_cand[1, best]
    none$fold <- stats_by_cand[2, best]
    if (is.finite(none$fold) && none$fold >= min_fold) {
      none$planted_label <- sub("^[^:]*:", "", candidates[[best]])
    }
    return(none)
  })
  return(dplyr::bind_rows(rows))
}


#' @title Score Cross-Trait Links Against The Planted Correspondence
#' @description Evaluate one cross-trait comparison against
#' `simulate_trait_pair()` ground truth. The estimand is the planted programme
#' correspondence: which program in trait A is the same program as which in trait
#' B, and whether the two traits are pushed the same way or opposite ways.
#'
#' A discovered link is a true positive when both of its programs map to the
#' **same** planted shared label. Recall is over planted shared programs, so
#' recovering one shared program through several discovered links counts once --
#' EBMF factor splitting should not inflate it. Precision is over discovered
#' links, so a method that links everything is penalised.
#'
#' Sign accuracy is scored only on true positives: a link that is not real has no
#' meaningful direction. It is also reported separately for concordant and
#' antagonistic planted programs, because a direction call that says
#' "concordant" for every link scores well overall whenever most planted
#' programs are concordant -- the split is what shows it cannot see antagonism
#' at all. The `direction` column scored is whichever the caller supplies, so
#' the same function scores the concordance sign or [module_rg()]'s
#' `rg_direction`. A link with no direction call (`NA` or `"undetermined"`) is
#' left out of the sign accuracy rather than counted as wrong.
#' @param pairs The `pairs` table from `compare_program_pairs_loadings()`. When
#'   it has a `linked` column only the linked pairs are scored; otherwise every
#'   row is treated as a link.
#' @param mapping Result of `map_programs_to_planted()`.
#' @param ground_truth The `ground_truth` component of `simulate_trait_pair()`.
#' @return A one-row dataframe: `n_links`, `n_true_positive`,
#'   `n_false_positive`, `precision`, `recall`, `f1`, `sign_accuracy`,
#'   `sign_accuracy_concordant`, `sign_accuracy_antagonistic`,
#'   `n_planted_shared`, `n_recovered_shared`.
#' @export
evaluate_multitrait_simulation <- function(pairs,
                                           mapping,
                                           ground_truth) {
  planted <- ground_truth$correspondence
  n_planted <- nrow(planted)
  empty <- data.frame(
    n_links = 0L, n_true_positive = 0L, n_false_positive = 0L,
    precision = NA_real_, recall = if (n_planted > 0) 0 else NA_real_,
    f1 = NA_real_, sign_accuracy = NA_real_,
    sign_accuracy_concordant = NA_real_, sign_accuracy_antagonistic = NA_real_,
    n_planted_shared = n_planted, n_recovered_shared = 0L,
    stringsAsFactors = FALSE
  )
  if (is.null(pairs) || nrow(pairs) == 0) {
    return(empty)
  }
  if ("linked" %in% names(pairs)) {
    pairs <- pairs[!is.na(pairs$linked) & pairs$linked, , drop = FALSE]
  }
  if (nrow(pairs) == 0) {
    return(empty)
  }

  label_of <- stats::setNames(mapping$planted_label, mapping$program_id)
  la <- label_of[as.character(pairs$program_id_a)]
  lb <- label_of[as.character(pairs$program_id_b)]
  shared_labels <- planted$program_label

  is_tp <- !is.na(la) & !is.na(lb) & la == lb & la %in% shared_labels
  n_links <- nrow(pairs)
  n_tp <- sum(is_tp)

  recovered <- unique(la[is_tp])
  direction_of <- stats::setNames(planted$direction, planted$program_label)
  sign_ok <- NA_real_
  sign_conc <- NA_real_
  sign_anta <- NA_real_
  if (n_tp > 0 && "direction" %in% names(pairs)) {
    truth_dir <- unname(direction_of[la[is_tp]])
    called <- pairs$direction[is_tp]
    called[called %in% "undetermined"] <- NA_character_
    sign_rate <- function(keep) {
      if (!any(keep)) {
        return(NA_real_)
      }
      return(mean(called[keep] == truth_dir[keep], na.rm = TRUE))
    }
    sign_ok <- sign_rate(truth_dir %in% c("concordant", "antagonistic"))
    sign_conc <- sign_rate(truth_dir == "concordant")
    sign_anta <- sign_rate(truth_dir == "antagonistic")
  }

  precision <- if (n_links > 0) n_tp / n_links else NA_real_
  recall <- if (n_planted > 0) length(recovered) / n_planted else NA_real_
  f1 <- if (is.finite(precision) && is.finite(recall) &&
              (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }

  return(data.frame(
    n_links = as.integer(n_links),
    n_true_positive = as.integer(n_tp),
    n_false_positive = as.integer(n_links - n_tp),
    precision = precision,
    recall = recall,
    f1 = f1,
    sign_accuracy = sign_ok,
    sign_accuracy_concordant = sign_conc,
    sign_accuracy_antagonistic = sign_anta,
    n_planted_shared = as.integer(n_planted),
    n_recovered_shared = as.integer(length(recovered)),
    stringsAsFactors = FALSE
  ))
}
