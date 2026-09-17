#' @title Evaluate EBMF Programs Against Simulated Ground Truth
#' @description Score a `run_univariate_clustering()` EBMF result against the
#' ground truth from `simulate_trait()`. Reports adjusted Rand index (hand-rolled
#' contingency-table implementation), detected program count vs planted K, and
#' how many unstructured background SNPs were absorbed into predicted programs.
#' @param simulation Result of `simulate_trait()`.
#' @param predicted_memberships Dataframe of predicted SNP-program memberships,
#'   typically `summarise_ebmf_programs()$assigned` restricted to valid
#'   programs. Overlapping EBMF memberships are preserved for coverage,
#'   background absorption, and module recall. ARI uses each SNP's strongest
#'   assigned program as a hard projection.
#' @return A list with:
#'   \itemize{
#'     \item k_planted: number of modules in the simulation (0 for null)
#'     \item k_hat: number of predicted programs
#'     \item k_error: absolute difference `|k_hat - k_planted|`
#'     \item ari_all: ARI over all SNPs, background as its own truth class and
#'       unassigned SNPs as their own predicted class
#'     \item ari_structured: ARI over module SNPs only (background excluded)
#'     \item background_absorbed: fraction of background SNPs placed into
#'       predicted programs (higher = more hallucinated structure). This floors
#'       at 0 whenever the lFSR/magnitude gate keeps background SNPs out, which
#'       it usually does; prefer `background_weight`.
#'     \item background_weight: share of the predicted programs' total loading
#'       mass sitting on background SNPs — the same idea without the floor
#'     \item v_structured: V-measure over module SNPs, the harmonic mean of
#'       homogeneity and completeness. Unlike ARI it distinguishes mixing two
#'       planted modules (a real error) from splitting one (which is correct
#'       when the planted truth is hierarchical)
#'     \item ari_driver_groups / v_driver_groups: the same scores against the
#'       finer planted partition (`ground_truth$driver_group_of_snp`), `NA`
#'       unless the simulation used `snp_driver_groups > 1`
#'     \item coverage: fraction of all SNPs placed into predicted programs
#'     \item module_recall: named numeric vector; for each planted module, the
#'       largest fraction of its SNPs contained in any single predicted program.
#'       Unlike ARI this does not require a one-to-one partition match, so it
#'       scores overlapping/merged programs fairly.
#'     \item mean_module_recall: mean of `module_recall` across planted modules
#'     \item confusion: truth-by-predicted contingency table over all SNPs
#'   }
#' @export
evaluate_univariate_simulation <- function(simulation, predicted_memberships) {
  if (is.null(simulation$ground_truth$module_of_snp)) {
    stop("simulation must be the result of simulate_trait()")
  }
  required <- c("snp_id", "program")
  if (!is.data.frame(predicted_memberships) ||
        !all(required %in% names(predicted_memberships))) {
    stop("predicted_memberships must contain snp_id and program columns")
  }

  truth <- simulation$ground_truth$module_of_snp
  pred_lab <- rep("unassigned", length(truth))
  names(pred_lab) <- names(truth)
  truth_lab <- ifelse(truth == 0L, "background", paste0("mod", truth))

  memberships <- predicted_memberships[
    !is.na(predicted_memberships$snp_id) &
      !is.na(predicted_memberships$program),
    ,
    drop = FALSE
  ]
  pred_ids <- unique(as.character(memberships$snp_id))
  pred_sets <- split(
    as.character(memberships$snp_id),
    as.character(memberships$program)
  )
  pred_sets <- lapply(pred_sets, unique)

  if (nrow(memberships) > 0) {
    if ("abs_loading" %in% names(memberships)) {
      loading <- memberships$abs_loading
      loading[is.na(loading)] <- -Inf
      memberships <- memberships[order(memberships$snp_id, -loading), , drop = FALSE]
    }
    winners <- memberships[!duplicated(memberships$snp_id), , drop = FALSE]
    common <- intersect(names(truth_lab), as.character(winners$snp_id))
    winner_idx <- match(common, as.character(winners$snp_id))
    pred_lab[common] <- paste0(
      "program", as.character(winners$program[winner_idx])
    )
  }
  reliable_ids <- sort(unique(as.character(memberships$program)))

  structured <- truth_lab != "background"

  recall_sets <- pred_sets
  module_memberships <- simulation$ground_truth$module_memberships
  if (is.null(module_memberships)) {
    module_ids <- sort(unique(truth[truth > 0]))
    module_memberships <- lapply(module_ids, function(m) {
      names(truth)[truth == m]
    })
  }
  module_ids <- seq_along(module_memberships)
  module_recall <- vapply(module_ids, function(m) {
    members <- module_memberships[[m]]
    if (length(recall_sets) == 0) {
      return(0)
    }
    max(vapply(recall_sets, function(p) {
      length(intersect(members, p)) / length(members)
    }, numeric(1)))
  }, numeric(1))
  module_recall <- stats::setNames(module_recall, sprintf("mod%d", module_ids))

  k_planted <- length(unique(truth[truth > 0]))
  k_hat <- length(reliable_ids)

  # ARI penalises splitting a planted module as heavily as merging two, but
  # when the generative model plants sub-structure (snp_driver_groups > 1) a
  # split is the method finding what is there. v_measure separates the two:
  # homogeneity falls only when a program mixes modules, completeness only when
  # a module is spread across programs.
  v_structured <- .v_measure(truth_lab[structured], pred_lab[structured])

  # Scored against the finer planted partition as well, when there is one.
  driver_truth <- simulation$ground_truth$driver_group_of_snp
  ari_driver_groups <- NA_real_
  v_driver_groups <- NA_real_
  if (!is.null(driver_truth)) {
    dg_lab <- ifelse(
      driver_truth == 0L, "background", paste0("grp", driver_truth)
    )
    names(dg_lab) <- names(driver_truth)
    common_dg <- intersect(names(dg_lab), names(pred_lab))
    dg_structured <- dg_lab[common_dg] != "background"
    if (any(dg_structured)) {
      ari_driver_groups <- .adjusted_rand_index(
        dg_lab[common_dg][dg_structured], pred_lab[common_dg][dg_structured]
      )
      v_driver_groups <- .v_measure(
        dg_lab[common_dg][dg_structured], pred_lab[common_dg][dg_structured]
      )
    }
  }

  # Weight-based background contamination. The membership-based
  # background_absorbed floors at 0 whenever the lFSR/magnitude gate keeps
  # background SNPs out, which it almost always does, so it carries no
  # information; this uses the share of each program's loading mass that sits
  # on background SNPs instead.
  background_weight <- NA_real_
  if ("abs_loading" %in% names(memberships) && any(!structured)) {
    bg_ids <- names(truth_lab)[!structured]
    w <- memberships$abs_loading
    w[!is.finite(w) | w < 0] <- 0
    total_w <- sum(w)
    background_weight <- if (total_w > 0) {
      sum(w[as.character(memberships$snp_id) %in% bg_ids]) / total_w
    } else {
      NA_real_
    }
  }

  return(list(
    k_planted = k_planted,
    k_hat = k_hat,
    k_error = abs(k_hat - k_planted),
    v_structured = v_structured,
    ari_driver_groups = ari_driver_groups,
    v_driver_groups = v_driver_groups,
    background_weight = background_weight,
    ari_all = .adjusted_rand_index(truth_lab, pred_lab),
    ari_structured = .adjusted_rand_index(
      truth_lab[structured],
      pred_lab[structured]
    ),
    background_absorbed = if (any(!structured)) {
      mean(names(truth_lab)[!structured] %in% pred_ids)
    } else {
      NA_real_
    },
    coverage = mean(names(truth_lab) %in% pred_ids),
    module_recall = module_recall,
    mean_module_recall = if (length(module_recall)) mean(module_recall) else NA_real_,
    confusion = table(truth = truth_lab, predicted = pred_lab)
  ))
}


# V-measure: harmonic mean of homogeneity (does a predicted program mix planted
# modules?) and completeness (is a planted module spread across programs?).
# Reported alongside ARI because the two failures are not equivalent here: with
# hierarchical planted truth, splitting a module is the method recovering
# sub-structure, while mixing modules is a genuine error.
.v_measure <- function(a, b, beta = 1) {
  t <- table(a, b)
  n <- sum(t)
  if (n < 2) {
    return(NA_real_)
  }
  entropy <- function(counts) {
    p <- counts[counts > 0] / sum(counts)
    return(-sum(p * log(p)))
  }
  h_class <- entropy(rowSums(t))
  h_cluster <- entropy(colSums(t))
  # H(class | cluster)
  h_class_given_cluster <- 0
  for (j in seq_len(ncol(t))) {
    col_total <- sum(t[, j])
    if (col_total > 0) {
      h_class_given_cluster <- h_class_given_cluster +
        (col_total / n) * entropy(t[, j])
    }
  }
  h_cluster_given_class <- 0
  for (i in seq_len(nrow(t))) {
    row_total <- sum(t[i, ])
    if (row_total > 0) {
      h_cluster_given_class <- h_cluster_given_class +
        (row_total / n) * entropy(t[i, ])
    }
  }
  homogeneity <- if (h_class == 0) 1 else 1 - h_class_given_cluster / h_class
  completeness <- if (h_cluster == 0) {
    1
  } else {
    1 - h_cluster_given_class / h_cluster
  }
  if (homogeneity + completeness == 0) {
    return(0)
  }
  return(
    (1 + beta) * homogeneity * completeness /
      (beta * homogeneity + completeness)
  )
}


.adjusted_rand_index <- function(a, b) {
  t <- table(a, b)
  n <- sum(t)
  if (n < 2) {
    return(NA_real_)
  }
  choose2 <- function(x) x * (x - 1) / 2
  sum_ij <- sum(choose2(t))
  sum_a <- sum(choose2(rowSums(t)))
  sum_b <- sum(choose2(colSums(t)))
  expected <- sum_a * sum_b / choose2(n)
  maximum <- (sum_a + sum_b) / 2
  denom <- maximum - expected
  if (denom == 0) {
    return(as.numeric(sum_ij == expected))
  }
  return((sum_ij - expected) / denom)
}
