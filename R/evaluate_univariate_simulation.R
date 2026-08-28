#' @title Evaluate Clustering Against Simulated Ground Truth
#' @description Score a `run_univariate_clustering()` result against the ground
#' truth from `simulate_trait()`. Reports adjusted Rand index (hand-rolled
#' contingency-table implementation), detected module count vs planted K, and
#' how many unstructured background SNPs were absorbed into reliable modules.
#' @param simulation Result of `simulate_trait()`.
#' @param clustering_result Result of `run_univariate_clustering()`. Reliable
#'   modules (`clusters_reliable`) are treated as the inferred partition; SNPs
#'   outside them are labelled `"unassigned"`.
#' @param predicted_memberships Optional dataframe of predicted SNP-program
#'   memberships, typically `summarise_ebmf_programs()$assigned` restricted to
#'   valid programs. This preserves overlapping EBMF memberships for coverage,
#'   background absorption, and module recall. ARI uses each SNP's strongest
#'   assigned program as a hard projection.
#' @return A list with:
#'   \itemize{
#'     \item k_planted: number of modules in the simulation (0 for null)
#'     \item k_hat: number of reliable modules found
#'     \item k_error: absolute difference `|k_hat - k_planted|`
#'     \item ari_all: ARI over all SNPs, background as its own truth class and
#'       unassigned SNPs as their own predicted class
#'     \item ari_structured: ARI over module SNPs only (background excluded)
#'     \item background_absorbed: fraction of background SNPs placed into
#'       reliable modules (higher = more hallucinated structure)
#'     \item coverage: fraction of all SNPs placed into reliable modules
#'     \item module_recall: named numeric vector; for each planted module, the
#'       largest fraction of its SNPs contained in any single reliable predicted
#'       module. Unlike ARI this does not require a one-to-one partition match,
#'       so it scores overlapping/merged programs (e.g. from matrix
#'       factorisation) fairly.
#'     \item mean_module_recall: mean of `module_recall` across planted modules
#'     \item confusion: truth-by-predicted contingency table over all SNPs
#'   }
#' @export
evaluate_univariate_simulation <- function(simulation, clustering_result,
                                           predicted_memberships = NULL) {
  if (is.null(simulation$ground_truth$module_of_snp)) {
    stop("simulation must be the result of simulate_trait()")
  }
  if (is.null(predicted_memberships) && is.null(clustering_result$module_quality)) {
    stop("clustering_result must be the result of run_univariate_clustering()")
  }

  truth <- simulation$ground_truth$module_of_snp
  pred_sets <- list()
  pred_ids <- character(0)

  truth_lab <- ifelse(truth == 0L, "background", paste0("mod", truth))
  pred_lab <- rep("unassigned", length(truth_lab))
  names(pred_lab) <- names(truth_lab)

  if (is.null(predicted_memberships)) {
    pred <- clustering_result$clusters_reliable
    common <- intersect(names(truth_lab), names(pred))
    pred_lab[common] <- as.character(pred[common])
    pred_ids <- unique(common)
    pred_sets <- split(common, as.character(pred[common]))
    reliable_ids <- clustering_result$module_quality$cluster[
      clustering_result$module_quality$reliable
    ]
  } else {
    required <- c("snp_id", "program")
    if (!is.data.frame(predicted_memberships) ||
        !all(required %in% names(predicted_memberships))) {
      stop("predicted_memberships must contain snp_id and program columns")
    }
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
  }

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

  return(list(
    k_planted = k_planted,
    k_hat = k_hat,
    k_error = abs(k_hat - k_planted),
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
