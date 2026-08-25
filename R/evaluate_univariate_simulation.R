#' @title Evaluate Clustering Against Simulated Ground Truth
#' @description Score a `run_univariate_clustering()` result against the ground
#' truth from `simulate_trait()`. Reports adjusted Rand index (hand-rolled
#' contingency-table implementation), detected module count vs planted K, and
#' how many unstructured background SNPs were absorbed into reliable modules.
#' @param simulation Result of `simulate_trait()`.
#' @param clustering_result Result of `run_univariate_clustering()`. Reliable
#'   modules (`clusters_reliable`) are treated as the inferred partition; SNPs
#'   outside them are labelled `"unassigned"`.
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
#'     \item confusion: truth-by-predicted contingency table over all SNPs
#'   }
#' @export
evaluate_univariate_simulation <- function(simulation, clustering_result) {
  if (is.null(simulation$ground_truth$module_of_snp)) {
    stop("simulation must be the result of simulate_trait()")
  }
  if (is.null(clustering_result$module_quality)) {
    stop("clustering_result must be the result of run_univariate_clustering()")
  }

  truth <- simulation$ground_truth$module_of_snp
  pred <- clustering_result$clusters_reliable

  truth_lab <- ifelse(truth == 0L, "background", paste0("mod", truth))
  pred_lab <- rep("unassigned", length(truth_lab))
  names(pred_lab) <- names(truth_lab)
  common <- intersect(names(truth_lab), names(pred))
  pred_lab[common] <- as.character(pred[common])

  structured <- truth_lab != "background"

  k_planted <- length(unique(truth[truth > 0]))
  reliable_ids <- clustering_result$module_quality$cluster[
    clustering_result$module_quality$reliable
  ]
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
      mean(pred_lab[!structured] != "unassigned")
    } else {
      NA_real_
    },
    coverage = mean(pred_lab != "unassigned"),
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
