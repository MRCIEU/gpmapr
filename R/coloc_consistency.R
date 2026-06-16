#' @title Pleiotropic Consistency Score (PCS)
#' @description Computes a Jaccard-style consistency score between two traits based on
#' colocalisation evidence at shared loci (LD blocks).
#'
#' PCS = n / (n + m1 + m2)
#' where:
#' - n = number of loci where both traits colocalize (share a coloc group)
#' - m1 = loci where only trait 1 has a causal signal
#' - m2 = loci where only trait 2 has a causal signal
#'
#' @param trait_id_1 Numeric trait ID for the first trait
#' @param trait_id_2 Numeric trait ID for the second trait
#' @param coloc_groups A dataframe of coloc_groups (e.g. from traits() or trait()).
#'   If NULL, will fetch via traits(c(trait_id_1, trait_id_2)).
#' @param p_threshold P-value threshold for considering a trait to have a signal at a locus.
#'   Defaults to 5e-8.
#' @param include_rare Logical, whether to also count rare result loci. Defaults to FALSE.
#' @return A list with:
#'   \itemize{
#'     \item pcs: the Pleiotropic Consistency Score (0 to 1)
#'     \item n_shared: number of shared loci
#'     \item m1: loci unique to trait 1
#'     \item m2: loci unique to trait 2
#'     \item trait_1_loci: total loci for trait 1
#'     \item trait_2_loci: total loci for trait 2
#'   }
#' @export
pleiotropic_consistency_score <- function(trait_id_1,
                                          trait_id_2,
                                          coloc_groups = NULL,
                                          p_threshold = 5e-8,
                                          include_rare = FALSE) {
  if (is.null(coloc_groups)) {
    result <- traits(c(trait_id_1, trait_id_2), include_associations = FALSE)
    coloc_groups <- result$coloc_groups
  }

  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    return(list(pcs = 0, n_shared = 0L, m1 = 0L, m2 = 0L,
                trait_1_loci = 0L, trait_2_loci = 0L))
  }

  cg <- coloc_groups |>
    dplyr::filter(min_p <= p_threshold)

  t1_loci <- cg |>
    dplyr::filter(trait_id == trait_id_1) |>
    dplyr::distinct(coloc_group_id) |>
    dplyr::pull(coloc_group_id)

  t2_loci <- cg |>
    dplyr::filter(trait_id == trait_id_2) |>
    dplyr::distinct(coloc_group_id) |>
    dplyr::pull(coloc_group_id)

  n_shared <- length(intersect(t1_loci, t2_loci))
  m1 <- length(setdiff(t1_loci, t2_loci))
  m2 <- length(setdiff(t2_loci, t1_loci))

  denom <- n_shared + m1 + m2
  pcs <- if (denom == 0) 0 else n_shared / denom

  return(list(
    pcs = pcs,
    n_shared = n_shared,
    m1 = m1,
    m2 = m2,
    trait_1_loci = length(t1_loci),
    trait_2_loci = length(t2_loci)
  ))
}


#' @title Pairwise Pleiotropic Consistency Score Matrix
#' @description Compute PCS for all pairs of traits in a set.
#' @param trait_ids A numeric vector of trait IDs
#' @param coloc_groups A dataframe of coloc_groups covering all requested traits.
#'   If NULL, will fetch via traits(trait_ids).
#' @param p_threshold P-value threshold. Defaults to 5e-8.
#' @return A symmetric matrix of PCS values with trait IDs as row/column names.
#' @export
pcs_matrix <- function(trait_ids, coloc_groups = NULL, p_threshold = 5e-8) {
  if (length(trait_ids) < 2) stop("At least 2 trait_ids required")

  if (is.null(coloc_groups)) {
    if (length(trait_ids) > 10) {
      coloc_groups <- dplyr::bind_rows(lapply(trait_ids, function(tid) {
        tryCatch(trait(tid)$coloc_groups, error = function(e) return(NULL))
      }))
    } else {
      coloc_groups <- traits(trait_ids)$coloc_groups
    }
  }

  n <- length(trait_ids)
  mat <- matrix(0, nrow = n, ncol = n)
  rownames(mat) <- colnames(mat) <- as.character(trait_ids)
  diag(mat) <- 1

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      res <- pleiotropic_consistency_score(
        trait_ids[i], trait_ids[j],
        coloc_groups = coloc_groups,
        p_threshold = p_threshold
      )
      mat[i, j] <- res$pcs
      mat[j, i] <- res$pcs
    }
  }

  return(mat)
}


#' @title Net Directionality Score
#' @description Computes the net directionality of effect between two traits at shared
#' colocalised loci. Uses association betas to determine concordance.
#'
#' NDS = (n_pos - n_neg) / (n_pos + n_neg)
#' where:
#' - n_pos = number of shared loci with concordant effect directions
#' - n_neg = number of shared loci with discordant effect directions
#'
#' @param trait_id_1 Numeric trait ID for the first trait
#' @param trait_id_2 Numeric trait ID for the second trait
#' @param coloc_groups A dataframe of coloc_groups with beta columns (from include_associations = TRUE).
#'   If NULL, will fetch via traits(c(trait_id_1, trait_id_2), include_associations = TRUE).
#' @param p_threshold P-value threshold. Defaults to 5e-8.
#' @return A list with:
#'   \itemize{
#'     \item nds: the Net Directionality Score (-1 to 1)
#'     \item n_pos: number of concordant loci
#'     \item n_neg: number of discordant loci
#'     \item n_total: total shared loci with available betas
#'   }
#' @export
net_directionality_score <- function(trait_id_1,
                                     trait_id_2,
                                     coloc_groups = NULL,
                                     p_threshold = 5e-8) {
  if (is.null(coloc_groups)) {
    result <- traits(c(trait_id_1, trait_id_2), include_associations = TRUE)
    coloc_groups <- result$coloc_groups
  }

  if (is.null(coloc_groups) || nrow(coloc_groups) == 0 || !"beta" %in% names(coloc_groups)) {
    return(list(nds = NA_real_, n_pos = 0L, n_neg = 0L, n_total = 0L))
  }

  cg <- coloc_groups |>
    dplyr::filter(min_p <= p_threshold, !is.na(beta))

  t1 <- cg |>
    dplyr::filter(trait_id == trait_id_1) |>
    dplyr::select("coloc_group_id", beta_1 = "beta")

  t2 <- cg |>
    dplyr::filter(trait_id == trait_id_2) |>
    dplyr::select("coloc_group_id", beta_2 = "beta")

  shared <- dplyr::inner_join(t1, t2, by = "coloc_group_id")

  if (nrow(shared) == 0) {
    return(list(nds = NA_real_, n_pos = 0L, n_neg = 0L, n_total = 0L))
  }

  shared <- shared |>
    dplyr::mutate(concordant = sign(beta_1) == sign(beta_2))

  n_pos <- sum(shared$concordant)
  n_neg <- sum(!shared$concordant)
  n_total <- n_pos + n_neg

  nds <- if (n_total == 0) NA_real_ else (n_pos - n_neg) / n_total

  return(list(
    nds = nds,
    n_pos = n_pos,
    n_neg = n_neg,
    n_total = n_total
  ))
}


#' @title Pairwise Net Directionality Score Matrix
#' @description Compute NDS for all pairs of traits in a set.
#' @param trait_ids A numeric vector of trait IDs
#' @param coloc_groups A dataframe of coloc_groups with betas (include_associations = TRUE).
#'   If NULL, will fetch via traits(trait_ids, include_associations = TRUE).
#' @param p_threshold P-value threshold. Defaults to 5e-8.
#' @return A matrix of NDS values with trait IDs as row/column names.
#'   Diagonal is 1. Values range from -1 (all discordant) to 1 (all concordant).
#' @export
nds_matrix <- function(trait_ids, coloc_groups = NULL, p_threshold = 5e-8) {
  if (length(trait_ids) < 2) stop("At least 2 trait_ids required")

  if (is.null(coloc_groups)) {
    if (length(trait_ids) > 10) {
      coloc_groups <- dplyr::bind_rows(lapply(trait_ids, function(tid) {
        tryCatch(trait(tid, include_associations = TRUE)$coloc_groups, error = function(e) return(NULL))
      }))
    } else {
      coloc_groups <- traits(trait_ids, include_associations = TRUE)$coloc_groups
    }
  }

  n <- length(trait_ids)
  mat <- matrix(NA_real_, nrow = n, ncol = n)
  rownames(mat) <- colnames(mat) <- as.character(trait_ids)
  diag(mat) <- 1

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      res <- net_directionality_score(
        trait_ids[i], trait_ids[j],
        coloc_groups = coloc_groups,
        p_threshold = p_threshold
      )
      mat[i, j] <- res$nds
      mat[j, i] <- res$nds
    }
  }

  return(mat)
}
