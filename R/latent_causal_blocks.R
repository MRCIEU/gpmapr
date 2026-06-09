#' @title Build Locus Z-Score Matrix
#' @description Constructs a loci x traits matrix of z-scores from colocalisation results,
#' suitable for Latent Profile Analysis. Each row is a locus (coloc_group_id),
#' each column is a trait, and values are z-scores (beta/se).
#'
#' Loci where a trait has no signal are set to 0 (no evidence of effect),
#' preserving directionality for traits with signals.
#'
#' @param trait_ids A numeric vector of trait IDs
#' @param coloc_groups A dataframe of coloc_groups with beta and se columns
#'   (from include_associations = TRUE). If NULL, will fetch via traits().
#' @param p_threshold P-value threshold for including a locus. Defaults to 5e-8.
#' @param min_traits Minimum number of traits that must have a signal at a locus
#'   for it to be included. Defaults to 2.
#' @return A list with:
#'   \itemize{
#'     \item z_matrix: a numeric matrix (loci x traits) of z-scores
#'     \item locus_info: a dataframe mapping row indices to coloc_group_id, chr, bp, display_snp
#'     \item trait_info: a dataframe mapping column indices to trait_id, trait_name
#'   }
#' @export
build_locus_zscore_matrix <- function(trait_ids,
                                      coloc_groups = NULL,
                                      p_threshold = 5e-8,
                                      min_traits = 2) {
  if (length(trait_ids) < 2) stop("At least 2 trait_ids required")

  if (is.null(coloc_groups)) {
    if (length(trait_ids) > 10) {
      coloc_groups <- dplyr::bind_rows(lapply(trait_ids, function(tid) {
        tryCatch(trait(tid, include_associations = TRUE)$coloc_groups, error = function(e) NULL)
      }))
    } else {
      coloc_groups <- traits(trait_ids, include_associations = TRUE)$coloc_groups
    }
  }

  if (is.null(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("No coloc_groups data available")
  }

  if (!all(c("beta", "se") %in% names(coloc_groups))) {
    stop("coloc_groups must include beta and se columns (use include_associations = TRUE)")
  }

  cg <- coloc_groups |>
    dplyr::filter(
      .data$trait_id %in% trait_ids,
      .data$min_p <= p_threshold,
      !is.na(.data$beta),
      !is.na(.data$se),
      .data$se > 0
    ) |>
    dplyr::mutate(z = .data$beta / .data$se)

  # For each locus-trait pair, take the z-score with the smallest p-value
  cg_summary <- cg |>
    dplyr::group_by(.data$coloc_group_id, .data$trait_id) |>
    dplyr::slice_min(.data$min_p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # Filter loci present in at least min_traits
  locus_trait_counts <- cg_summary |>
    dplyr::group_by(.data$coloc_group_id) |>
    dplyr::summarise(n_traits = dplyr::n_distinct(.data$trait_id), .groups = "drop") |>
    dplyr::filter(.data$n_traits >= min_traits)

  cg_summary <- cg_summary |>
    dplyr::filter(.data$coloc_group_id %in% locus_trait_counts$coloc_group_id)

  if (nrow(cg_summary) == 0) {
    stop("No loci with signals in >= min_traits traits after filtering")
  }

  # Build wide matrix
  z_wide <- cg_summary |>
    dplyr::select("coloc_group_id", "trait_id", "z") |>
    tidyr::pivot_wider(
      names_from = "trait_id",
      values_from = "z",
      values_fill = 0
    )

  locus_ids <- z_wide$coloc_group_id
  z_matrix <- as.matrix(z_wide[, -1])
  rownames(z_matrix) <- as.character(locus_ids)

  # Ensure all requested traits are columns (fill missing with 0)
  missing_traits <- setdiff(as.character(trait_ids), colnames(z_matrix))
  if (length(missing_traits) > 0) {
    zero_cols <- matrix(0, nrow = nrow(z_matrix), ncol = length(missing_traits))
    colnames(zero_cols) <- missing_traits
    z_matrix <- cbind(z_matrix, zero_cols)
  }
  z_matrix <- z_matrix[, as.character(trait_ids), drop = FALSE]

  # Build locus metadata
  locus_info <- cg_summary |>
    dplyr::distinct(.data$coloc_group_id, .keep_all = TRUE) |>
    dplyr::select(dplyr::any_of(c("coloc_group_id", "chr", "bp", "display_snp", "ld_block"))) |>
    dplyr::filter(.data$coloc_group_id %in% locus_ids)

  # Build trait metadata
  trait_info <- cg_summary |>
    dplyr::distinct(.data$trait_id, .data$trait_name)

  list(
    z_matrix = z_matrix,
    locus_info = locus_info,
    trait_info = trait_info
  )
}


#' @title Latent Profile Analysis on Colocalised Loci
#' @description Runs Latent Profile Analysis (LPA) on a z-score matrix of colocalised loci
#' to identify latent causal blocks — groups of loci with similar multi-trait effect profiles.
#'
#' Uses the tidyLPA package (must be installed).
#'
#' @param z_matrix A numeric matrix (loci x traits) of z-scores, as returned by
#'   build_locus_zscore_matrix()$z_matrix
#' @param n_profiles Integer or vector of integers specifying the number of profiles to fit.
#'   Defaults to 2:6.
#' @param models Integer vector of tidyLPA model types to try. Defaults to c(1, 2, 3).
#'   Model 1: equal variances, covariances fixed to 0.
#'   Model 2: varying variances, covariances fixed to 0.
#'   Model 3: equal variances, equal covariances.
#' @return A list with:
#'   \itemize{
#'     \item fit: the tidyLPA fit object (for model comparison)
#'     \item best_model: the best-fitting model summary
#'     \item assignments: a dataframe with coloc_group_id and assigned profile class
#'   }
#' @export
latent_causal_blocks <- function(z_matrix,
                                 n_profiles = 2:6,
                                 models = c(1, 2, 3)) {
  if (!requireNamespace("tidyLPA", quietly = TRUE)) {
    stop("Package 'tidyLPA' is required. Install with: install.packages('tidyLPA')")
  }

  z_df <- as.data.frame(z_matrix)
  z_df$locus_id <- rownames(z_matrix)

  fit <- tidyLPA::estimate_profiles(
    z_df[, colnames(z_matrix), drop = FALSE],
    n_profiles = n_profiles,
    models = models
  )

  best <- tidyLPA::get_fit(fit)
  best_idx <- which.min(best$BIC)
  best_model_info <- best[best_idx, ]

  best_fit <- tidyLPA::get_data(fit)
  # get_data returns a tibble with Class column for the best-fitting model
  # Filter to the best model
  best_fit_filtered <- best_fit |>
    dplyr::filter(
      .data$classes_number == best_model_info$Classes,
      .data$model_number == best_model_info$Model
    )

  assignments <- data.frame(
    locus_id = z_df$locus_id,
    profile = best_fit_filtered$Class,
    stringsAsFactors = FALSE
  )

  list(
    fit = fit,
    best_model = best_model_info,
    assignments = assignments
  )
}


#' @title Summarise Latent Causal Blocks
#' @description Summarise the biological profile of each latent class by computing
#' mean z-scores per trait within each profile.
#' @param z_matrix A numeric matrix (loci x traits) of z-scores
#' @param assignments A dataframe with locus_id and profile columns
#'   (as returned by latent_causal_blocks()$assignments)
#' @param trait_info Optional dataframe with trait_id and trait_name columns
#'   (as returned by build_locus_zscore_matrix()$trait_info)
#' @return A dataframe with columns: profile, trait, mean_z, sd_z, n_loci
#' @export
summarise_latent_blocks <- function(z_matrix, assignments, trait_info = NULL) {
  z_df <- as.data.frame(z_matrix)
  z_df$locus_id <- rownames(z_matrix)

  merged <- dplyr::left_join(z_df, assignments, by = "locus_id")

  long <- merged |>
    tidyr::pivot_longer(
      cols = colnames(z_matrix),
      names_to = "trait_id",
      values_to = "z"
    ) |>
    dplyr::group_by(.data$profile, .data$trait_id) |>
    dplyr::summarise(
      mean_z = mean(.data$z, na.rm = TRUE),
      sd_z = stats::sd(.data$z, na.rm = TRUE),
      n_loci = dplyr::n(),
      .groups = "drop"
    )

  if (!is.null(trait_info)) {
    long <- long |>
      dplyr::left_join(
        trait_info |> dplyr::mutate(trait_id = as.character(.data$trait_id)),
        by = "trait_id"
      )
  }

  long
}
