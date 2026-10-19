#' @title Run Univariate SNP Program Discovery
#' @description Run the full univariate analysis pipeline (pleiotropy matrix ->
#' sparse/ubiquitous trait filter -> orientation -> compression -> SNP cosine
#' similarity -> Empirical Bayes Matrix Factorization programs) from a single
#' entry point. The function takes a trait object as returned by
#' `trait(trait_id, include_associations = TRUE)` — or a simulated object with
#' the same shape — so no API calls are made.
#'
#' The defaults are the shared analysis settings in
#' `vignettes/_ebmf_settings.yml`. To run with a settings list, and to follow
#' with [summarise_ebmf_programs()], use [run_ebmf_pipeline()].
#'
#' The simulated `coloc_groups` dataframe must include at least: `coloc_group_id`,
#' `variant_id`, `display_snp`, `chr`, `bp`, `trait_id`, `trait_name`, `min_p`,
#' `beta`, and `se`. Rows where `trait_id == target_trait_id` define the target
#' SNPs (columns, named by `variant_id`); all other rows become background
#' trait profiles. Optional annotation columns (`gene_id`, `gene`, `tissue`,
#' `trait_category`) enrich downstream interpretation functions.
#' @param trait_object A trait result list (from `trait()` or a simulation) with
#'   a `coloc_groups` dataframe.
#' @param target_trait_id Target trait id. Defaults to `trait_object$trait$id`.
#' @param include_trans If `TRUE` (default), rows flagged `cis_trans ==
#'   "trans"` are kept; `FALSE` removes them before matrix construction
#'   (requires a `cis_trans` column; ignored otherwise).
#' @param trait_subset Which trait rows back the trait x SNP matrix:
#'   `"phenotypic"` (default) keeps only phenotype traits (rows whose
#'   `feature_type` is `"phenotypic"` or missing) plus the target trait,
#'   dropping molecular / gene-annotated rows, while `"all"` uses every
#'   available background trait. Molecular traits are much sparser, so
#'   clustering on phenotype traits only can recover structure that sparse
#'   molecular rows dilute. `"phenotypic"` requires `trait_info` with a
#'   `feature_type` column (present for `trait()` and simulated objects).
#' @param min_snp_signals Minimum non-`NA` SNPs per background trait row.
#' @param max_snp_fraction Maximum non-`NA` fraction per background trait row.
#' @param compress_method Effect compression passed to `compress_effect_matrix()`.
#' @param compress_scale Asinh scale passed to `compress_effect_matrix()`.
#' @param ebmf_greedy_Kmax Maximum number of EBMF factors.
#' @param ebmf_lfsr_threshold lFSR threshold for EBMF program membership.
#' @param ebmf_magnitude_threshold Minimum absolute posterior mean loading for
#'   EBMF program membership.
#' @param ebmf_drop_global If `TRUE` (default), drop at most one global
#'   mega-factor (via `identify_ebmf_global_factors()`) before extracting EBMF
#'   programs. This guards against the near-null "one giant cluster" failure
#'   mode.
#' @param ebmf_prior Prior family for EBMF loadings/factors:
#'   `"point_laplace"` (default; heavier tails, often better when true effects
#'   are large relative to noise) or `"point_normal"`.
#' @param ebmf_backfit If `TRUE` (default), cyclically refit all factors after
#'   the greedy phase. `FALSE` is faster but leaves greedy-order artifacts.
#' @param ebmf_missing_se Standard error of an unobserved cell (no colocalising
#'   association), which EBMF fits as an observed 0. See `run_ebmf()`.
#' @return A list with:
#'   \itemize{
#'     \item x_matrix: raw traits x SNPs z-score matrix after trait filtering
#'     \item x_star: oriented, compressed traits x SNPs matrix used for similarity
#'     \item trait_matrix: SNP x trait version of `x_star`
#'     \item s_matrix: SNP-by-SNP cosine similarity matrix, computed without the
#'       target row (which is positive at every SNP after orientation). It is
#'       not used by the EBMF fit; it backs `summarise_ebmf_programs()`.
#'     \item overlap_matrix, eligible_matrix: joint-observation counts / eligibility
#'     \item cluster_membership: SNP x program logical membership matrix from EBMF
#'     \item cluster_details: EBMF fit details (including the `flash` object)
#'     \item beta_matrix, se_matrix: trait x SNP beta / SE matrices when available
#'     \item ebmf_input: the features x SNPs matrix actually passed to flashier
#'     \item trait_info, snp_info: row / column metadata
#'     \item dropped_trait_ids: background traits removed by the feature-type
#'       (`trait_subset`) and/or sparse/ubiquitous filter
#'     \item coloc_groups: coloc groups actually used (after trans filtering)
#'     \item parameters: list of settings used
#'   }
#' @export
run_univariate_clustering <- function(trait_object,
                                      target_trait_id = NULL,
                                      include_trans = TRUE,
                                      trait_subset = c("phenotypic", "all"),
                                      min_snp_signals = 2L,
                                      max_snp_fraction = 0.8,
                                      compress_method = c("asinh", "none"),
                                      compress_scale = 5,
                                      ebmf_greedy_Kmax = 100L,
                                      ebmf_lfsr_threshold = 0.02,
                                      ebmf_magnitude_threshold = 0.5,
                                      ebmf_drop_global = TRUE,
                                      ebmf_prior = c("point_laplace", "point_normal"),
                                      ebmf_backfit = TRUE,
                                      ebmf_missing_se = 1) {
  trait_subset <- match.arg(trait_subset)
  compress_method <- match.arg(compress_method)
  ebmf_prior <- match.arg(ebmf_prior)

  if (is.null(trait_object)) {
    stop("trait_object is required")
  }
  coloc_groups <- trait_object$coloc_groups
  if (!is.data.frame(coloc_groups) || nrow(coloc_groups) == 0) {
    stop("trait_object$coloc_groups must be a non-empty dataframe")
  }

  if (is.null(target_trait_id)) {
    if (!is.null(trait_object$trait) && !is.null(trait_object$trait$id)) {
      target_trait_id <- trait_object$trait$id
    } else {
      stop("target_trait_id is required (could not infer from trait_object$trait$id)")
    }
  }
  target_id <- as.character(target_trait_id)

  if (!include_trans && "cis_trans" %in% names(coloc_groups)) {
    coloc_groups <- coloc_groups |>
      dplyr::filter(is.na(cis_trans) | tolower(cis_trans) != "trans")
    if (nrow(coloc_groups) == 0) {
      stop("No coloc_groups rows left after removing trans markers")
    }
  }

  locus_data <- .prepare_pleiotropy_locus_data(
    trait_id = target_trait_id,
    coloc_groups = coloc_groups
  )
  pleiotropy <- .finalize_pleiotropy_from_locus(locus_data, target_id)

  X <- pleiotropy$x_matrix
  if (!target_id %in% rownames(X)) {
    stop("target trait row missing from the built pleiotropy matrix")
  }

  subset_dropped <- character(0)
  if (trait_subset == "phenotypic") {
    info <- pleiotropy$trait_info
    if (is.null(info) || !"feature_type" %in% names(info)) {
      stop("trait_subset = 'phenotypic' requires trait_info with a feature_type column")
    }
    keep_ids <- as.character(info$trait_id[
      as.character(info$trait_id) == target_id |
        is.na(info$feature_type) |
        info$feature_type == "phenotypic"
    ])
    subset_dropped <- setdiff(rownames(X), keep_ids)
    if (length(subset_dropped) > 0) {
      X <- X[!rownames(X) %in% subset_dropped, , drop = FALSE]
      for (mat in c("beta_matrix", "se_matrix")) {
        if (!is.null(pleiotropy[[mat]])) {
          pleiotropy[[mat]] <- pleiotropy[[mat]][
            !rownames(pleiotropy[[mat]]) %in% subset_dropped, , drop = FALSE
          ]
        }
      }
      pleiotropy$trait_info <- info |>
        dplyr::filter(!as.character(trait_id) %in% subset_dropped)
    }
  }

  trait_snp_counts <- rowSums(!is.na(X))
  trait_snp_frac <- trait_snp_counts / ncol(X)
  sparse_trait_ids <- names(trait_snp_counts)[trait_snp_counts < min_snp_signals]
  ubiquitous_trait_ids <- names(trait_snp_frac)[trait_snp_frac > max_snp_fraction]
  drop_trait_ids <- setdiff(
    union(sparse_trait_ids, ubiquitous_trait_ids),
    target_id
  )
  if (length(drop_trait_ids) > 0) {
    X <- X[!rownames(X) %in% drop_trait_ids, , drop = FALSE]
    for (mat in c("beta_matrix", "se_matrix")) {
      if (!is.null(pleiotropy[[mat]])) {
        pleiotropy[[mat]] <- pleiotropy[[mat]][
          !rownames(pleiotropy[[mat]]) %in% drop_trait_ids, , drop = FALSE
        ]
      }
    }
    pleiotropy$trait_info <- pleiotropy$trait_info |>
      dplyr::filter(!as.character(trait_id) %in% drop_trait_ids)
  }

  oriented <- orient_pleiotropy_matrix(X, target_trait_id = target_trait_id)
  X_star <- compress_effect_matrix(
    oriented$x_matrix,
    method = compress_method,
    asinh_scale = compress_scale
  )
  trait_matrix <- t(X_star)

  # S is not used by the EBMF fit; it backs program validation and reporting.
  # The target row is left out: after orientation it is positive at every SNP,
  # so it adds the same similarity to every pair and hides the structure the
  # validation is looking for. SNPs observed only for the target then have an
  # empty profile, which is expected, so the zero-norm warning is muffled.
  similarity <- withCallingHandlers(
    snp_similarity_matrix(X_star[rownames(X_star) != target_id, , drop = FALSE]),
    warning = function(w) {
      if (grepl("zero norm", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  ebmf_x_input <- X_star

  ebmf_fit <- .cluster_snp_profiles_ebmf(
    ebmf_x_input,
    greedy_Kmax = ebmf_greedy_Kmax,
    lfsr_threshold = ebmf_lfsr_threshold,
    magnitude_threshold = ebmf_magnitude_threshold,
    drop_global = ebmf_drop_global,
    prior = ebmf_prior,
    backfit = ebmf_backfit,
    missing_se = ebmf_missing_se
  )

  return(list(
    x_matrix = X,
    x_star = X_star,
    trait_matrix = trait_matrix,
    s_matrix = similarity$s_matrix,
    overlap_matrix = similarity$overlap_matrix,
    eligible_matrix = similarity$eligible_matrix,
    cluster_membership = ebmf_fit$details$membership,
    cluster_details = ebmf_fit$details,
    beta_matrix = if (!is.null(pleiotropy$beta_matrix)) {
      pleiotropy$beta_matrix[rownames(X), , drop = FALSE]
    } else {
      NULL
    },
    se_matrix = if (!is.null(pleiotropy$se_matrix)) {
      pleiotropy$se_matrix[rownames(X), , drop = FALSE]
    } else {
      NULL
    },
    ebmf_input = ebmf_x_input,
    trait_info = pleiotropy$trait_info,
    snp_info = pleiotropy$snp_info,
    dropped_trait_ids = union(drop_trait_ids, subset_dropped),
    coloc_groups = coloc_groups,
    parameters = list(
      method = "ebmf",
      target_trait_id = target_id,
      include_trans = include_trans,
      trait_subset = trait_subset,
      min_snp_signals = min_snp_signals,
      max_snp_fraction = max_snp_fraction,
      compress_method = compress_method,
      compress_scale = compress_scale,
      ebmf_greedy_Kmax = ebmf_greedy_Kmax,
      ebmf_lfsr_threshold = ebmf_lfsr_threshold,
      ebmf_magnitude_threshold = ebmf_magnitude_threshold,
      ebmf_drop_global = ebmf_drop_global,
      ebmf_prior = ebmf_prior,
      ebmf_backfit = ebmf_backfit,
      ebmf_missing_se = ebmf_missing_se
    )
  ))
}


#' @title Run the Univariate EBMF Pipeline From One Settings List
#' @description Run [run_univariate_clustering()] and then
#' [summarise_ebmf_programs()] with one list of analysis settings, keyed as in
#' `vignettes/_ebmf_settings.yml`, typically from [resolve_ebmf_settings()].
#' Settings not supplied take the two functions' defaults, which match that
#' file.
#'
#' Keys passed to [run_univariate_clustering()] under the same name:
#' `include_trans`, `trait_subset`, `min_snp_signals`, `max_snp_fraction`,
#' `compress_method`, `compress_scale`, `ebmf_greedy_Kmax`,
#' `ebmf_lfsr_threshold`, `ebmf_magnitude_threshold`, `ebmf_drop_global`,
#' `ebmf_prior`, `ebmf_backfit`, `ebmf_missing_se`. Keys passed to
#' [summarise_ebmf_programs()]: `min_module_size`, `similarity_n_perm`,
#' `high_similarity_q`, `min_effective_traits` and `frac_traits` under the same
#' name, and `similarity_threshold` (as `edge_threshold`), `similarity_q`
#' (`similarity_q_threshold`), `high_connectedness`
#' (`high_min_connectedness`), `high_similarity_z` (`high_min_similarity_z`)
#' and `stability_top_n` (`top_n`). Any other key is an error.
#' @param trait_object A trait result list (from `trait()` or a simulation)
#'   with a `coloc_groups` dataframe.
#' @param settings Named list of analysis settings (see Description). `NULL`
#'   entries are ignored.
#' @param target_trait_id Target trait id. Defaults to `trait_object$trait$id`.
#' @param n_rep,seed,cores,verbose Passed to [summarise_ebmf_programs()].
#' @return A list with `fit` (the [run_univariate_clustering()] result) and
#'   `summary` (the [summarise_ebmf_programs()] result).
#' @export
run_ebmf_pipeline <- function(trait_object,
                              settings = list(),
                              target_trait_id = NULL,
                              n_rep = 10,
                              seed = 1,
                              cores = 1,
                              verbose = TRUE) {
  args <- .ebmf_pipeline_args(settings)
  fit <- do.call(
    run_univariate_clustering,
    c(
      list(trait_object = trait_object, target_trait_id = target_trait_id),
      args$clustering
    )
  )
  summary <- do.call(
    summarise_ebmf_programs,
    c(
      list(
        clustering_result = fit, n_rep = n_rep, seed = seed, cores = cores,
        verbose = verbose
      ),
      args$validation
    )
  )
  return(list(fit = fit, summary = summary))
}


# Settings keys (as in vignettes/_ebmf_settings.yml) and the argument each one
# sets, for run_univariate_clustering() and summarise_ebmf_programs().
.ebmf_settings_map <- function() {
  clustering <- c(
    "include_trans", "trait_subset", "min_snp_signals", "max_snp_fraction",
    "compress_method", "compress_scale", "ebmf_greedy_Kmax",
    "ebmf_lfsr_threshold", "ebmf_magnitude_threshold", "ebmf_drop_global",
    "ebmf_prior", "ebmf_backfit", "ebmf_missing_se"
  )
  return(list(
    clustering = stats::setNames(clustering, clustering),
    validation = c(
      min_module_size = "min_module_size",
      similarity_threshold = "edge_threshold",
      similarity_n_perm = "similarity_n_perm",
      similarity_q = "similarity_q_threshold",
      high_similarity_q = "high_similarity_q",
      high_connectedness = "high_min_connectedness",
      high_similarity_z = "high_min_similarity_z",
      min_effective_traits = "min_effective_traits",
      frac_traits = "frac_traits",
      stability_top_n = "top_n"
    )
  ))
}


# Split a settings list into the argument lists of the two pipeline steps.
.ebmf_pipeline_args <- function(settings) {
  if (is.null(settings)) {
    settings <- list()
  }
  if (!is.list(settings) || (length(settings) > 0 && is.null(names(settings)))) {
    stop("settings must be a named list")
  }
  settings <- settings[!vapply(settings, is.null, logical(1))]
  map <- .ebmf_settings_map()
  unknown <- setdiff(names(settings), c(names(map$clustering), names(map$validation)))
  if (length(unknown) > 0) {
    stop("Unknown EBMF setting(s): ", paste(unknown, collapse = ", "))
  }
  pick <- function(m) {
    keys <- intersect(names(m), names(settings))
    return(stats::setNames(settings[keys], unname(m[keys])))
  }
  return(list(clustering = pick(map$clustering), validation = pick(map$validation)))
}


.assert_ebmf_result <- function(clustering_result) {
  params <- clustering_result$parameters
  if (is.null(params) || !identical(params$method, "ebmf")) {
    stop("clustering_result must come from run_univariate_clustering()")
  }
  return(invisible(TRUE))
}
