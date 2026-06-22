#' @title Build EBMF Matrix with Annotations
#' @description Construct a flat traits x SNPs matrix and annotation vectors
#' for Empirical Bayes Matrix Factorization (EBMF). Rows represent traits (uncollapsed),
#' columns represent variants. Cell values are signed GWAS effect sizes (\code{beta});
#' matching standard errors are returned in \code{se_matrix} for use as the
#' \code{S} argument in \code{run_ebmf()}.
#' Rows and columns that are entirely missing or zero are removed before the
#' matrix is returned.
#'
#' Annotations encode biological metadata for each row and column. Row labels
#' depend on \code{label_scheme}; SNP columns are labelled by coloc group only.
#' \itemize{
#'   \item Trait row: \code{trait_category}
#'   \item Molecular QTL row (\code{pathway_gene_tissue}): \code{PathwayOrGene_Tissue}
#'   \item Molecular QTL row (\code{pathway_gene}): \code{PathwayOrGene}
#'   \item Molecular QTL row (\code{tissue}): \code{Tissue}
#'   \item SNP column: \code{ColocGroup\{ID\}}
#' }
#' @param label_scheme Row and SNP annotation scheme: \code{"pathway_gene_tissue"},
#'   \code{"pathway_gene"}, or \code{"tissue"}. Defaults to \code{"pathway_gene_tissue"}.
#' @inheritParams build_pleiotropy_matrix
#' @param pathway_source Optional pathway source for enrichment annotation.
#' @param pathway_p_value_threshold FDR threshold for pathway enrichment annotation.
#'   Defaults to 0.05.
#' @param minimum_count_in_network Minimum gene overlap per pathway term.
#' @return A list with:
#'   \itemize{
#'     \item beta_matrix: traits x SNPs matrix of signed effect sizes
#'     \item se_matrix: traits x SNPs matrix of standard errors (for \code{S})
#'     \item trait_annotations: named character vector of row annotations
#'     \item snp_annotations: named character vector of column annotations
#'     \item trait_info: dataframe of trait metadata
#'     \item snp_info: dataframe of SNP metadata
#'     \item target_trait_id: target trait ID
#'     \item pathway_enrichment: pathway enrichment results (if computed)
#'   }
#' @export
build_ebmf_matrix <- function(trait_id,
                              coloc_groups = NULL,
                              p_threshold = NULL,
                              snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                              label_scheme = c(
                                "pathway_gene_tissue", "pathway_gene", "tissue"
                              ),
                              pathway_source = NULL,
                              pathway_p_value_threshold = 0.05,
                              minimum_count_in_network = NULL) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }

  snp_key <- match.arg(snp_key)
  label_scheme <- match.arg(label_scheme)

  if (is.null(coloc_groups)) {
    coloc_groups <- trait(trait_id, include_associations = TRUE)$coloc_groups
  }

  assoc <- .build_ebmf_association_matrices(
    trait_id = trait_id,
    coloc_groups = coloc_groups,
    p_threshold = p_threshold,
    snp_key = snp_key
  )

  filtered <- .filter_ebmf_matrices(
    beta_matrix = assoc$beta_matrix,
    se_matrix = assoc$se_matrix
  )

  trait_annot <- .build_ebmf_trait_annotations(
    trait_ids = rownames(filtered$beta_matrix),
    coloc_groups = coloc_groups,
    label_scheme = label_scheme,
    pathway_source = pathway_source,
    pathway_p_value_threshold = pathway_p_value_threshold,
    minimum_count_in_network = minimum_count_in_network
  )

  snp_annot <- .build_ebmf_snp_annotations(
    snp_ids = colnames(filtered$beta_matrix),
    snp_info = assoc$snp_info
  )

  trait_info <- assoc$trait_info |>
    dplyr::filter(as.character(trait_id) %in% rownames(filtered$beta_matrix))

  snp_info <- assoc$snp_info |>
    dplyr::filter(as.character(snp_id) %in% colnames(filtered$beta_matrix))

  return(list(
    beta_matrix = filtered$beta_matrix,
    se_matrix = filtered$se_matrix,
    trait_annotations = trait_annot$annotations,
    snp_annotations = snp_annot,
    trait_info = trait_info,
    snp_info = snp_info,
    target_trait_id = trait_id,
    label_scheme = label_scheme,
    pathway_enrichment = trait_annot$pathway_enrichment
  ))
}


#' @title Run EBMF Factorization
#' @description Fit an Empirical Bayes Matrix Factorization model using
#' \pkg{flashier}. The model decomposes the traits x SNPs effect matrix
#' \eqn{Y = LF' + E} into latent factors (\eqn{F}, SNP programs) and
#' loadings (\eqn{L}, trait contributions), with sparse priors estimated
#' from the data.
#'
#' Known per-observation standard errors are passed to flashier via \code{S}
#' (typically \code{build_ebmf_matrix()$se_matrix}). Residual variance structure
#' is controlled by \code{var_type}.
#'
#' @param beta_matrix Traits x SNPs matrix of signed effect sizes, typically
#'   \code{build_ebmf_matrix()$beta_matrix}. Rows and columns that are entirely
#'   missing or zero are removed before fitting. \code{NA} entries are treated
#'   as missing data.
#' @param se_matrix Traits x SNPs matrix of standard errors aligned with
#'   \code{beta_matrix}. Passed to flashier as \code{S}.
#' @param greedy_Kmax Maximum number of factors to add greedily. Factors are
#'   only added while they improve the variational lower bound. Defaults to 50.
#' @param backfit Logical; if \code{TRUE} (default), all factors are cyclically
#'   updated after the greedy phase until convergence.
#' @param ebnm_fn EBNM function for priors on loadings and factors. Defaults to
#'   \code{ebnm::ebnm_point_normal} (point-normal prior). Pass a list of two
#'   functions to use different priors for loadings (traits) and factors (SNPs).
#' @param var_type Residual variance structure. Defaults to \code{c(1, 2)} for
#'   Kronecker (two-way) estimation. Use \code{0} for constant variance,
#'   \code{1} for per-row, or \code{2} for per-column.
#' @param verbose Verbosity: 0 = silent, 1 = summary, 2 = ELBO updates,
#'   3 = per-iteration.
#' @return A \code{flash} object from \pkg{flashier}. Key elements:
#'   \itemize{
#'     \item n_factors: number of discovered factors (programs)
#'     \item L_pm, L_lfsr: posterior means and lFSR for trait loadings (n x K)
#'     \item F_pm, F_lfsr: posterior means and lFSR for SNP factors (p x K)
#'     \item pve: proportion of variance explained per factor
#'     \item elbo: variational lower bound
#'   }
#' @export
run_ebmf <- function(beta_matrix,
                     se_matrix,
                     greedy_Kmax = 50L,
                     backfit = TRUE,
                     ebnm_fn = NULL,
                     var_type = c(1L, 2L),
                     verbose = 1L) {
  if (!requireNamespace("flashier", quietly = TRUE)) {
    stop("Package 'flashier' is required for EBMF clustering", call. = FALSE)
  }
  if (!requireNamespace("ebnm", quietly = TRUE)) {
    stop("Package 'ebnm' is required for EBMF clustering", call. = FALSE)
  }

  if (!is.matrix(beta_matrix)) {
    stop("beta_matrix must be a matrix")
  }
  if (!is.matrix(se_matrix)) {
    stop("se_matrix must be a matrix")
  }
  if (!identical(dim(beta_matrix), dim(se_matrix))) {
    stop("beta_matrix and se_matrix must have the same dimensions")
  }

  if (is.null(ebnm_fn)) {
    ebnm_fn <- ebnm::ebnm_point_normal
  }

  filtered <- .filter_ebmf_matrices(
    beta_matrix = beta_matrix,
    se_matrix = se_matrix
  )

  fit <- flashier::flash(
    data = filtered$beta_matrix,
    S = filtered$se_matrix,
    ebnm_fn = ebnm_fn,
    var_type = var_type,
    greedy_Kmax = greedy_Kmax,
    backfit = backfit,
    nullcheck = TRUE,
    verbose = verbose
  )

  return(fit)
}


#' @title Run EBMF Comparison Grid
#' @description Fit a grid of EBMF models varying label scheme, EBNM prior, and
#' residual variance structure. Each label scheme is built via
#' \code{build_ebmf_matrix()} once and reused across matching runs.
#' @inheritParams build_ebmf_matrix
#' @inheritParams run_ebmf
#' @param label_schemes Character vector of row-labelling schemes to compare.
#' @param ebnm_fns Character vector of EBNM priors: \code{"point_normal"} and/or
#'   \code{"point_laplace"}.
#' @param var_types List of \code{var_type} values passed to \code{run_ebmf()}.
#'   Defaults to per-row (\code{1L}) and Kronecker (\code{c(1L, 2L)}).
#' @param lfsr_threshold Passed to \code{extract_ebmf_clusters()}.
#' @param magnitude_threshold Passed to \code{extract_ebmf_clusters()}.
#' @param save_path Optional path to save results as \code{.rds}.
#' @return A list with \code{summary} (one-row-per-run dataframe) and
#'   \code{results} (named list of per-run outputs).
#' @export
run_ebmf_comparison <- function(trait_id,
                                coloc_groups = NULL,
                                p_threshold = NULL,
                                snp_key = c("variant_id", "display_snp", "coloc_group_id"),
                                label_schemes = c(
                                  "pathway_gene_tissue", "pathway_gene", "tissue"
                                ),
                                pathway_source = NULL,
                                pathway_p_value_threshold = 0.05,
                                minimum_count_in_network = NULL,
                                ebnm_fns = c("point_normal", "point_laplace"),
                                var_types = list(c(1L), c(1L, 2L)),
                                greedy_Kmax = 20L,
                                backfit = TRUE,
                                lfsr_threshold = 0.05,
                                magnitude_threshold = 0.10,
                                save_path = NULL,
                                verbose = 0L) {
  if (missing(trait_id) || is.null(trait_id)) {
    stop("trait_id is required")
  }

  snp_key <- match.arg(snp_key)
  label_schemes <- match.arg(
    label_schemes,
    c("pathway_gene_tissue", "pathway_gene", "tissue"),
    several.ok = TRUE
  )

  configs <- expand.grid(
    label_scheme = label_schemes,
    ebnm_fn = ebnm_fns,
    var_type_idx = seq_along(var_types),
    stringsAsFactors = FALSE
  )
  configs$var_type <- vapply(
    configs$var_type_idx,
    function(i) return(paste(var_types[[i]], collapse = ",")),
    character(1)
  )
  configs$run_id <- paste0(
    "labels_", configs$label_scheme,
    "__ebnm_", configs$ebnm_fn,
    "__var_", configs$var_type
  )

  n_runs <- nrow(configs)
  results <- vector("list", n_runs)
  names(results) <- configs$run_id
  summary_rows <- vector("list", n_runs)
  ebmf_data_by_scheme <- list()

  for (i in seq_len(n_runs)) {
    cfg <- configs[i, , drop = FALSE]
    run_id <- cfg$run_id
    label_scheme <- cfg$label_scheme
    message("EBMF run ", i, "/", n_runs, ": ", run_id)

    run_result <- tryCatch({
      if (is.null(ebmf_data_by_scheme[[label_scheme]])) {
        ebmf_data_by_scheme[[label_scheme]] <- build_ebmf_matrix(
          trait_id = trait_id,
          coloc_groups = coloc_groups,
          p_threshold = p_threshold,
          snp_key = snp_key,
          label_scheme = label_scheme,
          pathway_source = pathway_source,
          pathway_p_value_threshold = pathway_p_value_threshold,
          minimum_count_in_network = minimum_count_in_network
        )
      }
      ebmf_data <- ebmf_data_by_scheme[[label_scheme]]

      flash_fit <- run_ebmf(
        beta_matrix = ebmf_data$beta_matrix,
        se_matrix = ebmf_data$se_matrix,
        greedy_Kmax = greedy_Kmax,
        backfit = backfit,
        ebnm_fn = .resolve_ebnm_fn(cfg$ebnm_fn),
        var_type = var_types[[cfg$var_type_idx]],
        verbose = verbose
      )

      clusters <- extract_ebmf_clusters(
        flash_fit,
        lfsr_threshold = lfsr_threshold,
        magnitude_threshold = magnitude_threshold
      )

      drivers <- summarise_ebmf_program_drivers(
        flash_fit = flash_fit,
        trait_info = ebmf_data$trait_info
      )

      list(
        label_scheme = label_scheme,
        ebnm_fn = cfg$ebnm_fn,
        var_type = var_types[[cfg$var_type_idx]],
        ebmf_data = ebmf_data,
        flash_fit = flash_fit,
        ebmf_clusters = clusters,
        program_drivers = drivers
      )
    }, error = function(e) {
      warning("Run '", run_id, "' failed: ", conditionMessage(e), call. = FALSE)
      return(list(
        label_scheme = label_scheme,
        ebnm_fn = cfg$ebnm_fn,
        var_type = var_types[[cfg$var_type_idx]],
        error = conditionMessage(e)
      ))
    })

    results[[run_id]] <- run_result
    summary_rows[[i]] <- .ebmf_run_summary_row(run_id, run_result)
  }

  summary_df <- dplyr::bind_rows(summary_rows)
  output <- list(summary = summary_df, results = results)

  if (!is.null(save_path)) {
    saveRDS(output, save_path)
    message("Saved to ", save_path)
  }

  return(output)
}


#' @title Summarise EBMF Comparison Results
#' @description Rebuild or return the summary table from a comparison result.
#' @param comparison Output list from \code{run_ebmf_comparison()}.
#' @return A dataframe with one row per run: \code{run_id}, \code{label_scheme},
#'   \code{ebnm_fn}, \code{var_type}, \code{n_factors}, \code{n_programs_with_snps},
#'   \code{n_assigned}, \code{n_multi_program}, \code{program_sizes},
#'   \code{status}, and \code{error}.
#' @export
summarise_ebmf_comparison <- function(comparison) {
  if (!is.null(comparison$summary)) {
    return(comparison$summary)
  }
  summary_rows <- lapply(names(comparison$results), function(run_id) {
    return(.ebmf_run_summary_row(run_id, comparison$results[[run_id]]))
  })
  return(dplyr::bind_rows(summary_rows))
}


#' @title Summarise EBMF Program Trait Drivers
#' @description Rank trait loadings for each discovered EBMF program.
#' @param flash_fit A \code{flash} object from \code{run_ebmf()}.
#' @param trait_info Optional trait metadata with \code{trait_id} and \code{trait_name}.
#' @param lfsr_threshold lFSR threshold for trait loadings. Defaults to \code{0.05}.
#' @param magnitude_threshold Minimum absolute loading. Defaults to \code{0.05}.
#' @return A dataframe of significant trait loadings per program.
#' @export
summarise_ebmf_program_drivers <- function(flash_fit,
                                           trait_info = NULL,
                                           lfsr_threshold = 0.05,
                                           magnitude_threshold = 0.05) {
  if (!inherits(flash_fit, "flash")) {
    stop("flash_fit must be a flash object from run_ebmf()")
  }

  K <- flash_fit$n_factors
  if (K == 0) {
    return(data.frame(
      program = integer(0),
      trait_id = character(0),
      loading = numeric(0),
      lfsr = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  L_pm <- flash_fit$L_pm
  L_lfsr <- flash_fit$L_lfsr
  trait_ids <- rownames(L_pm)
  if (is.null(trait_ids)) {
    trait_ids <- as.character(seq_len(nrow(L_pm)))
  }

  summaries <- lapply(seq_len(K), function(k) {
    sig <- !is.na(L_lfsr[, k]) &
      L_lfsr[, k] < lfsr_threshold &
      abs(L_pm[, k]) > magnitude_threshold
    if (sum(sig) == 0) {
      return(NULL)
    }
    out <- data.frame(
      program = k,
      trait_id = trait_ids[sig],
      loading = L_pm[sig, k],
      lfsr = L_lfsr[sig, k],
      stringsAsFactors = FALSE
    ) |>
      dplyr::arrange(dplyr::desc(abs(loading)))
    return(out)
  })

  result <- dplyr::bind_rows(summaries)

  if (!is.null(trait_info) && nrow(result) > 0) {
    trait_info <- trait_info |>
      dplyr::mutate(trait_id = as.character(trait_id))
    result <- result |>
      dplyr::left_join(trait_info, by = "trait_id")
  }

  return(result)
}


#' @title Extract EBMF SNP Clusters
#' @description Extract overlapping SNP clusters from a fitted \code{flash}
#' object. For each discovered factor (program), SNPs are selected if their
#' local false sign rate (lFSR) is below \code{lfsr_threshold} \emph{and}
#' their absolute posterior mean factor loading exceeds
#' \code{magnitude_threshold}. SNPs may belong to multiple programs
#' (overlapping clusters).
#'
#' @param flash_fit A \code{flash} object from \code{run_ebmf()}.
#' @param lfsr_threshold lFSR threshold for including a SNP in a program.
#'   Defaults to 0.05.
#' @param magnitude_threshold Minimum absolute posterior mean factor loading.
#'   Defaults to 0.10.
#' @return A list with:
#'   \itemize{
#'     \item clusters: dataframe with columns \code{program}, \code{snp_id},
#'       \code{loading}, \code{lfsr}
#'     \item membership: logical matrix (SNPs x programs)
#'     \item n_programs: number of discovered programs
#'     \item n_assigned: number of SNPs assigned to at least one program
#'     \item n_multi_program: number of SNPs assigned to more than one program
#'     \item lfsr_threshold: threshold used
#'     \item magnitude_threshold: threshold used
#'   }
#' @export
extract_ebmf_clusters <- function(flash_fit,
                                  lfsr_threshold = 0.05,
                                  magnitude_threshold = 0.10) {
  if (!inherits(flash_fit, "flash")) {
    stop("flash_fit must be a flash object from run_ebmf()")
  }

  K <- flash_fit$n_factors
  if (K == 0) {
    return(list(
      clusters = data.frame(
        program = integer(0), snp_id = character(0),
        loading = numeric(0), lfsr = numeric(0),
        stringsAsFactors = FALSE
      ),
      membership = matrix(FALSE, nrow = 0, ncol = 0),
      n_programs = 0L,
      n_assigned = 0L,
      n_multi_program = 0L,
      lfsr_threshold = lfsr_threshold,
      magnitude_threshold = magnitude_threshold
    ))
  }

  F_pm <- flash_fit$F_pm
  F_lfsr <- flash_fit$F_lfsr

  snp_ids <- rownames(F_pm)
  if (is.null(snp_ids)) {
    snp_ids <- as.character(seq_len(nrow(F_pm)))
  }

  membership <- matrix(
    FALSE,
    nrow = nrow(F_pm),
    ncol = K,
    dimnames = list(snp_ids, paste0("program_", seq_len(K)))
  )

  cluster_rows <- vector("list", K)
  for (k in seq_len(K)) {
    if (is.na(lfsr_threshold)) {
      pass_lfsr <- TRUE
    } else {
      pass_lfsr <- !is.na(F_lfsr[, k]) & F_lfsr[, k] < lfsr_threshold
    }
    if (is.na(magnitude_threshold)) {
      pass_mag <- TRUE
    } else {
      pass_mag <- !is.na(F_pm[, k]) & abs(F_pm[, k]) > magnitude_threshold
    }
    selected <- pass_lfsr & pass_mag
    membership[, k] <- selected

    idx <- which(selected)
    if (length(idx) > 0) {
      cluster_rows[[k]] <- data.frame(
        program = k,
        snp_id = snp_ids[idx],
        loading = F_pm[idx, k],
        lfsr = F_lfsr[idx, k],
        stringsAsFactors = FALSE
      )
    }
  }

  clusters <- dplyr::bind_rows(cluster_rows)
  row_sums <- rowSums(membership)

  return(list(
    clusters = clusters,
    membership = membership,
    n_programs = K,
    n_assigned = sum(row_sums > 0),
    n_multi_program = sum(row_sums > 1),
    lfsr_threshold = lfsr_threshold,
    magnitude_threshold = magnitude_threshold
  ))
}


.ebmf_slice_has_signal <- function(x) {
  return(any(!is.na(x) & x != 0))
}


.ebmf_run_summary_row <- function(run_id, run_result) {
  label_scheme <- if (is.null(run_result$label_scheme)) {
    NA_character_
  } else {
    as.character(run_result$label_scheme)
  }

  if (!is.null(run_result$error)) {
    ebnm_fn <- if (is.null(run_result$ebnm_fn)) NA_character_ else as.character(run_result$ebnm_fn)
    var_type <- if (is.null(run_result$var_type)) NA_character_ else paste(run_result$var_type, collapse = ",")
    return(data.frame(
      run_id = run_id,
      label_scheme = label_scheme,
      ebnm_fn = ebnm_fn,
      var_type = var_type,
      n_factors = NA_integer_,
      n_programs_with_snps = NA_integer_,
      n_assigned = NA_integer_,
      n_multi_program = NA_integer_,
      program_sizes = NA_character_,
      status = "failed",
      error = run_result$error,
      stringsAsFactors = FALSE
    ))
  }

  clusters <- run_result$ebmf_clusters
  program_sizes <- integer(0)
  if (!is.null(clusters$membership) && ncol(clusters$membership) > 0) {
    program_sizes <- sort(colSums(clusters$membership), decreasing = TRUE)
  }

  return(data.frame(
    run_id = run_id,
    label_scheme = label_scheme,
    ebnm_fn = as.character(run_result$ebnm_fn),
    var_type = paste(run_result$var_type, collapse = ","),
    n_factors = run_result$flash_fit$n_factors,
    n_programs_with_snps = sum(program_sizes > 0),
    n_assigned = clusters$n_assigned,
    n_multi_program = clusters$n_multi_program,
    program_sizes = paste(program_sizes, collapse = ", "),
    status = "ok",
    error = NA_character_,
    stringsAsFactors = FALSE
  ))
}


.filter_ebmf_matrices <- function(beta_matrix, se_matrix) {
  if (!is.matrix(beta_matrix)) {
    stop("beta_matrix must be a matrix")
  }
  if (!is.matrix(se_matrix)) {
    stop("se_matrix must be a matrix")
  }

  row_keep <- apply(beta_matrix, 1, .ebmf_slice_has_signal)
  col_keep <- apply(beta_matrix, 2, .ebmf_slice_has_signal)

  removed_rows <- sum(!row_keep)
  removed_cols <- sum(!col_keep)
  if (removed_rows > 0 || removed_cols > 0) {
    warning(
      "Removed ", removed_rows, " all-zero/all-missing row(s) and ",
      removed_cols, " all-zero/all-missing column(s) for EBMF",
      call. = FALSE
    )
  }

  beta_matrix <- beta_matrix[row_keep, col_keep, drop = FALSE]
  se_matrix <- se_matrix[row_keep, col_keep, drop = FALSE]

  if (nrow(beta_matrix) == 0 || ncol(beta_matrix) == 0) {
    stop("No rows or columns remain after removing all-zero/all-missing slices")
  }

  return(list(
    beta_matrix = beta_matrix,
    se_matrix = se_matrix
  ))
}


.build_ebmf_association_matrices <- function(trait_id,
                                             coloc_groups,
                                             p_threshold,
                                             snp_key) {
  locus_data <- .prepare_pleiotropy_locus_data(
    trait_id = trait_id,
    coloc_groups = coloc_groups,
    p_threshold = p_threshold,
    snp_key = snp_key
  )

  assoc_long <- locus_data$cg |>
    dplyr::group_by(trait_id, trait_name, snp_id) |>
    dplyr::slice_min(min_p, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(trait_id, trait_name, snp_id, beta, se)

  beta_wide <- assoc_long |>
    tidyr::pivot_wider(
      id_cols = c(trait_id, trait_name),
      names_from = snp_id,
      values_from = beta
    )
  se_wide <- assoc_long |>
    tidyr::pivot_wider(
      id_cols = c(trait_id, trait_name),
      names_from = snp_id,
      values_from = se
    )

  snp_ids <- setdiff(names(beta_wide), c("trait_id", "trait_name"))
  beta_matrix <- as.matrix(beta_wide[, snp_ids, drop = FALSE])
  se_matrix <- as.matrix(se_wide[, snp_ids, drop = FALSE])
  rownames(beta_matrix) <- as.character(beta_wide$trait_id)
  rownames(se_matrix) <- as.character(se_wide$trait_id)

  trait_info <- beta_wide |>
    dplyr::select(trait_id, trait_name) |>
    dplyr::distinct()

  snp_info <- locus_data$target_snps |>
    dplyr::filter(snp_id %in% snp_ids) |>
    dplyr::distinct()

  return(list(
    beta_matrix = beta_matrix,
    se_matrix = se_matrix,
    trait_info = trait_info,
    snp_info = snp_info
  ))
}


.build_ebmf_trait_annotations <- function(trait_ids,
                                          coloc_groups,
                                          label_scheme,
                                          pathway_source = NULL,
                                          pathway_p_value_threshold = 0.05,
                                          minimum_count_in_network = NULL) {
  trait_meta <- coloc_groups |>
    dplyr::filter(as.character(trait_id) %in% trait_ids) |>
    dplyr::group_by(trait_id) |>
    dplyr::summarise(
      tissue = dplyr::first(tissue[!is.na(tissue) & tissue != ""]),
      gene = dplyr::first(gene[!is.na(gene) & gene != ""]),
      gene_id = dplyr::first(gene_id[!is.na(gene_id)]),
      trait_category = dplyr::first(
        trait_category[!is.na(trait_category) & trait_category != ""]
      ),
      is_qtl = any(!is.na(gene_id)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      trait_id = as.character(trait_id),
      tissue = dplyr::if_else(is.na(tissue), "None", tissue),
      gene = dplyr::if_else(is.na(gene), "None", gene),
      trait_category = dplyr::if_else(is.na(trait_category), "Unknown", trait_category)
    )

  gene_ids <- unique(stats::na.omit(trait_meta$gene_id[trait_meta$is_qtl]))
  feature_map <- data.frame(
    gene_id = integer(0),
    feature_name = character(0),
    stringsAsFactors = FALSE
  )
  pathway_enrichment <- NULL

  if (length(gene_ids) > 0 &&
      label_scheme %in% c("pathway_gene_tissue", "pathway_gene")) {
    tryCatch({
      map_result <- build_pathway_feature_map(
        genes = gene_ids,
        source = pathway_source,
        p_value_threshold = pathway_p_value_threshold,
        minimum_count_in_network = minimum_count_in_network
      )
      pathway_enrichment <- map_result$pathway_enrichment
      if (nrow(map_result$feature_map) > 0) {
        feature_map <- map_result$feature_map |>
          dplyr::group_by(gene_id) |>
          dplyr::summarise(
            feature_name = dplyr::first(feature_name),
            .groups = "drop"
          )
      }
    }, error = function(e) {
      warning("Pathway enrichment failed: ", conditionMessage(e), call. = FALSE)
      return(NULL)
    })
  }

  if (nrow(feature_map) > 0) {
    trait_meta <- trait_meta |>
      dplyr::left_join(feature_map, by = "gene_id")
  } else {
    trait_meta$feature_name <- NA_character_
  }

  trait_meta <- trait_meta |>
    dplyr::mutate(
      pathway_or_gene = dplyr::if_else(
        !is_qtl,
        NA_character_,
        dplyr::if_else(
          !is.na(feature_name),
          feature_name,
          dplyr::if_else(gene != "None", gene, paste0("gene:", gene_id))
        )
      )
    )

  meta_lookup <- stats::setNames(
    vapply(trait_meta$trait_id, function(tid) {
      row <- trait_meta[trait_meta$trait_id == tid, , drop = FALSE]
      if (!row$is_qtl) {
        return(row$trait_category)
      }
      if (label_scheme == "pathway_gene_tissue") {
        return(paste0(row$pathway_or_gene, "_", row$tissue))
      }
      if (label_scheme == "pathway_gene") {
        return(row$pathway_or_gene)
      }
      return(row$tissue)
    }, character(1)),
    trait_meta$trait_id
  )

  annotations <- vapply(trait_ids, function(tid) {
    if (tid %in% names(meta_lookup)) {
      return(meta_lookup[[tid]])
    }
    return("Unknown")
  }, character(1))
  names(annotations) <- trait_ids

  return(list(
    annotations = annotations,
    pathway_enrichment = pathway_enrichment
  ))
}


.build_ebmf_snp_annotations <- function(snp_ids, snp_info) {
  snp_meta <- snp_info |>
    dplyr::transmute(
      snp_id = as.character(snp_id),
      coloc_group_id = coloc_group_id
    )

  meta_lookup <- stats::setNames(
    paste0("ColocGroup", snp_meta$coloc_group_id),
    snp_meta$snp_id
  )

  annotations <- vapply(snp_ids, function(sid) {
    if (sid %in% names(meta_lookup)) {
      return(meta_lookup[[sid]])
    }
    return("ColocGroupNA")
  }, character(1))
  names(annotations) <- snp_ids

  return(annotations)
}


.resolve_ebnm_fn <- function(ebnm_fn) {
  if (is.function(ebnm_fn)) {
    return(ebnm_fn)
  }
  if (!is.character(ebnm_fn) || length(ebnm_fn) != 1L) {
    stop("ebnm_fn must be a function or character string")
  }

  resolved <- switch(
    ebnm_fn,
    point_normal = ebnm::ebnm_point_normal,
    point_laplace = ebnm::ebnm_point_laplace,
    stop("Unknown ebnm_fn: ", ebnm_fn, call. = FALSE)
  )

  return(resolved)
}
