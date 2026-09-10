#' @title Extract EBMF Program Loadings For Cross-Trait Comparison
#' @description Pull the SNP loadings and high-confidence membership of every
#' (optionally valid-only) EBMF program out of a `run_univariate_clustering()`
#' result so that programs discovered independently for different target traits
#' can be compared. Each program's SNP loading vector is sign-oriented to the
#' target trait before comparison: the factor sign in EBMF is not identifiable,
#' so `F_pm` is reflected so the target trait's loading (`L_pm`) is positive.
#' When the target trait does not load on a program (or `L_pm` is unavailable),
#' the program is alternatively anchored on the sign of its largest-|loading|
#' SNP. This makes a positive cross-trait correlation interpretable as a
#' concordant (target-aligned) architecture. High-confidence membership follows
#' the same lFSR / magnitude gate used elsewhere in the pipeline.
#' @param clustering_result Result of [run_univariate_clustering()] with
#'   `cluster_type = "ebmf"`.
#' @param program_summary Optional result of [summarise_ebmf_programs()], used
#'   only to restrict to programs with `status == "valid"` when
#'   `valid_only = TRUE`.
#' @param valid_only If `TRUE` (default) and `program_summary` is supplied, keep
#'   only programs whose `status == "valid"`. Otherwise all fitted programs are
#'   returned.
#' @return A list with:
#'   \itemize{
#'     \item loadings: one row per SNP x program with `program_id`, `trait_id`,
#'       `trait_name`, `program`, `snp_id`, `loading` (sign-oriented),
#'       `abs_loading`, `lfsr`, `high_confidence`
#'     \item programs: one row per program with `program_id`, `trait_id`,
#'       `trait_name`, `program`, `n_snps`, `n_high_confidence`, `sign`
#'   }
#' @export
extract_program_loadings <- function(clustering_result,
                                     program_summary = NULL,
                                     valid_only = TRUE) {
  params <- clustering_result$parameters
  if (is.null(params) || !identical(params$cluster_type, "ebmf")) {
    stop("clustering_result must come from cluster_type = 'ebmf'")
  }

  posterior <- ebmf_posterior_table(clustering_result)
  if (nrow(posterior) == 0) {
    return(list(
      loadings = .empty_program_loadings(),
      programs = .empty_program_meta()
    ))
  }

  fit <- clustering_result$cluster_details$flash_fit
  trait_id <- as.character(params$target_trait_id)
  trait_name <- .program_trait_name(clustering_result, trait_id)
  programs <- sort(unique(posterior$program))

  if (valid_only && !is.null(program_summary) &&
        !is.null(program_summary$programs) && nrow(program_summary$programs) > 0) {
    programs <- program_summary$programs$program[
      program_summary$programs$status == "valid"
    ]
  }
  posterior <- posterior[posterior$program %in% programs, , drop = FALSE]
  if (nrow(posterior) == 0) {
    return(list(
      loadings = .empty_program_loadings(),
      programs = .empty_program_meta()
    ))
  }

  signs <- .program_orientation_signs(fit, trait_id, sort(unique(posterior$program)))
  program_signs <- signs[as.character(posterior$program)]

  loadings <- data.frame(
    program_id = paste0(trait_id, ":", posterior$program),
    trait_id = trait_id,
    trait_name = trait_name,
    program = as.integer(posterior$program),
    snp_id = as.character(posterior$snp_id),
    loading = posterior$loading * program_signs,
    abs_loading = posterior$abs_loading,
    lfsr = posterior$lfsr,
    high_confidence = .program_high_confidence(
      posterior$lfsr,
      posterior$abs_loading,
      params$ebmf_lfsr_threshold,
      params$ebmf_magnitude_threshold
    ),
    stringsAsFactors = FALSE
  )

  programs_meta <- loadings |>
    dplyr::group_by(program_id, trait_id, trait_name, program) |>
    dplyr::summarise(
      n_snps = dplyr::n(),
      n_high_confidence = sum(high_confidence),
      .groups = "drop"
    ) |>
    dplyr::mutate(sign = signs[as.character(program)]) |>
    dplyr::select(
      program_id, trait_id, trait_name, program, n_snps, n_high_confidence, sign
    )

  return(list(loadings = loadings, programs = programs_meta))
}


#' @title Compare EBMF Programs Across Traits
#' @description Compare every cross-trait pair of EBMF programs. Program SNP
#' loading vectors are aligned on the union of their `variant_id`s (absent
#' SNPs treated as 0), entries that are effectively zero for both programs are
#' dropped, and the primary similarity is the Pearson correlation of the aligned
#' vectors. A supporting Jaccard overlap is computed on the high-confidence
#' membership sets. Because factor signs are arbitrary across independent fits
#' (see [extract_program_loadings()]), a positive `r` is a concordant
#' architecture and a negative `r` an antagonistic one.
#'
#' A matched null distribution for `r` is generated for each pair by permuting
#' one program's loadings across the aligned SNP slots, which preserves the
#' program's size and loading distribution while breaking the biological
#' correspondence between the two programs. The two-sided empirical p-value
#' `(1 + #{|r_null| >= |r_obs|}) / (1 + n_perm)` is BH-corrected across all
#' cross-trait pairs; a pair is `shared` when its FDR is at or below
#' `fdr_threshold`.
#' @param program_loadings Either a single extraction result from
#'   [extract_program_loadings()], a list of such results, or a combined long
#'   data.frame with columns `program_id`, `trait_id`, `program`, `snp_id`,
#'   `loading`, and `high_confidence` (plus optional `trait_name`).
#' @param n_perm Number of within-pair loading-label permutations. Defaults to
#'   `1000`.
#' @param fdr_threshold BH FDR at or below which a pair is called shared.
#'   Defaults to `0.05`.
#' @param zero_tol Absolute loading at or below which a cell is treated as
#'   effectively zero. Entries that are effectively zero for both programs are
#'   dropped before correlation. Defaults to `1e-8`.
#' @param seed RNG seed for the permutations.
#' @return A list with:
#'   \itemize{
#'     \item pairs: one row per cross-trait program pair with `program_id_a`,
#'       `program_id_b`, `trait_id_a`, `trait_id_b`, `program_a`, `program_b`,
#'       `n_snps_union`, `n_snps_used`, `r`, `jaccard`,
#'       `n_high_conf_shared`, `n_high_conf_union`, `p_emp`, `fdr`, `shared`,
#'       and `direction` (`concordant` / `antagonistic`)
#'     \item settings: settings used
#'   }
#' @export
compare_program_pairs <- function(program_loadings,
                                  n_perm = 1000L,
                                  fdr_threshold = 0.05,
                                  zero_tol = 1e-8,
                                  seed = 1) {
  loadings <- .normalise_program_loadings(program_loadings)
  keys <- .program_pair_keys(loadings)

  settings <- list(
    n_perm = as.integer(n_perm),
    fdr_threshold = fdr_threshold,
    zero_tol = zero_tol,
    seed = seed
  )

  if (nrow(keys) == 0) {
    return(list(pairs = .empty_pairs(), settings = settings))
  }

  set.seed(seed)
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    return(.compare_one_pair(loadings, keys[i, , drop = FALSE], n_perm, zero_tol))
  })
  out <- dplyr::bind_rows(rows)

  out$fdr <- NA_real_
  ok <- is.finite(out$p_emp)
  if (any(ok)) {
    out$fdr[ok] <- stats::p.adjust(out$p_emp[ok], method = "BH")
  }
  out$shared <- is.finite(out$fdr) & out$fdr <= fdr_threshold
  out$direction <- dplyr::case_when(
    !is.finite(out$r) ~ NA_character_,
    out$r > 0 ~ "concordant",
    out$r < 0 ~ "antagonistic",
    TRUE ~ NA_character_
  )

  return(list(pairs = out, settings = settings))
}


#' @title Build Multi-Trait Program Families
#' @description Turn a significant program-correspondence graph into multi-trait
#' program families. Programs are nodes and significant pairs are undirected
#' edges (signed by their correlation). Families are the weak connected
#' components of that graph, so a family is any chain of significantly
#' corresponding programs across traits rather than a forced one-to-one match.
#' Isolated programs form their own size-1 family.
#' @param pairs Either a result from [compare_program_pairs()] or its `pairs`
#'   data.frame.
#' @param program_info Optional data.frame with `program_id` and `trait_id`
#'   (one row per program). Defaults to the program ids appearing in `pairs`.
#' @param fdr_threshold FDR at or below which a pair is treated as an edge.
#'   Defaults to `0.05`.
#' @return A list with:
#'   \itemize{
#'     \item nodes: one row per program with `program_id`, `trait_id`,
#'       `trait_name` (when available), and `family`
#'     \item families: one row per family with `family`, `n_programs`,
#'       `n_traits`, `traits`, `n_edges`, `mean_r`, `median_r`, `min_r`,
#'       `max_r`, `n_concordant`, `n_antagonistic`, `mean_jaccard`,
#'       and `median_jaccard`
#'     \item graph: an `igraph` object (or `NULL` when there are no edges)
#'     \item edges: the significant pair table used as edges
#'   }
#' @export
build_program_families <- function(pairs,
                                   program_info = NULL,
                                   fdr_threshold = 0.05) {
  if (is.list(pairs) && !is.data.frame(pairs) && !is.null(pairs$pairs)) {
    pairs <- pairs$pairs
  }
  if (!is.data.frame(pairs)) {
    stop("pairs must be a compare_program_pairs() result or its pairs data.frame")
  }

  shared <- pairs[is.finite(pairs$fdr) & pairs$fdr <= fdr_threshold, , drop = FALSE]
  endpoints <- .normalise_program_info(pairs, NULL)
  program_info <- .merge_program_info(endpoints, program_info)

  if (nrow(shared) == 0 || nrow(program_info) == 0) {
    nodes <- program_info |>
      dplyr::mutate(family = paste0("family_", dplyr::row_number())) |>
      dplyr::select(program_id, trait_id, trait_name, family)
    return(list(
      nodes = nodes,
      families = .empty_families(),
      graph = NULL,
      edges = shared
    ))
  }

  graph <- igraph::graph_from_data_frame(
    shared[, c("program_id_a", "program_id_b"), drop = FALSE],
    directed = FALSE,
    vertices = data.frame(
      name = program_info$program_id,
      stringsAsFactors = FALSE
    )
  )
  graph <- igraph::set_edge_attr(graph, "r", value = shared$r)
  graph <- igraph::set_edge_attr(graph, "weight", value = abs(shared$r))
  components <- igraph::components(graph, mode = "weak")

  nodes <- data.frame(
    program_id = igraph::V(graph)$name,
    family = paste0("family_", components$membership),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      program_info |> dplyr::select(program_id, trait_id, trait_name),
      by = "program_id"
    ) |>
    dplyr::select(program_id, trait_id, trait_name, family)

  families <- .family_metrics(nodes, shared)
  return(list(
    nodes = nodes,
    families = families,
    graph = graph,
    edges = shared
  ))
}


.program_trait_name <- function(clustering_result, trait_id) {
  info <- clustering_result$trait_info
  if (is.null(info) || !all(c("trait_id", "trait_name") %in% names(info))) {
    return(NA_character_)
  }
  hit <- info$trait_name[as.character(info$trait_id) == trait_id]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  return(as.character(hit[1]))
}


.program_orientation_signs <- function(fit, target_trait_id, programs) {
  l_pm <- if (is.null(fit)) NULL else fit$L_pm
  f_pm <- if (is.null(fit)) NULL else fit$F_pm
  signs <- vapply(programs, function(k) {
    sign_k <- NA_real_
    if (!is.null(l_pm) && target_trait_id %in% rownames(l_pm) && k <= ncol(l_pm)) {
      value <- l_pm[target_trait_id, k]
      if (is.finite(value) && value != 0) {
        sign_k <- sign(value)
      }
    }
    if ((is.na(sign_k) || sign_k == 0) && !is.null(f_pm) && k <= ncol(f_pm)) {
      values <- f_pm[, k]
      if (any(is.finite(values) & values != 0)) {
        sign_k <- sign(values[which.max(abs(values))])
      }
    }
    if (is.na(sign_k) || sign_k == 0) {
      sign_k <- 1
    }
    return(sign_k)
  }, numeric(1))
  names(signs) <- as.character(programs)
  return(signs)
}


.program_high_confidence <- function(lfsr, abs_loading, lfsr_threshold, magnitude_threshold) {
  support <- is.finite(abs_loading) & abs_loading > 0
  lfsr_ok <- if (is.na(lfsr_threshold)) {
    TRUE
  } else {
    !is.na(lfsr) & lfsr < lfsr_threshold
  }
  magnitude_ok <- if (is.na(magnitude_threshold)) {
    TRUE
  } else {
    support & abs_loading > magnitude_threshold
  }
  return(support & lfsr_ok & magnitude_ok)
}


.normalise_program_loadings <- function(program_loadings) {
  if (is.data.frame(program_loadings)) {
    loadings <- program_loadings
  } else if (is.list(program_loadings) && "loadings" %in% names(program_loadings) &&
               is.data.frame(program_loadings$loadings)) {
    loadings <- program_loadings$loadings
  } else if (is.list(program_loadings)) {
    loadings <- dplyr::bind_rows(lapply(program_loadings, function(x) {
      if (is.data.frame(x)) {
        return(x)
      }
      return(x$loadings)
    }))
  } else {
    stop("program_loadings must be a data.frame or a list of extraction results")
  }

  required <- c("program_id", "trait_id", "program", "snp_id", "loading", "high_confidence")
  missing <- setdiff(required, names(loadings))
  if (length(missing) > 0) {
    stop("program_loadings is missing columns: ", paste(missing, collapse = ", "))
  }
  loadings$program_id <- as.character(loadings$program_id)
  loadings$trait_id <- as.character(loadings$trait_id)
  loadings$snp_id <- as.character(loadings$snp_id)
  loadings$high_confidence <- as.logical(loadings$high_confidence)
  if (!"trait_name" %in% names(loadings)) {
    loadings$trait_name <- NA_character_
  }
  return(loadings)
}


.normalise_program_info <- function(pairs, program_info) {
  if (!is.null(program_info)) {
    info <- as.data.frame(program_info, stringsAsFactors = FALSE)
    if (!"trait_name" %in% names(info)) {
      info$trait_name <- NA_character_
    }
    return(dplyr::distinct(info, program_id, trait_id, trait_name))
  }
  if (nrow(pairs) == 0) {
    return(data.frame(
      program_id = character(0), trait_id = character(0),
      trait_name = character(0), stringsAsFactors = FALSE
    ))
  }
  from_a <- data.frame(
    program_id = pairs$program_id_a, trait_id = pairs$trait_id_a,
    trait_name = NA_character_, stringsAsFactors = FALSE
  )
  from_b <- data.frame(
    program_id = pairs$program_id_b, trait_id = pairs$trait_id_b,
    trait_name = NA_character_, stringsAsFactors = FALSE
  )
  return(dplyr::distinct(
    dplyr::bind_rows(from_a, from_b),
    program_id, trait_id, trait_name
  ))
}


.merge_program_info <- function(endpoints, program_info) {
  if (is.null(program_info)) {
    return(endpoints)
  }
  info <- as.data.frame(program_info, stringsAsFactors = FALSE)
  if (!"trait_name" %in% names(info)) {
    info$trait_name <- NA_character_
  }
  info <- dplyr::distinct(info, program_id, trait_id, trait_name)
  extra <- info[!info$program_id %in% endpoints$program_id, , drop = FALSE]
  out <- dplyr::bind_rows(endpoints, extra)
  names_from_info <- dplyr::select(info, program_id, trait_name_info = trait_name)
  return(out |>
    dplyr::left_join(names_from_info, by = "program_id") |>
    dplyr::mutate(
      trait_name = dplyr::coalesce(trait_name, trait_name_info)
    ) |>
    dplyr::select(-trait_name_info))
}


.program_pair_keys <- function(loadings) {
  programs <- dplyr::distinct(loadings, program_id, trait_id, program)
  if (nrow(programs) < 2) {
    return(.empty_pair_keys())
  }
  idx <- utils::combn(nrow(programs), 2)
  keys <- data.frame(
    program_id_a = programs$program_id[idx[1, ]],
    program_id_b = programs$program_id[idx[2, ]],
    trait_id_a = programs$trait_id[idx[1, ]],
    trait_id_b = programs$trait_id[idx[2, ]],
    program_a = programs$program[idx[1, ]],
    program_b = programs$program[idx[2, ]],
    stringsAsFactors = FALSE
  )
  return(keys[keys$trait_id_a != keys$trait_id_b, , drop = FALSE])
}


.compare_one_pair <- function(loadings, key, n_perm, zero_tol) {
  a <- loadings[loadings$program_id == key$program_id_a, , drop = FALSE]
  b <- loadings[loadings$program_id == key$program_id_b, , drop = FALSE]
  union_snps <- union(a$snp_id, b$snp_id)
  va <- a$loading[match(union_snps, a$snp_id)]
  vb <- b$loading[match(union_snps, b$snp_id)]
  va[is.na(va)] <- 0
  vb[is.na(vb)] <- 0
  keep <- abs(va) > zero_tol | abs(vb) > zero_tol
  va <- va[keep]
  vb <- vb[keep]

  r <- NA_real_
  p_emp <- NA_real_
  if (length(va) >= 3 && stats::sd(va) > 0 && stats::sd(vb) > 0) {
    r <- stats::cor(va, vb)
    null_r <- vapply(seq_len(n_perm), function(i) {
      return(abs(stats::cor(va, sample(vb))))
    }, numeric(1))
    p_emp <- (1 + sum(null_r >= abs(r))) / (1 + n_perm)
  }

  inter <- length(intersect(a$snp_id[a$high_confidence], b$snp_id[b$high_confidence]))
  uni <- length(union(a$snp_id[a$high_confidence], b$snp_id[b$high_confidence]))

  return(data.frame(
    program_id_a = key$program_id_a,
    program_id_b = key$program_id_b,
    trait_id_a = key$trait_id_a,
    trait_id_b = key$trait_id_b,
    program_a = as.integer(key$program_a),
    program_b = as.integer(key$program_b),
    n_snps_union = length(union_snps),
    n_snps_used = length(va),
    r = r,
    jaccard = if (uni > 0) inter / uni else NA_real_,
    n_high_conf_shared = inter,
    n_high_conf_union = uni,
    p_emp = p_emp,
    stringsAsFactors = FALSE
  ))
}


.family_metrics <- function(nodes, shared) {
  families <- sort(unique(nodes$family))
  rows <- lapply(families, function(fam) {
    pids <- nodes$program_id[nodes$family == fam]
    edges <- shared[
      shared$program_id_a %in% pids & shared$program_id_b %in% pids,
      ,
      drop = FALSE
    ]
    traits <- sort(unique(nodes$trait_id[nodes$family == fam]))
    return(data.frame(
      family = fam,
      n_programs = length(pids),
      n_traits = length(traits),
      traits = paste(traits, collapse = ", "),
      n_edges = nrow(edges),
      mean_r = .safe_stat(edges$r, mean),
      median_r = .safe_stat(edges$r, stats::median),
      min_r = .safe_stat(edges$r, min),
      max_r = .safe_stat(edges$r, max),
      n_concordant = sum(edges$r > 0, na.rm = TRUE),
      n_antagonistic = sum(edges$r < 0, na.rm = TRUE),
      mean_jaccard = .safe_stat(edges$jaccard, mean),
      median_jaccard = .safe_stat(edges$jaccard, stats::median),
      stringsAsFactors = FALSE
    ))
  })
  out <- dplyr::bind_rows(rows)
  out <- dplyr::arrange(
    out,
    dplyr::desc(n_programs), dplyr::desc(n_edges), family
  )
  return(out)
}


.safe_stat <- function(x, fn) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  return(as.numeric(fn(x)))
}


.empty_program_loadings <- function() {
  return(data.frame(
    program_id = character(0),
    trait_id = character(0),
    trait_name = character(0),
    program = integer(0),
    snp_id = character(0),
    loading = numeric(0),
    abs_loading = numeric(0),
    lfsr = numeric(0),
    high_confidence = logical(0),
    stringsAsFactors = FALSE
  ))
}


.empty_program_meta <- function() {
  return(data.frame(
    program_id = character(0),
    trait_id = character(0),
    trait_name = character(0),
    program = integer(0),
    n_snps = integer(0),
    n_high_confidence = integer(0),
    sign = numeric(0),
    stringsAsFactors = FALSE
  ))
}


.empty_pair_keys <- function() {
  return(data.frame(
    program_id_a = character(0),
    program_id_b = character(0),
    trait_id_a = character(0),
    trait_id_b = character(0),
    program_a = integer(0),
    program_b = integer(0),
    stringsAsFactors = FALSE
  ))
}


.empty_pairs <- function() {
  return(data.frame(
    program_id_a = character(0),
    program_id_b = character(0),
    trait_id_a = character(0),
    trait_id_b = character(0),
    program_a = integer(0),
    program_b = integer(0),
    n_snps_union = integer(0),
    n_snps_used = integer(0),
    r = numeric(0),
    jaccard = numeric(0),
    n_high_conf_shared = integer(0),
    n_high_conf_union = integer(0),
    p_emp = numeric(0),
    fdr = numeric(0),
    shared = logical(0),
    direction = character(0),
    stringsAsFactors = FALSE
  ))
}


.empty_families <- function() {
  return(data.frame(
    family = character(0),
    n_programs = integer(0),
    n_traits = integer(0),
    traits = character(0),
    n_edges = integer(0),
    mean_r = numeric(0),
    median_r = numeric(0),
    min_r = numeric(0),
    max_r = numeric(0),
    n_concordant = integer(0),
    n_antagonistic = integer(0),
    mean_jaccard = numeric(0),
    median_jaccard = numeric(0),
    stringsAsFactors = FALSE
  ))
}
