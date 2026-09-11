#' @title Extract EBMF Program Loci And Pleiotropic Profiles
#' @description Pull the locus loadings and pleiotropic profiles of every
#' (optionally valid-only) EBMF program out of a `run_univariate_clustering()`
#' result so that programs discovered independently for different target traits
#' can be compared. Two program dimensions are extracted:
#'
#' \describe{
#'   \item{`loadings`}{one row per candidate SNP (coloc group) x program. These
#'   are the EBMF `F_pm` SNP loadings, joined to the coloc-group identifiers in
#'   `clustering_result$snp_info` so programs can be compared at the locus level.}
#'   \item{`profiles`}{one row per background feature (trait) x program, from
#'   `L_pm`. These describe the *pleiotropic profile* of the program (which
#'   background traits it loads on).}
#' }
#'
#' Each program is sign-oriented to its target trait before comparison: the
#' factor sign in EBMF is not identifiable, so both `F_pm` and `L_pm` are
#' reflected by the sign of the target trait's `L_pm` loading (making that
#' loading positive) while preserving the rank-1 product. When the target trait
#' does not load on a program (or `L_pm` is unavailable), the program is
#' alternatively anchored on the sign of its largest-|loading| SNP. This makes a
#' positive cross-trait comparison interpretable as a concordant
#' (target-aligned) architecture. High-confidence membership follows the same
#' lFSR / magnitude gate used elsewhere in the pipeline.
#' @param clustering_result Result of [run_univariate_clustering()] (EBMF).
#' @param program_summary Optional result of [summarise_ebmf_programs()], used
#'   only to restrict to programs with `status == "valid"` when
#'   `valid_only = TRUE`.
#' @param valid_only If `TRUE` (default) and `program_summary` is supplied, keep
#'   only programs whose `status == "valid"`. Otherwise all fitted programs are
#'   returned.
#' @return A list with:
#'   \itemize{
#'     \item loadings: one row per locus x program with `program_id`, `trait_id`,
#'       `trait_name`, `program`, `snp_id`, `coloc_group_id`, `chr`, `bp`,
#'       `loading` (sign-oriented), `abs_loading`, `lfsr`, `high_confidence`
#'     \item profiles: one row per feature x program with `program_id`,
#'       `trait_id`, `trait_name`, `program`, `feature_trait_id`,
#'       `feature_trait_name`, `loading` (sign-oriented), `lfsr`
#'     \item programs: one row per program with `program_id`, `trait_id`,
#'       `trait_name`, `program`, `n_snps`, `n_high_confidence`, `n_features`,
#'       `sign`
#'   }
#' @export
extract_program_loadings <- function(clustering_result,
                                     program_summary = NULL,
                                     valid_only = TRUE) {
  params <- clustering_result$parameters

  posterior <- ebmf_posterior_table(clustering_result)
  if (nrow(posterior) == 0) {
    return(.empty_program_extraction())
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
    return(.empty_program_extraction())
  }

  programs <- sort(unique(posterior$program))
  signs <- .program_orientation_signs(fit, trait_id, programs)
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
  loadings <- .attach_locus_meta(
    loadings, .program_locus_meta(clustering_result)
  )

  profiles <- .program_profiles(
    fit, trait_id, trait_name, programs, signs, clustering_result$trait_info
  )

  programs_meta <- loadings |>
    dplyr::group_by(program_id, trait_id, trait_name, program) |>
    dplyr::summarise(
      n_snps = dplyr::n(),
      n_high_confidence = sum(high_confidence),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      profiles |>
        dplyr::count(program_id, name = "n_features"),
      by = "program_id"
    ) |>
    dplyr::mutate(
      n_features = dplyr::coalesce(n_features, 0L),
      sign = signs[as.character(program)]
    ) |>
    dplyr::select(
      program_id, trait_id, trait_name, program,
      n_snps, n_high_confidence, n_features, sign
    )

  return(list(
    loadings = loadings,
    profiles = profiles,
    programs = programs_meta
  ))
}


.empty_program_extraction <- function() {
  return(list(
    loadings = .empty_program_loadings(),
    profiles = .empty_program_profiles(),
    programs = .empty_program_meta()
  ))
}


.program_locus_meta <- function(clustering_result) {
  info <- clustering_result$snp_info
  if (is.null(info) || !is.data.frame(info) || !"snp_id" %in% names(info)) {
    return(NULL)
  }
  if (!"coloc_group_id" %in% names(info)) {
    info$coloc_group_id <- NA_character_
  }
  if (!"chr" %in% names(info)) {
    info$chr <- NA_character_
  }
  if (!"bp" %in% names(info)) {
    info$bp <- NA_real_
  }
  return(
    info |>
      dplyr::transmute(
        snp_id = as.character(snp_id),
        coloc_group_id = as.character(coloc_group_id),
        chr = as.character(chr),
        bp = as.numeric(bp)
      ) |>
      dplyr::distinct(snp_id, .keep_all = TRUE)
  )
}


.attach_locus_meta <- function(loadings, locus_meta) {
  if (is.null(locus_meta) || nrow(locus_meta) == 0) {
    loadings$coloc_group_id <- NA_character_
    loadings$chr <- NA_character_
    loadings$bp <- NA_real_
    return(loadings)
  }
  out <- dplyr::left_join(loadings, locus_meta, by = "snp_id")
  out$coloc_group_id <- as.character(out$coloc_group_id)
  return(out)
}


.program_profiles <- function(fit, trait_id, trait_name, programs, signs,
                              trait_info = NULL) {
  l_pm <- if (is.null(fit)) NULL else fit$L_pm
  l_lfsr <- if (is.null(fit)) NULL else fit$L_lfsr
  if (is.null(l_pm) || length(programs) == 0) {
    return(.empty_program_profiles())
  }
  feature_ids <- rownames(l_pm)
  if (is.null(feature_ids)) {
    feature_ids <- as.character(seq_len(nrow(l_pm)))
  }
  feature_ids <- as.character(feature_ids)
  feature_names <- feature_ids
  if (!is.null(trait_info) &&
        all(c("trait_id", "trait_name") %in% names(trait_info))) {
    feature_names <- as.character(
      trait_info$trait_name[match(feature_ids, as.character(trait_info$trait_id))]
    )
    feature_names[is.na(feature_names)] <- feature_ids[is.na(feature_names)]
  }

  rows <- lapply(programs, function(k) {
    if (k > ncol(l_pm)) {
      return(NULL)
    }
    lfsr_k <- if (!is.null(l_lfsr) && k <= ncol(l_lfsr)) {
      l_lfsr[, k]
    } else {
      rep(NA_real_, length(feature_ids))
    }
    return(data.frame(
      program_id = paste0(trait_id, ":", k),
      trait_id = trait_id,
      trait_name = trait_name,
      program = as.integer(k),
      feature_trait_id = feature_ids,
      feature_trait_name = feature_names,
      loading = l_pm[, k] * signs[as.character(k)],
      lfsr = lfsr_k,
      stringsAsFactors = FALSE
    ))
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    return(.empty_program_profiles())
  }
  return(out)
}


#' @title Compare EBMF Programs Across Traits
#' @description Compare every cross-trait pair of EBMF programs with a
#' program-specific, calibrated correspondence test. Two channels are used:
#'
#' \describe{
#'   \item{Locus channel}{each program's locus membership is its high-confidence
#'   set (or top-`top_k` loci) of `coloc_group_id`s; pairs must share at least
#'   `min_shared_loci` loci. The overlap is tested against a hypergeometric null
#'   (a random membership set of the same size drawn from that trait's candidate
#'   loci).}
#'   \item{Profile channel}{pleiotropic profiles (`L`) are aligned on shared
#'   background `feature_trait_id`s and compared with a confidence-weighted
#'   Pearson correlation. When `partial_out_common = TRUE`, each trait's common
#'   (mean) profile is regressed out first, so correlation reflects
#'   program-specific rather than trait-wide structure.}
#' }
#'
#' Significance is calibrated **per program**, not per pair: for each program
#' the best locus-matched partner is taken (`max |profile_r|`) and compared
#' against a permutation null built by shuffling feature labels within programs
#' and recomputing the per-program maximum. This accounts for searching over
#' partners. Program-level empirical p-values are BH-corrected
#' (`q`) and a program is `significant` when `q <= fdr_threshold`. A pair is a
#' `shared` edge when it is the mutual best match of two significant programs.
#' `direction` (concordant / antagonistic) is the sign of the deciding profile
#' correlation. Both `F` and `L` are sign-oriented to each target trait by
#' [extract_program_loadings()].
#' @param program_data An extraction result from [extract_program_loadings()],
#'   or a list of such results (one per target trait).
#' @param n_perm Number of permutation replicates for the per-program null.
#'   Defaults to `1000`.
#' @param fdr_threshold BH FDR at or below which a program is significant.
#'   Defaults to `0.05`.
#' @param seed RNG seed for the permutations.
#' @param min_shared_loci Minimum shared `coloc_group_id`s for a pair to be a
#'   locus-matched candidate. Defaults to `3`.
#' @param min_shared_features Minimum shared background features for the profile
#'   channel to be testable. Defaults to `3`.
#' @param exclude_target_traits If `TRUE` (default), the two target traits are
#'   excluded from the shared feature axis before correlating profiles.
#' @param locus_membership `"high_confidence"` (default; the lFSR / magnitude
#'   gate) or `"top_k"` (top `top_k` loci by confidence-weighted loading).
#' @param top_k Number of loci per program when `locus_membership = "top_k"`.
#'   Defaults to `25`.
#' @param partial_out_common If `TRUE` (default), regress out each trait's mean
#'   profile before correlating.
#' @return A list with:
#'   \itemize{
#'     \item pairs: one row per cross-trait program pair (with a `candidate`
#'       flag for locus-matched pairs) carrying locus overlap, `p_locus`,
#'       `profile_r` (raw), `profile_r_resid` (after partialling), `p_profile`,
#'       `q_profile`, `direction`, and `shared`.
#'     \item matches: one row per program carrying `best_partner`,
#'       `best_profile_r`, `best_locus_jaccard`, `emp_p`, `q`, and
#'       `significant`.
#'     \item settings: settings used
#'   }
#' @export
compare_program_pairs <- function(program_data,
                                  n_perm = 1000L,
                                  fdr_threshold = 0.05,
                                  seed = 1,
                                  min_shared_loci = 3L,
                                  min_shared_features = 3L,
                                  exclude_target_traits = TRUE,
                                  locus_membership = c("high_confidence", "top_k"),
                                  top_k = 25L,
                                  partial_out_common = TRUE) {
  locus_membership <- match.arg(locus_membership)
  pd <- .normalise_program_data(program_data)
  loadings <- pd$loadings
  profiles <- pd$profiles
  keys <- .program_pair_keys(loadings)

  settings <- list(
    n_perm = as.integer(n_perm),
    fdr_threshold = fdr_threshold,
    seed = seed,
    min_shared_loci = as.integer(min_shared_loci),
    min_shared_features = as.integer(min_shared_features),
    exclude_target_traits = exclude_target_traits,
    locus_membership = locus_membership,
    top_k = as.integer(top_k),
    partial_out_common = partial_out_common
  )

  if (nrow(keys) == 0) {
    return(list(
      pairs = .empty_pairs(),
      matches = .empty_matches(),
      settings = settings
    ))
  }

  progs <- dplyr::distinct(loadings, program_id, trait_id, trait_name, program)
  trait_ids_by_program <- stats::setNames(
    as.character(progs$trait_id), as.character(progs$program_id)
  )
  program_index <- stats::setNames(
    seq_len(nrow(progs)), as.character(progs$program_id)
  )

  locus_sets <- .program_locus_sets(loadings, locus_membership, top_k)
  universe_by_trait <- .locus_universe_by_trait(loadings)

  out <- dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(i) {
    return(.locus_pair_stats(keys[i, , drop = FALSE], locus_sets, universe_by_trait))
  }))
  out$candidate <- !is.na(out$n_loci_shared) & out$n_loci_shared >= min_shared_loci
  out$n_features_shared <- 0L
  out$n_features_union <- 0L
  out$profile_r <- NA_real_
  out$profile_r_resid <- NA_real_
  out$p_profile <- NA_real_

  max_null <- matrix(
    -Inf,
    nrow = as.integer(n_perm),
    ncol = nrow(progs),
    dimnames = list(NULL, as.character(progs$program_id))
  )
  best <- data.frame(
    program_id = as.character(progs$program_id),
    best_partner = NA_character_,
    best_abs_r = NA_real_,
    best_row = NA_integer_,
    stringsAsFactors = FALSE
  )

  trait_ids <- sort(unique(as.character(progs$trait_id)))
  set.seed(seed)
  if (length(trait_ids) >= 2) {
    trait_pairs <- utils::combn(trait_ids, 2, simplify = FALSE)
    for (tp in trait_pairs) {
      res <- .profile_trait_pair(
        out = out,
        profiles = profiles,
        trait_ids_by_program = trait_ids_by_program,
        program_index = program_index,
        t_a = tp[[1]],
        t_b = tp[[2]],
        n_perm = as.integer(n_perm),
        min_shared_features = as.integer(min_shared_features),
        exclude_target_traits = exclude_target_traits,
        partial_out_common = partial_out_common,
        max_null = max_null,
        best = best
      )
      out <- res$out
      max_null <- res$max_null
      best <- res$best
    }
  }

  # Per-program calibrated significance.
  obs_max <- stats::setNames(
    best$best_abs_r, as.character(best$program_id)
  )
  matched <- which(is.finite(obs_max))
  emp_p <- stats::setNames(rep(NA_real_, nrow(progs)), as.character(progs$program_id))
  for (pid in names(matched)) {
    idx <- program_index[[pid]]
    null_vals <- max_null[, idx]
    null_vals <- null_vals[is.finite(null_vals)]
    if (length(null_vals) > 0) {
      emp_p[[pid]] <- (1 + sum(null_vals >= obs_max[[pid]])) /
        (1 + length(null_vals))
    }
  }
  q <- stats::setNames(rep(NA_real_, nrow(progs)), as.character(progs$program_id))
  ok <- is.finite(emp_p)
  if (any(ok)) {
    q[ok] <- stats::p.adjust(emp_p[ok], method = "BH")
  }

  matches <- progs |>
    dplyr::mutate(
      best_partner = best$best_partner[match(
        as.character(program_id), as.character(best$program_id)
      )],
      best_profile_r = best$best_abs_r[match(
        as.character(program_id), as.character(best$program_id)
      )],
      best_locus_jaccard = out$jaccard_loci[best$best_row[match(
        as.character(program_id), as.character(best$program_id)
      )]],
      best_n_loci_shared = out$n_loci_shared[best$best_row[match(
        as.character(program_id), as.character(best$program_id)
      )]],
      emp_p = emp_p[as.character(program_id)],
      q = q[as.character(program_id)]
    )
  matches$significant <- is.finite(matches$q) & matches$q <= fdr_threshold

  out$q_profile <- NA_real_
  ok_pair <- is.finite(out$p_profile)
  if (any(ok_pair)) {
    out$q_profile[ok_pair] <- stats::p.adjust(out$p_profile[ok_pair], method = "BH")
  }
  out$direction <- dplyr::case_when(
    !is.finite(out$profile_r_resid) ~ NA_character_,
    out$profile_r_resid > 0 ~ "concordant",
    out$profile_r_resid < 0 ~ "antagonistic",
    TRUE ~ NA_character_
  )
  out$shared <- .mutual_best_mask(out, matches, trait_ids_by_program)

  return(list(pairs = out, matches = matches, settings = settings))
}


.program_locus_sets <- function(loadings, membership = "high_confidence",
                                top_k = 25L) {
  programs <- unique(as.character(loadings$program_id))
  sets <- lapply(programs, function(pid) {
    sub <- loadings[
      as.character(loadings$program_id) == pid & !is.na(loadings$coloc_group_id),
      ,
      drop = FALSE
    ]
    if (nrow(sub) == 0) {
      return(character(0))
    }
    if (identical(membership, "high_confidence")) {
      keep <- as.logical(sub$high_confidence)
    } else {
      w <- .program_weight_factor(sub$lfsr, sub$loading)
      ord <- order(-w, na.last = TRUE)
      keep <- rep(FALSE, nrow(sub))
      keep[utils::head(ord, top_k)] <- TRUE
    }
    return(unique(as.character(sub$coloc_group_id[keep])))
  })
  names(sets) <- programs
  return(sets)
}


.locus_universe_by_trait <- function(loadings) {
  keep <- !is.na(loadings$coloc_group_id)
  if (!any(keep)) {
    return(list())
  }
  dt <- data.frame(
    trait_id = as.character(loadings$trait_id[keep]),
    coloc_group_id = as.character(loadings$coloc_group_id[keep]),
    stringsAsFactors = FALSE
  )
  return(tapply(dt$coloc_group_id, dt$trait_id, function(x) unique(x)))
}


.locus_pair_stats <- function(key, locus_sets, universe_by_trait) {
  set_a <- locus_sets[[key$program_id_a]]
  set_b <- locus_sets[[key$program_id_b]]
  universe_b <- universe_by_trait[[key$trait_id_b]]
  if (is.null(set_a)) set_a <- character(0)
  if (is.null(set_b)) set_b <- character(0)
  if (is.null(universe_b)) universe_b <- character(0)

  shared <- intersect(set_a, set_b)
  union_size <- length(union(set_a, set_b))
  n_a <- length(set_a)
  n_b <- length(set_b)
  n_shared <- length(shared)
  k <- length(intersect(set_a, universe_b))
  N <- length(universe_b)
  p_locus <- if (N > 0 && n_b > 0) {
    stats::phyper(n_shared - 1, k, N - k, n_b, lower.tail = FALSE)
  } else {
    NA_real_
  }

  return(data.frame(
    program_id_a = key$program_id_a,
    program_id_b = key$program_id_b,
    trait_id_a = key$trait_id_a,
    trait_id_b = key$trait_id_b,
    program_a = as.integer(key$program_a),
    program_b = as.integer(key$program_b),
    n_loci_a = n_a,
    n_loci_b = n_b,
    n_loci_shared = n_shared,
    n_loci_union = union_size,
    jaccard_loci = if (union_size > 0) n_shared / union_size else NA_real_,
    p_locus = p_locus,
    stringsAsFactors = FALSE
  ))
}


.align_profile_matrix <- function(profiles, program_ids, features) {
  loading <- matrix(
    NA_real_,
    nrow = length(features), ncol = length(program_ids),
    dimnames = list(features, program_ids)
  )
  conf <- loading
  for (j in seq_along(program_ids)) {
    sub <- profiles[as.character(profiles$program_id) == program_ids[j], , drop = FALSE]
    idx <- match(features, as.character(sub$feature_trait_id))
    loading[, j] <- sub$loading[idx]
    conf[, j] <- .confidence_factor(sub$lfsr[idx])
  }
  return(list(loading = loading, conf = conf))
}


.partial_out_common <- function(mat) {
  if (ncol(mat) < 2) {
    return(mat)
  }
  common <- rowMeans(mat, na.rm = TRUE)
  if (!is.finite(stats::sd(common, na.rm = TRUE)) ||
        stats::sd(common, na.rm = TRUE) == 0) {
    return(mat)
  }
  design <- cbind(1, common)
  resid <- mat
  for (j in seq_len(ncol(mat))) {
    y <- mat[, j]
    if (all(is.finite(y)) && all(is.finite(design))) {
      resid[, j] <- stats::lm.fit(design, y)$residuals
    }
  }
  return(resid)
}


.profile_trait_pair <- function(out, profiles, trait_ids_by_program, program_index,
                                t_a, t_b, n_perm, min_shared_features,
                                exclude_target_traits, partial_out_common,
                                max_null, best) {
  progs_a <- names(trait_ids_by_program)[trait_ids_by_program == t_a]
  progs_b <- names(trait_ids_by_program)[trait_ids_by_program == t_b]
  features_a <- unique(as.character(
    profiles$feature_trait_id[as.character(profiles$program_id) %in% progs_a]
  ))
  features_b <- unique(as.character(
    profiles$feature_trait_id[as.character(profiles$program_id) %in% progs_b]
  ))
  shared_features <- intersect(features_a, features_b)
  if (exclude_target_traits) {
    shared_features <- setdiff(shared_features, c(t_a, t_b))
  }

  rows <- which(
    as.character(out$trait_id_a) == t_a & as.character(out$trait_id_b) == t_b
  )
  out$n_features_shared[rows] <- length(shared_features)
  out$n_features_union[rows] <- length(union(features_a, features_b))

  cand_rows <- rows[out$candidate[rows]]
  if (length(shared_features) < min_shared_features || length(cand_rows) == 0) {
    return(list(out = out, max_null = max_null, best = best))
  }

  mat_a <- .align_profile_matrix(profiles, progs_a, shared_features)
  mat_b <- .align_profile_matrix(profiles, progs_b, shared_features)
  raw_a <- mat_a
  raw_b <- mat_b
  if (partial_out_common) {
    mat_a$loading <- .partial_out_common(mat_a$loading)
    mat_b$loading <- .partial_out_common(mat_b$loading)
  }

  a_idx <- match(as.character(out$program_id_a[cand_rows]), progs_a)
  b_idx <- match(as.character(out$program_id_b[cand_rows]), progs_b)

  pair_cor <- function(ma, mb, ca, cb) {
    return(vapply(seq_along(a_idx), function(k) {
      return(.weighted_pearson(
        ma[, a_idx[k]], mb[, b_idx[k]], ca[, a_idx[k]] * cb[, b_idx[k]]
      ))
    }, numeric(1)))
  }

  obs_dec <- pair_cor(mat_a$loading, mat_b$loading, mat_a$conf, mat_b$conf)
  obs_raw <- pair_cor(raw_a$loading, raw_b$loading, raw_a$conf, raw_b$conf)
  out$profile_r[cand_rows] <- obs_raw
  out$profile_r_resid[cand_rows] <- obs_dec

  # Track each program's best partner by |deciding correlation|.
  for (k in seq_along(cand_rows)) {
    if (!is.finite(obs_dec[k])) next
    for (pid in c(out$program_id_a[cand_rows[k]], out$program_id_b[cand_rows[k]])) {
      pos <- match(pid, as.character(best$program_id))
      if (is.na(best$best_abs_r[pos]) || abs(obs_dec[k]) > best$best_abs_r[pos]) {
        best$best_abs_r[pos] <- abs(obs_dec[k])
        best$best_partner[pos] <- setdiff(
          c(out$program_id_a[cand_rows[k]], out$program_id_b[cand_rows[k]]), pid
        )[1]
        best$best_row[pos] <- cand_rows[k]
      }
    }
  }

  n_ge <- integer(length(cand_rows))
  for (p in seq_len(n_perm)) {
    dec_perm <- mat_b$loading
    conf_perm <- mat_b$conf
    for (j in seq_len(ncol(mat_b$loading))) {
      idx <- sample.int(nrow(mat_b$loading))
      dec_perm[, j] <- mat_b$loading[idx, j]
      conf_perm[, j] <- mat_b$conf[idx, j]
    }
    rp <- pair_cor(mat_a$loading, dec_perm, mat_a$conf, conf_perm)
    ge <- is.finite(rp) & is.finite(obs_dec) & abs(rp) >= abs(obs_dec)
    n_ge <- n_ge + as.integer(ge)

    for (pa in unique(a_idx)) {
      vals <- abs(rp[a_idx == pa])
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) next
      slot <- program_index[[progs_a[pa]]]
      max_null[p, slot] <- max(max_null[p, slot], max(vals))
    }
    for (pb in unique(b_idx)) {
      vals <- abs(rp[b_idx == pb])
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) next
      slot <- program_index[[progs_b[pb]]]
      max_null[p, slot] <- max(max_null[p, slot], max(vals))
    }
  }
  out$p_profile[cand_rows] <- (1 + n_ge) / (1 + n_perm)

  return(list(out = out, max_null = max_null, best = best))
}


.mutual_best_mask <- function(out, matches, trait_ids_by_program) {
  n <- nrow(out)
  shared <- rep(FALSE, n)
  if (n == 0) {
    return(shared)
  }
  sig <- stats::setNames(matches$significant, as.character(matches$program_id))
  partner <- stats::setNames(
    matches$best_partner, as.character(matches$program_id)
  )
  for (i in seq_len(n)) {
    a <- as.character(out$program_id_a[i])
    b <- as.character(out$program_id_b[i])
    if (isTRUE(sig[[a]]) && isTRUE(sig[[b]]) &&
          identical(partner[[a]], b) && identical(partner[[b]], a)) {
      shared[i] <- TRUE
    }
  }
  return(shared)
}


.normalise_program_data <- function(program_data) {
  if (is.data.frame(program_data)) {
    data_list <- list(list(
      loadings = program_data,
      profiles = .empty_program_profiles()
    ))
  } else if (is.list(program_data) &&
               "loadings" %in% names(program_data) &&
               is.data.frame(program_data$loadings)) {
    data_list <- list(program_data)
  } else if (is.list(program_data)) {
    data_list <- program_data
  } else {
    stop("program_data must be an extract_program_loadings() result or a list of them")
  }

  loadings <- dplyr::bind_rows(lapply(data_list, function(x) {
    if (is.data.frame(x)) {
      return(x)
    }
    return(x$loadings)
  }))
  profiles <- dplyr::bind_rows(lapply(data_list, function(x) {
    if (is.data.frame(x) || is.null(x$profiles)) {
      return(.empty_program_profiles())
    }
    return(x$profiles)
  }))

  return(list(
    loadings = .normalise_program_loadings(loadings),
    profiles = .normalise_program_profiles(profiles)
  ))
}


#' @title Build Multi-Trait Program Families
#' @description Turn a calibrated program-correspondence result into multi-trait
#' program families. Programs are nodes and correspondence edges are drawn from
#' [compare_program_pairs()].
#'
#' By default (`rule = "reciprocal_best_match"`) an edge is used only when it is
#' the mutual best match of two significant programs, i.e. the calibrated
#' `shared` flag, so families are one-to-one correspondences rather than
#' transitive chains. `rule = "components"` instead uses every locus-matched pair
#' whose pair-level `q_profile <= fdr_threshold`, reproducing the (transitive)
#' connected-components behaviour for comparison. Isolated programs form their
#' own size-1 family.
#' @param result Either a result from [compare_program_pairs()] or its `pairs`
#'   data.frame.
#' @param program_info Optional data.frame with `program_id` and `trait_id`
#'   (one row per program). Defaults to the programs appearing in `result`.
#' @param fdr_threshold FDR at or below which a `"components"` edge is kept.
#'   Defaults to `0.05`.
#' @param rule Edge rule: `"reciprocal_best_match"` (default) or `"components"`.
#' @return A list with:
#'   \itemize{
#'     \item nodes: one row per program with `program_id`, `trait_id`,
#'       `trait_name` (when available), and `family`
#'     \item families: one row per family with `family`, `n_programs`,
#'       `n_traits`, `traits`, `n_edges`, `mean_locus_concordance`,
#'       `median_locus_concordance`, `mean_profile_r`, `n_concordant`,
#'       `n_antagonistic`, and `n_shared_loci`
#'     \item graph: an `igraph` object (or `NULL` when there are no edges)
#'     \item edges: the pair table used as edges (after the chosen rule)
#'   }
#' @export
build_program_families <- function(result,
                                   program_info = NULL,
                                   fdr_threshold = 0.05,
                                   rule = c("reciprocal_best_match",
                                            "components")) {
  rule <- match.arg(rule)
  matches <- NULL
  if (is.list(result) && !is.data.frame(result) && !is.null(result$pairs)) {
    matches <- result$matches
    pairs <- result$pairs
  } else {
    pairs <- result
  }
  if (!is.data.frame(pairs) ||
        !all(c("jaccard_loci", "profile_r_resid", "shared") %in% names(pairs))) {
    stop("result must come from compare_program_pairs()")
  }

  if (rule == "reciprocal_best_match") {
    edges <- pairs[pairs$shared, , drop = FALSE]
  } else {
    edges <- pairs[
      pairs$candidate & is.finite(pairs$q_profile) &
        pairs$q_profile <= fdr_threshold,
      ,
      drop = FALSE
    ]
  }
  if (nrow(edges) > 0) {
    edges$locus_concordance <- edges$jaccard_loci
    edges$profile_r <- edges$profile_r_resid
  } else {
    edges$locus_concordance <- numeric(0)
    edges$profile_r <- numeric(0)
  }

  if (!is.null(matches)) {
    endpoints <- matches |>
      dplyr::distinct(program_id, trait_id, trait_name)
  } else {
    endpoints <- .normalise_program_info(edges, NULL)
  }
  program_info <- .merge_program_info(endpoints, program_info)

  if (nrow(edges) == 0 || nrow(program_info) == 0) {
    nodes <- program_info |>
      dplyr::mutate(family = paste0("family_", dplyr::row_number())) |>
      dplyr::select(program_id, trait_id, trait_name, family)
    return(list(
      nodes = nodes,
      families = .empty_families(),
      graph = NULL,
      edges = edges
    ))
  }

  graph <- igraph::graph_from_data_frame(
    edges[, c("program_id_a", "program_id_b"), drop = FALSE],
    directed = FALSE,
    vertices = data.frame(
      name = program_info$program_id,
      stringsAsFactors = FALSE
    )
  )
  graph <- igraph::set_edge_attr(
    graph, "locus_concordance", value = edges$locus_concordance
  )
  graph <- igraph::set_edge_attr(
    graph, "profile_r", value = edges$profile_r
  )
  graph <- igraph::set_edge_attr(
    graph, "weight", value = abs(edges$profile_r)
  )
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

  families <- .family_metrics(nodes, edges)
  return(list(
    nodes = nodes,
    families = families,
    graph = graph,
    edges = edges
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
  if (!"coloc_group_id" %in% names(loadings)) {
    loadings$coloc_group_id <- NA_character_
  } else {
    loadings$coloc_group_id <- as.character(loadings$coloc_group_id)
  }
  if (!"lfsr" %in% names(loadings)) {
    # Missing lFSR falls back to magnitude-only weighting.
    loadings$lfsr <- NA_real_
  } else {
    loadings$lfsr <- as.numeric(loadings$lfsr)
  }
  if (!"trait_name" %in% names(loadings)) {
    loadings$trait_name <- NA_character_
  }
  return(loadings)
}


.normalise_program_profiles <- function(profiles) {
  if (is.null(profiles) || !is.data.frame(profiles) || nrow(profiles) == 0) {
    return(.empty_program_profiles())
  }
  required <- c("program_id", "trait_id", "program", "feature_trait_id", "loading")
  missing <- setdiff(required, names(profiles))
  if (length(missing) > 0) {
    stop("profiles is missing columns: ", paste(missing, collapse = ", "))
  }
  profiles$program_id <- as.character(profiles$program_id)
  profiles$trait_id <- as.character(profiles$trait_id)
  profiles$feature_trait_id <- as.character(profiles$feature_trait_id)
  if (!"lfsr" %in% names(profiles)) {
    profiles$lfsr <- NA_real_
  } else {
    profiles$lfsr <- as.numeric(profiles$lfsr)
  }
  if (!"trait_name" %in% names(profiles)) {
    profiles$trait_name <- NA_character_
  }
  if (!"feature_trait_name" %in% names(profiles)) {
    profiles$feature_trait_name <- profiles$feature_trait_id
  }
  return(profiles)
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


# Confidence weight for one program at one SNP: (1 - lFSR) x |loading|.
# Missing lFSR is treated as full confidence (weight falls back to magnitude).
.program_weight_factor <- function(lfsr, loading) {
  return(.confidence_factor(lfsr) * abs(loading))
}


# Clamp (1 - lFSR) into [0, 1]; missing lFSR is treated as full confidence.
.confidence_factor <- function(lfsr) {
  conf <- 1 - lfsr
  conf[is.na(conf)] <- 1
  return(pmin(pmax(conf, 0), 1))
}


# Weighted Pearson correlation with per-observation weights w. Returns NA when
# the weighted variance is degenerate (e.g. all weight on one program's SNPs).
.weighted_pearson <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w)
  if (!all(ok)) {
    x <- x[ok]
    y <- y[ok]
    w <- w[ok]
  }
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    return(NA_real_)
  }
  xbar <- sum(w * x) / sw
  ybar <- sum(w * y) / sw
  dx <- x - xbar
  dy <- y - ybar
  vx <- sum(w * dx * dx)
  vy <- sum(w * dy * dy)
  if (!is.finite(vx) || !is.finite(vy) || vx <= 0 || vy <= 0) {
    return(NA_real_)
  }
  return(sum(w * dx * dy) / sqrt(vx * vy))
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
    profile_r <- edges$profile_r
    traits <- sort(unique(nodes$trait_id[nodes$family == fam]))
    return(data.frame(
      family = fam,
      n_programs = length(pids),
      n_traits = length(traits),
      traits = paste(traits, collapse = ", "),
      n_edges = nrow(edges),
      mean_locus_concordance = .safe_stat(edges$locus_concordance, mean),
      median_locus_concordance = .safe_stat(edges$locus_concordance, stats::median),
      mean_profile_r = .safe_stat(profile_r, mean),
      n_concordant = sum(profile_r > 0, na.rm = TRUE),
      n_antagonistic = sum(profile_r < 0, na.rm = TRUE),
      n_shared_loci = sum(edges$n_loci_shared, na.rm = TRUE),
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
    coloc_group_id = character(0),
    chr = character(0),
    bp = numeric(0),
    loading = numeric(0),
    abs_loading = numeric(0),
    lfsr = numeric(0),
    high_confidence = logical(0),
    stringsAsFactors = FALSE
  ))
}


.empty_program_profiles <- function() {
  return(data.frame(
    program_id = character(0),
    trait_id = character(0),
    trait_name = character(0),
    program = integer(0),
    feature_trait_id = character(0),
    feature_trait_name = character(0),
    loading = numeric(0),
    lfsr = numeric(0),
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
    n_features = integer(0),
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
    n_loci_a = integer(0),
    n_loci_b = integer(0),
    n_loci_shared = integer(0),
    n_loci_union = integer(0),
    jaccard_loci = numeric(0),
    p_locus = numeric(0),
    candidate = logical(0),
    n_features_shared = integer(0),
    n_features_union = integer(0),
    profile_r = numeric(0),
    profile_r_resid = numeric(0),
    p_profile = numeric(0),
    q_profile = numeric(0),
    direction = character(0),
    shared = logical(0),
    stringsAsFactors = FALSE
  ))
}


.empty_matches <- function() {
  return(data.frame(
    program_id = character(0),
    trait_id = character(0),
    trait_name = character(0),
    program = integer(0),
    best_partner = character(0),
    best_profile_r = numeric(0),
    best_locus_jaccard = numeric(0),
    best_n_loci_shared = integer(0),
    emp_p = numeric(0),
    q = numeric(0),
    significant = logical(0),
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
    mean_locus_concordance = numeric(0),
    median_locus_concordance = numeric(0),
    mean_profile_r = numeric(0),
    n_concordant = integer(0),
    n_antagonistic = integer(0),
    n_shared_loci = integer(0),
    stringsAsFactors = FALSE
  ))
}
