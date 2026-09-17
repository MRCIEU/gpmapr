#' @title Compare EBMF Programs Across Traits By Loading Concordance
#' @description Decide whether two programs discovered independently in
#' different target traits describe the same biology, using their continuous
#' SNP loadings on the **locus** axis.
#'
#' For a pair of programs the locus axis is the intersection of the two traits'
#' `coloc_group_id` universes — **every** shared-universe locus, with no lFSR or
#' magnitude threshold. Each locus contributes in proportion to how strongly it
#' loads in both programs:
#' \deqn{\mathrm{concordance} = \sum_l w_l f_a(l) f_b(l), \quad
#' w_l = (1 - \mathrm{lfsr}_a(l))(1 - \mathrm{lfsr}_b(l))}
#'
#' The score is deliberately **covariance-like rather than correlation-like**: a
#' fully normalised (cosine) score is scale-free and so cannot distinguish a pair
#' co-loading on three loci from a pair co-loading on thirty. `cosine_loci` is
#' reported alongside as a diagnostic so the difference is visible.
#'
#' Significance is calibrated in two stages: a per-pair permutation null (locus
#' labels of the second program are shuffled, carrying its confidence weights,
#' so the null preserves each program's loading distribution and sparsity), then
#' a per-program max-statistic null over that program's partners, BH-corrected
#' across programs so that searching many partners does not inflate it. Links are reported in two tiers:
#' \describe{
#'   \item{`primary`}{the mutual best match of two significant programs — the
#'   strict, one-to-one correspondence, flagged by `shared = TRUE`.}
#'   \item{`candidate`}{any other pair whose own `q_concordance` clears
#'   `fdr_threshold`. These allow one program to correspond to several in the
#'   other trait, which EBMF factor splitting makes genuinely possible, at the
#'   cost of not being calibrated for the search over partners.}
#' }
#'
#' Because [extract_program_loadings()] sign-orients each program to its target
#' trait, `sign(concordance)` is interpretable: positive is a concordant
#' (target-aligned) architecture, negative an antagonistic one.
#'
#' IMPORTANT: loci are treated as independent observations. The candidate loci
#' are fine-mapped colocalising signals, so LD between them is largely handled
#' upstream, but no correction is made for sample overlap between the two target
#' traits' GWAS — for trait pairs drawn from the same cohort the concordance (and
#' [module_rg()]) will be biased away from the null.
#' @param program_data An extraction result from [extract_program_loadings()],
#'   or a list of such results (one per target trait).
#' @param n_perm Number of permutation replicates. Defaults to `1000`.
#' @param fdr_threshold BH FDR at or below which a program is significant.
#'   Defaults to `0.05`.
#' @param seed RNG seed for the permutations.
#' @param min_shared_loci Minimum shared high-confidence loci for a pair to be
#'   flagged in the `candidate` column. Defaults to `0`: this is a reporting
#'   flag only and does not gate the test.
#' @param coloc_groups Optional named list of coloc-group data.frames, keyed by
#'   trait id. When supplied, [module_rg()] is run on every linked pair (both
#'   tiers) and returned as `module_rg`.
#' @return A list with:
#'   \itemize{
#'     \item pairs: one row per cross-trait program pair with the locus-overlap
#'       columns (`n_loci_a`, `n_loci_b`, `n_loci_shared`, `jaccard_loci`,
#'       `p_locus`), the concordance columns (`n_loci_axis`,
#'       `locus_concordance`, `cosine_loci`, `n_eff_loci`, `concordance_z`,
#'       `p_concordance`, `q_concordance`), `direction`, and the two link
#'       columns `shared` and `link_tier`.
#'     \item matches: one row per program with `best_partner`,
#'       `best_concordance_z`, `emp_p`, `q` and `significant`.
#'     \item module_rg: cross-trait effect correlation per `shared` pair when
#'       `coloc_groups` is supplied, otherwise an empty table.
#'     \item settings: settings used.
#'   }
#' @export
compare_program_pairs_loadings <- function(program_data,
                                           n_perm = 1000L,
                                           fdr_threshold = 0.05,
                                           seed = 1,
                                           min_shared_loci = 0L,
                                           coloc_groups = NULL) {
  pd <- .normalise_program_data(program_data)
  loadings <- pd$loadings
  keys <- .program_pair_keys(loadings)

  settings <- list(
    n_perm = as.integer(n_perm),
    fdr_threshold = fdr_threshold,
    seed = seed,
    min_shared_loci = as.integer(min_shared_loci),
    statistic = "locus_concordance"
  )

  if (nrow(keys) == 0) {
    return(list(
      pairs = .empty_concordance_pairs(),
      matches = .empty_matches(),
      module_rg = .empty_module_rg(),
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

  # Locus-overlap columns (n_loci_shared, jaccard_loci, p_locus) describe how
  # much the two programs' high-confidence memberships overlap. They are
  # reported for context; the test itself uses the full shared locus axis.
  locus_sets <- .program_locus_sets(loadings, "high_confidence", 25L)
  universe_by_trait <- .locus_universe_by_trait(loadings)
  out <- dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(i) {
    return(.locus_pair_stats(keys[i, , drop = FALSE], locus_sets, universe_by_trait))
  }))
  out$candidate <- !is.na(out$n_loci_shared) & out$n_loci_shared >= min_shared_loci
  out$n_loci_axis <- 0L
  out$locus_concordance <- NA_real_
  out$cosine_loci <- NA_real_
  out$n_eff_loci <- NA_real_
  out$concordance_z <- NA_real_
  out$p_concordance <- NA_real_

  max_null <- matrix(
    -Inf,
    nrow = as.integer(n_perm), ncol = nrow(progs),
    dimnames = list(NULL, as.character(progs$program_id))
  )
  best <- data.frame(
    program_id = as.character(progs$program_id),
    best_partner = NA_character_,
    best_abs_z = NA_real_,
    best_row = NA_integer_,
    stringsAsFactors = FALSE
  )

  trait_ids <- sort(unique(as.character(progs$trait_id)))
  set.seed(seed)
  if (length(trait_ids) >= 2) {
    for (tp in utils::combn(trait_ids, 2, simplify = FALSE)) {
      res <- .concordance_trait_pair(
        out = out, loadings = loadings,
        trait_ids_by_program = trait_ids_by_program,
        program_index = program_index,
        t_a = tp[[1]], t_b = tp[[2]], n_perm = as.integer(n_perm),
        max_null = max_null, best = best
      )
      out <- res$out
      max_null <- res$max_null
      best <- res$best
    }
  }

  obs_max <- stats::setNames(best$best_abs_z, as.character(best$program_id))
  emp_p <- stats::setNames(
    rep(NA_real_, nrow(progs)), as.character(progs$program_id)
  )
  for (pid in names(obs_max)[is.finite(obs_max)]) {
    null_vals <- max_null[, program_index[[pid]]]
    null_vals <- null_vals[is.finite(null_vals)]
    if (length(null_vals) > 0) {
      emp_p[[pid]] <- (1 + sum(null_vals >= obs_max[[pid]])) / (1 + length(null_vals))
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
      best_concordance_z = best$best_abs_z[match(
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

  out$q_concordance <- NA_real_
  ok_pair <- is.finite(out$p_concordance)
  if (any(ok_pair)) {
    out$q_concordance[ok_pair] <- stats::p.adjust(
      out$p_concordance[ok_pair], method = "BH"
    )
  }
  out$direction <- dplyr::case_when(
    !is.finite(out$concordance_z) ~ NA_character_,
    out$concordance_z > 0 ~ "concordant",
    out$concordance_z < 0 ~ "antagonistic",
    TRUE ~ NA_character_
  )
  out$shared <- .mutual_best_mask(out, matches, trait_ids_by_program)
  out$link_tier <- dplyr::case_when(
    out$shared ~ "primary",
    is.finite(out$q_concordance) & out$q_concordance <= fdr_threshold ~ "candidate",
    TRUE ~ NA_character_
  )

  rg <- .empty_module_rg()
  if (!is.null(coloc_groups)) {
    rg <- module_rg(
      result = list(pairs = out),
      program_data = pd,
      coloc_groups = coloc_groups,
      seed = seed
    )
  }

  return(list(pairs = out, matches = matches, module_rg = rg, settings = settings))
}


#' @title Module-Restricted Cross-Trait Effect Correlation
#' @description For linked cross-trait program pairs, correlate the two target
#' traits' own SNP effects across the module's loci — a module-level analogue of
#' a local genetic correlation, where the stratification unit is a data-driven
#' program rather than a genomic window.
#'
#' Loci are the pair's shared `coloc_group_id` axis; each locus is weighted by
#' how strongly the module loads on it in both programs and by the precision of
#' the two effect estimates:
#' \deqn{w_l = (1 - \mathrm{lfsr}_a)(1 - \mathrm{lfsr}_b)
#' |f_a(l) f_b(l)| / (se_a^2 + se_b^2)}
#' so loci the module does not load on drop out continuously rather than by a
#' threshold. A positive `rg` means the alleles that raise trait A also raise
#' trait B across this module.
#'
#' Effect alleles are checked by requiring the two traits' `variant_id` to match
#' at a locus; mismatches are dropped and counted in `n_allele_mismatch`.
#'
#' `rg_direction` and the `direction` reported by
#' [compare_program_pairs_loadings()] answer different questions and need not
#' agree. `direction` is a statement about *latent factor structure* — whether
#' the two programs' loading patterns align once each has been sign-anchored to
#' its own target trait. `rg_direction` is a statement about *observed effects* —
#' whether the alleles that raise one trait raise or lower the other across the
#' module's loci. When the two disagree, `rg_direction` is the more directly
#' interpretable claim; read it together with the CI, which shows how well
#' determined it is.
#'
#' IMPORTANT: shared loci entered each trait's universe by clearing that trait's
#' own significance and fine-mapping bar, so the estimate is conditional on that
#' double selection and effect sizes are subject to winner's curse. No correction
#' is made for sample overlap between the two GWAS.
#' @param result A [compare_program_pairs_loadings()] result, or its `pairs`
#'   data.frame.
#' @param program_data An extraction result from [extract_program_loadings()] or
#'   a list of them — needed for the per-locus loading weights.
#' @param coloc_groups Named list of coloc-group data.frames keyed by trait id.
#'   Each trait's own effect is taken from its rows where `trait_id` equals that
#'   trait.
#' @param tiers Which link tiers to estimate for: any of `"primary"` and
#'   `"candidate"`. Defaults to both. Pass `NULL` to estimate every pair.
#' @param n_boot Percentile-bootstrap replicates over loci. Defaults to `1000`.
#' @param seed RNG seed for the bootstrap.
#' @return A data.frame with one row per pair: `program_id_a`, `program_id_b`,
#'   `trait_id_a`, `trait_id_b`, `link_tier`, `rg`, `rg_ci_lower`,
#'   `rg_ci_upper`, `n_loci_rg`, `n_allele_mismatch`, and `rg_direction`.
#' @export
module_rg <- function(result,
                      program_data,
                      coloc_groups,
                      tiers = c("primary", "candidate"),
                      n_boot = 1000L,
                      min_loci_rg = 10L,
                      seed = 1) {
  pairs <- if (is.data.frame(result)) result else result$pairs
  if (is.null(pairs) || nrow(pairs) == 0) {
    return(.empty_module_rg())
  }
  if (!is.null(tiers) && "link_tier" %in% names(pairs)) {
    pairs <- pairs[!is.na(pairs$link_tier) & pairs$link_tier %in% tiers, , drop = FALSE]
  }
  if (nrow(pairs) == 0) {
    return(.empty_module_rg())
  }
  if (!is.list(coloc_groups) || is.data.frame(coloc_groups)) {
    stop("coloc_groups must be a named list of coloc-group data.frames")
  }

  loadings <- .normalise_program_data(program_data)$loadings
  betas <- lapply(names(coloc_groups), function(tid) {
    return(.trait_locus_effects(coloc_groups[[tid]], tid))
  })
  names(betas) <- names(coloc_groups)

  set.seed(seed)
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    ta <- as.character(pairs$trait_id_a[i])
    tb <- as.character(pairs$trait_id_b[i])
    pa <- as.character(pairs$program_id_a[i])
    pb <- as.character(pairs$program_id_b[i])
    tier <- if ("link_tier" %in% names(pairs)) {
      as.character(pairs$link_tier[i])
    } else {
      NA_character_
    }
    empty_row <- data.frame(
      program_id_a = pa, program_id_b = pb, trait_id_a = ta, trait_id_b = tb,
      link_tier = tier, rg = NA_real_, rg_ci_lower = NA_real_, rg_ci_upper = NA_real_,
      n_loci_rg = 0L, n_allele_mismatch = 0L, rg_direction = NA_character_,
      stringsAsFactors = FALSE
    )
    if (is.null(betas[[ta]]) || is.null(betas[[tb]])) {
      return(empty_row)
    }

    # Restricted to the loci the two programs actually claim, not the whole
    # trait-pair intersection. Weighting the full axis by |loading| left every
    # pair with the same n_loci_rg -- the same trait-level rg lightly
    # reweighted -- including pairs whose programs share a single locus.
    loci <- .module_locus_axis(loadings, pa, pb, ta, tb)
    if (length(loci) == 0) {
      return(empty_row)
    }
    va <- .program_locus_vector(loadings, pa, loci)
    vb <- .program_locus_vector(loadings, pb, loci)

    ea <- betas[[ta]][match(loci, betas[[ta]]$coloc_group_id), , drop = FALSE]
    eb <- betas[[tb]][match(loci, betas[[tb]]$coloc_group_id), , drop = FALSE]
    usable <- is.finite(ea$beta) & is.finite(eb$beta) &
      is.finite(ea$se) & is.finite(eb$se) & ea$se > 0 & eb$se > 0
    mismatch <- usable & !is.na(ea$variant_id) & !is.na(eb$variant_id) &
      as.character(ea$variant_id) != as.character(eb$variant_id)
    n_mismatch <- sum(mismatch, na.rm = TRUE)
    usable <- usable & !mismatch

    if (sum(usable, na.rm = TRUE) < min_loci_rg) {
      empty_row$n_allele_mismatch <- as.integer(n_mismatch)
      empty_row$n_loci_rg <- as.integer(sum(usable, na.rm = TRUE))
      return(empty_row)
    }

    ba <- ea$beta[usable]
    bb <- eb$beta[usable]
    w <- va$weight[usable] * vb$weight[usable] *
      abs(va$loading[usable] * vb$loading[usable]) /
      (ea$se[usable]^2 + eb$se[usable]^2)
    w[!is.finite(w) | w < 0] <- 0

    rg <- .weighted_pearson(ba, bb, w)
    ci <- c(NA_real_, NA_real_)
    if (is.finite(rg) && n_boot > 0) {
      n <- length(ba)
      draws <- vapply(seq_len(n_boot), function(b) {
        idx <- sample.int(n, n, replace = TRUE)
        return(.weighted_pearson(ba[idx], bb[idx], w[idx]))
      }, numeric(1))
      draws <- draws[is.finite(draws)]
      if (length(draws) > 0) {
        ci <- unname(stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE))
      }
    }

    return(data.frame(
      program_id_a = pa, program_id_b = pb, trait_id_a = ta, trait_id_b = tb,
      link_tier = tier, rg = rg, rg_ci_lower = ci[1], rg_ci_upper = ci[2],
      n_loci_rg = as.integer(sum(usable, na.rm = TRUE)),
      n_allele_mismatch = as.integer(n_mismatch),
      rg_direction = if (!is.finite(rg)) {
        NA_character_
      } else if (rg > 0) {
        "concordant"
      } else {
        "antagonistic"
      },
      stringsAsFactors = FALSE
    ))
  })

  out <- dplyr::bind_rows(rows)
  rownames(out) <- NULL
  return(out)
}


# All loci both traits carry a loading for: the intersection of the two traits'
# coloc-group universes, with no lFSR / magnitude gate.
.shared_locus_axis <- function(loadings, trait_a, trait_b) {
  keep <- !is.na(loadings$coloc_group_id)
  a <- unique(as.character(loadings$coloc_group_id[
    keep & as.character(loadings$trait_id) == trait_a
  ]))
  b <- unique(as.character(loadings$coloc_group_id[
    keep & as.character(loadings$trait_id) == trait_b
  ]))
  return(sort(intersect(a, b)))
}


# Signed loading and confidence weight of one program at each locus in `loci`.
# When several SNPs map to one coloc group the SNP with the largest |loading| is
# used, so the locus is represented by its strongest evidence in that program.
.program_locus_vector <- function(loadings, program_id, loci) {
  loading <- stats::setNames(rep(0, length(loci)), loci)
  weight <- stats::setNames(rep(0, length(loci)), loci)
  s <- loadings[
    as.character(loadings$program_id) == program_id &
      !is.na(loadings$coloc_group_id),
    ,
    drop = FALSE
  ]
  s <- s[as.character(s$coloc_group_id) %in% loci, , drop = FALSE]
  if (nrow(s) > 0) {
    s <- s[order(as.character(s$coloc_group_id), -abs(s$loading)), , drop = FALSE]
    s <- s[!duplicated(as.character(s$coloc_group_id)), , drop = FALSE]
    i <- match(as.character(s$coloc_group_id), loci)
    loading[i] <- s$loading
    weight[i] <- .confidence_factor(s$lfsr)
  }
  return(list(loading = unname(loading), weight = unname(weight)))
}


# loci x programs loading and confidence matrices for one trait.
.program_locus_matrix <- function(loadings, trait_id, loci, program_ids) {
  fm <- matrix(
    0, nrow = length(loci), ncol = length(program_ids),
    dimnames = list(loci, program_ids)
  )
  wm <- fm
  for (j in seq_along(program_ids)) {
    v <- .program_locus_vector(loadings, program_ids[j], loci)
    fm[, j] <- v$loading
    wm[, j] <- v$weight
  }
  return(list(loading = fm, conf = wm))
}


# One trait pair: locus concordance and its permutation calibration.
.concordance_trait_pair <- function(out, loadings, trait_ids_by_program,
                                    program_index, t_a, t_b, n_perm,
                                    max_null, best) {
  progs_a <- names(trait_ids_by_program)[trait_ids_by_program == t_a]
  progs_b <- names(trait_ids_by_program)[trait_ids_by_program == t_b]
  rows <- which(
    as.character(out$trait_id_a) == t_a & as.character(out$trait_id_b) == t_b
  )
  if (length(rows) == 0 || length(progs_a) == 0 || length(progs_b) == 0) {
    return(list(out = out, max_null = max_null, best = best))
  }

  loci <- .shared_locus_axis(loadings, t_a, t_b)
  out$n_loci_axis[rows] <- length(loci)
  if (length(loci) < 2) {
    return(list(out = out, max_null = max_null, best = best))
  }

  ma <- .program_locus_matrix(loadings, t_a, loci, progs_a)
  mb <- .program_locus_matrix(loadings, t_b, loci, progs_b)
  aw <- ma$loading * ma$conf
  bw <- mb$loading * mb$conf

  obs <- t(aw) %*% bw
  den_a <- t(ma$loading^2 * ma$conf) %*% mb$conf
  den_b <- t(ma$conf) %*% (mb$loading^2 * mb$conf)
  cosine <- obs / sqrt(den_a * den_b)

  # Effective number of loci carrying the score (participation ratio of the
  # per-locus contributions): how much independent support the pair rests on.
  s1 <- t(abs(aw)) %*% abs(bw)
  s2 <- t(aw^2) %*% bw^2
  n_eff <- s1^2 / s2
  n_eff[!is.finite(n_eff)] <- NA_real_

  null_scores <- array(NA_real_, c(n_perm, length(progs_a), length(progs_b)))
  for (p in seq_len(n_perm)) {
    bp <- bw
    for (j in seq_len(ncol(bw))) {
      bp[, j] <- bw[sample.int(nrow(bw)), j]
    }
    null_scores[p, , ] <- t(aw) %*% bp
  }
  null_mean <- apply(null_scores, c(2, 3), mean)
  null_sd <- apply(null_scores, c(2, 3), stats::sd)
  null_sd[!is.finite(null_sd) | null_sd <= 0] <- NA_real_

  zobs <- (obs - null_mean) / null_sd
  pval <- matrix(
    NA_real_, nrow = length(progs_a), ncol = length(progs_b),
    dimnames = dimnames(obs)
  )
  for (ia in seq_along(progs_a)) {
    for (ib in seq_along(progs_b)) {
      nv <- null_scores[, ia, ib]
      nv <- nv[is.finite(nv)]
      if (length(nv) > 0 && is.finite(obs[ia, ib])) {
        pval[ia, ib] <- (1 + sum(abs(nv) >= abs(obs[ia, ib]))) / (1 + length(nv))
      }
    }
  }

  # Standardise each permutation draw on its own pair's null so the per-program
  # maximum is taken over comparable quantities.
  znull <- sweep(null_scores, c(2, 3), null_mean, "-")
  znull <- sweep(znull, c(2, 3), null_sd, "/")

  a_idx <- match(as.character(out$program_id_a[rows]), progs_a)
  b_idx <- match(as.character(out$program_id_b[rows]), progs_b)
  for (k in seq_along(rows)) {
    ia <- a_idx[k]
    ib <- b_idx[k]
    if (is.na(ia) || is.na(ib)) next
    out$locus_concordance[rows[k]] <- obs[ia, ib]
    out$cosine_loci[rows[k]] <- cosine[ia, ib]
    out$n_eff_loci[rows[k]] <- n_eff[ia, ib]
    out$concordance_z[rows[k]] <- zobs[ia, ib]
    out$p_concordance[rows[k]] <- pval[ia, ib]

    if (!is.finite(zobs[ia, ib])) next
    for (pid in c(out$program_id_a[rows[k]], out$program_id_b[rows[k]])) {
      pos <- match(pid, as.character(best$program_id))
      if (is.na(best$best_abs_z[pos]) || abs(zobs[ia, ib]) > best$best_abs_z[pos]) {
        best$best_abs_z[pos] <- abs(zobs[ia, ib])
        best$best_partner[pos] <- setdiff(
          c(out$program_id_a[rows[k]], out$program_id_b[rows[k]]), pid
        )[1]
        best$best_row[pos] <- rows[k]
      }
    }
  }

  for (ia in seq_along(progs_a)) {
    vals <- apply(abs(znull[, ia, , drop = FALSE]), 1, max, na.rm = TRUE)
    slot <- program_index[[progs_a[ia]]]
    vals[!is.finite(vals)] <- -Inf
    max_null[, slot] <- pmax(max_null[, slot], vals)
  }
  for (ib in seq_along(progs_b)) {
    vals <- apply(abs(znull[, , ib, drop = FALSE]), 1, max, na.rm = TRUE)
    slot <- program_index[[progs_b[ib]]]
    vals[!is.finite(vals)] <- -Inf
    max_null[, slot] <- pmax(max_null[, slot], vals)
  }

  return(list(out = out, max_null = max_null, best = best))
}


# One row per coloc group with that trait's own effect estimate.
.trait_locus_effects <- function(coloc_groups, trait_id) {
  required <- c("coloc_group_id", "trait_id", "beta", "se")
  if (!is.data.frame(coloc_groups) ||
        length(setdiff(required, names(coloc_groups))) > 0) {
    return(NULL)
  }
  out <- coloc_groups[
    !is.na(coloc_groups$trait_id) &
      as.character(coloc_groups$trait_id) == as.character(trait_id) &
      !is.na(coloc_groups$coloc_group_id),
    ,
    drop = FALSE
  ]
  if (nrow(out) == 0) {
    return(NULL)
  }
  variant <- if ("variant_id" %in% names(out)) {
    as.character(out$variant_id)
  } else {
    NA_character_
  }
  out <- data.frame(
    coloc_group_id = as.character(out$coloc_group_id),
    variant_id = variant,
    beta = as.numeric(out$beta),
    se = as.numeric(out$se),
    stringsAsFactors = FALSE
  )
  return(out[!duplicated(out$coloc_group_id), , drop = FALSE])
}


.empty_concordance_pairs <- function() {
  out <- .empty_pairs()
  out$n_loci_axis <- integer(0)
  out$link_tier <- character(0)
  out$locus_concordance <- numeric(0)
  out$cosine_loci <- numeric(0)
  out$n_eff_loci <- numeric(0)
  out$concordance_z <- numeric(0)
  out$p_concordance <- numeric(0)
  out$q_concordance <- numeric(0)
  return(out)
}


.empty_module_rg <- function() {
  return(data.frame(
    program_id_a = character(0),
    program_id_b = character(0),
    trait_id_a = character(0),
    trait_id_b = character(0),
    link_tier = character(0),
    rg = numeric(0),
    rg_ci_lower = numeric(0),
    rg_ci_upper = numeric(0),
    n_loci_rg = integer(0),
    n_allele_mismatch = integer(0),
    rg_direction = character(0),
    stringsAsFactors = FALSE
  ))
}


#' @title Shared Feature Universe Between Traits
#' @description Count the background studies that two traits' EBMF fits have in
#' common. This is the axis `compare_program_pairs_profiles()` scores on, and
#' it is deliberately reported before any comparison: the whole reason to match
#' programs on their trait profile rather than on shared loci is that the
#' feature axis should be far wider than the locus intersection. Where it is
#' not, profile matching has no more support than locus matching and the result
#' should be read accordingly.
#' @param program_data An `extract_program_loadings()` result or a list of them.
#' @return A dataframe with one row per trait pair: `trait_id_a`, `trait_id_b`,
#'   `n_features_a`, `n_features_b`, `n_features_shared`, `frac_shared`, and
#'   `n_loci_shared` for the locus axis alongside it.
#' @export
shared_feature_universe <- function(program_data) {
  pd <- .normalise_program_data(program_data)
  profiles <- pd$profiles
  loadings <- pd$loadings
  empty <- data.frame(
    trait_id_a = character(0), trait_id_b = character(0),
    n_features_a = integer(0), n_features_b = integer(0),
    n_features_shared = integer(0), frac_shared = numeric(0),
    n_loci_shared = integer(0), stringsAsFactors = FALSE
  )
  if (is.null(profiles) || nrow(profiles) == 0) {
    return(empty)
  }
  trait_ids <- sort(unique(as.character(profiles$trait_id)))
  if (length(trait_ids) < 2) {
    return(empty)
  }
  features_by_trait <- lapply(trait_ids, function(t) {
    unique(as.character(
      profiles$feature_trait_id[as.character(profiles$trait_id) == t]
    ))
  })
  names(features_by_trait) <- trait_ids

  idx <- utils::combn(length(trait_ids), 2)
  rows <- lapply(seq_len(ncol(idx)), function(i) {
    ta <- trait_ids[idx[1, i]]
    tb <- trait_ids[idx[2, i]]
    fa <- features_by_trait[[ta]]
    fb <- features_by_trait[[tb]]
    shared <- intersect(fa, fb)
    return(data.frame(
      trait_id_a = ta, trait_id_b = tb,
      n_features_a = length(fa), n_features_b = length(fb),
      n_features_shared = length(shared),
      frac_shared = if (length(union(fa, fb)) > 0) {
        length(shared) / length(union(fa, fb))
      } else {
        NA_real_
      },
      n_loci_shared = length(.shared_locus_axis(loadings, ta, tb)),
      stringsAsFactors = FALSE
    ))
  })
  return(dplyr::bind_rows(rows))
}


#' @title Compare EBMF Programs Across Traits On Their Trait Profiles
#' @description Score every cross-trait pair of programs on how similarly they
#' load across the background studies the two traits have in common, and
#' calibrate that against a permutation null.
#'
#' This is the sibling of `compare_program_pairs_loadings()`, which scores pairs
#' on shared *loci*. That statistic is an inner product over the loci both
#' traits carry, so two programs describing the same biology at different loci
#' score near zero by construction — and when two traits share few loci (BMI and
#' height share 73 of 3456) every locus-axis result rests on very little.
#' Matching on the trait profile drops the requirement that the programs act at
#' the same variants: two programs correspond when the same background studies
#' load on them, whatever loci carry that signal in each trait.
#'
#' Read `n_eff_features` alongside `concordance_z` exactly as `n_eff_loci` is
#' read for the locus axis: a large z resting on two or three effective
#' features is thin.
#' @param program_data An `extract_program_loadings()` result or a list of them.
#'   Requires the `profiles` component, which `extract_program_loadings()`
#'   populates from the EBMF fit's `L_pm`.
#' @param n_perm Permutations for the null. Defaults to `1000`.
#' @param fdr_threshold BH threshold defining the candidate tier. Defaults to
#'   `0.05`.
#' @param seed RNG seed.
#' @return A list with `pairs`, `matches` and `settings`, matching the shape of
#'   `compare_program_pairs_loadings()` but with feature-axis columns
#'   (`n_features_axis`, `profile_concordance`, `cosine_features`,
#'   `n_eff_features`).
#' @export
compare_program_pairs_profiles <- function(program_data,
                                           n_perm = 1000L,
                                           fdr_threshold = 0.05,
                                           seed = 1) {
  pd <- .normalise_program_data(program_data)
  profiles <- pd$profiles
  settings <- list(
    n_perm = as.integer(n_perm),
    fdr_threshold = fdr_threshold,
    seed = seed,
    statistic = "profile_concordance"
  )
  empty_pairs <- data.frame(
    program_id_a = character(0), program_id_b = character(0),
    trait_id_a = character(0), trait_id_b = character(0),
    n_features_axis = integer(0), profile_concordance = numeric(0),
    cosine_features = numeric(0), n_eff_features = numeric(0),
    concordance_z = numeric(0), p_concordance = numeric(0),
    q_concordance = numeric(0), shared = logical(0),
    link_tier = character(0), direction = character(0),
    stringsAsFactors = FALSE
  )
  if (is.null(profiles) || nrow(profiles) == 0) {
    return(list(pairs = empty_pairs, matches = .empty_matches(),
                settings = settings))
  }
  keys <- .program_pair_keys(profiles)
  if (nrow(keys) == 0) {
    return(list(pairs = empty_pairs, matches = .empty_matches(),
                settings = settings))
  }

  progs <- dplyr::distinct(profiles, program_id, trait_id, program)
  trait_ids_by_program <- stats::setNames(
    as.character(progs$trait_id), as.character(progs$program_id)
  )

  out <- keys
  out$n_features_axis <- 0L
  out$profile_concordance <- NA_real_
  out$cosine_features <- NA_real_
  out$n_eff_features <- NA_real_
  out$concordance_z <- NA_real_
  out$p_concordance <- NA_real_

  max_null <- matrix(
    -Inf, nrow = as.integer(n_perm), ncol = nrow(progs),
    dimnames = list(NULL, as.character(progs$program_id))
  )
  best <- data.frame(
    program_id = as.character(progs$program_id),
    best_partner = NA_character_, best_abs_z = NA_real_,
    stringsAsFactors = FALSE
  )

  trait_ids <- sort(unique(as.character(progs$trait_id)))
  set.seed(seed)
  for (i in seq_along(trait_ids)) {
    for (j in seq_along(trait_ids)) {
      if (j <= i) {
        next
      }
      res <- .profile_trait_pair(
        out, profiles, trait_ids_by_program, trait_ids[i], trait_ids[j],
        as.integer(n_perm), max_null, best
      )
      out <- res$out
      max_null <- res$max_null
      best <- res$best
    }
  }

  out$q_concordance <- stats::p.adjust(out$p_concordance, method = "BH")
  matches <- .profile_matches(best, max_null, fdr_threshold)
  out <- .profile_link_tiers(out, matches, fdr_threshold)
  return(list(pairs = out, matches = matches, settings = settings))
}


# One trait pair on the feature axis. Mirrors .concordance_trait_pair(): the
# statistic is the confidence-weighted inner product of the two programs'
# signed feature loadings, and the null permutes feature labels within each
# program of trait b.
.profile_trait_pair <- function(out, profiles, trait_ids_by_program,
                                t_a, t_b, n_perm, max_null, best) {
  progs_a <- names(trait_ids_by_program)[trait_ids_by_program == t_a]
  progs_b <- names(trait_ids_by_program)[trait_ids_by_program == t_b]
  rows <- which(
    as.character(out$trait_id_a) == t_a & as.character(out$trait_id_b) == t_b
  )
  if (length(rows) == 0 || length(progs_a) == 0 || length(progs_b) == 0) {
    return(list(out = out, max_null = max_null, best = best))
  }

  features <- .shared_feature_axis(profiles, t_a, t_b)
  out$n_features_axis[rows] <- length(features)
  if (length(features) < 2) {
    return(list(out = out, max_null = max_null, best = best))
  }

  ma <- .program_feature_matrix(profiles, features, progs_a)
  mb <- .program_feature_matrix(profiles, features, progs_b)
  aw <- ma$loading * ma$conf
  bw <- mb$loading * mb$conf

  obs <- t(aw) %*% bw
  den_a <- t(ma$loading^2 * ma$conf) %*% mb$conf
  den_b <- t(ma$conf) %*% (mb$loading^2 * mb$conf)
  cosine <- obs / sqrt(den_a * den_b)

  s1 <- t(abs(aw)) %*% abs(bw)
  s2 <- t(aw^2) %*% bw^2
  n_eff <- s1^2 / s2
  n_eff[!is.finite(n_eff)] <- NA_real_

  null_scores <- array(NA_real_, c(n_perm, length(progs_a), length(progs_b)))
  for (p in seq_len(n_perm)) {
    bp <- bw
    for (k in seq_len(ncol(bw))) {
      bp[, k] <- bw[sample.int(nrow(bw)), k]
    }
    null_scores[p, , ] <- t(aw) %*% bp
  }
  null_mean <- apply(null_scores, c(2, 3), mean)
  null_sd <- apply(null_scores, c(2, 3), stats::sd)
  null_sd[!is.finite(null_sd) | null_sd <= 0] <- NA_real_
  zobs <- (obs - null_mean) / null_sd

  pval <- matrix(NA_real_, nrow = length(progs_a), ncol = length(progs_b))
  for (ia in seq_along(progs_a)) {
    for (ib in seq_along(progs_b)) {
      nv <- null_scores[, ia, ib]
      nv <- nv[is.finite(nv)]
      if (length(nv) > 0 && is.finite(obs[ia, ib])) {
        pval[ia, ib] <- (1 + sum(abs(nv) >= abs(obs[ia, ib]))) / (1 + length(nv))
      }
    }
  }

  znull <- sweep(null_scores, c(2, 3), null_mean, "-")
  znull <- sweep(znull, c(2, 3), null_sd, "/")

  a_idx <- match(as.character(out$program_id_a[rows]), progs_a)
  b_idx <- match(as.character(out$program_id_b[rows]), progs_b)
  for (k in seq_along(rows)) {
    ia <- a_idx[k]
    ib <- b_idx[k]
    if (is.na(ia) || is.na(ib)) {
      next
    }
    out$profile_concordance[rows[k]] <- obs[ia, ib]
    out$cosine_features[rows[k]] <- cosine[ia, ib]
    out$n_eff_features[rows[k]] <- n_eff[ia, ib]
    out$concordance_z[rows[k]] <- zobs[ia, ib]
    out$p_concordance[rows[k]] <- pval[ia, ib]
  }

  # Max-statistic null per program, so a program's best partner is judged
  # against the best it would have found by chance across all candidates.
  for (ia in seq_along(progs_a)) {
    pid <- progs_a[ia]
    max_null[, pid] <- pmax(
      max_null[, pid],
      apply(abs(znull[, ia, , drop = FALSE]), 1, max, na.rm = TRUE)
    )
    obs_row <- abs(zobs[ia, ])
    if (any(is.finite(obs_row))) {
      w <- which.max(replace(obs_row, !is.finite(obs_row), -Inf))
      br <- which(best$program_id == pid)
      if (!is.finite(best$best_abs_z[br]) || obs_row[w] > best$best_abs_z[br]) {
        best$best_abs_z[br] <- obs_row[w]
        best$best_partner[br] <- progs_b[w]
      }
    }
  }
  for (ib in seq_along(progs_b)) {
    pid <- progs_b[ib]
    max_null[, pid] <- pmax(
      max_null[, pid],
      apply(abs(znull[, , ib, drop = FALSE]), 1, max, na.rm = TRUE)
    )
    obs_col <- abs(zobs[, ib])
    if (any(is.finite(obs_col))) {
      w <- which.max(replace(obs_col, !is.finite(obs_col), -Inf))
      br <- which(best$program_id == pid)
      if (!is.finite(best$best_abs_z[br]) || obs_col[w] > best$best_abs_z[br]) {
        best$best_abs_z[br] <- obs_col[w]
        best$best_partner[br] <- progs_a[w]
      }
    }
  }

  return(list(out = out, max_null = max_null, best = best))
}


# Features (background studies) both traits' fits carry.
.shared_feature_axis <- function(profiles, trait_a, trait_b) {
  keep <- !is.na(profiles$feature_trait_id)
  a <- unique(as.character(profiles$feature_trait_id[
    keep & as.character(profiles$trait_id) == trait_a
  ]))
  b <- unique(as.character(profiles$feature_trait_id[
    keep & as.character(profiles$trait_id) == trait_b
  ]))
  return(sort(intersect(a, b)))
}


# features x programs matrices of signed loading and confidence weight.
.program_feature_matrix <- function(profiles, features, program_ids) {
  loading <- matrix(
    0, nrow = length(features), ncol = length(program_ids),
    dimnames = list(features, program_ids)
  )
  conf <- loading
  for (j in seq_along(program_ids)) {
    s <- profiles[
      as.character(profiles$program_id) == program_ids[j] &
        !is.na(profiles$feature_trait_id),
      ,
      drop = FALSE
    ]
    if (nrow(s) == 0) {
      next
    }
    i <- match(as.character(s$feature_trait_id), features)
    ok <- !is.na(i)
    loading[i[ok], j] <- s$loading[ok]
    conf[i[ok], j] <- .confidence_factor(s$lfsr[ok])
  }
  return(list(loading = loading, conf = conf))
}


.profile_matches <- function(best, max_null, fdr_threshold) {
  emp_p <- vapply(best$program_id, function(pid) {
    nv <- max_null[, pid]
    nv <- nv[is.finite(nv)]
    obs <- best$best_abs_z[best$program_id == pid]
    if (length(nv) == 0 || !is.finite(obs)) {
      return(NA_real_)
    }
    return((1 + sum(nv >= obs)) / (1 + length(nv)))
  }, numeric(1))
  out <- data.frame(
    program_id = best$program_id,
    best_partner = best$best_partner,
    best_concordance_z = best$best_abs_z,
    emp_p = as.numeric(emp_p),
    stringsAsFactors = FALSE
  )
  out$q <- stats::p.adjust(out$emp_p, method = "BH")
  out$significant <- !is.na(out$q) & out$q < fdr_threshold
  return(out)
}


.profile_link_tiers <- function(out, matches, fdr_threshold) {
  sig <- stats::setNames(matches$significant, matches$program_id)
  partner <- stats::setNames(matches$best_partner, matches$program_id)
  a <- as.character(out$program_id_a)
  b <- as.character(out$program_id_b)
  reciprocal <- !is.na(partner[a]) & !is.na(partner[b]) &
    partner[a] == b & partner[b] == a
  both_sig <- !is.na(sig[a]) & !is.na(sig[b]) & sig[a] & sig[b]
  out$shared <- reciprocal & both_sig
  out$link_tier <- ifelse(
    out$shared, "primary",
    ifelse(!is.na(out$q_concordance) & out$q_concordance < fdr_threshold,
           "candidate", NA_character_)
  )
  out$direction <- ifelse(
    is.na(out$concordance_z), NA_character_,
    ifelse(out$concordance_z >= 0, "concordant", "antagonistic")
  )
  return(out[!is.na(out$link_tier) | !is.na(out$concordance_z), , drop = FALSE])
}


# Loci the two programs themselves claim, intersected with the loci both traits
# carry. This is what makes module_rg a statement about the module rather than
# about the trait pair: .shared_locus_axis() returns every locus the two traits
# have in common, so weighting it by |loading| still leaves every pair using the
# same loci and reports the same n_loci_rg for all of them.
.module_locus_axis <- function(loadings, program_a, program_b,
                               trait_a, trait_b) {
  shared <- .shared_locus_axis(loadings, trait_a, trait_b)
  if (length(shared) == 0) {
    return(character(0))
  }
  claimed <- function(pid) {
    s <- loadings[
      as.character(loadings$program_id) == pid &
        !is.na(loadings$coloc_group_id),
      ,
      drop = FALSE
    ]
    if (nrow(s) == 0) {
      return(character(0))
    }
    if ("high_confidence" %in% names(s)) {
      hc <- s$high_confidence
      hc[is.na(hc)] <- FALSE
      if (any(hc)) {
        s <- s[hc, , drop = FALSE]
      }
    }
    return(unique(as.character(s$coloc_group_id)))
  }
  own <- union(claimed(program_a), claimed(program_b))
  return(sort(intersect(shared, own)))
}
