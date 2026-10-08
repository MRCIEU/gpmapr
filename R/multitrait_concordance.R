#' @title Compare EBMF Programs Across Traits
#' @description Decide whether two programs discovered independently in
#' different target traits describe the same biology. Two questions are asked
#' of every cross-trait pair, on the two axes of the EBMF fit:
#'
#' \describe{
#'   \item{Same biology? -- profile congruence (`phi_traits`).}{Tucker's
#'   congruence coefficient (an uncentred correlation) between the two
#'   programs' **trait loadings** (`L_pm`, the background traits that drive
#'   them), over the background traits both fits carry. Each trait's loading
#'   is weighted by its confidence, \eqn{\tilde l = l\, w(\mathrm{lfsr})}, with
#'   no magnitude cutoff:
#'   \deqn{\phi = \frac{\sum_t \tilde l_a(t) \tilde l_b(t)}{\sqrt{\sum_t
#'   \tilde l_a(t)^2 \sum_t \tilde l_b(t)^2}}.}
#'   The comparable traits leave out both anchor (target) traits -- each is
#'   positive by construction in its own fit but an ordinary row in the other
#'   -- and any trait excluded from either fit as a proxy (`exclude_features`).
#'   `phi_traits` is signed: each trait's SNPs are oriented to its own risk
#'   allele, so two programs carrying the same biology in opposite directions
#'   have mirrored profiles and \eqn{\phi < 0}. `n_traits_shared` counts the
#'   comparable traits and `n_eff_traits` how many effectively carry the
#'   agreement (participation ratio of \eqn{|\tilde l_a \tilde l_b|}); a small
#'   value means a few traits dominate.}
#'   \item{Same SNPs? -- locus alignment (`alignment`) and comparability
#'   (`comparability`).}{The squared uncentred correlation of the two
#'   programs' **SNP loadings** (`F_pm`, weighted the same way) over the loci
#'   both traits carry, \eqn{A = (\sum_U \tilde f_a \tilde f_b)^2 / (\sum_U
#'   \tilde f_a^2 \sum_U \tilde f_b^2)}. Loci the other trait never measured
#'   cannot agree or disagree, so they are left out of `alignment` and
#'   reported instead as `comparability`, the share of each program's loading
#'   on the shared loci, multiplied: \eqn{C = (\sum_U \tilde f_a^2 / \sum_a
#'   \tilde f_a^2)(\sum_U \tilde f_b^2 / \sum_b \tilde f_b^2)}.}
#' }
#'
#' **Profile strength** (`profile_strength`): `strong` when
#' \eqn{|\phi| \ge} `phi_high`, `moderate` when `phi_low` \eqn{\le |\phi| <}
#' `phi_high`, `weak` otherwise.
#'
#' **Locus status** (`locus_status`): \describe{
#'   \item{`too_few_loci`}{the pair shares fewer than `min_shared_loci`
#'   high-confidence loci (`n_loci_shared`), or fewer than `min_eff_loci` loci
#'   effectively carry its locus score (`n_eff_loci`): the loci cannot settle
#'   whether the same SNPs carry the two programs.}
#'   \item{`agree`}{testable, `alignment >= min_alignment` and the permutation
#'   test passes (`q_concordance <= fdr_threshold`, below).}
#'   \item{`disagree`}{testable, but `alignment` or `q_concordance` fails.}
#' }
#'
#' **Tiers.** \describe{
#'   \item{`high`}{strong profile, loci agree: same biology, carried by the
#'   same SNPs.}
#'   \item{`medium`}{strong profile, too few loci to test; or moderate
#'   profile, loci agree.}
#'   \item{`low`}{strong profile, loci disagree (same background biology, but
#'   the shared SNPs do not carry it); moderate profile, too few loci; or weak
#'   profile, loci agree (shared loci, different biology, typical of large
#'   pleiotropic loci).}
#'   \item{`none`}{otherwise, including a moderate profile whose loci
#'   disagree.}
#' }
#' `linked` pairs are those of high or medium tier; [build_program_families()]
#' groups them into families. The default cutoffs follow Tucker's conventions
#' for factor similarity (0.85 "fair similarity", 0.95 "equal"); they are
#' judgements, not calibrated tests. [profile_contributions()] names the
#' background traits that carry each pair's \eqn{\phi}.
#'
#' **Locus permutation test.** For the locus axis, every pair is also scored by
#' its loading-weighted concordance on the shared loci either program claims,
#' \eqn{\sum_l \tilde f_a(l) \tilde f_b(l)}, against a two-sided permutation
#' null that shuffles the second program's loadings over the loci the two
#' traits share (`concordance_z`, `p_concordance`). `q_concordance` is the
#' Benjamini-Hochberg adjustment across all scored pairs (eligibility is not
#' independent of the statistic, so correcting over eligible pairs alone would
#' be anti-conservative). The sign of `concordance_z` is not a direction:
#' both programs are oriented to their own anchor, so a shared module loads
#' positively in both whichever way the traits move. Direction comes from
#' the sign of `phi_traits`, checked against [module_rg()]
#' (`direction_check`); `rg` never changes a tier.
#'
#' **Confidence weights.** \eqn{w(\mathrm{lfsr}) = (1 - 2\,\mathrm{lfsr})_+} by
#' default -- 0 for a coin-flip sign, 1 for a certain one -- set by
#' `options(gpmapr.lfsr_weight = "two_sided")`; `"one_sided"` gives
#' \eqn{1 - \mathrm{lfsr}}.
#'
#' IMPORTANT: loci are treated as independent observations, and no correction
#' is made for sample overlap between the two target traits' GWAS.
#' @param program_data An extraction result from [extract_program_loadings()],
#'   or a list of such results (one per target trait). Needs `profiles` for
#'   the profile congruence.
#' @param n_perm Number of permutation replicates for the locus test. Defaults
#'   to `1000`.
#' @param fdr_threshold BH FDR for the locus permutation test; a testable
#'   pair needs `q_concordance <= fdr_threshold` for its loci to `agree`.
#'   Defaults to `0.05`.
#' @param seed RNG seed for the permutations.
#' @param min_shared_loci Minimum shared high-confidence loci (`n_loci_shared`)
#'   for the loci to be testable; below it `locus_status` is `too_few_loci`.
#'   Defaults to `3`.
#' @param exclude_features Optional named list, keyed by trait id, of feature
#'   (background trait) ids excluded from that trait's fit as proxies; they
#'   are left out of every profile comparison involving either trait. The
#'   anchor traits themselves are always left out.
#' @param phi_high,phi_low Profile-congruence cutoffs for the tiers. Default
#'   `0.85` and `0.70`.
#' @param min_alignment Minimum locus `alignment` for a testable pair's loci
#'   to `agree`. Defaults to `0.5`.
#' @param min_eff_loci Minimum `n_eff_loci` for the loci to be testable; below
#'   it `locus_status` is `too_few_loci`. Defaults to `5`.
#' @param coloc_groups Optional named list of coloc-group data.frames, keyed by
#'   trait id. When supplied, [module_rg()] is run on every linked pair and
#'   returned as `module_rg`. Also required to build the LD locus bridge (see
#'   `min_ld_r`) for any pair involving an uploaded GWAS trait.
#' @param min_ld_r For a pair where at least one trait is an uploaded GWAS
#'   (its own coloc_group_id numbering is private to that upload, not shared
#'   with other traits -- see `.needs_locus_bridge()`), loci are matched first
#'   by the upload's own coloc record of the other trait
#'   (`existing_study_extraction_id`), then by exact `variant_id`, then by
#'   `ld_proxies()` LD looked up from both traits' leads. `min_ld_r` is the
#'   minimum `|r|` for that last step. GPMap's proxy table only holds pairs at
#'   `r^2 >= 0.8` (`|r|` about 0.89), so values below that have no effect.
#'   Defaults to `0.5`. Ignored for pairs where neither trait is an upload.
#' @return A list with:
#'   \itemize{
#'     \item pairs: one row per cross-trait program pair with the profile
#'       columns (`phi_traits`, `n_traits_shared`, `n_eff_traits`,
#'       `direction_profile`), the locus columns (`n_loci_a`, `n_loci_b`,
#'       `n_loci_shared`, `jaccard_loci`, `p_locus`, `n_loci_axis`,
#'       `locus_concordance`, `alignment`, `comparability`, `n_eff_loci`,
#'       `concordance_z`, `p_concordance`, `q_concordance`, `eligible`),
#'       and `profile_strength`, `locus_status`, `tier` and `linked`.
#'     \item programs: one row per program with `program_id`, `trait_id`,
#'       `trait_name`, `program` and `n_links`.
#'     \item module_rg: cross-trait effect correlation per linked pair when
#'       `coloc_groups` is supplied, otherwise an empty table.
#'     \item bridges: the per-trait-pair locus bridges built for uploaded GWAS
#'       traits (a named list, empty when none were needed), for reuse by
#'       [build_program_families()] and [family_module_rg()].
#'     \item settings: settings used, including `exclude_features`.
#'   }
#' @export
compare_program_pairs_loadings <- function(program_data,
                                           n_perm = 1000L,
                                           fdr_threshold = 0.05,
                                           seed = 1,
                                           min_shared_loci = 3L,
                                           exclude_features = NULL,
                                           phi_high = 0.85,
                                           phi_low = 0.70,
                                           min_alignment = 0.5,
                                           min_eff_loci = 5,
                                           coloc_groups = NULL,
                                           min_ld_r = 0.5,
                                           ld_proxies_fn = ld_proxies) {
  pd <- .normalise_program_data(program_data)
  loadings <- pd$loadings
  keys <- .program_pair_keys(loadings)

  settings <- list(
    n_perm = as.integer(n_perm),
    fdr_threshold = fdr_threshold,
    seed = seed,
    min_shared_loci = as.integer(min_shared_loci),
    phi_high = phi_high,
    phi_low = phi_low,
    min_alignment = min_alignment,
    min_eff_loci = min_eff_loci,
    exclude_features = exclude_features,
    lfsr_weight = .lfsr_weight_method(),
    statistic = "profile_congruence + locus_alignment"
  )

  progs <- dplyr::distinct(loadings, program_id, trait_id, trait_name, program)
  if (nrow(keys) == 0) {
    return(list(
      pairs = .empty_concordance_pairs(),
      programs = .program_link_counts(progs, .empty_concordance_pairs()),
      module_rg = .empty_module_rg(),
      bridges = list(),
      settings = settings
    ))
  }

  trait_ids_by_program <- stats::setNames(
    as.character(progs$trait_id), as.character(progs$program_id)
  )
  trait_ids <- sort(unique(as.character(progs$trait_id)))

  # One LD-based locus bridge per trait pair that needs it (an uploaded GWAS
  # on either side -- see .needs_locus_bridge()), built once here and reused
  # below so every raw coloc_group_id comparison for that pair -- eligibility
  # (.locus_pair_stats()), the locus scores (.concordance_trait_pair()) and
  # module_rg() -- agrees. Pairs that don't need it are absent from
  # `bridges`, which is a no-op everywhere it's read.
  bridges <- list()
  if (!is.null(coloc_groups) && length(trait_ids) >= 2) {
    for (tp in utils::combn(trait_ids, 2, simplify = FALSE)) {
      ta <- tp[[1]]
      tb <- tp[[2]]
      if (!.needs_locus_bridge(ta, tb)) {
        next
      }
      bridges[[.locus_bridge_key(ta, tb)]] <- .build_locus_bridge(
        coloc_groups[[ta]], coloc_groups[[tb]], ta, tb,
        min_r = min_ld_r, ld_proxies_fn = ld_proxies_fn
      )
    }
  }

  # Locus-overlap columns (n_loci_shared, jaccard_loci, p_locus) describe how
  # much the two programs' high-confidence memberships overlap; n_loci_shared
  # decides whether the locus evidence can count at all (`eligible`).
  locus_sets <- .program_locus_sets(loadings, "high_confidence", 25L)
  universe_by_trait <- .locus_universe_by_trait(loadings)
  out <- dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(i) {
    return(.locus_pair_stats(
      keys[i, , drop = FALSE], locus_sets, universe_by_trait, bridges = bridges
    ))
  }))
  out$eligible <- !is.na(out$n_loci_shared) & out$n_loci_shared >= min_shared_loci
  out$n_loci_axis <- 0L
  out$locus_concordance <- NA_real_
  out$alignment <- NA_real_
  out$comparability <- NA_real_
  out$n_eff_loci <- NA_real_
  out$concordance_z <- NA_real_
  out$p_concordance <- NA_real_

  energy <- .program_energy(loadings)
  set.seed(seed)
  # Drive the scoring loop from out's own (trait_id_a, trait_id_b) pairs, not
  # a freshly re-sorted trait_ids: .program_pair_keys() orients A/B by each
  # trait's row position in `progs` (first appearance in program_data), which
  # need not match alphabetical order. Re-deriving pairs via sort() + combn()
  # here could hand .concordance_trait_pair() the opposite orientation from
  # the one `out`'s rows actually use, silently finding no rows to score.
  trait_pairs <- unique(out[, c("trait_id_a", "trait_id_b")])
  for (i in seq_len(nrow(trait_pairs))) {
    out <- .concordance_trait_pair(
      out = out, loadings = loadings,
      trait_ids_by_program = trait_ids_by_program,
      energy = energy,
      t_a = as.character(trait_pairs$trait_id_a[i]),
      t_b = as.character(trait_pairs$trait_id_b[i]),
      n_perm = as.integer(n_perm), bridges = bridges
    )
  }

  # BH across every scored pair, not only the eligible ones: eligibility (enough
  # shared high-confidence loci) is not independent of the statistic, so the
  # p-values of eligible null pairs are not uniform and correcting over them
  # alone is anti-conservative.
  out$q_concordance <- NA_real_
  ok_pair <- is.finite(out$p_concordance)
  if (any(ok_pair)) {
    out$q_concordance[ok_pair] <- stats::p.adjust(
      out$p_concordance[ok_pair], method = "BH"
    )
  }

  out <- .profile_congruence(out, pd$profiles, exclude_features)
  out <- .assign_link_tiers(
    out,
    phi_high = phi_high, phi_low = phi_low,
    min_shared_loci = min_shared_loci, min_eff_loci = min_eff_loci,
    min_alignment = min_alignment, fdr_threshold = fdr_threshold
  )

  rg <- .empty_module_rg()
  if (!is.null(coloc_groups)) {
    rg <- module_rg(
      result = list(pairs = out),
      program_data = pd,
      coloc_groups = coloc_groups,
      seed = seed,
      bridges = bridges,
      min_ld_r = min_ld_r,
      ld_proxies_fn = ld_proxies_fn
    )
  }

  return(list(
    pairs = out,
    programs = .program_link_counts(progs, out),
    module_rg = rg,
    bridges = bridges,
    settings = settings
  ))
}


# Profile strength, locus status, tier and link of every pair.
#
# profile_strength: strong (|phi| >= phi_high), moderate (phi_low <= |phi| <
# phi_high) or weak.
# locus_status: too_few_loci when the pair shares fewer than min_shared_loci
# high-confidence loci or fewer than min_eff_loci loci effectively carry the
# locus score -- the loci cannot settle the question either way; otherwise
# agree when alignment >= min_alignment and q_concordance <= fdr_threshold,
# and disagree when the loci were testable and fail either.
.assign_link_tiers <- function(out, phi_high = 0.85, phi_low = 0.70,
                               min_shared_loci = 3L, min_eff_loci = 5,
                               min_alignment = 0.5, fdr_threshold = 0.05) {
  phi <- abs(out$phi_traits)
  phi[!is.finite(phi)] <- 0
  out$profile_strength <- dplyr::case_when(
    phi >= phi_high ~ "strong",
    phi >= phi_low ~ "moderate",
    TRUE ~ "weak"
  )

  testable <- is.finite(out$n_loci_shared) & out$n_loci_shared >= min_shared_loci &
    is.finite(out$n_eff_loci) & out$n_eff_loci >= min_eff_loci
  agrees <- is.finite(out$alignment) & out$alignment >= min_alignment &
    is.finite(out$q_concordance) & out$q_concordance <= fdr_threshold
  out$locus_status <- dplyr::case_when(
    !testable ~ "too_few_loci",
    agrees ~ "agree",
    TRUE ~ "disagree"
  )

  profile <- out$profile_strength
  locus <- out$locus_status
  out$tier <- dplyr::case_when(
    profile == "strong" & locus == "agree" ~ "high",
    profile == "strong" & locus == "too_few_loci" ~ "medium",
    profile == "moderate" & locus == "agree" ~ "medium",
    profile == "strong" & locus == "disagree" ~ "low",
    profile == "moderate" & locus == "too_few_loci" ~ "low",
    profile == "weak" & locus == "agree" ~ "low",
    TRUE ~ "none"
  )
  out$linked <- out$tier %in% c("high", "medium")
  return(out)
}


# Profile congruence of every pair: Tucker's congruence of the two programs'
# confidence-weighted trait loadings over the background traits both fits
# carry, leaving out both anchor traits and any trait excluded from either fit.
.profile_congruence <- function(out, profiles, exclude_features = NULL) {
  out$phi_traits <- NA_real_
  out$n_traits_shared <- 0L
  out$n_eff_traits <- NA_real_
  out$direction_profile <- NA_character_
  if (is.null(profiles) || nrow(profiles) == 0 || nrow(out) == 0) {
    return(out)
  }
  profiles <- profiles[!is.na(profiles$feature_trait_id), , drop = FALSE]
  profiles$weighted <- profiles$loading * .confidence_factor(profiles$lfsr)
  trait_pairs <- unique(out[, c("trait_id_a", "trait_id_b")])
  for (i in seq_len(nrow(trait_pairs))) {
    ta <- as.character(trait_pairs$trait_id_a[i])
    tb <- as.character(trait_pairs$trait_id_b[i])
    rows <- which(
      as.character(out$trait_id_a) == ta & as.character(out$trait_id_b) == tb
    )
    pa <- profiles[as.character(profiles$trait_id) == ta, , drop = FALSE]
    pb <- profiles[as.character(profiles$trait_id) == tb, , drop = FALSE]
    feats <- .comparable_features(pa, pb, ta, tb, exclude_features)
    out$n_traits_shared[rows] <- length(feats)
    if (length(feats) < 2) {
      next
    }
    ma <- .weighted_profile_matrix(pa, feats)
    mb <- .weighted_profile_matrix(pb, feats)
    num <- ma %*% t(mb)
    den <- sqrt(outer(rowSums(ma^2), rowSums(mb^2)))
    phi <- num / den
    phi[!is.finite(phi)] <- NA_real_
    s1 <- abs(ma) %*% t(abs(mb))
    s2 <- ma^2 %*% t(mb^2)
    n_eff <- s1^2 / s2
    n_eff[!is.finite(n_eff)] <- NA_real_
    ia <- match(as.character(out$program_id_a[rows]), rownames(ma))
    ib <- match(as.character(out$program_id_b[rows]), rownames(mb))
    ok <- !is.na(ia) & !is.na(ib)
    out$phi_traits[rows[ok]] <- phi[cbind(ia[ok], ib[ok])]
    out$n_eff_traits[rows[ok]] <- n_eff[cbind(ia[ok], ib[ok])]
  }
  out$direction_profile <- ifelse(
    is.na(out$phi_traits), NA_character_,
    ifelse(out$phi_traits >= 0, "concordant", "antagonistic")
  )
  return(out)
}


# The background traits two traits' profiles are compared over: those both
# fits carry, minus both anchors and anything excluded from either fit.
.comparable_features <- function(pa, pb, ta, tb, exclude_features = NULL) {
  drop <- c(ta, tb, as.character(exclude_features[[ta]]),
            as.character(exclude_features[[tb]]))
  return(setdiff(
    intersect(unique(pa$feature_trait_id), unique(pb$feature_trait_id)), drop
  ))
}


#' @title Background Traits Carrying a Profile Match
#' @description For each cross-trait program pair, list the background traits
#' that carry its profile congruence (`phi_traits`). Each comparable trait
#' \eqn{t} contributes \eqn{\tilde l_a(t) \tilde l_b(t)} to the numerator of
#' \eqn{\phi}, the product of the two programs' confidence-weighted trait
#' loadings. A trait's `share` is that contribution as a fraction of the total
#' absolute contribution, signed so that a positive share supports the sign
#' of \eqn{\phi} and a negative one works against it; the absolute shares of
#' all comparable traits sum to 1. `n_eff_traits` summarises the same shares
#' as one number; this function names the traits, so a match carried by
#' near-duplicate traits can be seen.
#'
#' The comparable traits are exactly those [compare_program_pairs_loadings()]
#' used: both fits carry them, both anchors are left out, and so is every
#' trait in `exclude_features` (read from `result$settings`).
#' @param result A [compare_program_pairs_loadings()] result.
#' @param program_data The extraction result(s) passed to
#'   [compare_program_pairs_loadings()]; their `profiles` are used.
#' @param pairs Optional subset of `result$pairs` to describe. Defaults to
#'   every pair.
#' @param top_n Number of traits to keep per pair. Defaults to `5`.
#' @return A list with:
#'   \itemize{
#'     \item traits: one row per pair and kept trait: `program_id_a`,
#'       `program_id_b`, `rank`, `feature_trait_id`, `feature_trait_name`,
#'       `contribution` and `share`.
#'     \item pairs: one row per pair: `program_id_a`, `program_id_b` and
#'       `top_traits`, the kept traits as one string
#'       (`"name 41%; name 22%; ..."`).
#'   }
#' @export
profile_contributions <- function(result, program_data, pairs = NULL, top_n = 5L) {
  if (is.null(pairs)) {
    pairs <- result$pairs
  }
  empty_traits <- data.frame(
    program_id_a = character(0), program_id_b = character(0),
    rank = integer(0), feature_trait_id = character(0),
    feature_trait_name = character(0), contribution = numeric(0),
    share = numeric(0), stringsAsFactors = FALSE
  )
  empty <- list(
    traits = empty_traits,
    pairs = data.frame(
      program_id_a = character(0), program_id_b = character(0),
      top_traits = character(0), stringsAsFactors = FALSE
    )
  )
  profiles <- .normalise_program_data(program_data)$profiles
  if (is.null(pairs) || nrow(pairs) == 0 || is.null(profiles) || nrow(profiles) == 0) {
    return(empty)
  }
  exclude_features <- result$settings$exclude_features
  profiles <- profiles[!is.na(profiles$feature_trait_id), , drop = FALSE]
  profiles$weighted <- profiles$loading * .confidence_factor(profiles$lfsr)
  feature_names <- stats::setNames(
    as.character(profiles$feature_trait_name), as.character(profiles$feature_trait_id)
  )
  feature_names <- feature_names[!duplicated(names(feature_names))]

  rows <- list()
  trait_pairs <- unique(pairs[, c("trait_id_a", "trait_id_b")])
  for (i in seq_len(nrow(trait_pairs))) {
    ta <- as.character(trait_pairs$trait_id_a[i])
    tb <- as.character(trait_pairs$trait_id_b[i])
    pa <- profiles[as.character(profiles$trait_id) == ta, , drop = FALSE]
    pb <- profiles[as.character(profiles$trait_id) == tb, , drop = FALSE]
    feats <- .comparable_features(pa, pb, ta, tb, exclude_features)
    if (length(feats) == 0) {
      next
    }
    ma <- .weighted_profile_matrix(pa, feats)
    mb <- .weighted_profile_matrix(pb, feats)
    sub <- pairs[as.character(pairs$trait_id_a) == ta &
                   as.character(pairs$trait_id_b) == tb, , drop = FALSE]
    for (k in seq_len(nrow(sub))) {
      pid_a <- as.character(sub$program_id_a[k])
      pid_b <- as.character(sub$program_id_b[k])
      if (!pid_a %in% rownames(ma) || !pid_b %in% rownames(mb)) {
        next
      }
      contrib <- ma[pid_a, ] * mb[pid_b, ]
      total <- sum(abs(contrib))
      if (!is.finite(total) || total == 0) {
        next
      }
      direction <- if (sum(contrib) < 0) -1 else 1
      keep <- utils::head(order(-abs(contrib)), top_n)
      keep <- keep[contrib[keep] != 0]
      rows[[length(rows) + 1]] <- data.frame(
        program_id_a = pid_a, program_id_b = pid_b,
        rank = seq_along(keep),
        feature_trait_id = feats[keep],
        feature_trait_name = unname(feature_names[feats[keep]]),
        contribution = unname(contrib[keep]),
        share = unname(direction * contrib[keep] / total),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(empty)
  }
  traits <- dplyr::bind_rows(rows)
  summary <- traits |>
    dplyr::group_by(program_id_a, program_id_b) |>
    dplyr::summarise(
      top_traits = paste0(
        feature_trait_name, " ", round(100 * share), "%",
        collapse = "; "
      ),
      .groups = "drop"
    ) |>
    as.data.frame(stringsAsFactors = FALSE)
  return(list(traits = traits, pairs = summary))
}


# programs x features matrix of confidence-weighted trait loadings (0 where a
# program has no loading on a feature).
.weighted_profile_matrix <- function(profiles, features) {
  pids <- unique(as.character(profiles$program_id))
  m <- matrix(0, nrow = length(pids), ncol = length(features),
              dimnames = list(pids, features))
  keep <- profiles$feature_trait_id %in% features
  p <- profiles[keep, , drop = FALSE]
  if (nrow(p) > 0) {
    m[cbind(match(as.character(p$program_id), pids),
            match(as.character(p$feature_trait_id), features))] <- p$weighted
  }
  m[!is.finite(m)] <- 0
  return(m)
}


# One row per program with how many links it takes part in.
.program_link_counts <- function(progs, pairs) {
  linked <- if ("linked" %in% names(pairs)) pairs[pairs$linked, , drop = FALSE] else pairs[0, ]
  ends <- c(as.character(linked$program_id_a), as.character(linked$program_id_b))
  counts <- table(factor(ends, levels = as.character(progs$program_id)))
  out <- as.data.frame(progs, stringsAsFactors = FALSE)
  out$program_id <- as.character(out$program_id)
  out$trait_id <- as.character(out$trait_id)
  out$n_links <- as.integer(counts[out$program_id])
  return(out)
}


#' @title Module-Restricted Cross-Trait Effect Correlation
#' @description For linked cross-trait program pairs, correlate the two target
#' traits' own SNP effects across the module's loci — a module-level analogue of
#' a local genetic correlation, where the stratification unit is a data-driven
#' program rather than a genomic window. This is the direction of a link: the
#' sign of [compare_program_pairs_loadings()]'s `concordance_z` is not (see its
#' documentation).
#'
#' Loci are the shared loci either program claims as a high-confidence member;
#' each locus is weighted by how strongly the module loads on it in both
#' programs and by the precision of the two effect estimates:
#' \deqn{w_l = w(\mathrm{lfsr}_a)\, w(\mathrm{lfsr}_b)
#' |f_a(l) f_b(l)| / (se_a^2 + se_b^2)}
#' with the confidence weight \eqn{w} of [compare_program_pairs_loadings()]
#' so loci the module loads on only weakly contribute little. `rg` is the
#' weighted correlation of the two traits' effects **through the origin**,
#' \eqn{\sum_l w_l \beta_A \beta_B / \sqrt{\sum_l w_l \beta_A^2 \sum_l w_l \beta_B^2}},
#' rather than a mean-centred Pearson correlation: each locus's effect allele is
#' arbitrary, and the uncentred form gives the same answer whichever allele a
#' locus is coded on. A positive `rg` means the alleles that raise trait A also
#' raise trait B across this module. Its confidence interval is a percentile
#' bootstrap over loci, and `rg_direction` is `"concordant"` or
#' `"antagonistic"` only when that interval excludes zero (`"undetermined"`
#' otherwise).
#'
#' Effect alleles must be the same in both traits for the product to mean
#' anything. Between two existing GPMap traits this is checked by requiring the
#' two traits' `variant_id` to match at a locus; mismatches are dropped and
#' counted in `n_allele_mismatch`. For a pair involving an uploaded GWAS, both
#' effects are read from the **upload's own coloc table** -- the upload's row
#' and the partner trait's row in the same upload coloc group, which sit at the
#' same variant and on the same allele. The partner trait's own table cannot be
#' used for this: GPMap reports an upload's associations and an existing
#' trait's associations on opposite alleles at about half of all variants (on
#' a T2D upload against GPMap's own T2D trait, every disagreement in sign was
#' matched by an effect-allele frequency of `1 - eaf`). Loci with no partner
#' row in the upload's table are dropped and counted in `n_allele_mismatch`.
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
#'   trait (from the upload's table for a pair involving an uploaded GWAS).
#' @param linked_only If `TRUE` (default), estimate only for pairs with
#'   `linked = TRUE`; `FALSE` estimates every pair supplied.
#' @param n_boot Percentile-bootstrap replicates over loci. Defaults to `1000`.
#' @param min_loci_rg Minimum usable loci for an estimate. Defaults to `10`.
#' @param seed RNG seed for the bootstrap.
#' @return A data.frame with one row per pair: `program_id_a`, `program_id_b`,
#'   `trait_id_a`, `trait_id_b`, `rg`, `rg_ci_lower`, `rg_ci_upper`,
#'   `n_loci_rg`, `n_allele_mismatch`, and `rg_direction`. When the pairs
#'   carry `direction_profile` (as [compare_program_pairs_loadings()]'s do),
#'   also `direction_check`: `"agree"` or `"conflict"` when `rg` calls a
#'   direction and it matches or contradicts the sign of `phi_traits`, and
#'   `"no call"` otherwise. `rg` is a check on a link's direction; it never
#'   changes a pair's tier.
#' @export
module_rg <- function(result,
                      program_data,
                      coloc_groups,
                      linked_only = TRUE,
                      n_boot = 1000L,
                      min_loci_rg = 10L,
                      seed = 1,
                      bridges = NULL,
                      min_ld_r = 0.5,
                      ld_proxies_fn = ld_proxies) {
  pairs <- if (is.data.frame(result)) result else result$pairs
  if (is.null(pairs) || nrow(pairs) == 0) {
    return(.empty_module_rg())
  }
  if (isTRUE(linked_only) && "linked" %in% names(pairs)) {
    pairs <- pairs[!is.na(pairs$linked) & pairs$linked, , drop = FALSE]
  }
  if (nrow(pairs) == 0) {
    return(.empty_module_rg())
  }
  if (!is.list(coloc_groups) || is.data.frame(coloc_groups)) {
    stop("coloc_groups must be a named list of coloc-group data.frames")
  }

  loadings <- .normalise_program_data(program_data)$loadings
  bridges <- .complete_bridges(
    unique(pairs[, c("trait_id_a", "trait_id_b")]), coloc_groups, bridges,
    min_ld_r = min_ld_r, ld_proxies_fn = ld_proxies_fn
  )
  effects <- .effects_cache(coloc_groups)

  set.seed(seed)
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    ta <- as.character(pairs$trait_id_a[i])
    tb <- as.character(pairs$trait_id_b[i])
    pa <- as.character(pairs$program_id_a[i])
    pb <- as.character(pairs$program_id_b[i])
    # Restricted to the loci the two programs actually claim, not the whole
    # trait-pair intersection. Weighting the full axis by |loading| left every
    # pair with the same n_loci_rg -- the same trait-level rg lightly
    # reweighted -- including pairs whose programs share a single locus.
    est <- .rg_for_programs(
      loadings, pa, pb, ta, tb, effects, bridges,
      n_boot = n_boot, min_loci_rg = min_loci_rg
    )
    return(data.frame(
      program_id_a = pa, program_id_b = pb, trait_id_a = ta, trait_id_b = tb,
      est, stringsAsFactors = FALSE
    ))
  })

  out <- dplyr::bind_rows(rows)
  rownames(out) <- NULL
  if ("direction_profile" %in% names(pairs)) {
    profile_dir <- pairs$direction_profile[match(
      paste(out$program_id_a, out$program_id_b),
      paste(pairs$program_id_a, pairs$program_id_b)
    )]
    out$direction_check <- .direction_check(profile_dir, out$rg_direction)
  }
  return(out)
}


# Does rg's direction agree with the profile's? "agree" / "conflict" when rg
# calls a direction (its CI excludes zero), "no call" otherwise.
.direction_check <- function(direction_profile, rg_direction) {
  called <- !is.na(rg_direction) & rg_direction %in% c("concordant", "antagonistic") &
    !is.na(direction_profile)
  return(dplyr::case_when(
    !called ~ "no call",
    rg_direction == direction_profile ~ "agree",
    TRUE ~ "conflict"
  ))
}


# Fill in any LD locus bridge a set of trait pairs needs (an uploaded trait on
# either side -- see .needs_locus_bridge()) and is not already in `bridges`.
# Callers that already built bridges (compare_program_pairs_loadings()) pass
# them in so this is not redone.
.complete_bridges <- function(trait_pairs, coloc_groups, bridges = NULL,
                              min_ld_r = 0.5, ld_proxies_fn = ld_proxies) {
  if (is.null(bridges)) {
    bridges <- list()
  }
  for (i in seq_len(nrow(trait_pairs))) {
    ta <- as.character(trait_pairs[[1]][i])
    tb <- as.character(trait_pairs[[2]][i])
    if (!.needs_locus_bridge(ta, tb)) {
      next
    }
    key <- .locus_bridge_key(ta, tb)
    if (!is.null(bridges[[key]])) {
      next
    }
    bridges[[key]] <- .build_locus_bridge(
      coloc_groups[[ta]], coloc_groups[[tb]], ta, tb,
      min_r = min_ld_r, ld_proxies_fn = ld_proxies_fn
    )
  }
  return(bridges)
}


# Memoised per-locus effects: effects(table_trait, trait) is `trait`'s own rows
# of coloc_groups[[table_trait]] (.trait_locus_effects()), computed once.
.effects_cache <- function(coloc_groups) {
  memo <- new.env(parent = emptyenv())
  return(function(table_trait, trait) {
    key <- paste(table_trait, trait, sep = "__")
    if (!exists(key, envir = memo, inherits = FALSE)) {
      cg <- coloc_groups[[as.character(table_trait)]]
      val <- if (is.null(cg)) NULL else .trait_locus_effects(cg, trait)
      assign(key, list(val), envir = memo)
    }
    return(get(key, envir = memo, inherits = FALSE)[[1]])
  })
}


# rg for one pair of program sets (one program each for a link, several for a
# family): the loci the programs claim on the shared axis, each trait's
# loading at a locus taken from whichever of its programs loads there most.
.rg_for_programs <- function(loadings, programs_a, programs_b, ta, tb, effects,
                             bridges = NULL, n_boot = 1000L, min_loci_rg = 10L) {
  loci <- .module_locus_axis(loadings, programs_a, programs_b, ta, tb, bridges = bridges)
  if (length(loci) == 0) {
    return(.empty_rg_estimate())
  }
  # loci is expressed in trait_b's space (see .shared_locus_axis()); trait_a's
  # programs and effects need the reverse translation, kept positional
  # (loci_a[i] <-> loci[i]) so everything stays aligned to loci's positions.
  loci_a <- .translate_locus_ids_positional(loci, tb, ta, bridges)
  va <- .programs_locus_vector(loadings, programs_a, loci, lookup_loci = loci_a)
  vb <- .programs_locus_vector(loadings, programs_b, loci)
  return(.rg_over_loci(
    loci_a, loci, ta, tb, va, vb, effects,
    n_boot = n_boot, min_loci_rg = min_loci_rg
  ))
}


# Allele-aligned effects of both traits at each locus, positionally:
# loci_a[i] (trait_a's own id) and loci_b[i] (trait_b's own id) are one locus.
# Returns beta/se for each side and `mismatch`, the loci with an effect in both
# traits that could not be put on a common allele.
.pair_locus_effects <- function(loci_a, loci_b, ta, tb, effects) {
  n <- length(loci_b)
  empty <- list(
    beta_a = rep(NA_real_, n), se_a = rep(NA_real_, n),
    beta_b = rep(NA_real_, n), se_b = rep(NA_real_, n),
    mismatch = rep(FALSE, n)
  )
  pick <- function(eff, ids) {
    return(eff[match(as.character(ids), eff$coloc_group_id), , drop = FALSE])
  }
  upload_a <- is_guid(ta) && !is_guid(tb)
  upload_b <- is_guid(tb) && !is_guid(ta)
  if (upload_a || upload_b) {
    # Both effects from the upload's own table, at the upload's coloc group:
    # the upload's row and the partner's row there share a variant and allele.
    upload <- if (upload_a) ta else tb
    partner <- if (upload_a) tb else ta
    upload_ids <- if (upload_a) loci_a else loci_b
    partner_ids <- if (upload_a) loci_b else loci_a
    own <- effects(upload, upload)
    in_upload <- effects(upload, partner)
    native <- effects(partner, partner)
    if (is.null(own)) {
      return(empty)
    }
    e_up <- pick(own, upload_ids)
    e_partner <- if (is.null(in_upload)) {
      e_up[0, ][rep(NA_integer_, n), , drop = FALSE]
    } else {
      pick(in_upload, upload_ids)
    }
    has_native <- if (is.null(native)) {
      rep(FALSE, n)
    } else {
      is.finite(pick(native, partner_ids)$beta)
    }
    mismatch <- is.finite(e_up$beta) & !is.finite(e_partner$beta) & has_native
    e_a <- if (upload_a) e_up else e_partner
    e_b <- if (upload_a) e_partner else e_up
    return(list(
      beta_a = e_a$beta, se_a = e_a$se, beta_b = e_b$beta, se_b = e_b$se,
      mismatch = mismatch
    ))
  }
  eff_a <- effects(ta, ta)
  eff_b <- effects(tb, tb)
  if (is.null(eff_a) || is.null(eff_b)) {
    return(empty)
  }
  e_a <- pick(eff_a, loci_a)
  e_b <- pick(eff_b, loci_b)
  differs <- !is.na(e_a$variant_id) & !is.na(e_b$variant_id) &
    as.character(e_a$variant_id) != as.character(e_b$variant_id)
  mismatch <- is.finite(e_a$beta) & is.finite(e_b$beta) & differs
  beta_a <- e_a$beta
  beta_b <- e_b$beta
  beta_a[differs] <- NA_real_
  beta_b[differs] <- NA_real_
  return(list(
    beta_a = beta_a, se_a = e_a$se, beta_b = beta_b, se_b = e_b$se,
    mismatch = mismatch
  ))
}


# The rg estimate over one positional locus axis, given each side's per-locus
# loading and confidence (va, vb from .programs_locus_vector()).
.rg_over_loci <- function(loci_a, loci_b, ta, tb, va, vb, effects,
                          n_boot = 1000L, min_loci_rg = 10L) {
  out <- .empty_rg_estimate()
  e <- .pair_locus_effects(loci_a, loci_b, ta, tb, effects)
  usable <- is.finite(e$beta_a) & is.finite(e$beta_b) &
    is.finite(e$se_a) & is.finite(e$se_b) & e$se_a > 0 & e$se_b > 0
  usable[is.na(usable)] <- FALSE
  out$n_allele_mismatch <- as.integer(sum(e$mismatch, na.rm = TRUE))
  out$n_loci_rg <- as.integer(sum(usable))
  if (sum(usable) < min_loci_rg) {
    return(out)
  }

  ba <- e$beta_a[usable]
  bb <- e$beta_b[usable]
  w <- va$weight[usable] * vb$weight[usable] *
    abs(va$loading[usable] * vb$loading[usable]) /
    (e$se_a[usable]^2 + e$se_b[usable]^2)
  w[!is.finite(w) | w < 0] <- 0

  rg <- .weighted_uncentred_cor(ba, bb, w)
  ci <- c(NA_real_, NA_real_)
  if (is.finite(rg) && n_boot > 0) {
    n <- length(ba)
    draws <- vapply(seq_len(n_boot), function(b) {
      idx <- sample.int(n, n, replace = TRUE)
      return(.weighted_uncentred_cor(ba[idx], bb[idx], w[idx]))
    }, numeric(1))
    draws <- draws[is.finite(draws)]
    if (length(draws) > 0) {
      ci <- unname(stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE))
    }
  }
  out$rg <- rg
  out$rg_ci_lower <- ci[1]
  out$rg_ci_upper <- ci[2]
  out$rg_direction <- .rg_direction(rg, ci[1], ci[2])
  return(out)
}


# Direction of an rg estimate: called only when its CI excludes zero.
.rg_direction <- function(rg, lower, upper) {
  if (!is.finite(rg)) {
    return(NA_character_)
  }
  if (is.finite(lower) && lower > 0) {
    return("concordant")
  }
  if (is.finite(upper) && upper < 0) {
    return("antagonistic")
  }
  return("undetermined")
}


.empty_rg_estimate <- function() {
  return(data.frame(
    rg = NA_real_, rg_ci_lower = NA_real_, rg_ci_upper = NA_real_,
    n_loci_rg = 0L, n_allele_mismatch = 0L, rg_direction = NA_character_,
    stringsAsFactors = FALSE
  ))
}


# Whether a cross-trait comparison needs LD-based locus bridging. An uploaded
# GWAS assigns its own coloc_group_id values independently of GPMap's shared
# numbering, so direct coloc_group_id equality finds almost nothing even when
# the same locus is genuinely shared (confirmed on real data: an uploaded T2D
# GWAS shared 1 coloc_group_id with an existing triglycerides trait by direct
# equality, ~40 by variant_id/LD). Existing (non-upload) trait pairs already
# share one coloc_group_id space -- confirmed on real data at 26-190 shared
# loci for several pairs -- and need no bridging.
.needs_locus_bridge <- function(trait_a, trait_b) {
  return(is_guid(trait_a) || is_guid(trait_b))
}


# Canonical, order-independent cache key for a trait pair.
.locus_bridge_key <- function(trait_a, trait_b) {
  return(paste(sort(c(as.character(trait_a), as.character(trait_b))), collapse = "__"))
}


# ld_proxies() rows (lead_variant_id, proxy_variant_id, signed r) for every id
# in variant_ids, as character ids. ld_proxies_fn is injectable so tests can
# supply a fake, network-free function. Chunked because the underlying API
# appends one query parameter per id (ld_proxies_by_variant_id_api()); a chunk
# that errors (rate limit, transient failure) is skipped rather than aborting
# the whole lookup.
.ld_proxy_table <- function(variant_ids, ld_proxies_fn = ld_proxies, chunk_size = 80L) {
  empty <- data.frame(
    lead_variant_id = character(0), proxy_variant_id = character(0),
    r = numeric(0), stringsAsFactors = FALSE
  )
  variant_ids <- unique(stats::na.omit(as.character(variant_ids)))
  if (length(variant_ids) == 0) {
    return(empty)
  }
  chunks <- split(variant_ids, ceiling(seq_along(variant_ids) / chunk_size))
  proxies <- dplyr::bind_rows(lapply(chunks, function(ids) {
    res <- tryCatch(
      ld_proxies_fn(variant_ids = as.numeric(ids)),
      error = function(e) NULL
    )
    lds <- if (is.list(res)) res$lds else NULL
    if (!is.data.frame(lds) || nrow(lds) == 0) {
      return(NULL)
    }
    return(data.frame(
      lead_variant_id = as.character(lds$lead_variant_id),
      proxy_variant_id = as.character(lds$proxy_variant_id),
      r = as.numeric(lds$r),
      stringsAsFactors = FALSE
    ))
  }))
  if (is.null(proxies) || nrow(proxies) == 0) {
    return(empty)
  }
  return(proxies[!duplicated(proxies[, c("lead_variant_id", "proxy_variant_id")]), , drop = FALSE])
}


# Signed LD r between each (variant_a[i], variant_b[i]) pair from a
# .ld_proxy_table(), looked up in either direction; 1 where the two are the
# same variant, NA where the table has no row for the pair.
.ld_pair_r <- function(variant_a, variant_b, proxies) {
  variant_a <- as.character(variant_a)
  variant_b <- as.character(variant_b)
  r <- rep(NA_real_, length(variant_a))
  if (nrow(proxies) > 0) {
    fwd <- paste(proxies$lead_variant_id, proxies$proxy_variant_id)
    rev <- paste(proxies$proxy_variant_id, proxies$lead_variant_id)
    key <- paste(variant_a, variant_b)
    r <- proxies$r[match(key, fwd)]
    r[is.na(r)] <- proxies$r[match(key[is.na(r)], rev)]
  }
  r[!is.na(variant_a) & variant_a == variant_b] <- 1
  return(r)
}


# Best LD match in target_variant_ids for each of query_variant_ids: r = 1 for
# an exact variant_id match (free, no API call), else the ld_proxies() proxy
# with the largest |r| at |r| >= min_r. r is returned signed: negative means the
# two variants' alleles are in opposite phase. Proxies are looked up from
# both sides (query and target ids), because GPMap only stores proxy lists for
# some variants -- on real data ~98% of Suzuki index SNVs and ~36% of an
# upload's own leads had none, while the same pairs were found from the other
# side. GPMap's proxy table itself stops at r^2 = 0.8 (|r| ~ 0.89), so min_r
# below that has no effect. Pass a precomputed .ld_proxy_table() as `proxies`
# to skip the lookup.
.ld_best_match <- function(query_variant_ids, target_variant_ids,
                           min_r = 0.5, ld_proxies_fn = ld_proxies,
                           chunk_size = 80L, proxies = NULL) {
  query_variant_ids <- unique(stats::na.omit(as.character(query_variant_ids)))
  target_variant_ids <- unique(stats::na.omit(as.character(target_variant_ids)))
  empty <- data.frame(
    query_variant_id = character(0), target_variant_id = character(0),
    r = numeric(0), stringsAsFactors = FALSE
  )
  if (length(query_variant_ids) == 0 || length(target_variant_ids) == 0) {
    return(empty)
  }

  exact <- intersect(query_variant_ids, target_variant_ids)
  exact_rows <- if (length(exact) > 0) {
    data.frame(
      query_variant_id = exact, target_variant_id = exact, r = 1,
      stringsAsFactors = FALSE
    )
  } else {
    empty
  }

  remaining <- setdiff(query_variant_ids, exact)
  if (length(remaining) == 0) {
    return(exact_rows)
  }
  if (is.null(proxies)) {
    proxies <- .ld_proxy_table(
      c(remaining, target_variant_ids), ld_proxies_fn = ld_proxies_fn,
      chunk_size = chunk_size
    )
  }
  if (nrow(proxies) == 0) {
    return(exact_rows)
  }
  fwd <- proxies[
    proxies$lead_variant_id %in% remaining &
      proxies$proxy_variant_id %in% target_variant_ids,
    ,
    drop = FALSE
  ]
  rev <- proxies[
    proxies$proxy_variant_id %in% remaining &
      proxies$lead_variant_id %in% target_variant_ids,
    ,
    drop = FALSE
  ]
  ld_rows <- data.frame(
    query_variant_id = c(fwd$lead_variant_id, rev$proxy_variant_id),
    target_variant_id = c(fwd$proxy_variant_id, rev$lead_variant_id),
    r = c(fwd$r, rev$r),
    stringsAsFactors = FALSE
  )
  ld_rows <- ld_rows[is.finite(ld_rows$r) & abs(ld_rows$r) >= min_r, , drop = FALSE]
  if (nrow(ld_rows) > 0) {
    ld_rows <- ld_rows[order(-abs(ld_rows$r)), , drop = FALSE]
    ld_rows <- ld_rows[!duplicated(ld_rows$query_variant_id), , drop = FALSE]
  } else {
    ld_rows <- empty
  }
  return(dplyr::bind_rows(exact_rows, ld_rows))
}


# A trait's own rows of its coloc-group table (trait_id == trait), with a
# coloc_group_id and variant_id.
.own_coloc_rows <- function(coloc_groups, trait) {
  return(coloc_groups[
    !is.na(coloc_groups$trait_id) &
      as.character(coloc_groups$trait_id) == as.character(trait) &
      !is.na(coloc_groups$coloc_group_id) & !is.na(coloc_groups$variant_id),
    ,
    drop = FALSE
  ])
}


# Upload <-> existing-trait locus matches recorded by the upload's own coloc
# run: each row of the upload's table for the existing trait carries that
# trait's existing_study_extraction_id, which is the study_extraction_id of a
# row in the existing trait's own table -- and so names its coloc_group_id.
# Exact and coloc-backed (38/38 shared loci recovered for a T2D upload vs
# triglycerides, against 34/38 by variant/LD matching alone). Returns
# upload_cg, existing_cg and each side's own lead variant_id; empty when either
# table lacks the columns.
.study_extraction_matches <- function(cg_upload, cg_existing, upload_trait, existing_trait) {
  empty <- data.frame(
    upload_cg = character(0), existing_cg = character(0),
    upload_variant = character(0), existing_variant = character(0),
    stringsAsFactors = FALSE
  )
  if (!"existing_study_extraction_id" %in% names(cg_upload) ||
        !"study_extraction_id" %in% names(cg_existing)) {
    return(empty)
  }
  partner <- cg_upload[
    !is.na(cg_upload$trait_id) &
      as.character(cg_upload$trait_id) == as.character(existing_trait) &
      !is.na(cg_upload$existing_study_extraction_id) &
      !is.na(cg_upload$coloc_group_id),
    ,
    drop = FALSE
  ]
  own_upload <- .own_coloc_rows(cg_upload, upload_trait)
  own_existing <- .own_coloc_rows(cg_existing, existing_trait)
  own_existing <- own_existing[!is.na(own_existing$study_extraction_id), , drop = FALSE]
  if (nrow(partner) == 0 || nrow(own_upload) == 0 || nrow(own_existing) == 0) {
    return(empty)
  }
  hit <- match(
    as.character(partner$existing_study_extraction_id),
    as.character(own_existing$study_extraction_id)
  )
  partner <- partner[!is.na(hit), , drop = FALSE]
  hit <- hit[!is.na(hit)]
  out <- data.frame(
    upload_cg = as.character(partner$coloc_group_id),
    existing_cg = as.character(own_existing$coloc_group_id[hit]),
    upload_variant = as.character(own_upload$variant_id[
      match(as.character(partner$coloc_group_id), as.character(own_upload$coloc_group_id))
    ]),
    existing_variant = as.character(own_existing$variant_id[hit]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$upload_variant), , drop = FALSE]
  return(out[!duplicated(out$upload_cg), , drop = FALSE])
}


# Per-trait-pair locus bridge: translates one trait's own coloc_group_id
# values into whatever the same (or LD-linked) locus is called in the other
# trait's own coloc_group_id space. Only meaningful when .needs_locus_bridge()
# is TRUE for the pair -- see its docs for why. Matches, in order:
#   1. upload <-> existing trait: the upload's own coloc record of that trait
#      (.study_extraction_matches());
#   2. the remaining loci by exact variant_id, then by LD (.ld_best_match()).
# `r` is the signed LD between the two traits' own lead variants (1 when they
# are the same variant; NA for a tier-1 match whose leads differ and have no
# proxy row) -- module_rg() uses it to align effect alleles.
.build_locus_bridge <- function(coloc_groups_a, coloc_groups_b,
                                trait_a, trait_b, min_r = 0.5,
                                ld_proxies_fn = ld_proxies) {
  empty <- data.frame(
    trait_a = character(0), trait_b = character(0),
    coloc_group_id_a = character(0), coloc_group_id_b = character(0),
    r = numeric(0), stringsAsFactors = FALSE
  )
  if (!is.data.frame(coloc_groups_a) || !is.data.frame(coloc_groups_b)) {
    return(empty)
  }
  own_a <- .own_coloc_rows(coloc_groups_a, trait_a)
  own_b <- .own_coloc_rows(coloc_groups_b, trait_b)
  if (nrow(own_a) == 0 || nrow(own_b) == 0) {
    return(empty)
  }

  study <- data.frame(
    cg_a = character(0), cg_b = character(0), v_a = character(0), v_b = character(0),
    stringsAsFactors = FALSE
  )
  if (is_guid(trait_a) && !is_guid(trait_b)) {
    m <- .study_extraction_matches(coloc_groups_a, coloc_groups_b, trait_a, trait_b)
    study <- data.frame(
      cg_a = m$upload_cg, cg_b = m$existing_cg,
      v_a = m$upload_variant, v_b = m$existing_variant, stringsAsFactors = FALSE
    )
  } else if (is_guid(trait_b) && !is_guid(trait_a)) {
    m <- .study_extraction_matches(coloc_groups_b, coloc_groups_a, trait_b, trait_a)
    study <- data.frame(
      cg_a = m$existing_cg, cg_b = m$upload_cg,
      v_a = m$existing_variant, v_b = m$upload_variant, stringsAsFactors = FALSE
    )
  }

  rest_a <- own_a[!as.character(own_a$coloc_group_id) %in% study$cg_a, , drop = FALSE]
  needs_ld <- unique(c(
    as.character(rest_a$variant_id[!as.character(rest_a$variant_id) %in% as.character(own_b$variant_id)]),
    study$v_a[study$v_a != study$v_b]
  ))
  proxies <- if (length(needs_ld) > 0) {
    .ld_proxy_table(
      c(needs_ld, as.character(own_b$variant_id), study$v_b[study$v_a != study$v_b]),
      ld_proxies_fn = ld_proxies_fn
    )
  } else {
    .ld_proxy_table(character(0))
  }
  study_rows <- data.frame(
    coloc_group_id_a = study$cg_a, coloc_group_id_b = study$cg_b,
    r = .ld_pair_r(study$v_a, study$v_b, proxies), stringsAsFactors = FALSE
  )

  matches <- .ld_best_match(
    rest_a$variant_id, own_b$variant_id,
    min_r = min_r, proxies = proxies
  )
  a_lookup <- stats::setNames(
    as.character(rest_a$coloc_group_id), as.character(rest_a$variant_id)
  )
  b_lookup <- stats::setNames(
    as.character(own_b$coloc_group_id), as.character(own_b$variant_id)
  )
  ld_match_rows <- data.frame(
    coloc_group_id_a = unname(a_lookup[matches$query_variant_id]),
    coloc_group_id_b = unname(b_lookup[matches$target_variant_id]),
    r = matches$r,
    stringsAsFactors = FALSE
  )
  rows <- dplyr::bind_rows(study_rows, ld_match_rows)
  if (nrow(rows) == 0) {
    return(empty)
  }
  return(data.frame(
    trait_a = as.character(trait_a),
    trait_b = as.character(trait_b),
    rows,
    stringsAsFactors = FALSE
  ))
}


# Positional translation of coloc_group_id values native to from_trait into
# to_trait's own space: out[i] corresponds to ids[i], NA where the pair needs
# bridging (.needs_locus_bridge()) but no bridge entry covers that id. Uses
# the bridge built for that (unordered) trait pair in `bridges` (a named list
# keyed by .locus_bridge_key()). bridges = NULL, or no bridge stored for this
# pair, is a no-op (ids returned unchanged) -- this is what keeps non-upload
# pairs unchanged, and is the right fallback when no bridge was built.
.translate_locus_ids_positional <- function(ids, from_trait, to_trait, bridges = NULL) {
  bridge <- .bridge_for(ids, from_trait, to_trait, bridges)
  if (is.null(bridge)) {
    return(ids)
  }
  lookup <- stats::setNames(bridge$to, bridge$from)
  return(unname(lookup[as.character(ids)]))
}


# The stored bridge for this trait pair, oriented as from/to columns; NULL when
# there is nothing to translate.
.bridge_for <- function(ids, from_trait, to_trait, bridges) {
  if (is.null(bridges) || length(ids) == 0) {
    return(NULL)
  }
  bridge <- bridges[[.locus_bridge_key(from_trait, to_trait)]]
  if (is.null(bridge) || nrow(bridge) == 0) {
    return(NULL)
  }
  forward <- as.character(from_trait) == bridge$trait_a[1]
  return(data.frame(
    from = if (forward) bridge$coloc_group_id_a else bridge$coloc_group_id_b,
    to = if (forward) bridge$coloc_group_id_b else bridge$coloc_group_id_a,
    r = if ("r" %in% names(bridge)) bridge$r else NA_real_,
    stringsAsFactors = FALSE
  ))
}


# Set version of .translate_locus_ids_positional(): translates and drops
# anything that didn't map (NA), for building/intersecting locus universes
# where position doesn't matter. See .translate_locus_ids_positional() for the
# position-preserving version needed when `ids` must stay aligned with
# another vector (e.g. module_rg()'s per-locus effect lookup).
.translate_locus_ids <- function(ids, from_trait, to_trait, bridges = NULL) {
  out <- .translate_locus_ids_positional(ids, from_trait, to_trait, bridges)
  return(unique(stats::na.omit(out)))
}


# All loci both traits carry a loading for: the intersection of the two traits'
# coloc-group universes, with no lFSR / magnitude gate. `bridges` (optional,
# named list keyed by .locus_bridge_key()) translates trait A's ids into
# trait B's space first when the pair needs LD bridging (.needs_locus_bridge());
# otherwise this is identical to plain coloc_group_id equality.
.shared_locus_axis <- function(loadings, trait_a, trait_b, bridges = NULL) {
  keep <- !is.na(loadings$coloc_group_id)
  a <- unique(as.character(loadings$coloc_group_id[
    keep & as.character(loadings$trait_id) == trait_a
  ]))
  b <- unique(as.character(loadings$coloc_group_id[
    keep & as.character(loadings$trait_id) == trait_b
  ]))
  a <- .translate_locus_ids(a, trait_a, trait_b, bridges)
  return(sort(intersect(a, b)))
}


# Each program's total confidence-weighted loading energy, sum of
# (loading x w(lfsr))^2 over all of its loci: the denominator of the
# comparability share. A coloc group is represented by its largest-|loading|
# SNP, as in .program_locus_vector(); a SNP with no coloc group is its own
# locus.
.program_energy <- function(loadings) {
  locus <- ifelse(
    is.na(loadings$coloc_group_id),
    paste0("snp:", loadings$snp_id),
    as.character(loadings$coloc_group_id)
  )
  d <- data.frame(
    program_id = as.character(loadings$program_id),
    locus = locus,
    loading = loadings$loading,
    weighted = loadings$loading * .confidence_factor(loadings$lfsr),
    stringsAsFactors = FALSE
  )
  d <- d[order(d$program_id, d$locus, -abs(d$loading)), , drop = FALSE]
  d <- d[!duplicated(d[, c("program_id", "locus")]), , drop = FALSE]
  return(tapply(d$weighted^2, d$program_id, sum))
}


# Signed loading and confidence weight of one program at each locus in `loci`.
# When several SNPs map to one coloc group the SNP with the largest |loading| is
# used, so the locus is represented by its strongest evidence in that program.
#' @param lookup_loci The coloc_group_id values to actually search for in this
#'   program's own (native) loadings, positionally aligned with `loci` --
#'   defaults to `loci` itself. Differs from `loci` only when the program's
#'   trait needs LD bridging to another trait's loci (see
#'   `.translate_locus_ids_positional()`): `loci` is the shared axis expressed
#'   in the *other* trait's space, while `lookup_loci[i]` is what locus `i`
#'   is called in *this* program's own trait.
#' @noRd
.program_locus_vector <- function(loadings, program_id, loci, lookup_loci = loci) {
  loading <- stats::setNames(rep(0, length(loci)), loci)
  weight <- stats::setNames(rep(0, length(loci)), loci)
  s <- loadings[
    as.character(loadings$program_id) == program_id &
      !is.na(loadings$coloc_group_id),
    ,
    drop = FALSE
  ]
  s <- s[as.character(s$coloc_group_id) %in% lookup_loci, , drop = FALSE]
  if (nrow(s) > 0) {
    s <- s[order(as.character(s$coloc_group_id), -abs(s$loading)), , drop = FALSE]
    s <- s[!duplicated(as.character(s$coloc_group_id)), , drop = FALSE]
    i <- match(as.character(s$coloc_group_id), lookup_loci)
    loading[i] <- s$loading
    weight[i] <- .confidence_factor(s$lfsr)
  }
  return(list(loading = unname(loading), weight = unname(weight)))
}


# Loading and confidence at each locus for a set of programs from one trait:
# at each locus, the program with the largest |loading| there. With a single
# program this is .program_locus_vector().
.programs_locus_vector <- function(loadings, program_ids, loci, lookup_loci = loci) {
  vs <- lapply(program_ids, function(pid) {
    return(.program_locus_vector(loadings, pid, loci, lookup_loci = lookup_loci))
  })
  if (length(vs) == 1) {
    return(vs[[1]])
  }
  loading <- do.call(cbind, lapply(vs, function(v) v$loading))
  weight <- do.call(cbind, lapply(vs, function(v) v$weight))
  pick <- max.col(abs(loading), ties.method = "first")
  idx <- cbind(seq_along(loci), pick)
  return(list(loading = loading[idx], weight = weight[idx]))
}


# loci x programs loading and confidence matrices for one trait.
.program_locus_matrix <- function(loadings, trait_id, loci, program_ids,
                                  lookup_loci = loci) {
  fm <- matrix(
    0, nrow = length(loci), ncol = length(program_ids),
    dimnames = list(loci, program_ids)
  )
  wm <- fm
  for (j in seq_along(program_ids)) {
    v <- .program_locus_vector(loadings, program_ids[j], loci, lookup_loci = lookup_loci)
    fm[, j] <- v$loading
    wm[, j] <- v$weight
  }
  return(list(loading = fm, conf = wm))
}


# One trait pair: locus concordance and its permutation calibration.
.concordance_trait_pair <- function(out, loadings, trait_ids_by_program,
                                    t_a, t_b, n_perm, energy, bridges = NULL) {
  progs_a <- names(trait_ids_by_program)[trait_ids_by_program == t_a]
  progs_b <- names(trait_ids_by_program)[trait_ids_by_program == t_b]
  rows <- which(
    as.character(out$trait_id_a) == t_a & as.character(out$trait_id_b) == t_b
  )
  if (length(rows) == 0 || length(progs_a) == 0 || length(progs_b) == 0) {
    return(out)
  }

  # loci is expressed in t_b's space (see .shared_locus_axis()); t_a's
  # programs need the positional reverse translation to look themselves up
  # in their own (native) loadings. A no-op when the pair needs no bridging.
  loci <- .shared_locus_axis(loadings, t_a, t_b, bridges = bridges)
  out$n_loci_axis[rows] <- length(loci)
  if (length(loci) < 2) {
    return(out)
  }
  lookup_a <- .translate_locus_ids_positional(loci, t_b, t_a, bridges)

  ma <- .program_locus_matrix(loadings, t_a, loci, progs_a, lookup_loci = lookup_a)
  mb <- .program_locus_matrix(loadings, t_b, loci, progs_b)
  aw <- ma$loading * ma$conf
  bw <- mb$loading * mb$conf
  claimed_a <- .program_claimed_matrix(loadings, loci, progs_a, lookup_loci = lookup_a)
  claimed_b <- .program_claimed_matrix(loadings, loci, progs_b)

  a_idx <- match(as.character(out$program_id_a[rows]), progs_a)
  b_idx <- match(as.character(out$program_id_b[rows]), progs_b)

  # Each pair is scored on its own axis: the shared loci at least one of the
  # two programs claims. Every program carries a small loading at every locus
  # (the target row loads on all of them), and across the whole shared axis
  # those near-zero tails can be faintly correlated between traits, which a
  # shuffle null standardises into a large z for loadings of no magnitude.
  # The null still shuffles the second program's loadings over the whole shared
  # axis before reading off the claimed loci: shuffling only within the claimed
  # loci would ask whether loading magnitudes correlate there, and two programs
  # claiming the same loci with similar loadings would then score nothing.
  for (k in seq_along(rows)) {
    ia <- a_idx[k]
    ib <- b_idx[k]
    if (is.na(ia) || is.na(ib)) next
    axis <- claimed_a[, ia] | claimed_b[, ib]
    out$n_loci_axis[rows[k]] <- sum(axis)
    if (sum(axis) < 2) next

    av <- aw[axis, ia]
    bv <- bw[axis, ib]
    obs <- sum(av * bv)
    # Effective number of loci carrying the score (participation ratio of the
    # per-locus contributions): how much independent support the pair rests on.
    n_eff <- sum(abs(av * bv))^2 / sum((av * bv)^2)

    perm <- replicate(n_perm, sample.int(nrow(bw)))[axis, , drop = FALSE]
    nv <- colSums(av * matrix(bw[perm, ib], nrow = sum(axis)))
    null_sd <- stats::sd(nv)

    out$locus_concordance[rows[k]] <- obs
    # Alignment: agreement of the two programs' loadings over every locus both
    # traits carry (no magnitude cutoff). Comparability: how much of each
    # program's loading lies on those loci at all.
    sa <- sum(aw[, ia]^2)
    sb <- sum(bw[, ib]^2)
    out$alignment[rows[k]] <- if (sa > 0 && sb > 0) {
      sum(aw[, ia] * bw[, ib])^2 / (sa * sb)
    } else {
      NA_real_
    }
    e_a <- energy[[progs_a[ia]]]
    e_b <- energy[[progs_b[ib]]]
    out$comparability[rows[k]] <- if (is.finite(e_a) && is.finite(e_b) &&
                                        e_a > 0 && e_b > 0) {
      min(sa / e_a, 1) * min(sb / e_b, 1)
    } else {
      NA_real_
    }
    out$n_eff_loci[rows[k]] <- if (is.finite(n_eff)) n_eff else NA_real_
    out$concordance_z[rows[k]] <- if (is.finite(null_sd) && null_sd > 0) {
      (obs - mean(nv)) / null_sd
    } else {
      NA_real_
    }
    # Two-sided: a pair whose loadings oppose at the shared loci links as
    # readily as one whose loadings agree (the sign is not a direction).
    out$p_concordance[rows[k]] <- (1 + sum(abs(nv) >= abs(obs))) / (1 + n_perm)
  }

  return(out)
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
  out$locus_concordance <- numeric(0)
  out$alignment <- numeric(0)
  out$comparability <- numeric(0)
  out$n_eff_loci <- numeric(0)
  out$concordance_z <- numeric(0)
  out$p_concordance <- numeric(0)
  out$q_concordance <- numeric(0)
  out$phi_traits <- numeric(0)
  out$n_traits_shared <- integer(0)
  out$n_eff_traits <- numeric(0)
  out$direction_profile <- character(0)
  out$profile_strength <- character(0)
  out$locus_status <- character(0)
  out$tier <- character(0)
  out$linked <- logical(0)
  return(out)
}


.empty_module_rg <- function() {
  return(data.frame(
    program_id_a = character(0),
    program_id_b = character(0),
    trait_id_a = character(0),
    trait_id_b = character(0),
    rg = numeric(0),
    rg_ci_lower = numeric(0),
    rg_ci_upper = numeric(0),
    n_loci_rg = integer(0),
    n_allele_mismatch = integer(0),
    rg_direction = character(0),
    stringsAsFactors = FALSE
  ))
}


# Loci the programs themselves claim, intersected with the loci both traits
# carry. This is what makes module_rg a statement about the module rather than
# about the trait pair: .shared_locus_axis() returns every locus the two traits
# have in common, so weighting it by |loading| still leaves every pair using the
# same loci and reports the same n_loci_rg for all of them. program_a and
# program_b may each be several programs of one trait (a family's members).
.module_locus_axis <- function(loadings, program_a, program_b,
                               trait_a, trait_b, bridges = NULL) {
  shared <- .shared_locus_axis(loadings, trait_a, trait_b, bridges = bridges)
  if (length(shared) == 0) {
    return(character(0))
  }
  claimed <- function(pids) {
    return(unique(unlist(lapply(pids, function(pid) {
      return(.program_claimed_loci(loadings, pid))
    }))))
  }
  # shared is expressed in trait_b's space (see .shared_locus_axis()), so
  # program_a's own claimed loci (native to trait_a) need the same
  # translation before the union; program_b's are already native to trait_b.
  own <- union(
    .translate_locus_ids(claimed(program_a), trait_a, trait_b, bridges),
    claimed(program_b)
  )
  return(sort(intersect(shared, own)))
}


# The loci a program claims: its high-confidence coloc groups, or all of its
# coloc groups when it has no high-confidence member.
.program_claimed_loci <- function(loadings, program_id) {
  s <- loadings[
    as.character(loadings$program_id) == program_id &
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


# loci x programs logical matrix: does each program claim each locus?
.program_claimed_matrix <- function(loadings, loci, program_ids,
                                    lookup_loci = loci) {
  m <- matrix(
    FALSE, nrow = length(loci), ncol = length(program_ids),
    dimnames = list(loci, program_ids)
  )
  for (j in seq_along(program_ids)) {
    m[, j] <- lookup_loci %in% .program_claimed_loci(loadings, program_ids[j])
  }
  return(m)
}
