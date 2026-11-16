#' @title Link SNP Programs Across Traits
#' @description Ask, for two or more traits, which of their programs describe
#' the same biology, which SNPs carry it, and whether those SNPs push the
#' traits the same way or opposite ways, from one call. Each trait is
#' decomposed independently ([decompose_trait()]), every cross-trait pair of
#' high- and medium-confidence programs is scored and tiered
#' ([compare_program_pairs_loadings()]), the direction of each link is checked
#' against the traits' own allele effects ([module_rg()]), and links are
#' grouped into multi-trait families ([build_program_families()]). By default
#' each family's SNPs are also enriched per trait for trait categories,
#' tissues and pathways ([enrich_snp_set()], [enrich_snp_set_pathways()]),
#' which play no part in forming the links.
#'
#' The result can be printed, plotted ([plot_families()] and friends),
#' queried ([links()], [families()], [family_loci()], [programs()],
#' [annotation()], [excluded()]), re-tiered at other cutoffs without
#' rescoring ([relink()]), or written up with [create_report()].
#'
#' **Settings.** When `traits` are given as ids, `settings` is used throughout.
#' When they are given as decompositions, their own fit and validation
#' settings are kept, and only the cross-trait settings of `settings`
#' (`n_perm`, the tier cutoffs, `mutual_best` and `seed`) apply.
#' Decompositions fit with different EBMF settings factorise different
#' matrices, so they cannot be compared and are an error.
#'
#' **Permutation floor.** A locus permutation p-value is at least
#' `1 / (n_perm + 1)`, and it is BH-corrected across every scored pair. The q
#' of pairs at that floor is returned as `q_floor`. When it is above half of
#' `fdr_threshold`, more replicates (`n_perm`) are needed for loci to `agree`,
#' and a warning says so.
#' @param traits A character vector of trait ids (GPMap ids or uploaded GWAS
#'   GUIDs), or a list whose elements are trait ids or [decompose_trait()]
#'   results.
#' @param settings Analysis settings from [gpmap_settings()].
#' @param proxies,remove_duplicates,annotate,duplicates,mappings,cores,verbose
#'   As in [decompose_trait()], applied to every trait decomposed here.
#'   `annotate` also controls the family enrichment.
#' @param proxy_trait_ids Named list (or vector), keyed by an uploaded GWAS
#'   GUID, of the GPMap trait id it measures (see [decompose_trait()]).
#' @param exclude_trait_ids Named list, keyed by trait id, of further trait
#'   ids to remove as that anchor's proxies. For a trait given as a
#'   decomposition, they are only left out of its profile comparisons.
#' @return An object of class `gpmap_correspondence`, a list with:
#'   \itemize{
#'     \item traits: one row per trait (`trait_id`, `trait_name`,
#'       `n_programs`, `n_loci`, `n_high_conf_loci`)
#'     \item trait_pairs: one row per trait pair with the loci the two traits
#'       share, program pairs by tier and links by direction
#'     \item pairs: every scored program pair, with its [module_rg()]
#'       estimate and the background traits carrying its profile match
#'       (`top_traits`, [profile_contributions()]); [links()] filters it
#'     \item families, nodes, family_loci, graph, modularity, n_unlinked: see
#'       [build_program_families()]
#'     \item loci: every program locus on GPMap's shared locus numbering
#'     \item family_enrichment: per family and trait, or `NULL`
#'     \item q_floor: q of pairs at the permutation floor (`NA` when none)
#'     \item scored, bridges: the scoring result reused by [relink()]
#'     \item decompositions: the per-trait [decompose_trait()] results
#'     \item settings, call, gpmapr_version, created
#'   }
#' @examplesIf gpmap_available()
#' \donttest{
#' x <- trait_correspondence(c("1631", "1630"), cores = 1)
#' x
#' links(x)
#' }
#' @concept trait decomposition
#' @export
trait_correspondence <- function(traits,
                                 settings = gpmap_settings(),
                                 proxies = TRUE,
                                 proxy_trait_ids = NULL,
                                 exclude_trait_ids = NULL,
                                 remove_duplicates = TRUE,
                                 duplicates = NULL,
                                 annotate = TRUE,
                                 mappings = NULL,
                                 cores = .default_cores(),
                                 verbose = TRUE) {
  call <- match.call()
  settings <- .as_gpmap_settings(settings)
  if (is.atomic(traits)) {
    traits <- as.list(as.character(traits))
  }
  if (!is.list(traits) || length(traits) < 2) {
    stop("At least two traits are required", call. = FALSE)
  }
  ids <- vapply(traits, function(t) {
    if (inherits(t, "gpmap_decomposition")) {
      return(as.character(t$trait$id))
    }
    if (length(t) != 1 || !is.atomic(t)) {
      stop("Each trait must be a trait id or a decompose_trait() result", call. = FALSE)
    }
    return(as.character(t))
  }, character(1))
  if (anyDuplicated(ids) > 0) {
    stop("Each trait may appear only once: ", paste(ids[duplicated(ids)], collapse = ", "), call. = FALSE)
  }
  proxy_trait_ids <- as.list(proxy_trait_ids)
  exclude_trait_ids <- as.list(exclude_trait_ids)

  to_fit <- !vapply(traits, inherits, logical(1), what = "gpmap_decomposition")
  if (any(to_fit) && (isTRUE(proxies) || isTRUE(remove_duplicates)) && is.null(duplicates)) {
    if (verbose) message("Fetching GPMap duplicate and proxy traits")
    duplicates <- trait_duplicates()
  }
  decomps <- lapply(seq_along(traits), function(i) {
    if (!to_fit[[i]]) {
      d <- traits[[i]]
      extra <- setdiff(as.character(exclude_trait_ids[[ids[[i]]]]), d$exclude_trait_ids)
      if (length(extra) > 0) {
        message(
          "Trait ", ids[[i]], " was already decomposed: exclude_trait_ids only leave ",
          "its extra traits out of profile comparisons"
        )
        d$exclude_trait_ids <- c(d$exclude_trait_ids, extra)
      }
      return(d)
    }
    return(decompose_trait(
      ids[[i]],
      settings = settings,
      proxies = proxies,
      proxy_trait_id = proxy_trait_ids[[ids[[i]]]],
      exclude_trait_ids = exclude_trait_ids[[ids[[i]]]],
      remove_duplicates = remove_duplicates,
      keep_trait_ids = ids,
      duplicates = duplicates,
      annotate = annotate,
      mappings = mappings,
      cores = cores,
      verbose = verbose
    ))
  })
  names(decomps) <- ids
  settings <- .correspondence_settings(decomps, settings)

  if (verbose) message("Scoring cross-trait program pairs")
  scored <- compare_program_pairs_loadings(
    lapply(decomps, function(d) return(d$extraction)),
    n_perm = settings$n_perm,
    fdr_threshold = settings$fdr_threshold,
    seed = settings$seed,
    min_shared_loci = settings$min_shared_loci,
    exclude_features = .exclude_features(decomps),
    phi_high = settings$phi_high,
    phi_low = settings$phi_low,
    min_alignment = settings$min_alignment,
    min_eff_loci = settings$min_eff_loci,
    coloc_groups = lapply(decomps, function(d) return(d$coloc_groups))
  )
  out <- .assemble_correspondence(scored, decomps, settings, annotate = annotate, rescore_rg = FALSE)
  out$call <- call
  return(out)
}


#' @title Re-tier a Trait Correspondence at Other Cutoffs
#' @description Apply other link cutoffs, or `mutual_best`, to a
#' [trait_correspondence()] result without refitting or rescoring: the stored
#' pair scores are re-tiered ([assign_link_tiers()]), the direction of every
#' link is re-estimated ([module_rg()]), and families, counts and (if the
#' result was annotated) family enrichment are rebuilt. With the original
#' cutoffs it returns the same result.
#'
#' Only `phi_high`, `phi_low`, `min_shared_loci`, `min_eff_loci`,
#' `min_alignment`, `fdr_threshold` and `mutual_best` can change here.
#' `n_perm` and `seed` need rescoring
#' (`trait_correspondence(x$decompositions, settings = ...)`), and the fit and
#' validation settings a new decomposition ([decompose_trait()]).
#' @param x A [trait_correspondence()] result.
#' @param settings Named list of the settings to change, or a
#'   [gpmap_settings()] object.
#' @return A `gpmap_correspondence` object.
#' @concept trait decomposition
#' @export
relink <- function(x, settings = list()) {
  if (!inherits(x, "gpmap_correspondence")) {
    stop("x must come from trait_correspondence()", call. = FALSE)
  }
  allowed <- c(
    "phi_high", "phi_low", "min_shared_loci", "min_eff_loci", "min_alignment",
    "fdr_threshold", "mutual_best"
  )
  overrides <- if (inherits(settings, "gpmap_settings")) {
    changed <- names(settings)[!vapply(names(settings), function(k) {
      return(identical(settings[[k]], x$settings[[k]]))
    }, logical(1))]
    unclass(settings)[changed]
  } else {
    .flatten_settings_args(list(settings))
  }
  invalid <- setdiff(names(overrides), allowed)
  if (length(invalid) > 0) {
    rescore <- intersect(invalid, c("n_perm", "seed"))
    refit <- setdiff(invalid, rescore)
    hints <- c(
      if (length(rescore) > 0) {
        paste0(
          paste(rescore, collapse = ", "), " need rescoring: ",
          "trait_correspondence(x$decompositions, settings = ...)"
        )
      },
      if (length(refit) > 0) {
        paste0(
          paste(refit, collapse = ", "), " need new decompositions: decompose_trait(...)"
        )
      }
    )
    stop("relink() cannot change ", paste(invalid, collapse = ", "), ". ",
      paste(hints, collapse = "; "),
      call. = FALSE
    )
  }
  new_settings <- gpmap_settings(unclass(x$settings), overrides)
  out <- .assemble_correspondence(
    x$scored, x$decompositions, new_settings,
    annotate = !is.null(x$family_enrichment), rescore_rg = TRUE
  )
  out$call <- x$call
  return(out)
}


#' @export
print.gpmap_correspondence <- function(x, ...) {
  cat("<gpmap_correspondence> ", nrow(x$traits), " traits\n", sep = "")
  for (i in seq_len(nrow(x$traits))) {
    cat(sprintf(
      "  %-40s %3d programs, %4d loci\n",
      .truncate(x$traits$trait_name[i], 40), as.integer(x$traits$n_programs[i]),
      as.integer(x$traits$n_loci[i])
    ))
  }
  cat("  Links (high + medium tier) per trait pair:\n")
  tp <- x$trait_pairs
  for (i in seq_len(nrow(tp))) {
    cat(sprintf(
      "    %-50s %2d high, %2d medium (%d concordant, %d antagonistic); %d shared loci\n",
      .truncate(paste(tp$trait_name_a[i], "~", tp$trait_name_b[i]), 50),
      as.integer(tp$n_high[i]), as.integer(tp$n_medium[i]),
      as.integer(tp$n_concordant[i]), as.integer(tp$n_antagonistic[i]),
      as.integer(tp$n_loci_shared[i])
    ))
  }
  if (nrow(x$families) > 0) {
    cat(sprintf(
      "  %d families from %d linked programs (%d unlinked); modularity %.2f\n",
      nrow(x$families), nrow(x$nodes), as.integer(x$n_unlinked), x$modularity
    ))
  } else {
    cat("  No links, so no families.\n")
  }
  if (is.finite(x$q_floor)) {
    cat(sprintf(
      "  Pairs at the permutation floor (n_perm = %d) have q = %.3g\n",
      as.integer(x$settings$n_perm), x$q_floor
    ))
  }
  return(invisible(x))
}


#' @export
summary.gpmap_correspondence <- function(object, ...) {
  return(object$trait_pairs)
}


.truncate <- function(x, n) {
  x <- as.character(x)
  return(ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x))
}


# Settings of a correspondence: the decompositions' fit and validation
# settings (which must agree on the fit), and the link settings of `settings`.
.correspondence_settings <- function(decomps, settings) {
  fit_keys <- unname(.ebmf_settings_map()$clustering)
  validation_keys <- c(names(.ebmf_settings_map()$validation), "n_rep")
  first <- decomps[[1]]$settings
  for (d in decomps[-1]) {
    differ <- fit_keys[!vapply(fit_keys, function(k) {
      return(identical(d$settings[[k]], first[[k]]))
    }, logical(1))]
    if (length(differ) > 0) {
      stop(
        "Traits ", decomps[[1]]$trait$id, " and ", d$trait$id, " were fit with different ",
        "EBMF settings (", paste(differ, collapse = ", "), "), so their programs are not ",
        "comparable",
        call. = FALSE
      )
    }
    vdiff <- validation_keys[!vapply(validation_keys, function(k) {
      return(identical(d$settings[[k]], first[[k]]))
    }, logical(1))]
    if (length(vdiff) > 0) {
      warning(
        "Traits ", decomps[[1]]$trait$id, " and ", d$trait$id, " were validated with ",
        "different settings (", paste(vdiff, collapse = ", "), ")",
        call. = FALSE
      )
    }
  }
  link_keys <- c(
    "n_perm", "phi_high", "phi_low", "min_shared_loci", "min_eff_loci",
    "min_alignment", "fdr_threshold", "mutual_best", "seed"
  )
  out <- unclass(first)
  out[link_keys] <- unclass(settings)[link_keys]
  return(gpmap_settings(out))
}


.exclude_features <- function(decomps) {
  out <- lapply(decomps, function(d) return(as.character(d$exclude_trait_ids)))
  out <- out[lengths(out) > 0]
  return(if (length(out) == 0) NULL else out)
}


# Everything downstream of the pair scores: tiers (and rg) when re-tiering,
# profile contributions, families, locus tables, counts and family enrichment.
.assemble_correspondence <- function(scored, decomps, settings, annotate, rescore_rg) {
  program_data <- lapply(decomps, function(d) return(d$extraction))
  coloc_by_trait <- lapply(decomps, function(d) return(d$coloc_groups))
  trait_names <- vapply(decomps, function(d) return(d$trait$name), character(1))
  cutoffs <- .settings_for(settings, "link")

  core <- scored
  if (rescore_rg) {
    core$pairs <- do.call(assign_link_tiers, c(list(out = scored$pairs), cutoffs))
    core$module_rg <- if (any(core$pairs$linked)) {
      module_rg(
        result = list(pairs = core$pairs),
        program_data = program_data,
        coloc_groups = coloc_by_trait,
        seed = settings$seed,
        bridges = scored$bridges
      )
    } else {
      .empty_module_rg()
    }
  }

  pairs <- core$pairs
  rg_cols <- c(
    "program_id_a", "program_id_b", "rg", "rg_ci_lower", "rg_ci_upper",
    "n_loci_rg", "n_allele_mismatch", "rg_direction", "direction_check"
  )
  rg <- core$module_rg
  for (col in setdiff(rg_cols, names(rg))) {
    rg[[col]] <- rep(NA, nrow(rg))
  }
  pairs <- dplyr::left_join(pairs, rg[, rg_cols, drop = FALSE], by = c("program_id_a", "program_id_b"))
  scored_pairs <- pairs[pairs$tier != "none", , drop = FALSE]
  if (nrow(scored_pairs) > 0) {
    pairs <- dplyr::left_join(
      pairs,
      profile_contributions(core, program_data, pairs = scored_pairs),
      by = c("program_id_a", "program_id_b")
    )
  } else {
    pairs$top_traits <- rep(NA_character_, nrow(pairs))
  }
  pairs$trait_name_a <- unname(trait_names[as.character(pairs$trait_id_a)])
  pairs$trait_name_b <- unname(trait_names[as.character(pairs$trait_id_b)])

  families <- build_program_families(
    core,
    seed = settings$seed,
    program_data = program_data,
    mutual_best = settings$mutual_best
  )

  loci <- .canonical_program_loci(decomps, scored$bridges)
  traits <- .correspondence_traits(decomps, loci)
  trait_pairs <- .correspondence_trait_pairs(decomps, loci, pairs)
  q_floor <- .permutation_q_floor(pairs, settings$n_perm)
  if (is.finite(q_floor) && q_floor > settings$fdr_threshold / 2) {
    warning(
      "Program pairs at the permutation floor (p = 1/", settings$n_perm + 1, ") have q = ",
      signif(q_floor, 3), " across ", sum(is.finite(pairs$q_concordance)), " scored pairs, ",
      "close to or above fdr_threshold = ", settings$fdr_threshold, "; increase n_perm.",
      call. = FALSE
    )
  }

  family_enrichment <- NULL
  if (isTRUE(annotate)) {
    family_enrichment <- .family_enrichment(families$nodes, decomps)
  }

  out <- list(
    traits = traits,
    trait_pairs = trait_pairs,
    pairs = pairs,
    families = families$families,
    nodes = families$nodes,
    family_loci = families$family_loci,
    edges = families$edges,
    graph = families$graph,
    modularity = families$modularity,
    n_unlinked = families$n_unlinked,
    loci = loci,
    family_enrichment = family_enrichment,
    q_floor = q_floor,
    scored = scored,
    bridges = scored$bridges,
    decompositions = decomps,
    settings = settings,
    call = NULL,
    gpmapr_version = as.character(utils::packageVersion("gpmapr")),
    created = Sys.time()
  )
  return(structure(out, class = "gpmap_correspondence"))
}


# Every program locus of every trait, on GPMap's shared coloc_group_id
# numbering (an upload's loci through its locus bridges).
.canonical_program_loci <- function(decomps, bridges) {
  rows <- lapply(decomps, function(d) {
    l <- d$loadings
    if (is.null(l) || nrow(l) == 0) {
      return(NULL)
    }
    l <- l[!is.na(l$coloc_group_id), , drop = FALSE]
    out <- dplyr::distinct(
      data.frame(
        trait_id = as.character(l$trait_id), trait_name = as.character(l$trait_name),
        coloc_group_id = as.character(l$coloc_group_id), high_confidence = l$high_confidence,
        stringsAsFactors = FALSE
      )
    )
    out$locus <- canonical_locus_ids(out$coloc_group_id, d$trait$id, bridges)
    return(out)
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    return(data.frame(
      trait_id = character(0), trait_name = character(0), coloc_group_id = character(0),
      high_confidence = logical(0), locus = character(0), stringsAsFactors = FALSE
    ))
  }
  return(out)
}


.correspondence_traits <- function(decomps, loci) {
  return(dplyr::bind_rows(lapply(decomps, function(d) {
    l <- loci[loci$trait_id == d$trait$id, , drop = FALSE]
    return(data.frame(
      trait_id = d$trait$id,
      trait_name = d$trait$name,
      n_programs = length(unique(d$loadings$program_id)),
      n_loci = length(unique(l$locus)),
      n_high_conf_loci = length(unique(l$locus[l$high_confidence])),
      stringsAsFactors = FALSE
    ))
  })))
}


.correspondence_trait_pairs <- function(decomps, loci, pairs) {
  ids <- names(decomps)
  loci_by_trait <- split(loci$locus, loci$trait_id)
  tier_levels <- c("high", "medium", "low", "none")
  return(dplyr::bind_rows(lapply(utils::combn(ids, 2, simplify = FALSE), function(tp) {
    in_pair <- (as.character(pairs$trait_id_a) == tp[1] & as.character(pairs$trait_id_b) == tp[2]) |
      (as.character(pairs$trait_id_a) == tp[2] & as.character(pairs$trait_id_b) == tp[1])
    p <- pairs[in_pair, , drop = FALSE]
    tiers <- table(factor(p$tier, levels = tier_levels))
    linked <- p[p$linked, , drop = FALSE]
    check <- if ("direction_check" %in% names(linked)) linked$direction_check else character(0)
    return(data.frame(
      trait_id_a = tp[1],
      trait_id_b = tp[2],
      trait_name_a = decomps[[tp[1]]]$trait$name,
      trait_name_b = decomps[[tp[2]]]$trait$name,
      n_loci_shared = length(intersect(loci_by_trait[[tp[1]]], loci_by_trait[[tp[2]]])),
      n_pairs = nrow(p),
      n_high = as.integer(tiers[["high"]]),
      n_medium = as.integer(tiers[["medium"]]),
      n_low = as.integer(tiers[["low"]]),
      n_none = as.integer(tiers[["none"]]),
      n_concordant = sum(linked$direction_profile == "concordant", na.rm = TRUE),
      n_antagonistic = sum(linked$direction_profile == "antagonistic", na.rm = TRUE),
      n_rg_agree = sum(check == "agree", na.rm = TRUE),
      n_rg_conflict = sum(check == "conflict", na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  })))
}


# q of the pairs whose permutation p-value sits at its floor, 1 / (n_perm + 1).
.permutation_q_floor <- function(pairs, n_perm) {
  floor_p <- 1 / (n_perm + 1)
  at_floor <- is.finite(pairs$p_concordance) & abs(pairs$p_concordance - floor_p) < 1e-12
  if (!any(at_floor)) {
    return(NA_real_)
  }
  return(min(pairs$q_concordance[at_floor], na.rm = TRUE))
}


# Per family and trait: the high-confidence SNPs of the family's programs in
# that trait, enriched against the trait's own SNPs for trait categories and
# tissues, and for pathways through their colocalizing QTL genes.
.family_enrichment <- function(nodes, decomps) {
  empty <- list(
    summary = data.frame(
      family = character(0), trait_id = character(0), trait_name = character(0),
      programs = character(0), n_snps = integer(0), status = character(0),
      stringsAsFactors = FALSE
    ),
    categories = NULL, tissues = NULL, pathways = NULL
  )
  if (is.null(nodes) || nrow(nodes) == 0) {
    return(empty)
  }
  groups <- split(nodes, list(nodes$family, nodes$trait_id), drop = TRUE)
  results <- lapply(groups, function(m) {
    tid <- as.character(m$trait_id[1])
    d <- decomps[[tid]]
    l <- d$loadings
    snps <- unique(l$snp_id[l$program_id %in% m$program_id & l$high_confidence])
    baseline <- unique(l$snp_id)
    info <- data.frame(
      family = m$family[1], trait_id = tid, trait_name = d$trait$name,
      programs = paste(sub("^.*:", "P", sort(m$program_id)), collapse = ", "),
      n_snps = length(snps), status = "ok", stringsAsFactors = FALSE
    )
    if (length(snps) < 5) {
      info$status <- "too few high-confidence SNPs"
      return(list(summary = info))
    }
    tag <- function(df) {
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      return(data.frame(family = info$family, trait_id = tid, df, stringsAsFactors = FALSE))
    }
    pw <- tryCatch(
      enrich_snp_set_pathways(snps, d$coloc_groups),
      error = function(e) return(NULL)
    )
    if (is.null(pw)) {
      info$status <- "pathways unavailable"
    }
    return(list(
      summary = info,
      categories = tag(enrich_snp_set(snps, baseline, d$coloc_groups, annotation = "trait_category")),
      tissues = tag(enrich_snp_set(snps, baseline, d$coloc_groups, annotation = "tissue")),
      pathways = if (is.null(pw)) NULL else tag(pw$pathways)
    ))
  })
  pick <- function(part) {
    out <- dplyr::bind_rows(lapply(results, function(r) return(r[[part]])))
    return(if (nrow(out) == 0) NULL else out)
  }
  summary <- pick("summary")
  family_num <- suppressWarnings(as.integer(sub("^family_", "", summary$family)))
  summary <- summary[order(family_num, summary$trait_name), , drop = FALSE]
  rownames(summary) <- NULL
  return(list(
    summary = summary,
    categories = pick("categories"),
    tissues = pick("tissues"),
    pathways = pick("pathways")
  ))
}
