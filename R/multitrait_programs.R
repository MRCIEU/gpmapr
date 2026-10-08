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
#' alternatively anchored on the sign of its largest-|loading| SNP. Because each
#' trait's matrix is also oriented to its own risk allele, two programs that
#' share a module have positive SNP loadings at its loci whether the two traits
#' move together or apart, so the sign of a cross-trait SNP-loading comparison
#' is not a direction (see [compare_program_pairs_loadings()] and
#' [module_rg()]). High-confidence membership follows the same lFSR / magnitude
#' gate used elsewhere in the pipeline.
#' @param clustering_result Result of [run_univariate_clustering()] (EBMF).
#' @param program_summary Optional result of [summarise_ebmf_programs()], used
#'   only to restrict to programs of high or medium confidence when
#'   `valid_only = TRUE`.
#' @param valid_only If `TRUE` (default) and `program_summary` is supplied, keep
#'   only programs of high or medium confidence (not `"low"` or
#'   `"single_trait"`). Otherwise all fitted programs are returned.
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
      as.character(program_summary$programs$confidence_tier) %in% c("high", "medium")
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


.locus_pair_stats <- function(key, locus_sets, universe_by_trait, bridges = NULL) {
  set_a <- locus_sets[[key$program_id_a]]
  set_b <- locus_sets[[key$program_id_b]]
  universe_b <- universe_by_trait[[key$trait_id_b]]
  if (is.null(set_a)) set_a <- character(0)
  if (is.null(set_b)) set_b <- character(0)
  if (is.null(universe_b)) universe_b <- character(0)
  n_a <- length(set_a)
  n_b <- length(set_b)

  # set_a is native to trait_id_a; translate into trait_id_b's space before
  # comparing against set_b/universe_b, which are already native to trait_id_b
  # (a no-op unless this pair needs LD bridging -- see .needs_locus_bridge()).
  # n_a/n_b stay as the untranslated native sizes: a program's own size
  # shouldn't depend on which partner it's being compared to.
  set_a <- .translate_locus_ids(set_a, key$trait_id_a, key$trait_id_b, bridges)

  shared <- intersect(set_a, set_b)
  union_size <- length(union(set_a, set_b))
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
#' @description Turn a program-correspondence result into multi-trait program
#' families. Programs are nodes and the `linked` pairs from
#' [compare_program_pairs_loadings()] are edges; programs with no link are left
#' out (their count is returned as `n_unlinked`).
#'
#' By default (`rule = "communities"`) families are the communities of the
#' link graph found by Louvain modularity optimisation
#' ([igraph::cluster_louvain()]), with each edge weighted by its profile
#' congruence `|phi_traits|`. Plain connected components (`rule = "components"`)
#' chain through any shared partner: on real data a few dozen links between
#' four traits joined most programs into one component, while their
#' communities were well separated. Louvain is randomised, so `seed` fixes it.
#'
#' With `mutual_best = TRUE` a link is kept only when each of its two programs
#' is the other's closest partner (highest `|phi_traits|`) among the linked
#' programs of the other's trait, so a program keeps at most one partner per
#' other trait and a family tends to hold one program per trait.
#'
#' When `program_data` is supplied, each family is also described by the loci
#' its programs claim as high-confidence members, on one locus id space across
#' traits (an uploaded GWAS trait's own `coloc_group_id`s are translated through
#' `result$bridges`, see [compare_program_pairs_loadings()]).
#' @param result A result from [compare_program_pairs_loadings()].
#' @param program_info Optional data.frame with `program_id`, `trait_id` and
#'   `trait_name` (one row per program), used to fill in trait names.
#' @param rule Family rule: `"communities"` (default) or `"components"`.
#' @param seed RNG seed for the community detection.
#' @param program_data Optional extraction result(s) from
#'   [extract_program_loadings()], needed for the family locus sets.
#' @param mutual_best If `TRUE`, keep only links that are each program's
#'   highest-`|phi_traits|` link to the other program's trait. Defaults to
#'   `FALSE`.
#' @return A list with:
#'   \itemize{
#'     \item nodes: one row per linked program with `program_id`, `trait_id`,
#'       `trait_name` and `family`
#'     \item families: one row per family with `family`, `n_programs`,
#'       `n_traits`, `traits`, `n_edges`, `mean_phi`, `n_tier_high`,
#'       `n_tier_medium`,
#'       `n_rg_concordant`, `n_rg_antagonistic`, `n_rg_undetermined` (edges by
#'       [module_rg()] direction; links with no estimate count in none), and
#'       `n_core_loci` (claimed by every member) and `n_union_loci` (by any
#'       member)
#'     \item family_loci: one row per family x locus with `family`, `locus`,
#'       `n_programs`, `n_traits`, `programs`, `traits` and `snp_ids`
#'     \item edges: the linked pairs, with the [module_rg()] columns joined on
#'     \item graph: an `igraph` object (or `NULL` when there are no links)
#'     \item modularity: modularity of the family partition
#'     \item n_unlinked: number of programs with no link
#'   }
#' @export
build_program_families <- function(result,
                                   program_info = NULL,
                                   rule = c("communities", "components"),
                                   seed = 1,
                                   program_data = NULL,
                                   mutual_best = FALSE) {
  rule <- match.arg(rule)
  if (!is.list(result) || is.data.frame(result) || is.null(result$pairs) ||
        !"linked" %in% names(result$pairs)) {
    stop("result must come from compare_program_pairs_loadings()")
  }
  pairs <- result$pairs
  edges <- pairs[!is.na(pairs$linked) & pairs$linked, , drop = FALSE]
  if (isTRUE(mutual_best) && nrow(edges) > 0) {
    edges <- edges[.mutual_best_links(edges), , drop = FALSE]
  }
  if (nrow(edges) > 0 && !is.null(result$module_rg) && nrow(result$module_rg) > 0) {
    rg_cols <- c(
      "program_id_a", "program_id_b", "rg", "rg_ci_lower", "rg_ci_upper",
      "n_loci_rg", "rg_direction"
    )
    edges <- dplyr::left_join(
      edges, result$module_rg[, rg_cols, drop = FALSE],
      by = c("program_id_a", "program_id_b")
    )
  }
  for (col in c("rg", "rg_ci_lower", "rg_ci_upper")) {
    if (!col %in% names(edges)) edges[[col]] <- rep(NA_real_, nrow(edges))
  }
  if (!"n_loci_rg" %in% names(edges)) edges$n_loci_rg <- rep(NA_integer_, nrow(edges))
  if (!"rg_direction" %in% names(edges)) {
    edges$rg_direction <- rep(NA_character_, nrow(edges))
  }

  endpoints <- if (!is.null(result$programs)) {
    dplyr::distinct(result$programs, program_id, trait_id, trait_name)
  } else {
    .normalise_program_info(pairs, NULL)
  }
  info <- .merge_program_info(endpoints, program_info)
  linked_ids <- unique(c(edges$program_id_a, edges$program_id_b))
  n_unlinked <- sum(!info$program_id %in% linked_ids)

  if (nrow(edges) == 0) {
    return(list(
      nodes = data.frame(
        program_id = character(0), trait_id = character(0),
        trait_name = character(0), family = character(0),
        stringsAsFactors = FALSE
      ),
      families = .empty_families(),
      family_loci = .empty_family_loci(),
      edges = edges,
      graph = NULL,
      modularity = NA_real_,
      n_unlinked = n_unlinked
    ))
  }

  graph <- igraph::graph_from_data_frame(
    edges[, c("program_id_a", "program_id_b"), drop = FALSE],
    directed = FALSE,
    vertices = data.frame(name = linked_ids, stringsAsFactors = FALSE)
  )
  for (col in c("phi_traits", "alignment", "tier", "profile_strength",
                "locus_status", "direction_profile",
                "concordance_z", "rg", "rg_ci_lower", "rg_ci_upper", "rg_direction")) {
    graph <- igraph::set_edge_attr(graph, col, value = edges[[col]])
  }
  graph <- igraph::set_edge_attr(graph, "weight", value = abs(edges$phi_traits))

  membership <- if (rule == "communities") {
    set.seed(seed)
    igraph::membership(igraph::cluster_louvain(graph))
  } else {
    igraph::components(graph, mode = "weak")$membership
  }
  membership <- as.integer(membership)
  # Number families largest first, so family_1 is the biggest.
  sizes <- table(membership)
  rank <- stats::setNames(
    seq_along(sizes),
    names(sizes)[order(-as.integer(sizes), as.integer(names(sizes)))]
  )
  family <- paste0("family_", rank[as.character(membership)])
  modularity <- igraph::modularity(graph, membership)

  nodes <- data.frame(
    program_id = igraph::V(graph)$name,
    family = family,
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      info |> dplyr::select(program_id, trait_id, trait_name),
      by = "program_id"
    ) |>
    dplyr::select(program_id, trait_id, trait_name, family)
  graph <- igraph::set_vertex_attr(graph, "family", value = nodes$family)
  graph <- igraph::set_vertex_attr(graph, "trait_id", value = nodes$trait_id)

  family_loci <- .empty_family_loci()
  if (!is.null(program_data)) {
    loadings <- .normalise_program_data(program_data)$loadings
    family_loci <- .family_loci(nodes, loadings, result$bridges)
  }
  families <- .family_metrics(nodes, edges, family_loci)
  return(list(
    nodes = nodes[order(nodes$family, nodes$trait_id, nodes$program_id), , drop = FALSE],
    families = families,
    family_loci = family_loci,
    edges = edges,
    graph = graph,
    modularity = modularity,
    n_unlinked = n_unlinked
  ))
}


#' @title Pooled Effect Correlation Per Family
#' @description [module_rg()] for whole program families: for every family
#' from [build_program_families()] and every pair of traits it spans, the two
#' traits' own effects are correlated across the union of the loci that the
#' family's programs in those two traits claim (on the loci both traits
#' carry). Each trait's loading at a locus is taken from whichever of its
#' family programs loads there most. Pooling gives a direction to trait pairs
#' whose individual links rest on too few loci for [module_rg()].
#' @param families A [build_program_families()] result.
#' @param program_data Extraction result(s) from [extract_program_loadings()].
#' @param coloc_groups Named list of coloc-group data.frames keyed by trait id,
#'   as for [module_rg()].
#' @param bridges Locus bridges, normally `compare_program_pairs_loadings()$bridges`.
#'   Any missing bridge is built.
#' @param n_boot,min_loci_rg,seed As in [module_rg()].
#' @return A data.frame with one row per family x trait pair: `family`,
#'   `trait_id_a`, `trait_id_b`, `n_programs_a`, `n_programs_b`, `rg`,
#'   `rg_ci_lower`, `rg_ci_upper`, `n_loci_rg`, `n_allele_mismatch` and
#'   `rg_direction`.
#' @export
family_module_rg <- function(families,
                             program_data,
                             coloc_groups,
                             bridges = NULL,
                             n_boot = 1000L,
                             min_loci_rg = 10L,
                             seed = 1,
                             ld_proxies_fn = ld_proxies) {
  nodes <- families$nodes
  if (is.null(nodes) || nrow(nodes) == 0) {
    return(.empty_family_rg())
  }
  loadings <- .normalise_program_data(program_data)$loadings
  combos <- dplyr::bind_rows(lapply(split(nodes, nodes$family), function(fam) {
    traits <- sort(unique(as.character(fam$trait_id)))
    if (length(traits) < 2) {
      return(NULL)
    }
    tp <- utils::combn(traits, 2)
    return(data.frame(
      family = fam$family[1], trait_id_a = tp[1, ], trait_id_b = tp[2, ],
      stringsAsFactors = FALSE
    ))
  }))
  if (nrow(combos) == 0) {
    return(.empty_family_rg())
  }
  bridges <- .complete_bridges(
    unique(combos[, c("trait_id_a", "trait_id_b")]), coloc_groups, bridges,
    ld_proxies_fn = ld_proxies_fn
  )
  effects <- .effects_cache(coloc_groups)

  set.seed(seed)
  rows <- lapply(seq_len(nrow(combos)), function(i) {
    fam <- nodes[nodes$family == combos$family[i], , drop = FALSE]
    ta <- combos$trait_id_a[i]
    tb <- combos$trait_id_b[i]
    pa <- fam$program_id[fam$trait_id == ta]
    pb <- fam$program_id[fam$trait_id == tb]
    est <- .rg_for_programs(
      loadings, pa, pb, ta, tb, effects, bridges,
      n_boot = n_boot, min_loci_rg = min_loci_rg
    )
    return(data.frame(
      combos[i, , drop = FALSE],
      n_programs_a = length(pa), n_programs_b = length(pb),
      est, stringsAsFactors = FALSE
    ))
  })
  out <- dplyr::bind_rows(rows)
  rownames(out) <- NULL
  return(out)
}


# Links that are the closest (highest |phi_traits|) link of both of their
# programs to the other program's trait.
.mutual_best_links <- function(edges) {
  r2 <- abs(edges$phi_traits)
  r2[!is.finite(r2)] <- -Inf
  is_best <- function(program, partner_trait) {
    key <- paste(program, partner_trait)
    best <- tapply(r2, key, max)
    return(r2 == best[key])
  }
  return(
    is_best(edges$program_id_a, edges$trait_id_b) &
      is_best(edges$program_id_b, edges$trait_id_a)
  )
}


# A locus id that is comparable across traits: GPMap's shared coloc_group_id.
# An uploaded GWAS trait's own ids are translated through any bridge to an
# existing trait; an upload locus no bridge covers keeps a trait-prefixed id,
# so it counts towards a union but never matches another trait's locus.
.canonical_locus_ids <- function(ids, trait_id, bridges = NULL) {
  ids <- as.character(ids)
  trait_id <- as.character(trait_id)
  if (!is_guid(trait_id)) {
    return(ids)
  }
  out <- rep(NA_character_, length(ids))
  for (bridge in bridges) {
    if (is.null(bridge) || nrow(bridge) == 0) next
    partners <- setdiff(c(bridge$trait_a[1], bridge$trait_b[1]), trait_id)
    if (!trait_id %in% c(bridge$trait_a[1], bridge$trait_b[1]) ||
          length(partners) != 1 || is_guid(partners)) {
      next
    }
    todo <- is.na(out)
    out[todo] <- .translate_locus_ids_positional(ids[todo], trait_id, partners, bridges)
  }
  out[is.na(out)] <- paste0(trait_id, ":", ids[is.na(out)])
  return(out)
}


# Every family's claimed loci on the canonical locus id space.
.family_loci <- function(nodes, loadings, bridges = NULL) {
  rows <- lapply(seq_len(nrow(nodes)), function(i) {
    pid <- nodes$program_id[i]
    tid <- as.character(nodes$trait_id[i])
    loci <- .program_claimed_loci(loadings, pid)
    if (length(loci) == 0) {
      return(NULL)
    }
    s <- loadings[
      as.character(loadings$program_id) == pid &
        as.character(loadings$coloc_group_id) %in% loci,
      ,
      drop = FALSE
    ]
    return(data.frame(
      family = nodes$family[i],
      program_id = pid,
      trait_id = tid,
      locus = .canonical_locus_ids(s$coloc_group_id, tid, bridges),
      snp_id = as.character(s$snp_id),
      stringsAsFactors = FALSE
    ))
  })
  claims <- dplyr::bind_rows(rows)
  if (nrow(claims) == 0) {
    return(.empty_family_loci())
  }
  out <- claims |>
    dplyr::group_by(family, locus) |>
    dplyr::summarise(
      n_programs = dplyr::n_distinct(program_id),
      n_traits = dplyr::n_distinct(trait_id),
      programs = paste(sort(unique(program_id)), collapse = ", "),
      traits = paste(sort(unique(trait_id)), collapse = ", "),
      snp_ids = paste(sort(unique(snp_id)), collapse = ", "),
      .groups = "drop"
    ) |>
    dplyr::arrange(family, dplyr::desc(n_programs), locus)
  return(as.data.frame(out, stringsAsFactors = FALSE))
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


# Confidence weight of a loading from its lFSR, in [0, 1]. The default,
# "two_sided", is (1 - 2 lfsr)_+: 0 for a coin-flip sign (lfsr = 0.5), 1 for a
# certain one. "one_sided" is 1 - lfsr, which still gives a coin flip half
# weight. Set with options(gpmapr.lfsr_weight = ...). Missing lFSR is treated
# as full confidence.
.confidence_factor <- function(lfsr) {
  conf <- if (identical(.lfsr_weight_method(), "one_sided")) {
    1 - lfsr
  } else {
    1 - 2 * lfsr
  }
  conf[is.na(conf)] <- 1
  return(pmin(pmax(conf, 0), 1))
}


.lfsr_weight_method <- function() {
  method <- getOption("gpmapr.lfsr_weight", "two_sided")
  if (!method %in% c("two_sided", "one_sided")) {
    stop("options(gpmapr.lfsr_weight) must be \"two_sided\" or \"one_sided\"")
  }
  return(method)
}


# Weighted correlation through the origin, sum(w x y) / sqrt(sum(w x^2)
# sum(w y^2)): unlike a mean-centred Pearson correlation it does not change
# when the coding allele of any one observation is flipped (x and y both
# negated). NA when either weighted sum of squares is degenerate.
.weighted_uncentred_cor <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w)
  x <- x[ok]
  y <- y[ok]
  w <- w[ok]
  vx <- sum(w * x * x)
  vy <- sum(w * y * y)
  if (length(x) == 0 || !is.finite(vx) || !is.finite(vy) || vx <= 0 || vy <= 0) {
    return(NA_real_)
  }
  return(sum(w * x * y) / sqrt(vx * vy))
}


.family_metrics <- function(nodes, edges, family_loci = .empty_family_loci()) {
  families <- sort(unique(nodes$family))
  rows <- lapply(families, function(fam) {
    members <- nodes[nodes$family == fam, , drop = FALSE]
    pids <- members$program_id
    fam_edges <- edges[
      edges$program_id_a %in% pids & edges$program_id_b %in% pids,
      ,
      drop = FALSE
    ]
    loci <- family_loci[family_loci$family == fam, , drop = FALSE]
    has_loci <- nrow(family_loci) > 0
    trait_names <- if ("trait_name" %in% names(members)) {
      dplyr::coalesce(members$trait_name, members$trait_id)
    } else {
      members$trait_id
    }
    return(data.frame(
      family = fam,
      n_programs = length(pids),
      n_traits = dplyr::n_distinct(members$trait_id),
      traits = paste(sort(unique(trait_names)), collapse = ", "),
      n_edges = nrow(fam_edges),
      mean_phi = .safe_stat(fam_edges$phi_traits, mean),
      n_tier_high = sum(fam_edges$tier == "high", na.rm = TRUE),
      n_tier_medium = sum(fam_edges$tier == "medium", na.rm = TRUE),
      n_rg_concordant = sum(fam_edges$rg_direction == "concordant", na.rm = TRUE),
      n_rg_antagonistic = sum(fam_edges$rg_direction == "antagonistic", na.rm = TRUE),
      n_rg_undetermined = sum(fam_edges$rg_direction == "undetermined", na.rm = TRUE),
      n_core_loci = if (has_loci) sum(loci$n_programs == length(pids)) else NA_integer_,
      n_union_loci = if (has_loci) nrow(loci) else NA_integer_,
      stringsAsFactors = FALSE
    ))
  })
  out <- dplyr::bind_rows(rows)
  out <- out[order(as.integer(sub("^family_", "", out$family))), , drop = FALSE]
  rownames(out) <- NULL
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
    eligible = logical(0),
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
    mean_phi = numeric(0),
    n_tier_high = integer(0),
    n_tier_medium = integer(0),
    n_rg_concordant = integer(0),
    n_rg_antagonistic = integer(0),
    n_rg_undetermined = integer(0),
    n_core_loci = integer(0),
    n_union_loci = integer(0),
    stringsAsFactors = FALSE
  ))
}


.empty_family_loci <- function() {
  return(data.frame(
    family = character(0),
    locus = character(0),
    n_programs = integer(0),
    n_traits = integer(0),
    programs = character(0),
    traits = character(0),
    snp_ids = character(0),
    stringsAsFactors = FALSE
  ))
}


.empty_family_rg <- function() {
  return(data.frame(
    family = character(0),
    trait_id_a = character(0),
    trait_id_b = character(0),
    n_programs_a = integer(0),
    n_programs_b = integer(0),
    rg = numeric(0),
    rg_ci_lower = numeric(0),
    rg_ci_upper = numeric(0),
    n_loci_rg = integer(0),
    n_allele_mismatch = integer(0),
    rg_direction = character(0),
    stringsAsFactors = FALSE
  ))
}
