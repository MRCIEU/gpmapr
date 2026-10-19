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
#'       `trait_name`, `program`, `snp_id`, `coloc_group_id`, `loading`
#'       (sign-oriented), `abs_loading`, `lfsr`, `high_confidence`
#'     \item profiles: one row per feature x program with `program_id`,
#'       `trait_id`, `trait_name`, `program`, `feature_trait_id`,
#'       `feature_trait_name`, `loading` (sign-oriented), `lfsr`
#'     \item programs: one row per program with `program_id`, `trait_id`,
#'       `trait_name` and `program`
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
    dplyr::distinct(program_id, trait_id, trait_name, program) |>
    dplyr::arrange(program)

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
  return(
    info |>
      dplyr::transmute(
        snp_id = as.character(snp_id),
        coloc_group_id = as.character(coloc_group_id)
      ) |>
      dplyr::distinct(snp_id, .keep_all = TRUE)
  )
}


.attach_locus_meta <- function(loadings, locus_meta) {
  if (is.null(locus_meta) || nrow(locus_meta) == 0) {
    loadings$coloc_group_id <- NA_character_
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


# Each program's high-confidence loci (coloc_group_ids native to its trait).
.program_locus_sets <- function(loadings) {
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
    keep <- as.logical(sub$high_confidence)
    return(unique(as.character(sub$coloc_group_id[keep])))
  })
  names(sets) <- programs
  return(sets)
}


# How many high-confidence loci two programs share. set_a is native to
# trait_id_a; it is translated into trait_id_b's space before comparing against
# set_b (a no-op unless this pair needs LD bridging -- see
# .needs_locus_bridge()).
.locus_pair_stats <- function(key, locus_sets, bridges = NULL) {
  set_a <- locus_sets[[key$program_id_a]]
  set_b <- locus_sets[[key$program_id_b]]
  if (is.null(set_a)) set_a <- character(0)
  if (is.null(set_b)) set_b <- character(0)
  set_a <- .translate_locus_ids(set_a, key$trait_id_a, key$trait_id_b, bridges)

  return(data.frame(
    program_id_a = key$program_id_a,
    program_id_b = key$program_id_b,
    trait_id_a = key$trait_id_a,
    trait_id_b = key$trait_id_b,
    n_loci_shared = length(intersect(set_a, set_b)),
    stringsAsFactors = FALSE
  ))
}


.normalise_program_data <- function(program_data) {
  if (!is.list(program_data) || is.data.frame(program_data)) {
    stop("program_data must be an extract_program_loadings() result or a list of them")
  }
  data_list <- if ("loadings" %in% names(program_data) &&
                     is.data.frame(program_data$loadings)) {
    list(program_data)
  } else {
    program_data
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
#' Families are the communities of the link graph found by Louvain modularity
#' optimisation ([igraph::cluster_louvain()]), with each edge weighted by its
#' profile congruence `|phi_traits|`. Plain connected components chain through
#' any shared partner: on real data a few dozen links between four traits
#' joined most programs into one component, while their communities were well
#' separated. Louvain is randomised, so `seed` fixes it.
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
#'       whose edges carry `phi_traits`, `direction_profile`, `rg_direction`
#'       and `weight`, and whose vertices carry `family` and `trait_id`
#'     \item modularity: modularity of the family partition
#'     \item n_unlinked: number of programs with no link
#'   }
#' @export
build_program_families <- function(result,
                                   seed = 1,
                                   program_data = NULL,
                                   mutual_best = FALSE) {
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

  info <- dplyr::distinct(result$programs, program_id, trait_id, trait_name)
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
  for (col in c("phi_traits", "direction_profile", "rg_direction")) {
    graph <- igraph::set_edge_attr(graph, col, value = edges[[col]])
  }
  graph <- igraph::set_edge_attr(graph, "weight", value = abs(edges$phi_traits))

  set.seed(seed)
  membership <- as.integer(igraph::membership(igraph::cluster_louvain(graph)))
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


#' @title Canonical Locus Ids Across Traits
#' @description Put one trait's locus ids on an id space that is comparable
#' across traits: GPMap's shared `coloc_group_id`. An uploaded GWAS trait's own
#' ids are translated through any bridge to an existing trait; an upload locus
#' no bridge covers keeps a trait-prefixed id (`"<trait_id>:<id>"`), so it
#' counts towards a union but never matches another trait's locus. Ids of an
#' existing (non-upload) trait are returned unchanged.
#' @param ids The trait's own `coloc_group_id`s.
#' @param trait_id The trait the ids belong to.
#' @param bridges Locus bridges, normally
#'   `compare_program_pairs_loadings()$bridges`.
#' @return A character vector of canonical locus ids, one per element of `ids`.
#' @export
canonical_locus_ids <- function(ids, trait_id, bridges = NULL) {
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
      locus = canonical_locus_ids(s$coloc_group_id, tid, bridges),
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


.normalise_program_loadings <- function(loadings) {
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
    stringsAsFactors = FALSE
  )
  return(keys[keys$trait_id_a != keys$trait_id_b, , drop = FALSE])
}


# Confidence weight of a loading from its lFSR, in [0, 1]: (1 - 2 lfsr)_+, so 0
# for a coin-flip sign (lfsr = 0.5) and 1 for a certain one. Missing lFSR is
# treated as full confidence.
.confidence_factor <- function(lfsr) {
  conf <- 1 - 2 * lfsr
  conf[is.na(conf)] <- 1
  return(pmin(pmax(conf, 0), 1))
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
    stringsAsFactors = FALSE
  ))
}


.empty_pair_keys <- function() {
  return(data.frame(
    program_id_a = character(0),
    program_id_b = character(0),
    trait_id_a = character(0),
    trait_id_b = character(0),
    stringsAsFactors = FALSE
  ))
}


.empty_pairs <- function() {
  return(data.frame(
    program_id_a = character(0),
    program_id_b = character(0),
    trait_id_a = character(0),
    trait_id_b = character(0),
    n_loci_shared = integer(0),
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
