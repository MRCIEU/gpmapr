#' @title Decompose a Trait Into SNP Programs
#' @description Run the single-trait decomposition from one call: remove the
#' trait's proxy and duplicate traits, build its pleiotropy matrix, factorise
#' it with EBMF into overlapping SNP programs, give each program a confidence
#' tier, and (by default) annotate the programs. The result can be printed,
#' plotted ([plot_programs()] and friends), queried ([programs()],
#' [annotation()], [excluded()], [program_cards()]), passed to
#' [trait_correspondence()], or written up with [create_report()].
#'
#' **Which traits are removed.** Three kinds, all recorded in [excluded()]:
#' \describe{
#'   \item{Proxies of the anchor (`proxies = TRUE`)}{the anchor's own
#'   duplicate group and proxy group in GPMap ([anchor_proxies()]). An uploaded
#'   GWAS is not in GPMap's duplicate table: set `proxy_trait_id` to the GPMap
#'   trait it measures, whose groups (and which itself) are then removed.}
#'   \item{Hand additions (`exclude_trait_ids`)}{any other trait ids to remove
#'   as proxies of the anchor.}
#'   \item{Background duplicates (`remove_duplicates = TRUE`)}{every other
#'   repeated measurement in the background keeps one copy
#'   ([remove_duplicate_traits()]).}
#' }
#' The first two are also left out of every cross-trait profile comparison
#' involving this anchor in [trait_correspondence()].
#'
#' **Annotation** (`annotate = TRUE`) adds module-specific traits
#' ([summarise_module_specific_traits()]), loading-based enrichment of trait
#' categories, tissues and pathways ([enrich_program_loadings_trait_categories()],
#' [enrich_program_loadings_tissues()], [enrich_program_loadings_pathways()]),
#' program distinctness ([program_distinctness()]) and the colocalizing QTL
#' genes at each program's SNPs ([genes_at_snps()]). Only the pathway mappings
#' need the API; if they cannot be fetched, pathway enrichment is skipped.
#' @param trait A trait id (GPMap id or uploaded GWAS GUID), fetched with
#'   `trait(id, include_associations = TRUE)`, or a trait object of that shape
#'   (from [trait()] or [simulate_trait()]).
#' @param settings Analysis settings from [gpmap_settings()] (or a named list
#'   of overrides).
#' @param proxies If `TRUE` (default), remove the anchor's proxy traits found
#'   by [anchor_proxies()].
#' @param proxy_trait_id For an uploaded GWAS, the GPMap trait id it measures,
#'   used to find its proxies. Ignored when `proxies = FALSE`.
#' @param exclude_trait_ids Further trait ids to remove as proxies of the
#'   anchor.
#' @param remove_duplicates If `TRUE` (default), keep one copy of each
#'   duplicated background trait.
#' @param keep_trait_ids Trait ids [remove_duplicate_traits()] must not remove
#'   (the anchor is always kept). [trait_correspondence()] passes every anchor.
#' @param duplicates Optional [trait_duplicates()] table, fetched when needed
#'   and not supplied.
#' @param annotate If `TRUE` (default), annotate the programs.
#' @param mappings Optional [pathway_mappings()] result for pathway
#'   enrichment, fetched when needed and not supplied.
#' @param cores Cores for the replication refits. Defaults to half the
#'   machine's physical cores; results do not depend on it.
#' @param verbose Print progress messages.
#' @return An object of class `gpmap_decomposition`, a list with:
#'   \itemize{
#'     \item trait: `id` and `name`
#'     \item programs: one row per program of at least `min_module_size` SNPs,
#'       ordered by [order_programs_by_confidence()], with `program_id`,
#'       `trait_id` and `trait_name`
#'     \item memberships: the gated SNP x program memberships
#'     \item loadings, profiles: SNP and trait loadings of the high- and
#'       medium-confidence programs ([extract_program_loadings()]), and
#'       `extraction`, that result itself
#'     \item annotation: `specific_traits`, `categories`, `tissues`,
#'       `pathways`, `distinctness` and `genes`, or `NULL`
#'     \item excluded: the removed traits ([excluded()])
#'     \item exclude_trait_ids: the anchor's proxies (sources
#'       `duplicate_of_anchor`, `proxy_group`, `anchor_trait` and `user`)
#'     \item fit, summary: the [run_ebmf_pipeline()] results
#'     \item coloc_groups: the trait's coloc groups before any removal
#'     \item settings, call, gpmapr_version, created
#'   }
#' @examplesIf gpmap_available()
#' \donttest{
#' d <- decompose_trait("1631", cores = 1)
#' d
#' programs(d)
#' }
#' @concept trait decomposition
#' @export
decompose_trait <- function(trait,
                            settings = gpmap_settings(),
                            proxies = TRUE,
                            proxy_trait_id = NULL,
                            exclude_trait_ids = NULL,
                            remove_duplicates = TRUE,
                            keep_trait_ids = NULL,
                            duplicates = NULL,
                            annotate = TRUE,
                            mappings = NULL,
                            cores = .default_cores(),
                            verbose = TRUE) {
  call <- match.call()
  settings <- .as_gpmap_settings(settings)
  input <- .resolve_trait_input(trait)
  obj <- input$object
  target_id <- input$id
  target_name <- .trait_object_name(obj, target_id)
  raw_coloc_groups <- obj$coloc_groups
  if (is.null(raw_coloc_groups) || nrow(raw_coloc_groups) == 0) {
    stop("The trait has no coloc_groups to decompose", call. = FALSE)
  }

  if ((isTRUE(proxies) || isTRUE(remove_duplicates)) && is.null(duplicates)) {
    if (verbose) message("Fetching GPMap duplicate and proxy traits")
    duplicates <- trait_duplicates()
  }
  excluded <- .anchor_exclusions(
    target_id, raw_coloc_groups, duplicates,
    proxies = proxies, proxy_trait_id = proxy_trait_id,
    exclude_trait_ids = exclude_trait_ids
  )
  anchor_exclude <- excluded$trait_id
  cg <- raw_coloc_groups[!as.character(raw_coloc_groups$trait_id) %in% anchor_exclude, , drop = FALSE]

  if (isTRUE(remove_duplicates)) {
    before <- cg
    cg <- remove_duplicate_traits(
      before, duplicates,
      keep_trait_ids = unique(c(target_id, as.character(keep_trait_ids)))
    )
    excluded <- dplyr::bind_rows(
      excluded, .background_duplicate_rows(before, cg, target_id)
    )
  }
  obj$coloc_groups <- cg

  if (verbose) message("Fitting EBMF programs for ", target_name)
  pipeline <- run_ebmf_pipeline(
    obj,
    settings = .settings_for(settings, "pipeline"),
    target_trait_id = target_id,
    n_rep = settings$n_rep,
    seed = settings$seed,
    cores = cores,
    verbose = verbose
  )
  fit <- pipeline$fit
  summary <- pipeline$summary

  extraction <- extract_program_loadings(fit, summary, valid_only = TRUE)
  for (part in c("loadings", "profiles", "programs")) {
    if (nrow(extraction[[part]]) > 0) {
      extraction[[part]]$trait_id <- target_id
      extraction[[part]]$trait_name <- target_name
    }
  }

  programs <- summary$programs
  programs$program_id <- if (nrow(programs) > 0) paste0(target_id, ":", programs$program) else character(0)
  programs$trait_id <- rep(target_id, nrow(programs))
  programs$trait_name <- rep(target_name, nrow(programs))

  memberships <- summary$assigned
  if (!is.null(memberships) && nrow(memberships) > 0) {
    memberships <- memberships[memberships$program %in% programs$program, , drop = FALSE]
    rownames(memberships) <- NULL
  }

  annotation <- NULL
  if (isTRUE(annotate)) {
    if (verbose) message("Annotating programs")
    annotation <- .annotate_decomposition(fit, programs, memberships, target_id, settings, mappings)
  }

  out <- list(
    trait = list(id = target_id, name = target_name),
    programs = programs,
    memberships = memberships,
    loadings = extraction$loadings,
    profiles = extraction$profiles,
    extraction = extraction,
    annotation = annotation,
    excluded = excluded,
    exclude_trait_ids = anchor_exclude,
    fit = fit,
    summary = summary,
    coloc_groups = raw_coloc_groups,
    settings = settings,
    call = call,
    gpmapr_version = as.character(utils::packageVersion("gpmapr")),
    created = Sys.time()
  )
  return(structure(out, class = "gpmap_decomposition"))
}


#' @export
print.gpmap_decomposition <- function(x, ...) {
  p <- x$programs
  tiers <- c("high", "medium", "single_trait", "low")
  counts <- table(factor(as.character(p$confidence_tier), levels = tiers))
  cat("<gpmap_decomposition> ", x$trait$name, " (", x$trait$id, ")\n", sep = "")
  cat(
    "  ", ncol(x$fit$x_matrix), " SNPs x ", nrow(x$fit$x_matrix), " background traits; ",
    nrow(p), " programs: ",
    paste0(counts, " ", names(counts), collapse = ", "), "\n",
    sep = ""
  )
  if (nrow(x$excluded) > 0) {
    by_source <- table(x$excluded$source)
    cat(
      "  Removed traits: ",
      paste0(as.integer(by_source), " ", names(by_source), collapse = ", "), "\n",
      sep = ""
    )
  }
  kept <- p[as.character(p$confidence_tier) %in% c("high", "medium"), , drop = FALSE]
  if (nrow(kept) > 0) {
    top <- .top_specific_traits(x, n = 3)
    cat("  High and medium programs:\n")
    shown <- utils::head(kept, 10)
    for (i in seq_len(nrow(shown))) {
      traits <- top[as.character(shown$program[i])]
      cat(sprintf(
        "    P%-4s %-7s %4d SNPs  %s\n", shown$program[i], shown$confidence_tier[i],
        as.integer(shown$n_snps_filtered[i]),
        if (length(traits) == 0 || is.na(traits)) "" else traits
      ))
    }
    if (nrow(kept) > nrow(shown)) {
      cat("    ... and ", nrow(kept) - nrow(shown), " more\n", sep = "")
    }
  }
  if (is.null(x$annotation)) {
    cat("  Not annotated (annotate = FALSE).\n")
  }
  return(invisible(x))
}


#' @export
summary.gpmap_decomposition <- function(object, ...) {
  return(object$programs)
}


# A trait id (fetched) or a trait object, with the id the analysis keys it
# by: an upload's GUID, otherwise its trait id.
.resolve_trait_input <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (is.null(x$coloc_groups)) {
      stop("A trait object needs a coloc_groups table", call. = FALSE)
    }
    id <- if (!is.null(x$trait$guid)) x$trait$guid else x$trait$id
    if (is.null(id) || length(id) != 1) {
      stop("The trait object has no trait$id", call. = FALSE)
    }
    return(list(object = x, id = as.character(id)))
  }
  if (length(x) != 1 || is.na(x)) {
    stop("trait must be one trait id or a trait object", call. = FALSE)
  }
  id <- as.character(x)
  return(list(object = trait(id, include_associations = TRUE), id = id))
}


.trait_object_name <- function(obj, id) {
  for (nm in c(obj$trait$trait_name, obj$trait$name)) {
    if (!is.null(nm) && length(nm) == 1 && !is.na(nm) && nzchar(nm)) {
      return(as.character(nm))
    }
  }
  return(as.character(id))
}


# The anchor's proxies (from GPMap's groups) and hand additions, with how
# closely each tracks the anchor at its own loci. The anchor itself is never
# listed.
.anchor_exclusions <- function(target_id, coloc_groups, duplicates, proxies = TRUE,
                               proxy_trait_id = NULL, exclude_trait_ids = NULL) {
  table <- .empty_proxy_table()
  if (isTRUE(proxies)) {
    lookup <- if (is.null(proxy_trait_id)) target_id else as.character(proxy_trait_id)
    if (is.null(proxy_trait_id) && is_guid(target_id)) {
      warning(
        "Trait ", target_id, " is an uploaded GWAS, which has no GPMap proxy traits; ",
        "set proxy_trait_id to the GPMap trait it measures to remove its proxies.",
        call. = FALSE
      )
    }
    table <- .proxy_table(lookup, duplicates, include_self = lookup != target_id)
  }
  user <- setdiff(unique(as.character(exclude_trait_ids)), c(table$trait_id, target_id))
  if (length(user) > 0) {
    table <- dplyr::bind_rows(table, data.frame(
      trait_id = user, trait_name = .trait_names_in(coloc_groups, user),
      source = "user", rule = NA_character_, stringsAsFactors = FALSE
    ))
  }
  table <- table[table$trait_id != target_id, , drop = FALSE]
  table$anchor_id <- rep(target_id, nrow(table))
  table <- .add_proxy_evidence(table, coloc_groups, target_id)
  table$n_coloc_groups <- rep(NA_integer_, nrow(table))
  rownames(table) <- NULL
  return(table)
}


# Background traits remove_duplicate_traits() dropped, with the number of the
# anchor's coloc groups each was observed in.
.background_duplicate_rows <- function(before, after, target_id) {
  gone <- before[!as.character(before$trait_id) %in% as.character(after$trait_id), , drop = FALSE]
  if (nrow(gone) == 0) {
    return(NULL)
  }
  out <- gone |>
    dplyr::mutate(trait_id = as.character(trait_id)) |>
    dplyr::group_by(trait_id) |>
    dplyr::summarise(
      trait_name = as.character(trait_name[1]),
      n_coloc_groups = dplyr::n_distinct(coloc_group_id),
      .groups = "drop"
    ) |>
    dplyr::mutate(anchor_id = target_id, source = "background_duplicate", rule = NA_character_)
  return(as.data.frame(out, stringsAsFactors = FALSE))
}


.trait_names_in <- function(coloc_groups, ids) {
  if (is.null(coloc_groups) || !"trait_name" %in% names(coloc_groups)) {
    return(rep(NA_character_, length(ids)))
  }
  hit <- match(ids, as.character(coloc_groups$trait_id))
  return(as.character(coloc_groups$trait_name[hit]))
}


# Module-specific traits, loading-based enrichment, distinctness and genes of
# a fit's programs.
.annotate_decomposition <- function(fit, programs, memberships, target_id, settings, mappings) {
  out <- list(
    specific_traits = NULL, categories = NULL, tissues = NULL, pathways = NULL,
    pathways_status = "not run", distinctness = NULL, genes = .empty_program_genes()
  )
  if (nrow(programs) == 0) {
    return(out)
  }
  cg <- fit$coloc_groups
  groups <- data.frame(
    snp_id = as.character(memberships$snp_id), group = as.character(memberships$program),
    stringsAsFactors = FALSE
  )
  if (nrow(groups) > 0) {
    out$specific_traits <- summarise_module_specific_traits(
      trait_matrix = fit$trait_matrix,
      groups = groups,
      trait_info = fit$trait_info,
      coloc_groups = cg,
      exclude_trait_ids = target_id,
      min_group_size = settings$min_module_size
    )
  }
  out$categories <- enrich_program_loadings_trait_categories(
    clustering_result = fit, coloc_groups = cg, snp_key = "variant_id"
  )
  out$tissues <- enrich_program_loadings_tissues(
    clustering_result = fit, coloc_groups = cg, snp_key = "variant_id"
  )
  if (is.null(mappings)) {
    mappings <- tryCatch(pathway_mappings(), error = function(e) {
      message("Pathway mappings unavailable (", conditionMessage(e), "); skipping pathway enrichment")
      return(NULL)
    })
  }
  if (is.null(mappings)) {
    out$pathways_status <- "unavailable"
  } else {
    out$pathways <- enrich_program_loadings_pathways(
      clustering_result = fit, coloc_groups = cg, mappings = mappings,
      snp_key = "variant_id", sources = c("KEGG", "Reactome")
    )
    out$pathways_status <- "ok"
  }
  out$distinctness <- program_distinctness(fit, programs$program)
  if (nrow(groups) > 0) {
    genes <- genes_at_snps(unique(groups$snp_id), coloc_groups = cg, snp_key = "variant_id")
    if (nrow(genes) > 0) {
      out$genes <- memberships |>
        dplyr::transmute(program = as.character(program), snp_id = as.character(snp_id), abs_loading) |>
        dplyr::inner_join(
          genes |> dplyr::transmute(snp_id = as.character(snp_id), gene_id, gene) |> dplyr::distinct(),
          by = "snp_id", relationship = "many-to-many"
        ) |>
        dplyr::arrange(program, dplyr::desc(abs_loading)) |>
        as.data.frame(stringsAsFactors = FALSE)
    }
  }
  return(out)
}


.empty_program_genes <- function() {
  return(data.frame(
    program = character(0), snp_id = character(0), abs_loading = numeric(0),
    gene_id = integer(0), gene = character(0), stringsAsFactors = FALSE
  ))
}
