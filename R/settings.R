#' @title Analysis Settings for Trait Decomposition and Correspondence
#' @description Build the one list of settings that [decompose_trait()],
#' [trait_correspondence()] and [relink()] read. It holds every setting that
#' changes the answer; arguments of those functions only say which data go in
#' and how to run (`cores`, `verbose`).
#'
#' The settings fall into four groups:
#' \describe{
#'   \item{EBMF fit}{`include_trans`, `trait_subset`, `min_snp_signals`,
#'   `max_snp_fraction`, `compress_method`, `compress_scale`,
#'   `ebmf_greedy_Kmax`, `ebmf_lfsr_threshold`, `ebmf_magnitude_threshold`,
#'   `ebmf_drop_global`, `ebmf_prior`, `ebmf_backfit`, `ebmf_missing_se`: see
#'   [run_univariate_clustering()].}
#'   \item{Program validation}{`min_module_size`, `similarity_threshold`,
#'   `similarity_n_perm`, `similarity_q`, `high_similarity_q`,
#'   `high_connectedness`, `high_similarity_z`, `min_effective_traits`,
#'   `frac_traits`, `stability_top_n`, and `n_rep`, the number of
#'   trait-subsample refits behind the reported `replication` column (it does
#'   not change a program's tier): see [summarise_ebmf_programs()].}
#'   \item{Cross-trait links}{`n_perm`, the permutation replicates of the locus
#'   test, and the tier cutoffs `phi_high`, `phi_low`, `min_shared_loci`,
#'   `min_eff_loci`, `min_alignment` and `fdr_threshold`: see
#'   [compare_program_pairs_loadings()]. Every permutation p-value is at least
#'   `1 / (n_perm + 1)` and is corrected across every scored pair, so the more
#'   pairs are scored, the more replicates a link needs: [trait_correspondence()]
#'   warns when the floor comes close to `fdr_threshold`.}
#'   \item{Families}{`mutual_best`. By default (`FALSE`) a program may link to
#'   several programs of another trait, because EBMF can split one biological
#'   module into two programs within a trait. With `TRUE`, a link is kept only
#'   when each of its programs is the other's closest match (highest
#'   \eqn{|\phi|}) among the linked programs of the other trait, so each program
#'   keeps at most one partner per other trait and a family tends to hold one
#'   program per trait; the weaker half of a split module loses its link. See
#'   [build_program_families()].}
#' }
#' `seed` seeds every random step (similarity test, replication refits,
#' locus permutations, rg bootstrap and community detection).
#'
#' The defaults are the defaults of the functions these settings are passed to,
#' which are the shared analysis settings used throughout the package's
#' investigations.
#' @param ... Settings to override, as `name = value` pairs, or one or more
#'   named lists of them (for example from `resolve_ebmf_settings()`); later
#'   values win. `NULL` values are ignored. Unknown names are an error.
#' @return A named list of class `gpmap_settings` with every setting.
#' @examples
#' gpmap_settings()
#' gpmap_settings(n_perm = 3000, phi_high = 0.9)
#' @concept trait decomposition
#' @export
gpmap_settings <- function(...) {
  overrides <- .flatten_settings_args(list(...))
  defaults <- .gpmap_settings_defaults()
  unknown <- setdiff(names(overrides), names(defaults))
  if (length(unknown) > 0) {
    stop("Unknown setting(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  out <- utils::modifyList(defaults, overrides)
  for (nm in names(overrides)) {
    value <- out[[nm]]
    if (length(value) != 1 || is.na(value)) {
      stop("Setting `", nm, "` must be a single non-missing value", call. = FALSE)
    }
    if (is.numeric(defaults[[nm]]) && !is.numeric(value)) {
      stop("Setting `", nm, "` must be numeric", call. = FALSE)
    }
    if (is.logical(defaults[[nm]]) && !is.logical(value)) {
      stop("Setting `", nm, "` must be TRUE or FALSE", call. = FALSE)
    }
  }
  out <- out[names(defaults)]
  return(structure(out, class = c("gpmap_settings", "list")))
}


#' @export
print.gpmap_settings <- function(x, ...) {
  groups <- .gpmap_settings_groups()
  cat("<gpmap_settings>\n")
  for (g in names(groups)) {
    keys <- groups[[g]]
    cat("  ", g, ":\n", sep = "")
    for (k in keys) {
      cat("    ", format(k, width = 26), format(x[[k]]), "\n", sep = "")
    }
  }
  return(invisible(x))
}


# Settings keys by group, in display order.
.gpmap_settings_groups <- function() {
  map <- .ebmf_settings_map()
  return(list(
    "EBMF fit" = unname(map$clustering),
    "Program validation" = c(names(map$validation), "n_rep"),
    "Cross-trait links" = c(
      "n_perm", "phi_high", "phi_low", "min_shared_loci", "min_eff_loci",
      "min_alignment", "fdr_threshold"
    ),
    "Families" = "mutual_best",
    "Random seed" = "seed"
  ))
}


# Defaults read from the formals of the functions each setting is passed to,
# so they cannot drift from them. A match.arg() choice vector contributes its
# first element.
.gpmap_settings_defaults <- function() {
  first <- function(f, arg) {
    value <- eval(formals(f)[[arg]])
    return(if (is.character(value)) value[[1]] else value)
  }
  map <- .ebmf_settings_map()
  out <- list()
  for (k in names(map$clustering)) {
    out[[k]] <- first(run_univariate_clustering, map$clustering[[k]])
  }
  for (k in names(map$validation)) {
    out[[k]] <- first(summarise_ebmf_programs, map$validation[[k]])
  }
  out$n_rep <- first(summarise_ebmf_programs, "n_rep")
  for (k in c(
    "n_perm", "phi_high", "phi_low", "min_shared_loci", "min_eff_loci",
    "min_alignment", "fdr_threshold"
  )) {
    out[[k]] <- first(compare_program_pairs_loadings, k)
  }
  out$mutual_best <- first(build_program_families, "mutual_best")
  out$seed <- first(summarise_ebmf_programs, "seed")
  return(out[unlist(.gpmap_settings_groups(), use.names = FALSE)])
}


# Merge `...` of gpmap_settings(): unnamed list arguments are spliced in
# order, named arguments override them.
.flatten_settings_args <- function(args) {
  out <- list()
  if (length(args) == 0) {
    return(out)
  }
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }
  for (i in seq_along(args)) {
    if (!nzchar(arg_names[[i]])) {
      if (is.null(args[[i]])) {
        next
      }
      if (!is.list(args[[i]]) || (length(args[[i]]) > 0 && is.null(names(args[[i]])))) {
        stop("Unnamed arguments to gpmap_settings() must be named lists", call. = FALSE)
      }
      piece <- unclass(args[[i]])
    } else {
      piece <- args[i]
    }
    piece <- piece[!vapply(piece, is.null, logical(1))]
    out[names(piece)] <- piece
  }
  return(out)
}


# Coerce a settings argument (a gpmap_settings object or a plain list of
# overrides) to a complete gpmap_settings object.
.as_gpmap_settings <- function(settings) {
  if (inherits(settings, "gpmap_settings")) {
    return(settings)
  }
  if (is.null(settings)) {
    return(gpmap_settings())
  }
  if (!is.list(settings)) {
    stop("settings must come from gpmap_settings() or be a named list", call. = FALSE)
  }
  return(gpmap_settings(settings))
}


# The keys of one settings group, as named arguments of the step that takes
# them: "pipeline" (run_ebmf_pipeline()'s settings list), "link" (the tier
# cutoffs of compare_program_pairs_loadings() / assign_link_tiers()).
.settings_for <- function(settings, part = c("pipeline", "link")) {
  part <- match.arg(part)
  map <- .ebmf_settings_map()
  keys <- switch(part,
    pipeline = c(names(map$clustering), names(map$validation)),
    link = c(
      "phi_high", "phi_low", "min_shared_loci", "min_eff_loci",
      "min_alignment", "fdr_threshold"
    )
  )
  return(unclass(settings)[keys])
}


# Half the physical cores (at least 1). Cores only parallelise the n_rep
# replication refits, which seed themselves, so results do not depend on it.
.default_cores <- function(n = suppressWarnings(parallel::detectCores(logical = FALSE))) {
  if (length(n) != 1 || is.na(n)) {
    return(1L)
  }
  return(max(1L, as.integer(floor(n / 2))))
}
