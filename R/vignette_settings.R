#' @title Locate the Shared EBMF Settings File
#' @description Find `_ebmf_settings.yml`, the canonical analysis and
#' program-validation settings shared by every investigation vignette. Searches
#' the working directory first, then the usual vignette-relative locations, so
#' the same call works from a vignette, from the package root, and from
#' `tests/testthat`.
#' @param path Optional explicit path to the settings file. When supplied it is
#'   used directly and must exist.
#' @return Path to the settings file.
#' @export
ebmf_settings_path <- function(path = NULL) {
  if (!is.null(path)) {
    if (!file.exists(path)) {
      stop("settings file not found: ", path)
    }
    return(path)
  }
  candidates <- c(
    "_ebmf_settings.yml",
    file.path("vignettes", "_ebmf_settings.yml"),
    file.path("..", "vignettes", "_ebmf_settings.yml"),
    file.path("..", "..", "vignettes", "_ebmf_settings.yml")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop(
      "could not locate _ebmf_settings.yml; looked in: ",
      paste(candidates, collapse = ", ")
    )
  }
  return(found[[1]])
}


#' @title Read the Shared EBMF Settings
#' @description Read `_ebmf_settings.yml` without applying any overrides.
#' @param path Optional explicit path, passed to `ebmf_settings_path()`.
#' @return Named list of settings.
#' @export
read_ebmf_settings <- function(path = NULL) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read the EBMF settings file")
  }
  settings <- yaml::read_yaml(ebmf_settings_path(path))
  if (!is.list(settings) || length(settings) == 0) {
    stop("_ebmf_settings.yml did not parse to a non-empty list")
  }
  return(settings)
}


#' @title Resolve Shared EBMF Settings Against Vignette Parameters
#' @description Merge a vignette's `params` over the shared settings in
#' `_ebmf_settings.yml`. Only the keys present in the settings file are
#' returned, so a vignette can pass its whole `params` list without its
#' run-specific entries (target trait, replicate counts, seeds) leaking into the
#' shared set. A parameter that is `NULL` or absent falls back to the shared
#' value, which is what makes the render-time overrides used by
#' `scripts/run_investigation_*_grid.sh` work.
#' @param overrides Named list of vignette parameters, typically `params`.
#'   Entries that are `NULL` are ignored.
#' @param path Optional explicit path, passed to `ebmf_settings_path()`.
#' @return Named list with one entry per key in the settings file.
#' @export
resolve_ebmf_settings <- function(overrides = list(), path = NULL) {
  defaults <- read_ebmf_settings(path)
  if (is.null(overrides)) {
    overrides <- list()
  }
  if (!is.list(overrides)) {
    stop("overrides must be a list")
  }
  supplied <- overrides[intersect(names(overrides), names(defaults))]
  supplied <- supplied[!vapply(supplied, is.null, logical(1))]
  supplied <- supplied[!vapply(supplied, function(x) {
    length(x) == 1 && is.character(x) && !is.na(x) && x == ""
  }, logical(1))]
  return(utils::modifyList(defaults, supplied))
}


#' @title Signature of a Resolved EBMF Settings List
#' @description A single deterministic string summarising the resolved shared
#' settings. Two rendered vignettes carrying the same signature were run on the
#' same analysis settings; this is the cross-study check that the simulation,
#' null-simulation, univariate and multi-trait studies have not drifted apart.
#' @param settings Result of `resolve_ebmf_settings()`.
#' @return Single character string of sorted `key=value` pairs.
#' @export
ebmf_settings_signature <- function(settings) {
  if (!is.list(settings) || is.null(names(settings))) {
    stop("settings must be a named list")
  }
  keys <- sort(names(settings))
  parts <- vapply(keys, function(k) {
    paste0(k, "=", paste(format(settings[[k]], trim = TRUE), collapse = ","))
  }, character(1))
  return(paste(parts, collapse = "; "))
}
