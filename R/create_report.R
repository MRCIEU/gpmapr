#' @title Write an HTML Report of a Result
#' @description Render a self-contained HTML report of an analysis result:
#' its settings, the traits removed, the tables and the standard figures, with
#' short notes on how to read each one. The report is built only from the
#' result, so it makes no API calls and renders quickly from a saved object.
#'
#' Methods exist for [decompose_trait()] results (programs, confidence,
#' distinctness, program cards and enrichment) and [trait_correspondence()]
#' results (links, directions, families and family enrichment, and with
#' `include_traits = TRUE` each trait's decomposition too).
#'
#' Rendering needs the 'rmarkdown' and 'knitr' packages and pandoc (bundled
#' with RStudio; otherwise see <https://pandoc.org/installing.html>).
#' @param x A result to report on.
#' @param output_file Path of the HTML file to write. Defaults to a temporary
#'   file.
#' @param ... Passed to the method.
#' @param open If `TRUE` (default in interactive sessions), open the report in
#'   a browser.
#' @param include_traits For a correspondence, also include each trait's
#'   decomposition. Defaults to `FALSE`.
#' @param quiet If `TRUE` (default), render without knitr's progress output.
#' @return The path of the written report, invisibly.
#' @examplesIf gpmap_available() && rmarkdown::pandoc_available()
#' \donttest{
#' d <- decompose_trait("1631", cores = 1)
#' create_report(d)
#' }
#' @concept trait decomposition
#' @export
create_report <- function(x, output_file = NULL, ..., open = interactive()) {
  UseMethod("create_report")
}


#' @rdname create_report
#' @export
create_report.default <- function(x, output_file = NULL, ..., open = interactive()) {
  supported <- sub("^create_report\\.", "", utils::methods("create_report"))
  supported <- setdiff(supported, "default")
  stop(
    "create_report() has no method for an object of class ", paste(class(x), collapse = "/"),
    ". Supported: ", paste(supported, collapse = ", "),
    call. = FALSE
  )
}


#' @rdname create_report
#' @export
create_report.gpmap_decomposition <- function(x, output_file = NULL, ..., open = interactive(),
                                              quiet = TRUE) {
  return(.render_report(
    "decomposition.Rmd",
    params = list(x = x),
    output_file = output_file,
    open = open,
    quiet = quiet
  ))
}


#' @rdname create_report
#' @export
create_report.gpmap_correspondence <- function(x, output_file = NULL, ..., open = interactive(),
                                               include_traits = FALSE, quiet = TRUE) {
  return(.render_report(
    "correspondence.Rmd",
    params = list(x = x, include_traits = isTRUE(include_traits)),
    output_file = output_file,
    open = open,
    quiet = quiet
  ))
}


# Render a template from inst/reports/ in a scratch copy of that directory,
# so nothing is written next to the installed package or into the working
# directory except output_file.
.render_report <- function(template, params, output_file = NULL, open = FALSE, quiet = TRUE) {
  for (pkg in c("rmarkdown", "knitr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("create_report() needs the '", pkg, "' package: install.packages(\"", pkg, "\")", call. = FALSE)
    }
  }
  if (!rmarkdown::pandoc_available()) {
    stop(
      "create_report() needs pandoc, which RStudio bundles; otherwise install it from ",
      "https://pandoc.org/installing.html",
      call. = FALSE
    )
  }
  source_dir <- system.file("reports", package = "gpmapr")
  if (!nzchar(source_dir) || !file.exists(file.path(source_dir, template))) {
    stop("Report template ", template, " not found in the gpmapr installation", call. = FALSE)
  }
  if (is.null(output_file)) {
    output_file <- tempfile(pattern = "gpmapr-report-", fileext = ".html")
  }
  output_file <- file.path(
    normalizePath(dirname(output_file), mustWork = TRUE), basename(output_file)
  )

  work_dir <- tempfile("gpmapr-report-")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  file.copy(list.files(source_dir, full.names = TRUE), work_dir, recursive = TRUE)

  rmarkdown::render(
    input = file.path(work_dir, template),
    output_file = basename(output_file),
    output_dir = dirname(output_file),
    intermediates_dir = work_dir,
    knit_root_dir = work_dir,
    params = params,
    envir = new.env(parent = globalenv()),
    quiet = quiet
  )
  if (isTRUE(open)) {
    utils::browseURL(output_file)
  }
  return(invisible(output_file))
}
