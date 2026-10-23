.onLoad <- function(libname, pkgname) {
  # No network access at load: the API version is checked on the first request
  # (see ensure_api_compatible()), so the package loads offline.
  if (is.null(getOption("gpmap_url"))) {
    options(gpmap_url = "https://gpmap.opengwas.io/api")
  }
}
