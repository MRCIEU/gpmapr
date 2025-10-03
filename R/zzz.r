.onLoad <- function(libname, pkgname) {
  select_api("production")
  api_version <- version_api()
  package_version <- utils::packageVersion("gpmapr")

  if (!as.character(package_version) %in% api_to_package_version[[api_version]]) {
    warning("WARNING: Version mismatch between package and API. \n\nPlease update the package to the latest version.\n")
  }
}
