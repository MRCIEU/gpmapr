# Live tests (test_live_api.R) call the real GPMap API. They are skipped on
# CRAN, which has no reliable network access, and when the API is unreachable,
# so an outage reads as skipped rather than failed. They run locally and on
# GitHub Actions, where they catch API changes that leave the recorded
# fixtures behind.
#
# They call production by default.  Set GPMAP_URL to run them against another
# API, e.g. a local one: GPMAP_URL=http://localhost:8000
if (nzchar(Sys.getenv("GPMAP_URL"))) {
  options(gpmap_url = Sys.getenv("GPMAP_URL"))
}

#
# They are also skipped when this gpmapr doesn't support the API's version (e.g.
# production before an API release this gpmapr needs is deployed); the mocked
# tests in test_version_gate.R check which versions are supported.
skip_if_no_api <- function() {
  testthat::skip_on_cran()
  health <- tryCatch(health_api(), error = function(e) NULL)
  if (is.null(health)) {
    testthat::skip("GPMap API is not reachable")
  }
  unsupported <- tryCatch(ensure_api_compatible(), error = function(e) {
    return(conditionMessage(e))
  })
  if (is.character(unsupported)) {
    testthat::skip(unsupported)
  }
}
