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

skip_if_no_api <- function() {
  testthat::skip_on_cran()
  if (!gpmap_available()) {
    testthat::skip("GPMap API is not reachable")
  }
}
