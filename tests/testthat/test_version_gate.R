library(testthat)

# The package/API version gate (ensure_api_compatible()), tested offline by
# mocking the version lookup.

reset_api_state <- function() {
  rm(list = ls(.api_state), envir = .api_state)
}

test_that("a supported API version passes and is cached for the session", {
  reset_api_state()
  on.exit(reset_api_state())
  lookups <- 0
  local_mocked_bindings(fetch_api_version = function(url) {
    lookups <<- lookups + 1
    return("1.0.1")
  })
  local_mocked_bindings(installed_package_version = function() "0.1.0")
  expect_true(ensure_api_compatible())
  expect_true(ensure_api_compatible())
  expect_equal(lookups, 1)
})

test_that("gpmapr 1.0.0 needs API 1.1.0, and older versions are told to update for 1.1.0", {
  reset_api_state()
  local_mocked_bindings(installed_package_version = function() "1.0.0")
  local_mocked_bindings(fetch_api_version = function(url) "1.1.0")
  expect_true(ensure_api_compatible())
  reset_api_state()
  local_mocked_bindings(fetch_api_version = function(url) "1.0.1")
  expect_error(ensure_api_compatible(), "gpmapr 1.0.0 does not support version 1.0.1")
  reset_api_state()
  on.exit(reset_api_state())
  local_mocked_bindings(fetch_api_version = function(url) "1.1.0")
  local_mocked_bindings(installed_package_version = function() "0.1.0")
  expect_error(ensure_api_compatible(), "gpmapr 0.1.0 does not support version 1.1.0")
})

test_that("an unsupported API version stops with upgrade instructions", {
  reset_api_state()
  on.exit(reset_api_state())
  local_mocked_bindings(fetch_api_version = function(url) "99.0.0")
  expect_error(ensure_api_compatible(), "does not support version 99.0.0.*install.packages")
})

test_that("an older gpmapr is told to update when the API has moved on", {
  reset_api_state()
  on.exit(reset_api_state())
  local_mocked_bindings(fetch_api_version = function(url) "1.0.1")
  local_mocked_bindings(installed_package_version = function() "0.0.0.9000")
  expect_error(ensure_api_compatible(), "gpmapr 0.0.0.9000 does not support version 1.0.1")
})

test_that("an unreachable version endpoint is not cached and does not block the request", {
  reset_api_state()
  on.exit(reset_api_state())
  local_mocked_bindings(fetch_api_version = function(url) stop("timeout"))
  expect_false(ensure_api_compatible())
  expect_null(.api_state$compatible)
})

test_that("each API URL is checked separately", {
  reset_api_state()
  on.exit(reset_api_state())
  old <- options(gpmap_url = "https://gpmap.opengwas.io/api")
  on.exit(options(old), add = TRUE)
  local_mocked_bindings(fetch_api_version = function(url) if (grepl("localhost", url)) "99.0.0" else "1.0.1")
  local_mocked_bindings(installed_package_version = function() "0.1.0")
  expect_true(ensure_api_compatible())
  options(gpmap_url = "http://localhost:8000")
  expect_error(ensure_api_compatible(), "99.0.0")
})

test_that("the version and health endpoints skip the gate", {
  local_mocked_bindings(ensure_api_compatible = function() stop("gate should not run"))
  local_mocked_bindings(GET = function(...) fake_api_response('{"version": "1.0.1", "status": "healthy"}'), .package = "httr")
  expect_equal(version_api(), "1.0.1")
  expect_equal(health_api()$status, "healthy")
})

test_that("loading the package makes no network call and keeps a chosen API", {
  local_mocked_bindings(GET = function(...) stop("no network at load"), .package = "httr")
  old <- options(gpmap_url = "http://localhost:8000")
  on.exit(options(old))
  expect_no_error(.onLoad(NULL, "gpmapr"))
  expect_equal(getOption("gpmap_url"), "http://localhost:8000")
  options(gpmap_url = NULL)
  .onLoad(NULL, "gpmapr")
  expect_equal(getOption("gpmap_url"), "https://gpmap.opengwas.io/api")
})

test_that("gpmap_available() is FALSE for an API version this gpmapr doesn't support", {
  reset_api_state()
  on.exit(reset_api_state())
  local_mocked_bindings(GET = function(...) fake_api_response('{"status": "healthy"}'), .package = "httr")
  local_mocked_bindings(installed_package_version = function() "1.0.0")
  local_mocked_bindings(fetch_api_version = function(url) "1.1.0")
  expect_true(gpmap_available())
  reset_api_state()
  local_mocked_bindings(fetch_api_version = function(url) "1.0.1")
  expect_false(gpmap_available())
})
