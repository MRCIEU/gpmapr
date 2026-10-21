# Mocked GPMap API for the offline tests (test_mocked_api.R).
#
# Recorded, trimmed API responses live in fixtures/, indexed by request URL
# relative to the API root (see data-raw/record_api_fixtures.R, which records
# them). local_mock_api() replaces api_get() / api_post() so that each request
# is answered from those fixtures, with no network access; the parsing in each
# *_api() function and the reshaping in each exported function still run.

fixture_path <- function(...) {
  return(testthat::test_path("fixtures", ...))
}

fake_api_response <- function(json, status = 200L, url = "https://gpmap.opengwas.io/api") {
  return(structure(
    list(
      url = url,
      status_code = as.integer(status),
      headers = list(`content-type` = "application/json"),
      content = charToRaw(enc2utf8(json))
    ),
    class = "response"
  ))
}

relative_api_url <- function(url) {
  return(sub("^https?://[^/]+(/api)?", "", url))
}

local_mock_api <- function(env = parent.frame()) {
  index <- jsonlite::fromJSON(fixture_path("index.json"))
  serve <- function(url, ...) {
    relative <- relative_api_url(url)
    file <- index[[relative]]
    if (is.null(file)) {
      stop("No recorded API fixture for ", relative, call. = FALSE)
    }
    json <- paste(readLines(fixture_path(file), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    return(fake_api_response(json, url = url))
  }
  testthat::local_mocked_bindings(api_get = serve, api_post = serve, .env = env)
}
