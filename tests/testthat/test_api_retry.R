library(testthat)

fake_response <- function(status, headers = list()) {
  return(structure(list(status_code = status, headers = headers), class = "response"))
}

test_that("request_with_backoff retries on 429 and returns the eventual success", {
  calls <- 0
  verb <- function(...) {
    calls <<- calls + 1
    if (calls <= 2) {
      return(fake_response(429L))
    }
    return(fake_response(200L))
  }

  result <- request_with_backoff(verb, base_delay = 0.001, max_delay = 0.002)

  expect_equal(httr::status_code(result), 200L)
  expect_equal(calls, 3)
})

test_that("request_with_backoff passes non-429 responses straight through", {
  calls <- 0
  verb <- function(...) {
    calls <<- calls + 1
    return(fake_response(404L))
  }

  result <- request_with_backoff(verb, base_delay = 0.001, max_delay = 0.002)

  expect_equal(httr::status_code(result), 404L)
  expect_equal(calls, 1)
})

test_that("request_with_backoff errors after exhausting retries on persistent 429s", {
  calls <- 0
  verb <- function(...) {
    calls <<- calls + 1
    return(fake_response(429L))
  }

  expect_error(
    request_with_backoff(verb, max_retries = 2, base_delay = 0.001, max_delay = 0.002),
    "rate limit exceeded"
  )
  expect_equal(calls, 3)
})

test_that("retry_delay_seconds honors a numeric Retry-After header", {
  response <- fake_response(429L, headers = list("retry-after" = "5"))
  expect_equal(retry_delay_seconds(0, response), 5)
})

test_that("retry_delay_seconds falls back to bounded exponential backoff without a header", {
  response <- fake_response(429L)
  delays <- vapply(0:4, function(attempt) {
    return(retry_delay_seconds(attempt, response, base_delay = 1, max_delay = 8))
  }, numeric(1))

  expect_true(all(delays >= 0))
  expect_true(all(delays <= 8))
})

test_that("api_get retries transparently on 429 via httr::GET", {
  calls <- 0
  fake_get <- function(...) {
    calls <<- calls + 1
    if (calls == 1) {
      return(fake_response(429L))
    }
    return(fake_response(200L))
  }
  testthat::local_mocked_bindings(GET = fake_get, .package = "httr")

  result <- api_get("http://example.com", base_delay = 0.001, max_delay = 0.002)
  expect_equal(httr::status_code(result), 200L)
  expect_equal(calls, 2)
})
