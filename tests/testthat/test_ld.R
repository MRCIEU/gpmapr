library(testthat)

test_that("ld_proxies() and ld_matrix() add the r-squared", {
  local_mock_api()
  proxies <- ld_proxies(c(5553693, 5553694))$lds
  expect_equal(proxies$rsq, round(proxies$r^2, 2))
  matrix <- ld_matrix(c(5553693, 5553688))$lds
  expect_gt(nrow(matrix), 0)
  expect_equal(matrix$rsq, round(matrix$r^2, 2))
})

test_that("ld_matrix() returns an empty table when no pairs are stored", {
  local_mocked_bindings(api_get = function(url, ...) fake_api_response('{"lds": []}'))
  result <- ld_matrix(c(5553693, 80717))
  expect_equal(nrow(result$lds), 0)
  expect_equal(names(result$lds), c("lead_variant_id", "proxy_variant_id", "ld_block_id", "r", "rsq"))
})

test_that("ld_matrix() still reports other API errors", {
  local_mocked_bindings(api_get = function(url, ...) {
    return(stop_for_api_status(fake_api_response('{"detail": "No SNPs found"}', status = 400L)))
  })
  expect_error(ld_matrix(c("not-a-snp", "nor-this")), "HTTP 400.*No SNPs found")
})
