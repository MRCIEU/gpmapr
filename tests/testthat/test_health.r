library(testthat)

test_that("health() returns expected output", {
    # Example: Assuming health() returns a status message
    result <- health()
    expect_type(result, "list")
    expect_true(result$status %in% c("healthy", "unhealthy"))
})
