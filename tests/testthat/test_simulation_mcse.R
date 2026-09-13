library(testthat)

test_that("mcse_prop matches the ADEMP worked example", {
  x <- c(rep(1, 1520), rep(0, 80))
  expect_equal(length(x), 1600)
  expect_equal(round(mcse_prop(x), 4), 0.0054)
})

test_that("n_sim_for_mcse reproduces ADEMP sample-size rules", {
  expect_equal(n_sim_for_mcse(0.95, 0.005), 1900)
  expect_equal(n_sim_for_mcse(0.5, 0.005), 10000)
  expect_error(n_sim_for_mcse(1.5, 0.005), "between 0 and 1")
  expect_error(n_sim_for_mcse(0.5, 0), "positive")
})

test_that("mcse_mean and mcse_bias use sd/sqrt(n)", {
  x <- c(1, 2, 3, 4)
  expect_equal(mcse_mean(x), stats::sd(x) / 2)
  expect_equal(mcse_bias(x, truth = 0), stats::sd(x) / 2)
  expect_equal(mcse_bias(x, truth = 100), stats::sd(x) / 2)
  expect_true(is.na(mcse_mean(1)))
})

test_that("mcse_prop handles all-zero and all-one vectors", {
  expect_equal(mcse_prop(rep(0, 10)), 0)
  expect_equal(mcse_prop(rep(1, 10)), 0)
  expect_true(is.na(mcse_prop(numeric(0))))
})

test_that("ademp_summarise reports estimate, mcse and interval by group", {
  df <- data.frame(
    version = rep(c("a", "b"), each = 4),
    ari = c(0.8, 0.9, 0.85, 0.95, 0.1, 0.2, 0.15, 0.25),
    pass = c(1, 1, 0, 1, 0, 0, 0, 1)
  )
  out <- ademp_summarise(df, "ari", group_cols = "version", type = "mean")
  expect_equal(nrow(out), 2)
  expect_equal(out$n, c(4L, 4L))
  expect_equal(out$estimate[out$version == "a"], mean(df$ari[df$version == "a"]))
  expect_equal(out$mcse[out$version == "a"], stats::sd(df$ari[df$version == "a"]) / 2)
  expect_true(all(out$conf_low < out$estimate))
  expect_true(all(out$conf_high > out$estimate))

  prop_out <- ademp_summarise(df, "pass", group_cols = "version", type = "prop")
  expect_equal(prop_out$estimate[prop_out$version == "a"], 0.75)
  expect_equal(prop_out$mcse[prop_out$version == "a"], sqrt(0.75 * 0.25 / 4))
})

test_that("ademp_summarise supports overall and bias summaries", {
  df <- data.frame(x = c(2, 4, 6))
  overall <- ademp_summarise(df, "x", type = "mean")
  expect_equal(overall$n, 3L)
  expect_equal(overall$estimate, 4)

  bias <- ademp_summarise(df, "x", type = "bias", truth = 1)
  expect_equal(bias$estimate, 3)
  expect_equal(bias$mcse, stats::sd(df$x) / sqrt(3))

  expect_error(ademp_summarise(df, "missing"), "single column")
  expect_error(ademp_summarise(df, "x", type = "bias"), "truth")
})
