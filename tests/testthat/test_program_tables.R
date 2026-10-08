test_that(".program_confidence_tier applies the similarity and connectedness cuts", {
  tier <- .program_confidence_tier(
    similarity_q = c(0.005, 0.005, 0.03, 0.05, 0.01, NA),
    connectedness = c(0.4, 0.39, 0.9, 0.9, 0.9, 0.9)
  )
  expect_true(is.ordered(tier))
  expect_equal(
    as.character(tier),
    c("high", "medium", "medium", "low", "medium", "low")
  )
})

test_that(".program_confidence_tier marks single-trait programs only when they pass the similarity test", {
  tier <- .program_confidence_tier(
    similarity_q = c(0.005, 0.005, 0.03, 0.5, 0.005, 0.03, NA),
    connectedness = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.1, 0.9),
    effective_n_traits = c(1.2, 1.5, NA, 1, 4, 1.1, 1),
    min_effective_traits = 1.5
  )
  # A single-trait program that fails the similarity test (or has no
  # estimate) is low; one that passes it at either level is single_trait.
  expect_equal(
    as.character(tier),
    c("single_trait", "high", "medium", "low", "high", "single_trait", "low")
  )
  expect_equal(levels(tier), c("low", "single_trait", "medium", "high"))
  # Without effective_n_traits the tier is the similarity tier alone.
  expect_equal(
    as.character(.program_confidence_tier(0.005, 0.9)),
    "high"
  )
})

make_program_table <- function() {
  return(data.frame(
    program = c(1L, 2L, 3L, 4L, 5L, 6L),
    confidence_tier = factor(
      c("medium", "high", "high", "low", "high", "medium"),
      levels = c("low", "medium", "high"), ordered = TRUE
    ),
    parent_program = c(NA, 1L, NA, 2L, NA, 99L),
    similarity_z = c(5, 9, 4, 1, 8, 7),
    n_snps_filtered = c(40, 10, 20, 3, 15, 6),
    stringsAsFactors = FALSE
  ))
}

test_that("order_programs_by_confidence puts children directly under their parent", {
  out <- order_programs_by_confidence(make_program_table())
  # Top level: high tier first by similarity_z (5, 3), then medium (6, whose
  # parent 99 is not in the table, then 1), with 1's child 2 and grandchild 4
  # directly under it.
  expect_equal(out$program, c(5L, 3L, 6L, 1L, 2L, 4L))
  expect_equal(out$depth, c(0L, 0L, 0L, 0L, 1L, 2L))
  expect_equal(out$program_label[1], "5")
  expect_equal(out$program_label[5], "  ↳ 2")
  expect_equal(out$program_label[6], "    ↳ 4")
  expect_equal(nrow(out), 6)
})

test_that("order_programs_by_confidence handles missing optional columns and empty input", {
  tab <- make_program_table()[, c("program", "confidence_tier")]
  out <- order_programs_by_confidence(tab)
  expect_equal(as.character(out$confidence_tier), c("high", "high", "high", "medium", "medium", "low"))
  expect_true(all(out$depth == 0L))
  empty <- order_programs_by_confidence(tab[0, ])
  expect_equal(nrow(empty), 0)
  expect_true(all(c("depth", "program_label") %in% names(empty)))
  expect_error(order_programs_by_confidence(data.frame(program = 1)), "confidence_tier")
})

test_that("order_programs_by_confidence sorts single_trait programs between medium and low", {
  tab <- data.frame(
    program = 1:4,
    confidence_tier = c("single_trait", "low", "high", "medium"),
    similarity_z = c(20, 1, 2, 3),
    stringsAsFactors = FALSE
  )
  out <- order_programs_by_confidence(tab)
  expect_equal(out$program, c(3L, 4L, 1L, 2L))
})

test_that("kable_confidence colours each body row by its tier", {
  df <- data.frame(
    program = c("1", "2", "3", "4"),
    confidence_tier = c("high", "medium", "low", "single_trait"),
    stringsAsFactors = FALSE
  )
  html <- as.character(kable_confidence(df, caption = "test"))
  expect_equal(lengths(regmatches(html, gregexpr("<tr style=", html, fixed = TRUE))), 4)
  expect_match(html, "background-color:#d9f2d9", fixed = TRUE)
  expect_match(html, "background-color:#ffe5b4", fixed = TRUE)
  expect_match(html, "background-color:#f8d0d0", fixed = TRUE)
  expect_match(html, "background-color:#e6e6e6", fixed = TRUE)
  expect_equal(attr(kable_confidence(df), "format"), "html")
  expect_error(kable_confidence(df, tier_col = "nope"), "tier_col")
})

test_that("kable_confidence can colour rows without showing the tier column", {
  df <- data.frame(
    program = c("1", "2"),
    confidence_tier = c("high", "medium"),
    stringsAsFactors = FALSE
  )
  html <- as.character(kable_confidence(df, drop_tier_col = TRUE))
  expect_equal(lengths(regmatches(html, gregexpr("<tr style=", html, fixed = TRUE))), 2)
  expect_match(html, "background-color:#d9f2d9", fixed = TRUE)
  expect_match(html, "background-color:#ffe5b4", fixed = TRUE)
  expect_false(grepl("confidence_tier", html, fixed = TRUE))
})
