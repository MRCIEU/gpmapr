test_that("remove_duplicate_traits() drops duplicates and keeps parents", {
  duplicates <- data.frame(
    trait_id = c(11L, 12L, 21L),
    trait_name = c("A copy 1", "A copy 2", "B copy"),
    parent_trait_id = c(1L, 1L, 2L),
    parent_trait_name = c("A", "A", "B")
  )
  data <- data.frame(trait_id = c(1, 11, 12, 2, 21, 3, 1), value = 1:7)

  result <- remove_duplicate_traits(data, duplicates)
  expect_equal(result$trait_id, c(1, 2, 3, 1))

  # a kept duplicate stays, and its parent goes instead
  result <- remove_duplicate_traits(data, duplicates, keep_trait_ids = "21")
  expect_equal(result$trait_id, c(1, 21, 3, 1))

  expect_error(remove_duplicate_traits(data.frame(id = 1), duplicates), "trait_id")
})

test_that("remove_duplicate_traits() reads duplicate groups with a representative", {
  # The same duplicates as above, in the duplicate-group shape: each group's
  # representative is the parent, and trait 3 is in a group of its own.
  duplicates <- data.frame(
    trait_id = c(1L, 11L, 12L, 2L, 21L, 3L),
    trait_name = c("A", "A copy 1", "A copy 2", "B", "B copy", "C"),
    duplicate_group = c(1L, 1L, 1L, 2L, 2L, 3L),
    is_representative = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE)
  )
  data <- data.frame(trait_id = c(1, 11, 12, 2, 21, 3, 1), value = 1:7)

  result <- remove_duplicate_traits(data, duplicates)
  expect_equal(result$trait_id, c(1, 2, 3, 1))
  result <- remove_duplicate_traits(data, duplicates, keep_trait_ids = "21")
  expect_equal(result$trait_id, c(1, 21, 3, 1))

  expect_error(
    remove_duplicate_traits(data, data.frame(trait_id = 1L, group = 1L)),
    "duplicate_group"
  )
})
