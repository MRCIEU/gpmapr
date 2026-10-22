library(testthat)

# Field names sent by API 1.0.1 are renamed to the API 1.1.0 names gpmapr returns
# (legacy_column_names and legacy_upload_column_names in R/common.R).

test_that("legacy upload coloc pair columns are renamed in place", {
  pairs <- data.frame(
    existing_study_extraction_id_a = c(NA, 46415),
    study_extraction_id_a = c(7293, NA),
    existing_study_extraction_id_b = c(1, NA),
    study_extraction_id_b = c(NA, 7296),
    h4 = c(0.9, 0.95)
  )
  result <- standardise_column_names(pairs)
  expect_equal(names(result), c(
    "existing_study_extraction_a_id", "upload_study_extraction_a_id",
    "existing_study_extraction_b_id", "upload_study_extraction_b_id", "h4"
  ))
  expect_equal(result$upload_study_extraction_a_id, c(7293, NA))
  expect_equal(result$existing_study_extraction_a_id, c(NA, 46415))
})

test_that("current names pass through, and win when both are sent", {
  current <- data.frame(upload_study_extraction_a_id = 1, ld_block_id = 9, pipeline_id = "ebi-a-1")
  expect_equal(standardise_column_names(current), current)

  both <- data.frame(study_extraction_id_a = 1, upload_study_extraction_a_id = 2)
  expect_equal(standardise_column_names(both), data.frame(upload_study_extraction_a_id = 2))
})

test_that("data frames nested in an api object are renamed, other elements untouched", {
  api_object <- list(
    trait = list(name = "My GWAS"),
    coloc_pairs = data.frame(study_extraction_id_a = 1, existing_study_extraction_id_b = 2)
  )
  result <- standardise_column_names(api_object)
  expect_equal(result$trait, api_object$trait)
  expect_equal(names(result$coloc_pairs), c("upload_study_extraction_a_id", "existing_study_extraction_b_id"))
})

test_that("trait tables get pipeline_id and variable_type, including nested study tables", {
  traits_table <- data.frame(id = 1, trait = "ebi-a-GCST000922", trait_name = "Openness", category = "continuous")
  traits_table$common_study <- data.frame(id = 10, category = "continuous")
  result <- standardise_column_names(traits_table)
  expect_equal(names(result), c("id", "pipeline_id", "trait_name", "variable_type", "common_study"))
  expect_equal(names(result$common_study), c("id", "variable_type"))
})

test_that("trait records get pipeline_id and variable_type, without renaming the api object's trait element", {
  api_object <- list(
    trait = list(id = 1, trait = "ebi-a-GCST000922", common_study = list(id = 10, category = "continuous")),
    coloc_groups = data.frame(study_extraction_id = 1)
  )
  result <- cleanup_api_object(api_object)
  expect_equal(names(result), c("trait", "coloc_groups"))
  expect_equal(result$trait$pipeline_id, "ebi-a-GCST000922")
  expect_null(result$trait$trait)
  expect_equal(result$trait$common_study$variable_type, "continuous")
})

test_that("study_extraction_id is renamed only in GWAS upload tables", {
  upload_groups <- data.frame(existing_study_extraction_id = c(NA, 5), study_extraction_id = c(7, NA))
  expect_equal(
    names(standardise_column_names(upload_groups)),
    c("existing_study_extraction_id", "upload_study_extraction_id")
  )

  coloc_groups <- data.frame(coloc_group_id = 1, study_extraction_id = 5)
  expect_equal(standardise_column_names(coloc_groups), coloc_groups)

  # traits() binds renamed upload rows onto normal coloc groups, which keep their study_extraction_id
  combined <- data.frame(
    study_extraction_id = c(5, NA), existing_study_extraction_id = c(NA, 6), upload_study_extraction_id = c(NA, 7)
  )
  expect_equal(standardise_column_names(combined), combined)
})

test_that("upload coloc groups say which side each row comes from", {
  groups <- data.frame(
    gwas_upload_id = 90, existing_study_extraction_id = c(NA, 5), upload_study_extraction_id = c(7, NA), h4 = 1
  )
  result <- add_upload_source(groups)
  expect_equal(result$source, c("upload", "existing"))
  expect_equal(names(result), c("gwas_upload_id", "existing_study_extraction_id", "upload_study_extraction_id",
                                "source", "h4"))
  normal <- data.frame(coloc_group_id = 1, study_extraction_id = 5)
  expect_equal(add_upload_source(normal), normal)
})

test_that("get_gwas() returns the current names from an API sending the old ones", {
  local_mock_api()
  result <- get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155")
  expect_true(all(c(
    "upload_study_extraction_a_id", "existing_study_extraction_a_id",
    "upload_study_extraction_b_id", "existing_study_extraction_b_id"
  ) %in% names(result$coloc_pairs)))
  expect_false(any(names(legacy_column_names) %in% names(result$coloc_pairs)))
  expect_true(all(c("upload_study_extraction_id", "existing_study_extraction_id", "source", "variable_type")
                  %in% names(result$coloc_groups)))
  expect_false(any(c("study_extraction_id", "category") %in% names(result$coloc_groups)))
  expect_equal(result$coloc_groups$source,
               ifelse(is.na(result$coloc_groups$upload_study_extraction_id), "existing", "upload"))
})

test_that("all_traits() and trait() return pipeline_id and variable_type", {
  local_mock_api()
  traits_table <- all_traits()
  expect_true(all(c("pipeline_id", "variable_type") %in% names(traits_table)))
  expect_false(any(c("trait", "category") %in% names(traits_table)))

  result <- trait(5020)
  expect_true(is.character(result$trait$pipeline_id))
  expect_equal(result$trait$source_url, create_source_url(result$trait$pipeline_id))
  expect_null(result$trait$trait)
})
