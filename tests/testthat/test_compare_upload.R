library(testthat)

# compare_upload() against the hand-written GWAS upload fixture (data-raw/build_upload_fixture.R): coloc group 1
# has the upload (AD GWAS, upload 90), allergic rhinitis (trait 2527), a compared upload (144) and type 1
# diabetes, whose p-value is smaller than all of theirs.

upload_with_associations <- function() {
  local_mock_api(parent.frame())
  return(get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155", include_associations = TRUE))
}

test_that("compare_upload() lines up the upload's own effect with the trait's in shared coloc groups", {
  comparison <- compare_upload(upload_with_associations(), 2527)
  groups <- comparison$coloc_groups
  expect_equal(nrow(groups), 1)
  expect_equal(groups$coloc_group_id, 1)
  expect_equal(groups$variant_id, 100)
  # The upload's own row, not type 1 diabetes, which has the smallest p-value in the group.
  expect_equal(groups$upload_study_extraction_id, 7223)
  expect_equal(groups$upload_beta, 0.5)
  expect_equal(groups$trait_name, "Allergic rhinitis")
  expect_equal(groups$trait_study_extraction_id, 32958)
  expect_equal(groups$trait_beta, 0.2)
})

test_that("compare_upload() returns the coloc pairs between the upload and the trait, upload first", {
  pairs <- compare_upload(upload_with_associations(), 2527)$coloc_pairs
  expect_equal(pairs$upload_study_extraction_id, c(7223, 7225))
  expect_equal(pairs$trait_study_extraction_id, c(32958, 40000))
  expect_equal(pairs$h4, c(0.9, 0.95))
  expect_equal(pairs$in_shared_coloc_group, c(TRUE, FALSE))
})

test_that("compare_upload() compares with another upload by its GUID", {
  comparison <- compare_upload(upload_with_associations(), "11111111-2222-3333-4444-555555555555")
  expect_equal(comparison$coloc_groups$trait_name, "Compared GWAS")
  expect_equal(comparison$coloc_groups$trait_beta, 0.45)
  expect_equal(comparison$coloc_pairs$trait_study_extraction_id, 8001)
})

test_that("compare_upload() skips rows whose association is a placeholder", {
  # Coloc group 2 pairs the upload with granulocyte percentage, whose association there is a placeholder.
  comparison <- compare_upload(upload_with_associations(), 21)
  expect_equal(nrow(comparison$coloc_groups), 0)
  expect_equal(comparison$coloc_pairs$trait_study_extraction_id, 7223)
})

test_that("compare_upload() returns empty tables for a trait the upload doesn't colocalise with", {
  comparison <- compare_upload(upload_with_associations(), 999999)
  expect_equal(nrow(comparison$coloc_groups), 0)
  expect_equal(nrow(comparison$coloc_pairs), 0)
})

test_that("compare_upload() explains what it needs", {
  local_mock_api()
  without_associations <- get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155")
  expect_error(compare_upload(without_associations, 2527), "include_associations = TRUE")
  expect_error(compare_upload(list(trait = list(id = 1)), 2527), "GWAS upload result")
  expect_error(compare_upload(upload_with_associations(), c(1, 2)), "single trait id")
  expect_error(
    compare_upload(upload_with_associations(), "99999999-2222-3333-4444-555555555555"),
    "not a GWAS upload this upload was compared with"
  )
})
