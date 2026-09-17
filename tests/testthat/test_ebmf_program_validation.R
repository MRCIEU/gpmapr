make_block_similarity <- function(n = 60, block = 1:15, rho = 0.9, seed = 1) {
  set.seed(seed)
  s <- matrix(0, n, n)
  s[block, block] <- rho
  noise <- matrix(stats::rnorm(n * n, 0, 0.02), n, n)
  noise <- (noise + t(noise)) / 2
  s <- s + noise
  diag(s) <- 1
  dimnames(s) <- list(paste0("snp", seq_len(n)), paste0("snp", seq_len(n)))
  return(s)
}

make_loadings <- function(n = 60, block = 1:15, seed = 1) {
  set.seed(seed)
  f <- matrix(stats::rnorm(n * 2, 0, 0.05), n, 2)
  f[block, 1] <- stats::rnorm(length(block), 1, 0.05)
  dimnames(f) <- list(paste0("snp", seq_len(n)), NULL)
  return(f)
}

test_that(".program_loading_coherence recovers a planted block", {
  s <- make_block_similarity()
  w <- rep(0, 60)
  w[1:15] <- 1

  on_block <- .program_loading_coherence(s, w, edge_threshold = 0.2)
  expect_gt(on_block$weighted_internal, 0.8)
  expect_equal(on_block$weighted_quorum, 1)
  expect_equal(on_block$n_eff, 15, tolerance = 1e-8)

  spread <- .program_loading_coherence(s, rep(1, 60), edge_threshold = 0.2)
  expect_lt(spread$weighted_internal, on_block$weighted_internal)
  expect_equal(spread$n_eff, 60, tolerance = 1e-8)
})

test_that(".program_loading_coherence is invariant to weight scaling", {
  s <- make_block_similarity()
  w <- abs(stats::rnorm(60))
  a <- .program_loading_coherence(s, w, edge_threshold = 0.2)
  b <- .program_loading_coherence(s, w * 1000, edge_threshold = 0.2)
  expect_equal(a$weighted_internal, b$weighted_internal)
  expect_equal(a$weighted_quorum, b$weighted_quorum)
  expect_equal(a$n_eff, b$n_eff)
})

test_that(".program_loading_coherence handles degenerate weights", {
  s <- make_block_similarity()
  expect_true(is.na(.program_loading_coherence(s, rep(0, 60), 0.2)$n_eff))
  point <- .program_loading_coherence(s, c(1, rep(0, 59)), 0.2)
  expect_equal(point$n_eff, 1)
  expect_true(is.na(point$weighted_internal))
})

test_that("calibrate_program_coherence separates a real program from noise", {
  s <- make_block_similarity()
  f <- make_loadings()
  res <- calibrate_program_coherence(s, f, n_perm = 499, edge_threshold = 0.2,
                                     seed = 42)

  expect_equal(nrow(res), 2)
  expect_equal(res$program, 1:2)
  # Program 1 loads on the planted block and sits at the permutation floor.
  # Program 2 is loading noise; its p is a draw from the null, so asserting it
  # clears any particular level would be flaky by construction -- the claim
  # worth testing is the separation, and the null graph test below covers the
  # false-positive rate properly.
  expect_equal(res$coherence_emp_p[1], 1 / 500)
  expect_lt(res$coherence_emp_p[1], res$coherence_emp_p[2])
  expect_lt(res$coherence_q[1], 0.05)
  expect_gt(res$weighted_internal[1], res$weighted_internal[2])
})

test_that("calibrate_program_coherence is uncalibrated-free on a null graph", {
  # A similarity matrix with no block structure: no loading vector should look
  # aligned with it, whatever its shape.
  set.seed(3)
  n <- 60
  s <- matrix(stats::rnorm(n * n, 0.35, 0.05), n, n)
  s <- (s + t(s)) / 2
  diag(s) <- 1
  dimnames(s) <- list(paste0("snp", seq_len(n)), paste0("snp", seq_len(n)))
  f <- matrix(stats::rnorm(n * 8), n, 8)
  rownames(f) <- paste0("snp", seq_len(n))

  res <- calibrate_program_coherence(s, f, n_perm = 499, seed = 7)
  expect_equal(sum(res$coherence_q < 0.05), 0)
  # p-values should not pile up at either end
  expect_gt(min(res$coherence_emp_p), 0.01)
})

test_that("calibrate_program_coherence batching matches a single chunk", {
  s <- make_block_similarity()
  f <- make_loadings()
  a <- calibrate_program_coherence(s, f, n_perm = 200, seed = 11,
                                   chunk_size = 200)
  b <- calibrate_program_coherence(s, f, n_perm = 200, seed = 11,
                                   chunk_size = 37)
  # Same seed and same total permutation count: the observed statistics are
  # identical and the p-values agree closely regardless of chunking.
  expect_equal(a$weighted_internal, b$weighted_internal)
  expect_equal(a$n_eff, b$n_eff)
  expect_equal(a$coherence_emp_p, b$coherence_emp_p, tolerance = 0.05)
})

test_that("calibrate_program_coherence handles empty and tiny inputs", {
  s <- make_block_similarity()
  expect_equal(nrow(calibrate_program_coherence(NULL, NULL)), 0)
  expect_equal(nrow(calibrate_program_coherence(s, matrix(0, 60, 0))), 0)

  tiny <- s[1:2, 1:2]
  f <- matrix(1, 2, 1, dimnames = list(rownames(tiny), NULL))
  expect_equal(nrow(calibrate_program_coherence(tiny, f)), 0)
})

test_that("resolve_ebmf_settings merges overrides over the shared file", {
  defaults <- read_ebmf_settings()
  expect_true(all(
    c("ebmf_prior", "ebmf_magnitude_threshold", "min_module_size",
      "coherence_q", "include_trans") %in% names(defaults)
  ))

  resolved <- resolve_ebmf_settings(list(
    ebmf_magnitude_threshold = 0.75,
    target_trait_id = 1992,
    ebmf_prior = NULL
  ))
  expect_equal(resolved$ebmf_magnitude_threshold, 0.75)
  # run-specific params do not leak into the shared set
  expect_false("target_trait_id" %in% names(resolved))
  # NULL falls back to the shared value
  expect_equal(resolved$ebmf_prior, defaults$ebmf_prior)
  expect_equal(sort(names(resolved)), sort(names(defaults)))
})

test_that("ebmf_settings_signature is stable and order-independent", {
  a <- resolve_ebmf_settings()
  b <- a[rev(names(a))]
  expect_equal(ebmf_settings_signature(a), ebmf_settings_signature(b))
  changed <- a
  changed$ebmf_magnitude_threshold <- 0.9
  expect_false(identical(
    ebmf_settings_signature(a), ebmf_settings_signature(changed)
  ))
})
