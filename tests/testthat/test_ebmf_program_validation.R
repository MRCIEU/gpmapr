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

assigned_of <- function(sets) {
  return(data.frame(
    snp_id = unlist(sets),
    program = rep(seq_along(sets), lengths(sets)),
    stringsAsFactors = FALSE
  ))
}

# Member SNP ids per program, named by program id.
members_of_sets <- function(sets) {
  return(stats::setNames(lapply(sets, as.character), seq_along(sets)))
}

test_that(".program_similarity_null separates a planted block from a random set", {
  s <- make_block_similarity()
  set.seed(5)
  members <- members_of_sets(list(
    paste0("snp", 1:6),
    paste0("snp", sample(16:60, 6))
  ))
  res <- .program_similarity_null(s, members, n_perm = 499, seed = 42)

  expect_equal(res$program, 1:2)
  # The planted block sits at the permutation floor; the random set is a draw
  # from the null, so only the separation is asserted for it.
  expect_equal(res$similarity_emp_p[1], 1 / 500)
  expect_lt(res$similarity_emp_p[1], res$similarity_emp_p[2])
  expect_gt(res$excess_similarity[1], 0.5)
  expect_gt(res$excess_connectedness[1], 0.5)
  expect_lt(abs(res$excess_similarity[2]), 0.2)
  expect_gt(res$similarity_z[1], 5)
})

test_that("excess metrics do not grow with program size", {
  # Two programs drawn from the same tight block, one small and one large, have
  # the same excess similarity; the larger one only gets a larger z.
  s <- make_block_similarity(block = 1:30)
  members <- members_of_sets(list(paste0("snp", 1:4), paste0("snp", 1:30)))
  res <- .program_similarity_null(s, members, n_perm = 499, seed = 1)
  expect_equal(res$excess_similarity[1], res$excess_similarity[2], tolerance = 0.05)
  expect_equal(res$excess_connectedness[1], res$excess_connectedness[2], tolerance = 0.05)
  expect_gt(res$similarity_z[2], res$similarity_z[1])
})

test_that("connectedness counts positive similarity only", {
  s <- make_block_similarity()
  s[1:3, 4:6] <- -0.9
  s[4:6, 1:3] <- -0.9
  coh <- .program_internal_coherence(s, paste0("snp", 1:6), edge_threshold = 0.2)
  # 6 of 15 pairs are within {1,2,3} or {4,5,6} (positive); 9 are between (negative).
  expect_equal(coh$connectedness, 6 / 15)
})

test_that(".program_similarity_null handles tiny programs and SNPs with no profile", {
  s <- make_block_similarity()
  # snp60 has no similarity to anything (e.g. observed only for the target row)
  s["snp60", ] <- 0
  s[, "snp60"] <- 0
  s["snp60", "snp60"] <- 1
  members <- members_of_sets(list("snp1", paste0("snp", 1:3)))
  res <- .program_similarity_null(s, members, n_perm = 99, seed = 1)
  expect_true(is.na(res$similarity_emp_p[1]))
  expect_true(is.na(res$excess_similarity[1]))
  expect_true(is.finite(res$similarity_emp_p[2]))
  expect_equal(nrow(.program_similarity_null(s, list())), 0)
})

test_that(".similarity_background summarises all SNP pairs", {
  s <- make_block_similarity()
  bg <- .similarity_background(s, NULL, edge_threshold = 0.2)
  vals <- s[upper.tri(s)]
  expect_equal(bg$mean_similarity, mean(vals))
  expect_equal(bg$frac_pairs_edge, mean(vals >= 0.2))
  expect_true(is.na(bg$median_shared_traits))
})

test_that(".program_parents links a nested program to the larger one", {
  snp_info <- data.frame(
    snp_id = paste0("s", 1:30),
    chr = rep(c("1", "2", "3"), each = 10),
    bp = rep(seq(1e6, 1e7, length.out = 10), 3),
    stringsAsFactors = FALSE
  )
  assigned <- assigned_of(list(
    paste0("s", 1:10),       # 1: chr1, 10 SNPs
    paste0("s", c(1, 2, 3)), # 2: nested inside 1
    paste0("s", 21:24)       # 3: chr3, disjoint
  ))
  par <- .program_parents(assigned, snp_info)
  expect_equal(par$program, 1:3)
  expect_equal(par$parent_program, c(NA, 1L, NA))
  expect_equal(par$parent_containment[2], 1)

  # Nearby but not shared SNPs still count within the window
  near <- snp_info
  near$bp[near$snp_id == "s21"] <- near$bp[near$snp_id == "s1"] + 1e5
  near$chr[near$snp_id == "s21"] <- "1"
  par_near <- .program_parents(assigned_of(list(paste0("s", 1:10), c("s21", "s2"))), near)
  expect_equal(par_near$parent_program, c(NA, 1L))

  # No positions -> no parents
  expect_true(all(is.na(.program_parents(assigned, NULL)$parent_program)))

  # 3 of 5 SNPs inside the larger program (60%) is below the 0.9 cut
  partial <- assigned_of(list(paste0("s", 1:10), c("s1", "s2", "s3", "s21", "s25")))
  expect_true(all(is.na(.program_parents(partial, snp_info)$parent_program)))
  expect_equal(
    .program_parents(partial, snp_info, min_containment = 0.6)$parent_program,
    c(NA, 1L)
  )
})

test_that(".program_effective_traits counts the traits a program loads on", {
  L <- matrix(
    c(5, 2, 0, 0,   # 1: target plus one trait
      0, 1, 1, 1,   # 2: three equal traits
      3, 0, 0, 0,   # 3: target only
      0, 2, 1, 0),  # 4: two unequal traits
    nrow = 4,
    dimnames = list(c("100", "2", "3", "4"), NULL)
  )
  fit <- list(L_pm = L)
  eff <- .program_effective_traits(fit, 1:5, target_trait_id = 100)
  expect_equal(eff[1:2], c(1, 3))
  expect_true(is.na(eff[3]))
  expect_equal(eff[4], 1 / ((4 / 5)^2 + (1 / 5)^2))
  expect_true(is.na(eff[5]))
  # The target row counts when it is not excluded.
  expect_equal(.program_effective_traits(fit, 1)[1], 1 / ((25 / 29)^2 + (4 / 29)^2))
  expect_true(all(is.na(.program_effective_traits(NULL, 1:2))))
})

test_that(".stability_best_overlap scores Jaccard, not one-sided containment", {
  ref <- paste0("s", 1:3)
  big <- paste0("s", 1:50)
  expect_equal(.stability_best_overlap(ref, list(big)), 3 / 50)
  expect_equal(.stability_best_overlap(ref, list(big, ref)), 1)
  expect_equal(.stability_best_overlap(ref, list()), 0)
  expect_equal(.stability_best_overlap(character(0), list(big)), 0)
})

test_that("resolve_ebmf_settings merges overrides over the shared file", {
  defaults <- read_ebmf_settings()
  expect_true(all(
    c("ebmf_prior", "ebmf_magnitude_threshold", "min_module_size",
      "similarity_q", "high_similarity_q", "high_connectedness",
      "min_effective_traits",
      "include_trans") %in% names(defaults)
  ))
  expect_false(any(c("coherence_q", "coherence_n_perm", "strength_q",
                     "min_excess_similarity", "min_excess_connectedness") %in%
                     names(defaults)))

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
