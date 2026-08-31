#' @title Simulate a Trait Object With Planted SNP Modules
#' @description Generate a synthetic trait object (same shape as
#' `trait(id, include_associations = TRUE)`) with an optional planted latent
#' module structure over SNPs, for testing the univariate clustering pipeline.
#' The returned object feeds directly into `run_univariate_clustering()`, and
#' the accompanying ground truth supports recovery scoring via
#' `evaluate_univariate_simulation()`.
#'
#' Structure: `n_coloc_groups` SNPs (one coloc group / locus each) define the
#' total number of SNPs (trait x SNP matrix columns). `K` contiguous module
#' regions are laid out left-to-right; adjacent regions share `overlap` SNPs
#' (shared SNPs have driver traits from both modules). SNPs not covered by any
#' module are unstructured background: sparse random pleiotropic profiles with
#' no shared latent structure.
#'
#' Effects: each module has `drivers_per_module` driver traits with
#' `effect_size`-scale effects on their module's SNPs. Driver support is thinned
#' by `p_structural_zero` (structural zeros inside true support). Cells outside
#' true support receive small noise effects with probability `p_spurious`.
#' Background traits activate on each SNP independently with probability
#' `p_active_background`, giving realistic marginal sparsity without structure.
#'
#' Trait composition: the total number of traits (rows) is
#' `1` (target) + `sum(drivers_per_module)` (drivers) + `n_background_traits`,
#' or is fixed directly via `n_traits`. A fraction `p_molecular` of the
#' non-target traits is designated molecular: these carry a gene annotation
#' (`gene_id` non-missing) and are observed much more sparsely, with per-SNP
#' probability `p_active_molecular` instead of `p_active_background`
#' (background) or `1 - p_structural_zero` (within-module drivers). This
#' emulates eQTL-like traits whose signal is confined to a small subset of
#' SNPs. The remaining traits are phenotype traits (no gene annotation) at
#' normal density.
#'
#' Annotations: `module_annotations` optionally plants genes / tissues /
#' trait categories on each module's driver traits so enrichment machinery can
#' be verified; unplanted components draw randomly from built-in pools, and
#' `annotation_noise` controls contamination of planted labels.
#' @param n_coloc_groups Integer number of SNPs (= coloc groups) for the target trait.
#' @param K Integer number of planted modules. If `NULL`, `K` is inferred from
#'   the length of `module_sizes`, or set to `0` when `module_sizes` is `NULL`.
#'   `K = 0` gives the pure null (all SNPs unstructured).
#' @param module_sizes Optional integer vector of length `K`: SNPs per module
#'   region (shared SNPs counted once per adjacent pair). Defaults to an equal
#'   split of the non-background SNPs.
#' @param n_background_snps Minimum number of trailing SNPs reserved as
#'   unstructured background.
#' @param overlap Either `"disjoint"` (no shared SNPs), `"moderate"`
#'   (~50 shared SNPs between adjacent modules), `"strong"` (~half of the
#'   smaller adjacent module shared), or an explicit non-negative integer
#'   applied to every adjacent pair. Overlap is capped at half the smaller
#'   module minus one.
#' @param driver_sharing Fraction (0–1) of each module's drivers whose support
#'   extends into the previous module, so adjacent modules share some driver
#'   traits. Zero (default) gives every module an exclusive driver set, which
#'   makes same-signature modules statistically exchangeable in trait space
#'   and lets factor models merge them without loss. Positive sharing creates
#'   chain-like overlap structure that penalises merging.
#' @param drivers_per_module Number of driver traits per module. Either a single
#'   positive integer applied to every module, or an integer vector of length
#'   `K` giving one value per module (so modules can carry different numbers of
#'   drivers). Total drivers = `sum(drivers_per_module)`.
#' @param n_background_traits Number of unstructured background traits.
#' @param n_traits Optional total number of traits (rows) in the trait x SNP
#'   matrix, including the target trait. When provided,
#'   `n_background_traits` is derived as
#'   `n_traits - 1 - sum(drivers_per_module)` and must be non-negative. When
#'   `NULL` (default), `n_background_traits` is used directly. Total traits =
#'   1 + sum(drivers_per_module) + n_background_traits.
#' @param effect_size Mean absolute z-score of driver effects inside their
#'   module. Either a single number applied to every module, or a numeric
#'   vector of length `K` giving one value per module. Varying effect sizes
#'   across modules breaks the exchangeability of same-signature modules so
#'   that methods which operate in trait space (e.g. EBMF) can separate them.
#' @param noise_sd Standard deviation of noise effects added to target,
#'   driver, and background trait observations.
#' @param sign_pattern Driver-sign regime: `"coherent"` (all drivers positive),
#'   `"flipped"` (signs alternate across drivers within a module, so
#'   anti-correlated profiles share a module), `"random"` (independent signs).
#' @param p_structural_zero Probability that a cell inside a module's true
#'   support is absent entirely (structural zero). Either a single number or a
#'   numeric vector of length `K`.
#' @param p_spurious Probability that a cell outside true support carries a
#'   small noise effect rather than being absent.
#' @param p_active_background Per-SNP activation probability for background
#'   trait profiles.
#' @param p_molecular Fraction (0–1) of background traits designated
#'   molecular, and — when `p_molecular_drivers = NULL` (default) — also the
#'   fraction of each module's drivers designated molecular. Molecular traits
#'   carry a gene annotation (`gene_id` non-missing) and are observed much more
#'   sparsely: drivers fill within-module support with probability
#'   `p_active_molecular` instead of `1 - p_structural_zero`, and background
#'   traits activate per SNP with probability `p_active_molecular` instead of
#'   `p_active_background`. Traits whose module annotation explicitly plants
#'   `genes` are always molecular. Default 0 = every non-target trait is a
#'   phenotype trait (no gene annotation, normal density).
#' @param p_molecular_drivers Optional fraction (0–1) of each module's driver
#'   traits that are molecular, decoupled from `p_molecular` (which then
#'   governs background traits only). Either a single number applied to every
#'   module, or a numeric vector of length `K`. Useful when `p_molecular` is
#'   high but the modules must keep a dense (phenotype) driver core — if every
#'   driver in a module is a sparse molecular trait the module becomes
#'   unrecoverable. `NULL` (default) makes drivers follow `p_molecular`.
#' @param p_active_molecular Per-SNP observation probability for molecular
#'   traits. Set much lower than `p_active_background` to emulate eQTL-like
#'   sparsity, where a molecular trait is observed at only a small subset of
#'   SNPs. Default 0.01.
#' @param snps_per_trait Optional range `c(min, max)`; trait rows exceeding the
#'   maximum have excess observations randomly removed (approximate caps;
#'   minima are not enforced).
#' @param traits_per_snp Optional range `c(min, max)`; SNP columns exceeding the
#'   maximum have excess observations randomly removed.
#' @param module_annotations Optional list of length `K` (shorter lists are
#'   padded with `NULL`). Each element is a list with any of `genes`,
#'   `tissues`, `trait_categories`; these are planted on the module's driver
#'   trait rows (`NULL` elements draw randomly). Planted genes make the driver
#'   traits molecular (non-missing `gene_id`).
#' @param annotation_noise Probability that any single planted annotation is
#'   replaced by a random draw from the corresponding pool.
#' @param log_se_sd Standard deviation (on the log scale) of per-cell standard
#'   errors, drawn as `se = exp(N(0, log_se_sd))`. Zero (default) gives
#'   `se = 1` everywhere, so beta equals the z-score. Positive values emulate
#'   heterogeneous GWAS sample sizes: observed betas keep their raw scale,
#'   z-scores become noisy rescalings of them, and methods given access to the
#'   true SEs can down-weight imprecise cells. This is what makes it possible
#'   to test whether SE-aware clustering (EBMF `ebmf_se_mode = "matrix"`)
#'   recovers structure better than unit-z clustering.
#' @param snp_driver_groups Number of sub-groups (`>= 1`) into which each
#'   module's SNPs are partitioned. Group g responds only to driver-subset g
#'   of the module, so with `snp_driver_groups > 1` SNPs within a module no
#'   longer share one uniform driver pattern: each module becomes several
#'   distinct sub-patterns instead of a single rank-1 axis. This breaks the
#'   exchangeability that lets factor models merge otherwise identical modules
#'   at no cost, and is closer to real biology where variants in a pathway tap
#'   different subsets of it.
#' @param seed Optional RNG seed; a random one is drawn and recorded when `NULL`.
#' @return A list with:
#'   \itemize{
#'     \item trait_object: simulated trait result usable by
#'       `run_univariate_clustering()`
#'     \item ground_truth: list with `module_of_snp` (named vector, `0` =
#'       background; shared SNPs take their earlier module), `multi_module_snps`
#'       (shared-SNP names), `module_memberships` (SNP index list, shared SNPs
#'       appear in both), `driver_traits` (per-module trait names),
#'       `molecular_traits` (trait names designated molecular), `planted_annotations`,
#'       `seed`, and `parameters`
#'   }
#' @export
simulate_trait <- function(n_coloc_groups = 100,
                           K = NULL,
                           module_sizes = NULL,
                           n_background_snps = 0,
                           overlap = "disjoint",
                           drivers_per_module = NULL,
                           driver_sharing = 0,
                           n_background_traits = 40L,
                           n_traits = NULL,
                           effect_size = 3,
                           noise_sd = 0.5,
                           sign_pattern = c("coherent", "flipped", "random"),
                           p_structural_zero = 0.4,
                           p_spurious = 0.05,
                           p_active_background = 0.08,
                           p_molecular = 0,
                           p_molecular_drivers = NULL,
                           p_active_molecular = 0.01,
                           snps_per_trait = NULL,
                           traits_per_snp = NULL,
                           module_annotations = NULL,
                           annotation_noise = 0.1,
                           log_se_sd = 0,
                           snp_driver_groups = 1,
                           seed = NULL) {
  sign_pattern <- match.arg(sign_pattern)

  if (!is.numeric(n_coloc_groups) || length(n_coloc_groups) != 1 ||
        n_coloc_groups < 2 || n_coloc_groups != floor(n_coloc_groups)) {
    stop("n_coloc_groups must be an integer >= 2")
  }
  if (is.null(K)) {
    K <- if (is.null(module_sizes)) 0L else length(module_sizes)
  }
  if (!is.numeric(K) || length(K) != 1 || !is.finite(K) || K < 0 ||
      K != floor(K)) {
    stop("K must be a non-negative integer")
  }
  n_snps <- as.integer(n_coloc_groups)
  K <- as.integer(K)

  if (K == 0) {
    drivers_per_module <- integer(0)
  } else if (is.null(drivers_per_module)) {
    drivers_per_module <- rep(10L, K)
  } else {
    if (length(drivers_per_module) == 1L) {
      drivers_per_module <- rep(drivers_per_module, K)
    }
    if (length(drivers_per_module) != K) {
      stop("drivers_per_module must have length 1 or ", K)
    }
    if (any(drivers_per_module < 1) ||
        any(drivers_per_module != floor(drivers_per_module))) {
      stop("drivers_per_module must be integers >= 1 (single value or one per module)")
    }
    drivers_per_module <- as.integer(drivers_per_module)
  }
  n_drivers <- if (K == 0) 0L else sum(drivers_per_module)
  if (!is.null(n_traits)) {
    if (!is.numeric(n_traits) || length(n_traits) != 1L || !is.finite(n_traits) ||
        n_traits < 1L + n_drivers || n_traits != floor(n_traits)) {
      stop("n_traits must be a single integer >= 1 + sum(drivers_per_module) (",
           1L + n_drivers, ")")
    }
    n_background_traits <- as.integer(n_traits) - 1L - n_drivers
  }
  if (!is.numeric(n_background_traits) || length(n_background_traits) != 1L ||
      !is.finite(n_background_traits) || n_background_traits < 0 ||
      n_background_traits != floor(n_background_traits)) {
    stop("n_background_traits must be a non-negative integer")
  }
  n_background_traits <- as.integer(n_background_traits)
  if (!is.numeric(n_background_snps) || length(n_background_snps) != 1L ||
      !is.finite(n_background_snps) || n_background_snps < 0 ||
      n_background_snps != floor(n_background_snps)) {
    stop("n_background_snps must be a non-negative integer")
  }
  n_background_snps <- as.integer(n_background_snps)
  if (!is.numeric(driver_sharing) || length(driver_sharing) != 1L ||
      !is.finite(driver_sharing) || driver_sharing < 0 || driver_sharing > 1) {
    stop("driver_sharing must be a number between 0 and 1")
  }
  if (!is.numeric(snp_driver_groups) || length(snp_driver_groups) != 1L ||
      !is.finite(snp_driver_groups) || snp_driver_groups < 1 ||
      snp_driver_groups != floor(snp_driver_groups) ||
      any(snp_driver_groups > drivers_per_module)) {
    stop("snp_driver_groups must be an integer between 1 and drivers_per_module")
  }
  snp_driver_groups <- as.integer(snp_driver_groups)
  if (!is.numeric(noise_sd) || length(noise_sd) != 1L ||
      !is.finite(noise_sd) || noise_sd < 0) {
    stop("noise_sd must be a non-negative finite number")
  }
  if (!is.numeric(log_se_sd) || length(log_se_sd) != 1L ||
      !is.finite(log_se_sd) || log_se_sd < 0) {
    stop("log_se_sd must be a non-negative finite number")
  }
  probabilities <- c(
    p_spurious, p_active_background, annotation_noise,
    p_molecular, p_active_molecular
  )
  if (any(!is.finite(probabilities)) || any(probabilities < 0) ||
      any(probabilities > 1)) {
    stop("p_spurious, p_active_background, p_molecular, p_active_molecular, and annotation_noise must be between 0 and 1")
  }
  if (!is.null(p_molecular_drivers)) {
    if (!is.numeric(p_molecular_drivers) || anyNA(p_molecular_drivers) ||
        any(!is.finite(p_molecular_drivers)) ||
        any(p_molecular_drivers < 0) || any(p_molecular_drivers > 1)) {
      stop("p_molecular_drivers must be NULL or a number between 0 and 1")
    }
  }

  effect_size <- .expand_per_module_parameter(
    effect_size, K, "effect_size", allow_empty = TRUE
  )
  p_structural_zero <- .expand_per_module_parameter(
    p_structural_zero, K, "p_structural_zero", allow_empty = TRUE
  )
  if (any(p_structural_zero > 1)) {
    stop("p_structural_zero must be between 0 and 1")
  }
  p_molecular_drivers <- if (is.null(p_molecular_drivers)) {
    NULL
  } else {
    .expand_per_module_parameter(
      p_molecular_drivers, K, "p_molecular_drivers", allow_empty = TRUE
    )
  }

  overlap_n <- .resolve_overlap(overlap, K)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("seed must be a single number")
    }
    seed <- as.integer(seed)
  } else {
    seed <- sample.int(.Machine$integer.max, 1L)
  }
  set.seed(seed)

  geometry <- .plant_module_regions(
    n_snps = n_snps,
    K = K,
    module_sizes = module_sizes,
    n_background_snps = n_background_snps,
    overlap_n = overlap_n
  )

  module_of <- rep(0L, n_snps)
  for (m in seq_len(K)) {
    unfilled <- module_of[geometry$regions[[m]]] == 0L
    module_of[geometry$regions[[m]][unfilled]] <- m
  }

  n_bg_traits <- n_background_traits
  trait_ids <- 1:(1 + n_drivers + n_bg_traits)
  target_tid <- 1L
  driver_tids <- if (n_drivers > 0) 2:(1 + n_drivers) else integer(0)
  bg_tids <- if (n_bg_traits > 0) seq_len(n_bg_traits) + 1L + n_drivers else integer(0)
  driver_module <- if (n_drivers > 0) {
    rep(seq_len(K), times = drivers_per_module)
  } else {
    integer(0)
  }
  is_molecular <- .assign_molecular_traits(
    n_drivers = n_drivers,
    n_bg_traits = n_bg_traits,
    p_molecular = p_molecular,
    p_molecular_drivers = p_molecular_drivers,
    driver_module = driver_module,
    drivers_per_module = drivers_per_module,
    module_annotations = module_annotations
  )

  M <- matrix(NA_real_, nrow = length(trait_ids), ncol = n_snps,
              dimnames = list(as.character(trait_ids), .sim_variant_ids(n_snps)))
  snp_module <- module_of

  target_signal <- ifelse(
    snp_module > 0,
    effect_size[pmax(snp_module, 1)] * ifelse(snp_module %% 2 == 1, 1, -1),
    0
  )
  M[as.character(target_tid), ] <- target_signal +
    stats::rnorm(n_snps, 0, noise_sd)

  driver_signs <- switch(sign_pattern,
    coherent = rep(1, n_drivers),
    flipped = unlist(lapply(seq_len(K), function(m) {
      ifelse(seq_len(drivers_per_module[m]) %% 2 == 1, 1, -1)
    })),
    random = sample(c(-1, 1), n_drivers, replace = TRUE)
  )
  if (n_drivers > 0) {
    n_shared_per_module <- floor(drivers_per_module * driver_sharing)
    n_groups <- max(1L, as.integer(snp_driver_groups))
    snp_groups <- vector("list", K)
    if (n_groups > 1) {
      for (m in seq_len(K)) {
        idx <- geometry$memberships[[m]]
        g <- cut(seq_along(idx), breaks = n_groups, labels = FALSE)
        snp_groups[[m]] <- stats::setNames(g, colnames(M)[idx])
      }
    }
    drivers_per_group <- drivers_per_module / n_groups
    for (d in seq_len(n_drivers)) {
      m <- driver_module[d]
      row <- as.character(driver_tids[d])
      in_module <- geometry$memberships[[m]]
      slot <- (d - 1) %% drivers_per_module[m]
      if (m > 1 && slot < n_shared_per_module[m]) {
        in_module <- union(in_module, geometry$memberships[[m - 1]])
      }
      if (n_groups > 1) {
        driver_group <- min(n_groups, floor(slot / drivers_per_group[m]) + 1)
        grp <- snp_groups[[m]][as.character(colnames(M)[in_module])]
        keep <- !is.na(grp) & grp == driver_group
        in_module <- in_module[keep]
        if (length(in_module) == 0) {
          in_module <- geometry$memberships[[m]][1]
        }
      }
      M[row, in_module] <- ifelse(
        stats::runif(length(in_module)) < if (is_molecular[d]) {
          p_active_molecular
        } else {
          1 - p_structural_zero[m]
        },
        effect_size[m] * driver_signs[d] +
          stats::rnorm(length(in_module), 0, noise_sd),
        NA_real_
      )
      off_module <- setdiff(seq_len(n_snps), in_module)
      spur <- stats::runif(length(off_module)) < p_spurious
      vals <- rep(NA_real_, length(off_module))
      vals[spur] <- stats::rnorm(sum(spur), 0, noise_sd)
      M[row, off_module] <- vals
    }
  }
  if (n_bg_traits > 0) {
    for (t in seq_along(bg_tids)) {
      row <- as.character(bg_tids[t])
      act_prob <- if (is_molecular[n_drivers + t]) {
        p_active_molecular
      } else {
        p_active_background
      }
      act <- stats::runif(n_snps) < act_prob
      vals <- rep(NA_real_, n_snps)
      vals[act] <- stats::rnorm(sum(act), 0, noise_sd)
      M[row, ] <- vals
    }
  }

  M <- .cap_matrix_density(M, snps_per_trait, traits_per_snp,
                           protected_row = as.character(target_tid))

  annotations <- .sim_trait_annotations(
    trait_ids = trait_ids,
    target_tid = target_tid,
    driver_tids = driver_tids,
    driver_module = driver_module,
    bg_tids = bg_tids,
    module_annotations = module_annotations,
    annotation_noise = annotation_noise,
    is_molecular = is_molecular,
    K = K
  )

  se_matrix <- matrix(
    exp(stats::rnorm(length(M), 0, log_se_sd)),
    nrow = nrow(M),
    ncol = ncol(M),
    dimnames = dimnames(M)
  )

  cg <- .sim_coloc_groups_from_matrix(M, annotations, se_matrix)
  trait_object <- list(
    trait = list(
      id = target_tid,
      trait = "SIMULATED",
      trait_name = "Simulated target trait",
      source_url = NA_character_
    ),
    coloc_groups = cg
  )

  memberships_named <- lapply(geometry$memberships, function(idx) colnames(M)[idx])
  multi_snps <- names(which(table(unlist(memberships_named)) > 1))

  driver_names <- lapply(seq_len(K), function(m) {
    return(annotations$trait_name[driver_tids[driver_module == m]])
  })

  molecular_trait_ids <- c(
    driver_tids[is_molecular[seq_len(n_drivers)]],
    bg_tids[is_molecular[n_drivers + seq_len(n_bg_traits)]]
  )
  molecular_traits <- annotations$trait_name[
    annotations$trait_id %in% molecular_trait_ids
  ]

  return(list(
    trait_object = trait_object,
    ground_truth = list(
      module_of_snp = stats::setNames(module_of, colnames(M)),
      multi_module_snps = multi_snps,
      module_memberships = memberships_named,
      driver_traits = stats::setNames(driver_names, sprintf("module_%d", seq_len(K))),
      molecular_traits = molecular_traits,
      planted_annotations = module_annotations,
      seed = seed,
      parameters = list(
        n_coloc_groups = n_snps,
        K = K,
        module_sizes = geometry$sizes,
        overlap_n = geometry$overlaps,
        n_background_snps_reserved = n_background_snps,
        n_background_snps_actual = sum(module_of == 0L),
        drivers_per_module = drivers_per_module,
        driver_sharing = driver_sharing,
        n_background_traits = n_bg_traits,
        n_traits = as.integer(1L + n_drivers + n_bg_traits),
        effect_size = effect_size,
        noise_sd = noise_sd,
        sign_pattern = sign_pattern,
        p_structural_zero = p_structural_zero,
        p_spurious = p_spurious,
        p_active_background = p_active_background,
        p_molecular = p_molecular,
        p_molecular_drivers = p_molecular_drivers,
        p_active_molecular = p_active_molecular,
        annotation_noise = annotation_noise,
        log_se_sd = log_se_sd,
        snp_driver_groups = snp_driver_groups
      )
    )
  ))
}


.expand_per_module_parameter <- function(x, K, name, allow_empty = FALSE) {
  if (!is.numeric(x) || length(x) == 0 || anyNA(x) || any(!is.finite(x))) {
    stop(name, " must be a non-NA numeric vector")
  }
  if (length(x) == 1L) {
    x <- rep(x, max(K, 1L))
  } else if (length(x) != K) {
    stop(name, " must have length 1 or ", K)
  }
  if (any(x < 0)) {
    stop(name, " must be non-negative")
  }
  return(x)
}


.resolve_overlap <- function(overlap, K) {
  if (K <= 1) {
    return(rep(0L, max(K - 1, 0)))
  }
  if (length(overlap) == 1 && is.character(overlap)) {
    o <- switch(match.arg(overlap, c("disjoint", "moderate", "strong")),
      disjoint = 0L,
      moderate = 50L,
      strong = NA_integer_
    )
    return(rep(o, K - 1L))
  }
  if (length(overlap) == 1 && is.numeric(overlap)) {
    if (overlap < 0 || overlap != floor(overlap)) {
      stop("overlap must be a named level or a non-negative integer")
    }
    return(rep(as.integer(overlap), K - 1L))
  }
  stop("overlap must be 'disjoint', 'moderate', 'strong', or a single integer")
}


.assign_molecular_traits <- function(n_drivers, n_bg_traits, p_molecular,
                                     p_molecular_drivers, driver_module,
                                     drivers_per_module, module_annotations) {
  n <- n_drivers + n_bg_traits
  is_molecular <- rep(FALSE, n)
  if (n == 0) {
    return(is_molecular)
  }
  # Planted module genes always make their drivers molecular.
  if (n_drivers > 0 && !is.null(module_annotations)) {
    for (d in seq_len(n_drivers)) {
      m <- driver_module[d]
      spec <- if (m <= length(module_annotations)) module_annotations[[m]] else NULL
      if (!is.null(spec) && !is.null(spec$genes) && length(spec$genes) > 0) {
        is_molecular[d] <- TRUE
      }
    }
  }
  # Driver molecular assignment, per module. When p_molecular_drivers is NULL
  # the drivers follow p_molecular (one fraction for every module); otherwise
  # each module uses its own p_molecular_drivers[m], so a module can keep a
  # dense phenotype driver core even when p_molecular (background) is high.
  if (n_drivers > 0) {
    driver_frac <- if (is.null(p_molecular_drivers)) {
      rep(p_molecular, length(drivers_per_module))
    } else {
      p_molecular_drivers
    }
    for (m in seq_along(drivers_per_module)) {
      idx <- which(driver_module == m)
      if (length(idx) == 0) {
        next
      }
      desired <- round(driver_frac[m] * length(idx))
      forced <- sum(is_molecular[idx])
      if (forced < desired) {
        candidates <- idx[!is_molecular[idx]]
        pick <- sample(candidates, min(desired - forced, length(candidates)))
        is_molecular[pick] <- TRUE
      }
    }
  }
  # Background molecular assignment.
  if (n_bg_traits > 0 && p_molecular > 0) {
    bg_idx <- n_drivers + seq_len(n_bg_traits)
    desired <- round(p_molecular * n_bg_traits)
    candidates <- bg_idx[!is_molecular[bg_idx]]
    pick <- sample(candidates, min(desired, length(candidates)))
    is_molecular[pick] <- TRUE
  }
  return(is_molecular)
}


.plant_module_regions <- function(n_snps, K, module_sizes, n_background_snps, overlap_n) {
  if (K == 0) {
    return(list(regions = list(), memberships = list(), sizes = integer(0), overlaps = integer(0)))
  }
  capacity <- n_snps - max(0L, as.integer(n_background_snps))
  if (capacity < K) {
    stop("not enough SNPs for ", K, " modules given n_background_snps")
  }
  if (is.null(module_sizes)) {
    base <- capacity %/% K
    rem <- capacity %% K
    sizes <- rep(base, K)
    sizes[seq_len(rem)] <- sizes[seq_len(rem)] + 1L
  } else {
    if (!is.numeric(module_sizes) || length(module_sizes) != K ||
          any(module_sizes < 2) || any(module_sizes != floor(module_sizes))) {
      stop("module_sizes must be a vector of ", K, " integers >= 2")
    }
    sizes <- as.integer(module_sizes)
  }

  overlaps <- integer(max(K - 1L, 0))
  for (i in seq_along(overlaps)) {
    cap <- floor(min(sizes[i], sizes[i + 1]) / 2) - 1L
    if (is.na(overlap_n[i])) {
      overlaps[i] <- max(cap, 0L)
    } else {
      overlaps[i] <- as.integer(min(overlap_n[i], max(cap, 0L)))
    }
  }
  consumed <- sum(sizes) - sum(overlaps)
  if (consumed > n_snps) {
    stop("modules need ", consumed, " SNPs but only ", n_snps, " available")
  }

  regions <- vector("list", K)
  memberships <- vector("list", K)
  start <- 1L
  for (m in seq_len(K)) {
    end <- start + sizes[m] - 1L
    regions[[m]] <- start:end
    memberships[[m]] <- start:end
    if (m < K) {
      o <- overlaps[m]
      if (o > 0) {
        shared <- (end - o + 1L):end
        memberships[[m + 1]] <- c(shared, regions[[m + 1]])
      }
      start <- end - o + 1L
    }
  }

  return(list(
    regions = regions,
    memberships = memberships,
    sizes = sizes,
    overlaps = overlaps
  ))
}


.sim_variant_ids <- function(n) {
  return(sprintf("1:%d:A:G", seq_len(n) * 1000L))
}


.cap_matrix_density <- function(M, snps_per_trait, traits_per_snp, protected_row = NULL) {
  if (!is.null(snps_per_trait)) {
    if (length(snps_per_trait) != 2 || any(snps_per_trait < 0)) {
      stop("snps_per_trait must be c(min, max)")
    }
    counts <- rowSums(!is.na(M))
    over <- counts > snps_per_trait[2]
    if (!is.null(protected_row)) {
      over[protected_row] <- FALSE
    }
    for (r in rownames(M)[over]) {
      drop <- sample(which(!is.na(M[r, ])), counts[r] - snps_per_trait[2])
      M[r, drop] <- NA_real_
    }
  }
  if (!is.null(traits_per_snp)) {
    if (length(traits_per_snp) != 2 || any(traits_per_snp < 0) ||
        any(!is.finite(traits_per_snp)) ||
        traits_per_snp[2] != floor(traits_per_snp[2])) {
      stop("traits_per_snp must be c(min, max)")
    }
    counts <- colSums(!is.na(M))
    over <- counts > traits_per_snp[2]
    for (c in colnames(M)[over]) {
      max_count <- as.integer(traits_per_snp[2])
      protected_idx <- if (!is.null(protected_row)) {
        match(protected_row, rownames(M))
      } else {
        NA_integer_
      }
      if (!is.na(protected_idx)) {
        max_count <- max(max_count, 1L)
      }
      drop_n <- counts[c] - max_count
      eligible <- which(!is.na(M[, c]))
      if (!is.na(protected_idx)) {
        eligible <- setdiff(eligible, protected_idx)
      }
      if (length(eligible) < drop_n) {
        stop("traits_per_snp cap would remove the protected target-trait observation")
      }
      drop <- sample(eligible, drop_n)
      M[drop, c] <- NA_real_
    }
  }
  return(M)
}


.sim_annotation_pools <- function() {
  return(list(
    tissues = c(
      "Liver", "Brain", "Blood", "Heart", "Muscle", "Lung",
      "Kidney", "Pancreas", "Spleen", "Adipose"
    ),
    trait_categories = c(
      "Blood", "Metabolic", "Cardiovascular", "Immune",
      "Digestive", "Neurological", "Renal", "Respiratory"
    ),
    genes = c(
      "TMPRSS6", "HFE", "TFR2", "HAMP", "SLC40A1", "APOE", "PCSK9", "LDLR",
      "LPL", "APOA5", "FTO", "MC4R", "LEPR", "PPARG", "TCF7L2", "KCNJ11",
      "ABCG8", "HMGCR", "SORT1", "CETP"
    )
  ))
}


.sim_trait_annotations <- function(trait_ids,
                                   target_tid,
                                   driver_tids,
                                   driver_module,
                                   bg_tids,
                                   module_annotations,
                                   annotation_noise,
                                   is_molecular,
                                   K) {
  pools <- .sim_annotation_pools()
  specs <- vector("list", K)
  if (!is.null(module_annotations)) {
    for (m in seq_len(K)) {
      specs[m] <- list(
        if (m <= length(module_annotations)) module_annotations[[m]] else NULL
      )
    }
  }

  draw <- function(pool, prob = 1) {
    if (stats::runif(1) > prob) {
      return(NA_character_)
    }
    return(sample(pool, 1))
  }
  maybe_contaminate <- function(value, pool, planted) {
    if (is.null(value)) {
      return(draw(pool, 0.5))
    }
    if (stats::runif(1) < annotation_noise) {
      return(sample(setdiff(pool, planted), 1))
    }
    return(value)
  }

  out <- data.frame(
    trait_id = as.integer(trait_ids),
    trait_name = NA_character_,
    tissue = NA_character_,
    gene = NA_character_,
    gene_id = NA_integer_,
    trait_category = NA_character_,
    stringsAsFactors = FALSE
  )
  out$trait_name[out$trait_id == target_tid] <- "Simulated target trait"

  gene_registry <- new.env(parent = emptyenv())
  next_gene_id <- local({
    counter <- 0L
    function(sym) {
      if (!exists(sym, envir = gene_registry, inherits = FALSE)) {
        counter <<- counter + 1L
        assign(sym, counter, envir = gene_registry, inherits = FALSE)
      }
      return(get(sym, envir = gene_registry, inherits = FALSE))
    }
  })
  register_gene <- function(sym) {
    if (is.na(sym)) {
      return(NA_integer_)
    }
    return(as.integer(next_gene_id(sym)))
  }

  if (length(driver_tids) > 0) {
    for (d in seq_along(driver_tids)) {
      m <- driver_module[d]
      spec <- specs[[m]]
      tid_row <- out$trait_id == driver_tids[d]
      out$trait_name[tid_row] <- sprintf("driver_mod%d_%03d", m, d)
      if (is_molecular[d]) {
        if (!is.null(spec) && !is.null(spec$genes) && length(spec$genes) > 0) {
          sym <- spec$genes[(d - 1) %% length(spec$genes) + 1]
        } else {
          sym <- draw(pools$genes, 1)
        }
        out$gene[tid_row] <- sym
        out$gene_id[tid_row] <- register_gene(sym)
      }
      if (!is.null(spec) && !is.null(spec$tissues) && length(spec$tissues) > 0) {
        out$tissue[tid_row] <- maybe_contaminate(
          spec$tissues[(d - 1) %% length(spec$tissues) + 1],
          pools$tissues,
          spec$tissues
        )
      } else {
        out$tissue[tid_row] <- draw(pools$tissues, 0.5)
      }
      if (!is.null(spec) && !is.null(spec$trait_categories) && length(spec$trait_categories) > 0) {
        out$trait_category[tid_row] <- maybe_contaminate(
          spec$trait_categories[(d - 1) %% length(spec$trait_categories) + 1],
          pools$trait_categories,
          spec$trait_categories
        )
      } else {
        out$trait_category[tid_row] <- draw(pools$trait_categories, 0.5)
      }
    }
  }

  if (length(bg_tids) > 0) {
    for (b in seq_along(bg_tids)) {
      tid_row <- out$trait_id == bg_tids[b]
      out$trait_name[tid_row] <- sprintf("background_trait_%03d", b)
      if (is_molecular[length(driver_tids) + b]) {
        sym <- draw(pools$genes, 1)
        out$gene[tid_row] <- sym
        out$gene_id[tid_row] <- register_gene(sym)
      }
      out$tissue[tid_row] <- draw(pools$tissues, 0.5)
      out$trait_category[tid_row] <- draw(pools$trait_categories, 0.5)
    }
  }

  return(out)
}


.sim_coloc_groups_from_matrix <- function(M, annotations, se_matrix = NULL) {
  observed <- which(!is.na(M), arr.ind = TRUE)
  if (nrow(observed) == 0) {
    stop("simulated matrix has no observations")
  }
  trait_idx <- as.integer(rownames(M))[observed[, "row"]]
  snp_idx <- observed[, "col"]
  beta <- M[observed]

  if (is.null(se_matrix)) {
    se <- rep(1, length(beta))
  } else {
    se <- se_matrix[observed]
  }
  z <- beta / se
  p <- pmax(2 * stats::pnorm(-abs(z)), 1e-300)

  anno <- annotations[match(trait_idx, annotations$trait_id), , drop = FALSE]

  cg <- data.frame(
    coloc_group_id = as.integer(snp_idx),
    variant_id = colnames(M)[snp_idx],
    display_snp = colnames(M)[snp_idx],
    chr = 1L,
    bp = as.integer(snp_idx) * 1000L,
    trait_id = anno$trait_id,
    trait_name = anno$trait_name,
    min_p = p,
    beta = beta,
    se = se,
    tissue = anno$tissue,
    gene = anno$gene,
    gene_id = anno$gene_id,
    trait_category = anno$trait_category,
    stringsAsFactors = FALSE
  )
  cg <- cg[order(cg$coloc_group_id, cg$trait_id), , drop = FALSE]
  rownames(cg) <- NULL
  return(cg)
}
