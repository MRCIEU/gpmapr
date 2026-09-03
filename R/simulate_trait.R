#' @title Simulate a Trait Object With Planted SNP Modules
#' @description Generate a synthetic trait object (same shape as
#' `trait(id, include_associations = TRUE)`) with an optional planted latent
#' module structure over SNPs, for testing the univariate clustering pipeline.
#' The returned object feeds directly into `run_univariate_clustering()`, and
#' the accompanying ground truth supports recovery scoring via
#' `evaluate_univariate_simulation()`.
#'
#' By default the generative parameters are calibrated to the messy marginals
#' of a real pleiotropy matrix: every observed cell is a significant hit
#' (`min_abs_z = 4.5`), driver effects are heavy-tailed across traits
#' (`effect_tail = 0.4`) with a minority of cells disagreeing with the target
#' orientation (`p_negative = 0.05`), background densities are heavy-tailed
#' per trait (`background_sparsity_sd = 1.2`), a proportion of pleiotropic
#' **hub** traits is planted automatically (`n_hub_traits = NULL`), and the
#' target row is dense and significant at every SNP (`target_pattern =
#' "dense"`). Simulations built this way are genuinely hard — recovery is
#' partial and noisy rather than trivially perfect, so the validation gates
#' and parameter trade-offs can actually be tested. Each realism knob can be
#' turned off individually (e.g. `min_abs_z = 0`, `n_hub_traits = 0L`).
#'
#' Structure: `n_coloc_groups` SNPs (one coloc group / locus each) define the
#' total number of SNPs (trait x SNP matrix columns). `K` module regions are
#' laid out left-to-right as **disjoint** (non-overlapping) SNP sets by default,
#' or, with `overlap = "nested"`, as a parent module (module 1) containing every
#' other module's SNPs as correlated subsets. SNPs not covered by any module
#' are unstructured background: sparse random pleiotropic profiles with no
#' shared latent structure.
#'
#' Effects: each module has `n_traits_per_module` driver traits with
#' `effect_size`-scale effects on their module's SNPs. A fraction
#' `trait_overlap` of each module's driver traits is **shared with the other
#' modules**: a shared trait drives every module that draws it, with per-module
#' effects that are correlated across modules (`trait_effect_corr`) but not
#' identical. This creates modules with distinct SNP sets but overlapping,
#' correlated trait signatures. Driver support is thinned by
#' `p_structural_zero` (structural zeros inside true support). Cells outside
#' true support receive small noise effects with probability `p_spurious`.
#' Background traits activate on each SNP independently with probability
#' `p_active_background`, giving realistic marginal sparsity without structure.
#'
#' Trait composition: the total number of traits (rows) is
#' `1` (target) + `n_driver_rows` (unique driver traits; fewer than
#' `sum(n_traits_per_module)` when `trait_overlap > 0` because shared traits
#' are counted once) + `n_background_traits`, or is fixed directly via
#' `n_traits`. A fraction `p_molecular` of the non-target traits is designated
#' molecular: these carry a gene annotation (`gene_id` non-missing) and are
#' observed much more sparsely, with per-SNP probability `p_active_molecular`
#' instead of `p_active_background` (background) or `1 - p_structural_zero`
#' (within-module drivers). This emulates eQTL-like traits whose signal is
#' confined to a small subset of SNPs. The remaining traits are phenotype
#' traits (no gene annotation) at normal density.
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
#'   region (disjoint regions, laid out left-to-right). Defaults to an equal
#'   split of the non-background SNPs.
#' @param n_background_snps Minimum number of trailing SNPs reserved as
#'   unstructured background.
#' @param overlap How the planted modules are arranged in SNP space:
#'   `"disjoint"` (default) gives every module its own contiguous, non-
#'   overlapping SNP set; `"nested"` makes the largest module (module 1) the
#'   parent and samples every other module's SNPs as a subset of it, so the
#'   smaller modules are contained inside the parent and may overlap one
#'   another; `"partial"` keeps disjoint cores but lets each module borrow a
#'   fraction (`overlap_fraction`) of the previous module's SNPs, so adjacent
#'   modules share boundary SNPs.
#' @param nested_effect_corr Correlation between a nested module's effects and
#'   the parent module's effects at SNPs shared by both, used when
#'   `overlap = "nested"`. Effects are generated as
#'   `z_child = corr*z_parent + sqrt(1-corr^2)*N(0,1)` per shared SNP, so a
#'   nested module is a correlated (not arbitrary) subset of the parent.
#'   Default 0.8.
#' @param n_traits_per_module Number of driver traits per module. Either a single
#'   positive integer applied to every module, or an integer vector of length
#'   `K` giving one value per module. Total driver positions =
#'   `sum(n_traits_per_module)`; the number of unique driver trait rows is
#'   smaller when `trait_overlap > 0`.
#' @param trait_overlap Fraction (0–1) of each module's driver traits that are
#'   shared with the other modules. `0` (default) gives every module an
#'   exclusive driver set (independent trait signatures); `0.5` shares half the
#'   drivers; `0.75` gives strongly overlapping signatures; `0.9` almost
#'   identical signatures. Shared traits drive every module that draws them.
#' @param trait_effect_corr Correlation between the effect sizes of a shared
#'   trait across the modules it drives. Effects are generated as
#'   `z1 ~ N(0,1)` and `z2 = corr*z1 + sqrt(1-corr^2)*N(0,1)`, so shared
#'   effects are correlated but identical only when `corr = 1`. Default 0.8.
#' @param driver_sharing Fraction (0–1) of each module's drivers whose support
#'   extends into the previous module's SNP set (a SNP-space mechanism, kept
#'   separate from `trait_overlap`). Zero (default) confines every driver's
#'   support to its own module's SNPs.
#' @param n_background_traits Number of unstructured background traits.
#' @param n_traits Optional total number of traits (rows) in the trait x SNP
#'   matrix, including the target trait. When provided,
#'   `n_background_traits` is derived as
#'   `n_traits - 1 - n_driver_rows` and must be non-negative. When
#'   `NULL` (default), `n_background_traits` is used directly.
#' @param effect_size Mean absolute z-score of driver effects inside their
#'   module. Either a single number applied to every module, or a numeric
#'   vector of length `K` giving one value per module. Varying effect sizes
#'   across modules breaks the exchangeability of same-signature modules so
#'   that methods which operate in trait space (e.g. EBMF) can separate them.
#'   Defaults to `6`, matching the median observed |z| of real coloc cells and
#'   keeping driver signal above the `min_abs_z` floor.
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
#' @param background_sparsity_sd Standard deviation (on the log scale) of
#'   per-trait multipliers on the background activation probabilities. Zero
#'   gives every background trait the same density; positive values
#'   draw heavy-tailed per-trait rates so most background traits are observed
#'   at only a couple of SNPs while a tail is much denser, matching the
#'   per-trait sparsity distribution of real pleiotropy matrices.
#'   Default 1.2.
#' @param min_abs_z Significance floor for observed cells. When positive
#'   (`4.5` by default, the approximate floor of real coloc z-scores), every
#'   *observed* cell is a significant hit: background and spurious cells are
#'   drawn as `sign * (min_abs_z + Exp)` and all other observed cells are
#'   floored at `min_abs_z`, so noise is encoded as absence (NA) rather than
#'   as small observed values — the way real colocalisation data behaves.
#'   Set to `0` to keep the classic behaviour where background and spurious
#'   cells carry small `N(0, noise_sd)` values.
#' @param effect_tail Standard deviation (on the log scale) of a per-*trait*
#'   lognormal multiplier applied to a driver trait's effects across all of its
#'   module's SNPs. Zero keeps every driver trait at (approximately)
#'   the same `effect_size`; positive values give the heavy-tailed
#'   across-trait heterogeneity of real modules (different driver traits have
#'   different typical effect sizes) while preserving the correlation of a
#'   trait's values across the SNPs it drives. Default 0.4.
#' @param n_hub_traits Number of background traits planted as pleiotropic
#'   **hubs**: dense rows (observed at a fraction of SNPs drawn from
#'   `hub_snp_fraction`) whose cells are significant and mostly positive in
#'   the oriented frame (sign flips follow `p_negative`), creating the global
#'   baseline similarity and mega-factor structure that competes with the
#'   planted modules — the way real pleiotropy matrices contain a core of
#'   highly pleiotropic traits hitting most loci with target-aligned effects.
#'   Hubs are drawn from the non-molecular background so they stay gene-free.
#'   `NULL` (default) auto-plants ~3% of background traits (clamped to the
#'   number of non-molecular background traits available); set to `0L` to
#'   disable.
#' @param hub_snp_fraction Range `c(min, max)` of the fraction of SNPs at
#'   which each hub trait is observed. Defaults to `c(0.3, 1)`.
#' @param target_pattern Target-trait row regime: `"dense"` (default) gives
#'   the target trait a significant, heavy-tailed positive z-score at **every**
#'   SNP (magnitude `effect_size[1] * lognormal(0, 0.5)`, floored at
#'   `min_abs_z` when set), matching a real target-trait row in the oriented
#'   frame — informative about nothing but its own strength; `"module"` gives
#'   the target trait elevated signal at module SNPs and noise elsewhere.
#' @param p_negative Probability that an individual driver cell within module
#'   support has its sign flipped, emulating the minority of cells that
#'   disagree with the target-trait orientation in real data. Default 0.05.
#'   Set to `NULL` to keep the `sign_pattern` regime untouched.
#' @param overlap_fraction Fraction of each module's SNPs (for `overlap =
#'   "partial"`) borrowed from the previous module's core, giving modules
#'   soft SNP boundaries and non-empty `multi_module_snps` without full
#'   nesting. Defaults to `0.2`; only used when `overlap = "partial"`.
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
#'       background; each SNP is labelled by its most specific module), 
#'       `multi_module_snps` (SNP names shared by more than one module; empty
#'       for `overlap = "disjoint"`), `module_memberships` (SNP index list),
#'       `driver_traits` (per-module trait names, including shared traits),
#'       `driver_overlap` (K x K matrix of the fraction of each module's driver
#'       traits shared with every other module), `molecular_traits` (trait
#'       names designated molecular), `planted_annotations`, `seed`, and
#'       `parameters`
#'   }
#' @export
simulate_trait <- function(n_coloc_groups = 100,
                           K = NULL,
                           module_sizes = NULL,
                           n_background_snps = 0,
                           overlap = "disjoint",
                           nested_effect_corr = 0.8,
                           n_traits_per_module = NULL,
                           trait_overlap = 0,
                           trait_effect_corr = 0.8,
                           driver_sharing = 0,
                           n_background_traits = 40L,
                           n_traits = NULL,
                           effect_size = 6,
                           noise_sd = 0.5,
                           sign_pattern = c("coherent", "flipped", "random"),
                           p_structural_zero = 0.4,
                           p_spurious = 0.05,
                           p_active_background = 0.08,
                           background_sparsity_sd = 1.2,
                           min_abs_z = 4.5,
                           effect_tail = 0.4,
                           n_hub_traits = NULL,
                           hub_snp_fraction = c(0.3, 1),
                           target_pattern = c("dense", "module"),
                           p_negative = 0.05,
                           overlap_fraction = 0.2,
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
  overlap <- match.arg(overlap, c("disjoint", "nested", "partial"))
  target_pattern <- match.arg(target_pattern)

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
    n_traits_per_module <- integer(0)
  } else if (is.null(n_traits_per_module)) {
    n_traits_per_module <- rep(10L, K)
  } else {
    if (length(n_traits_per_module) == 1L) {
      n_traits_per_module <- rep(n_traits_per_module, K)
    }
    if (length(n_traits_per_module) != K) {
      stop("n_traits_per_module must have length 1 or ", K)
    }
    if (any(n_traits_per_module < 1) ||
        any(n_traits_per_module != floor(n_traits_per_module))) {
      stop("n_traits_per_module must be integers >= 1 (single value or one per module)")
    }
    n_traits_per_module <- as.integer(n_traits_per_module)
  }
  if (!is.numeric(trait_overlap) || length(trait_overlap) != 1L ||
      !is.finite(trait_overlap) || trait_overlap < 0 || trait_overlap > 1) {
    stop("trait_overlap must be a number between 0 and 1")
  }
  if (!is.numeric(trait_effect_corr) || length(trait_effect_corr) != 1L ||
      !is.finite(trait_effect_corr) ||
      trait_effect_corr < -1 || trait_effect_corr > 1) {
    stop("trait_effect_corr must be a number between -1 and 1")
  }
  n_shared <- if (K > 0) as.integer(round(n_traits_per_module * trait_overlap)) else integer(0)
  n_specific <- if (K > 0) n_traits_per_module - n_shared else integer(0)
  n_shared_pool <- if (K > 0) max(n_shared) else 0L
  n_driver_rows <- if (K > 0) n_shared_pool + sum(n_specific) else 0L

  if (!is.null(n_traits)) {
    if (!is.numeric(n_traits) || length(n_traits) != 1L || !is.finite(n_traits) ||
        n_traits < 1L + n_driver_rows || n_traits != floor(n_traits)) {
      stop("n_traits must be a single integer >= 1 + n_driver_rows (",
           1L + n_driver_rows, ")")
    }
    n_background_traits <- as.integer(n_traits) - 1L - n_driver_rows
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
  if (!is.numeric(nested_effect_corr) || length(nested_effect_corr) != 1L ||
      !is.finite(nested_effect_corr) ||
      nested_effect_corr < -1 || nested_effect_corr > 1) {
    stop("nested_effect_corr must be a number between -1 and 1")
  }
  if (!is.numeric(driver_sharing) || length(driver_sharing) != 1L ||
      !is.finite(driver_sharing) || driver_sharing < 0 || driver_sharing > 1) {
    stop("driver_sharing must be a number between 0 and 1")
  }
  if (!is.numeric(snp_driver_groups) || length(snp_driver_groups) != 1L ||
      !is.finite(snp_driver_groups) || snp_driver_groups < 1 ||
      snp_driver_groups != floor(snp_driver_groups) ||
      any(snp_driver_groups > n_traits_per_module)) {
    stop("snp_driver_groups must be an integer between 1 and n_traits_per_module")
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
  scalar_nonneg <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
      stop(name, " must be a non-negative finite number")
    }
  }
  scalar_nonneg(min_abs_z, "min_abs_z")
  scalar_nonneg(effect_tail, "effect_tail")
  scalar_nonneg(background_sparsity_sd, "background_sparsity_sd")
  n_hub_traits_auto <- is.null(n_hub_traits)
  if (!n_hub_traits_auto) {
    if (!is.numeric(n_hub_traits) || length(n_hub_traits) != 1L ||
        !is.finite(n_hub_traits) || n_hub_traits < 0 ||
        n_hub_traits != floor(n_hub_traits)) {
      stop("n_hub_traits must be NULL or a non-negative integer")
    }
    n_hub_traits <- as.integer(n_hub_traits)
  }
  if (!is.numeric(hub_snp_fraction) || length(hub_snp_fraction) != 2L ||
      any(!is.finite(hub_snp_fraction)) || any(hub_snp_fraction < 0) ||
      any(hub_snp_fraction > 1) || hub_snp_fraction[1] > hub_snp_fraction[2]) {
    stop("hub_snp_fraction must be c(min, max) within [0, 1] with min <= max")
  }
  if (!is.null(p_negative)) {
    if (!is.numeric(p_negative) || length(p_negative) != 1L ||
        !is.finite(p_negative) || p_negative < 0 || p_negative > 1) {
      stop("p_negative must be NULL or a number between 0 and 1")
    }
  }
  if (!is.numeric(overlap_fraction) || length(overlap_fraction) != 1L ||
      !is.finite(overlap_fraction) || overlap_fraction < 0 ||
      overlap_fraction > 1) {
    stop("overlap_fraction must be a number between 0 and 1")
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
    overlap = overlap,
    overlap_fraction = overlap_fraction
  )

  module_of <- rep(0L, n_snps)
  fill_order <- order(
    vapply(geometry$memberships, length, integer(1)),
    seq_len(K)
  )
  for (m in fill_order) {
    unfilled <- module_of[geometry$memberships[[m]]] == 0L
    module_of[geometry$memberships[[m]][unfilled]] <- m
  }

  snp_z <- NULL
  if (identical(overlap, "nested") && K > 0) {
    z_parent <- rep(NA_real_, n_snps)
    parent_idx <- geometry$memberships[[1]]
    z_parent[parent_idx] <- stats::rnorm(length(parent_idx))
    snp_z <- vector("list", K)
    snp_z[[1]] <- z_parent
    for (m in 2:K) {
      idx <- geometry$memberships[[m]]
      snp_z[[m]] <- rep(NA_real_, n_snps)
      snp_z[[m]][idx] <- nested_effect_corr * z_parent[idx] +
        sqrt(1 - nested_effect_corr^2) * stats::rnorm(length(idx))
    }
  }

  n_bg_traits <- n_background_traits
  trait_ids <- 1:(1 + n_driver_rows + n_bg_traits)
  target_tid <- 1L
  driver_tids <- if (n_driver_rows > 0) 2:(1 + n_driver_rows) else integer(0)
  bg_tids <- if (n_bg_traits > 0) seq_len(n_bg_traits) + 1L + n_driver_rows else integer(0)

  roster <- .build_trait_roster(K, n_shared, n_specific)
  driver_module <- roster$module
  trait_first_module <- if (n_driver_rows > 0) {
    vapply(seq_len(n_driver_rows), function(r) {
      min(roster$module[roster$trait_row == r])
    }, integer(1))
  } else {
    integer(0)
  }

  is_molecular <- .assign_molecular_traits(
    n_driver_rows = n_driver_rows,
    n_bg_traits = n_bg_traits,
    p_molecular = p_molecular,
    p_molecular_drivers = p_molecular_drivers,
    n_specific = n_specific,
    roster = roster,
    module_annotations = module_annotations,
    K = K
  )

  M <- matrix(NA_real_, nrow = length(trait_ids), ncol = n_snps,
              dimnames = list(as.character(trait_ids), .sim_variant_ids(n_snps)))
  snp_module <- module_of

  if (identical(target_pattern, "module")) {
    target_signal <- ifelse(
      snp_module > 0,
      effect_size[pmax(snp_module, 1)] * ifelse(snp_module %% 2 == 1, 1, -1),
      0
    )
    M[as.character(target_tid), ] <- target_signal +
      stats::rnorm(n_snps, 0, noise_sd)
  } else {
    # Dense target: significant, heavy-tailed z at every SNP, positive in the
    # oriented frame (matching a real target-trait row after orientation).
    target_vals <- effect_size[1] * stats::rlnorm(n_snps, 0, 0.5) +
      stats::rnorm(n_snps, 0, noise_sd)
    if (min_abs_z > 0) {
      target_vals <- sign(target_vals) * pmax(min_abs_z, abs(target_vals))
    }
    M[as.character(target_tid), ] <- target_vals
  }

  trait_signs <- switch(sign_pattern,
    coherent = rep(1L, n_driver_rows),
    flipped = vapply(seq_len(n_driver_rows), function(r) {
      slot <- min(roster$slot[roster$trait_row == r])
      ifelse(slot %% 2 == 1, 1L, -1L)
    }, integer(1)),
    random = sample(c(-1L, 1L), n_driver_rows, replace = TRUE)
  )

  roster$effect_mult <- trait_signs[roster$trait_row]
  if (n_driver_rows > 0 && any(roster$is_shared)) {
    for (r in seq_len(n_driver_rows)) {
      pos <- which(roster$trait_row == r)
      if (length(pos) > 1) {
        z <- numeric(length(pos))
        z[1] <- stats::rnorm(1)
        for (k in 2:length(pos)) {
          z[k] <- trait_effect_corr * z[k - 1] +
            sqrt(1 - trait_effect_corr^2) * stats::rnorm(1)
        }
        roster$effect_mult[pos] <- trait_signs[r] * z
      }
    }
  }

  if (n_driver_rows > 0) {
    n_shared_per_module <- floor(n_traits_per_module * driver_sharing)
    n_groups <- max(1L, as.integer(snp_driver_groups))
    snp_groups <- vector("list", K)
    if (n_groups > 1) {
      for (m in seq_len(K)) {
        idx <- geometry$memberships[[m]]
        g <- cut(seq_along(idx), breaks = n_groups, labels = FALSE)
        snp_groups[[m]] <- stats::setNames(g, colnames(M)[idx])
      }
    }
    drivers_per_group <- n_traits_per_module / n_groups
    for (r in seq_len(n_driver_rows)) {
      row <- as.character(driver_tids[r])
      pos <- which(roster$trait_row == r)
      filled <- integer(0)
      for (k in seq_along(pos)) {
        m <- roster$module[pos[k]]
        slot0 <- roster$slot[pos[k]] - 1L
        in_module <- geometry$memberships[[m]]
        if (m > 1 && slot0 < n_shared_per_module[m]) {
          in_module <- union(in_module, geometry$memberships[[m - 1]])
        }
        if (n_groups > 1) {
          driver_group <- min(n_groups, floor(slot0 / drivers_per_group[m]) + 1)
          grp <- snp_groups[[m]][as.character(colnames(M)[in_module])]
          keep <- !is.na(grp) & grp == driver_group
          in_module <- in_module[keep]
          if (length(in_module) == 0) {
            in_module <- geometry$memberships[[m]][1]
          }
        }
        in_module <- setdiff(in_module, filled)
        if (length(in_module) > 0) {
          fill_prob <- if (is_molecular[r]) {
            p_active_molecular
          } else {
            1 - p_structural_zero[m]
          }
          effect_z <- if (!is.null(snp_z)) {
            zv <- z_parent[in_module]
            if (m > 1) {
              in_mod <- in_module %in% geometry$memberships[[m]]
              zv[in_mod] <- snp_z[[m]][in_module[in_mod]]
            }
            zv
          } else {
            rep(1, length(in_module))
          }
          n_fill <- length(in_module)
          active <- stats::runif(n_fill) < fill_prob
          # effect_tail is a per-TRAIT-ROW scale (shared across the module's
          # SNPs) so that a driver trait's values stay correlated across the
          # SNPs it drives — an independent per-cell multiplier would destroy
          # that shared signal. Per-cell noise_sd keeps within-module spread.
          trait_scale <- if (effect_tail > 0) {
            stats::rlnorm(1, 0, effect_tail)
          } else {
            1
          }
          vals <- effect_size[m] * roster$effect_mult[pos[k]] * effect_z *
            trait_scale + stats::rnorm(n_fill, 0, noise_sd)
          if (!is.null(p_negative)) {
            vals <- ifelse(stats::runif(n_fill) < p_negative, -vals, vals)
          }
          if (min_abs_z > 0) {
            vals <- sign(vals) * pmax(min_abs_z, abs(vals))
          }
          M[row, in_module] <- ifelse(active, vals, NA_real_)
          filled <- c(filled, in_module)
        }
      }
      off_module <- setdiff(seq_len(n_snps), filled)
      spur <- stats::runif(length(off_module)) < p_spurious
      vals <- rep(NA_real_, length(off_module))
      if (any(spur)) {
        if (min_abs_z > 0) {
          # Spurious associations are significant hits too, with random sign.
          vals[spur] <- sign(stats::rnorm(sum(spur))) *
            (min_abs_z + stats::rexp(sum(spur), rate = 1 / 3))
        } else {
          vals[spur] <- stats::rnorm(sum(spur), 0, noise_sd)
        }
      }
      M[row, off_module] <- vals
    }
  }
  if (n_bg_traits > 0) {
    bg_rate_mult <- if (background_sparsity_sd > 0) {
      stats::rlnorm(n_bg_traits, 0, background_sparsity_sd)
    } else {
      NULL
    }
    hub_idx <- integer(0)
    hub_scale <- NULL
    hub_latent <- NULL
    if (n_hub_traits_auto) {
      # Realism default: plant ~3% of background traits as pleiotropic hubs.
      n_hub_traits <- as.integer(round(0.03 * n_bg_traits))
    }
    if (n_hub_traits > 0) {
      bg_positions <- seq_len(n_bg_traits)
      hub_candidates <- bg_positions[
        !is_molecular[n_driver_rows + bg_positions]
      ]
      if (length(hub_candidates) < n_hub_traits) {
        if (n_hub_traits_auto) {
          n_hub_traits <- length(hub_candidates)
        } else {
          stop(
            "n_hub_traits (", n_hub_traits, ") exceeds the non-molecular ",
            "background traits available (", length(hub_candidates), ")"
          )
        }
      }
      if (n_hub_traits > 0) {
        hub_idx <- sample(hub_candidates, n_hub_traits)
      # Per-hub magnitude scale and a shared per-SNP strength latent: dense,
      # mostly-positive (in the oriented frame) rows that are also mutually
      # correlated through the shared latent, so they read as one global
      # factor competing with the modules — the way real pleiotropy matrices
      # contain a core of traits hitting most loci with target-aligned effects.
        hub_scale <- stats::rlnorm(n_hub_traits, 0, 0.25)
        # Wide shared per-SNP latent so hub rows are strongly mutually
        # correlated despite per-cell sign flips.
        hub_latent <- stats::rlnorm(n_snps, 0, 0.9)
      }
    }
    for (t in seq_along(bg_tids)) {
      row <- as.character(bg_tids[t])
      if (n_hub_traits > 0 && t %in% hub_idx) {
        h <- match(t, hub_idx)
        hub_frac <- stats::runif(
          1,
          min = hub_snp_fraction[1],
          max = hub_snp_fraction[2]
        )
        act <- stats::runif(n_snps) < hub_frac
        vals <- rep(NA_real_, n_snps)
        # Magnitude rides the shared per-SNP latent (hub rows correlate) on a
        # significant floor; the per-hub scale sets each trait's typical size.
        hub_mag <- hub_scale[h] * hub_latent * pmax(min_abs_z, 3)
        hub_sign <- if (!is.null(p_negative)) {
          ifelse(stats::runif(n_snps) < p_negative, -1, 1)
        } else {
          rep(1, n_snps)
        }
        vals[act] <- hub_sign[act] * hub_mag[act] +
          stats::rnorm(sum(act), 0, noise_sd)
        if (min_abs_z > 0) {
          vals[act] <- sign(vals[act]) * pmax(min_abs_z, abs(vals[act]))
        }
        M[row, ] <- vals
      } else {
        base_prob <- if (is_molecular[n_driver_rows + t]) {
          p_active_molecular
        } else {
          p_active_background
        }
        act_prob <- if (!is.null(bg_rate_mult)) {
          min(1, base_prob * bg_rate_mult[t])
        } else {
          base_prob
        }
        act <- stats::runif(n_snps) < act_prob
        vals <- rep(NA_real_, n_snps)
        if (min_abs_z > 0) {
          # Every observed background cell is a significant hit; noise is
          # absence (NA), never a small observed value.
          n_act <- sum(act)
          if (n_act > 0) {
            hit_sign <- sign(stats::rnorm(n_act))
            vals[act] <- hit_sign *
              (min_abs_z + stats::rexp(n_act, rate = 1 / 3))
          }
        } else {
          vals[act] <- stats::rnorm(sum(act), 0, noise_sd)
        }
        M[row, ] <- vals
      }
    }
  }

  M <- .cap_matrix_density(M, snps_per_trait, traits_per_snp,
                           protected_row = as.character(target_tid))

  annotations <- .sim_trait_annotations(
    trait_ids = trait_ids,
    target_tid = target_tid,
    driver_tids = driver_tids,
    trait_first_module = trait_first_module,
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
    rows <- unique(roster$trait_row[roster$module == m])
    return(annotations$trait_name[driver_tids[rows]])
  })

  molecular_trait_ids <- c(
    driver_tids[is_molecular[seq_len(n_driver_rows)]],
    bg_tids[is_molecular[n_driver_rows + seq_len(n_bg_traits)]]
  )
  molecular_traits <- annotations$trait_name[
    annotations$trait_id %in% molecular_trait_ids
  ]

  driver_overlap <- .driver_overlap_matrix(roster, K)

  return(list(
    trait_object = trait_object,
    ground_truth = list(
      module_of_snp = stats::setNames(module_of, colnames(M)),
      multi_module_snps = multi_snps,
      module_memberships = memberships_named,
      driver_traits = stats::setNames(driver_names, sprintf("module_%d", seq_len(K))),
      driver_overlap = driver_overlap,
      molecular_traits = molecular_traits,
      planted_annotations = module_annotations,
      seed = seed,
      parameters = list(
        n_coloc_groups = n_snps,
        K = K,
        module_sizes = geometry$sizes,
        n_background_snps_reserved = n_background_snps,
        n_background_snps_actual = sum(module_of == 0L),
        overlap = overlap,
        nested_effect_corr = nested_effect_corr,
        n_traits_per_module = n_traits_per_module,
        trait_overlap = trait_overlap,
        trait_effect_corr = trait_effect_corr,
        n_shared = n_shared,
        n_specific = n_specific,
        n_driver_rows = n_driver_rows,
        driver_sharing = driver_sharing,
        n_background_traits = n_bg_traits,
        n_traits = as.integer(1L + n_driver_rows + n_bg_traits),
        effect_size = effect_size,
        noise_sd = noise_sd,
        sign_pattern = sign_pattern,
        p_structural_zero = p_structural_zero,
        p_spurious = p_spurious,
        p_active_background = p_active_background,
        background_sparsity_sd = background_sparsity_sd,
        min_abs_z = min_abs_z,
        effect_tail = effect_tail,
        n_hub_traits = n_hub_traits,
        hub_snp_fraction = hub_snp_fraction,
        target_pattern = target_pattern,
        p_negative = p_negative,
        overlap_fraction = overlap_fraction,
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


.build_trait_roster <- function(K, n_shared, n_specific) {
  if (K == 0) {
    return(data.frame(
      trait_row = integer(0), module = integer(0), slot = integer(0),
      is_shared = logical(0), stringsAsFactors = FALSE
    ))
  }
  rows <- list()
  n_shared_pool <- max(n_shared)
  if (n_shared_pool > 0) {
    for (p in seq_len(n_shared_pool)) {
      for (m in which(n_shared >= p)) {
        rows[[length(rows) + 1L]] <- data.frame(
          trait_row = p, module = m, slot = p, is_shared = TRUE,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  row_counter <- n_shared_pool
  for (m in seq_len(K)) {
    if (n_specific[m] > 0) {
      for (j in seq_len(n_specific[m])) {
        row_counter <- row_counter + 1L
        rows[[length(rows) + 1L]] <- data.frame(
          trait_row = row_counter, module = m,
          slot = n_shared[m] + j, is_shared = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) {
    return(data.frame(
      trait_row = integer(0), module = integer(0), slot = integer(0),
      is_shared = logical(0), stringsAsFactors = FALSE
    ))
  }
  return(do.call(rbind, rows))
}


.assign_molecular_traits <- function(n_driver_rows, n_bg_traits, p_molecular,
                                     p_molecular_drivers, n_specific, roster,
                                     module_annotations, K) {
  n <- n_driver_rows + n_bg_traits
  is_molecular <- rep(FALSE, n)
  if (n == 0) {
    return(is_molecular)
  }
  frac <- if (is.null(p_molecular_drivers)) {
    rep(p_molecular, K)
  } else {
    p_molecular_drivers
  }
  if (n_driver_rows > 0) {
    first_mod <- vapply(seq_len(n_driver_rows), function(r) {
      min(roster$module[roster$trait_row == r])
    }, integer(1))
    # Planted module genes always make their drivers molecular.
    if (!is.null(module_annotations)) {
      for (r in seq_len(n_driver_rows)) {
        mods <- unique(roster$module[roster$trait_row == r])
        planted <- any(vapply(mods, function(m) {
          spec <- if (m <= length(module_annotations)) module_annotations[[m]] else NULL
          !is.null(spec) && !is.null(spec$genes) && length(spec$genes) > 0
        }, logical(1)))
        if (planted) {
          is_molecular[r] <- TRUE
        }
      }
    }
    # Shared traits: bernoulli draw by the first driving module's fraction.
    shared_rows <- unique(roster$trait_row[roster$is_shared])
    for (r in shared_rows) {
      if (!is_molecular[r]) {
        is_molecular[r] <- stats::runif(1) < frac[first_mod[r]]
      }
    }
    # Specific traits: exact per-module counts keep a dense phenotype core.
    for (m in seq_len(K)) {
      idx <- unique(roster$trait_row[!roster$is_shared & roster$module == m])
      if (length(idx) == 0) {
        next
      }
      desired <- round(frac[m] * length(idx))
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
    bg_idx <- n_driver_rows + seq_len(n_bg_traits)
    desired <- round(p_molecular * n_bg_traits)
    candidates <- bg_idx[!is_molecular[bg_idx]]
    pick <- sample(candidates, min(desired, length(candidates)))
    is_molecular[pick] <- TRUE
  }
  return(is_molecular)
}


.driver_overlap_matrix <- function(roster, K) {
  if (K == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  traits_by_mod <- lapply(seq_len(K), function(m) {
    unique(roster$trait_row[roster$module == m])
  })
  mat <- matrix(0, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      mat[i, j] <- length(intersect(traits_by_mod[[i]], traits_by_mod[[j]])) /
        length(traits_by_mod[[i]])
    }
  }
  return(mat)
}


.plant_module_regions <- function(n_snps, K, module_sizes, n_background_snps,
                                  overlap = "disjoint", overlap_fraction = 0.2) {
  if (K == 0) {
    return(list(regions = list(), memberships = list(), sizes = integer(0),
                overlap = overlap))
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

  if (identical(overlap, "nested")) {
    if (any(sizes[1] < sizes[-1])) {
      stop("for overlap = 'nested', module_sizes[1] must be the largest module (the parent)")
    }
    if (sizes[1] > capacity) {
      stop("modules need ", sizes[1], " SNPs for the parent module but only ",
           capacity, " available")
    }
    memberships <- vector("list", K)
    memberships[[1]] <- seq_len(sizes[1])
    if (K > 1) {
      for (m in 2:K) {
        memberships[[m]] <- sort(sample(memberships[[1]], sizes[m], replace = FALSE))
      }
    }
    return(list(
      regions = memberships,
      memberships = memberships,
      sizes = sizes,
      overlap = overlap
    ))
  }

    consumed <- sum(sizes)
  if (consumed > n_snps) {
    stop("modules need ", consumed, " SNPs but only ", n_snps, " available")
  }

  regions <- vector("list", K)
  start <- 1L
  for (m in seq_len(K)) {
    end <- start + sizes[m] - 1L
    regions[[m]] <- start:end
    start <- end + 1L
  }

  if (identical(overlap, "partial")) {
    # Disjoint cores, but each module m > 1 additionally claims a fraction of
    # the previous module's core, so adjacent modules share boundary SNPs.
    memberships <- vector("list", K)
    memberships[[1]] <- regions[[1]]
    for (m in seq_len(K)[-1]) {
      n_borrow <- as.integer(round(overlap_fraction * sizes[m]))
      borrowed <- integer(0)
      if (n_borrow > 0) {
        borrowed <- sample(
          regions[[m - 1]],
          min(n_borrow, length(regions[[m - 1]]))
        )
      }
      memberships[[m]] <- sort(union(regions[[m]], borrowed))
    }
    return(list(
      regions = regions,
      memberships = memberships,
      sizes = sizes,
      overlap = overlap
    ))
  }

  memberships <- regions
  return(list(
    regions = regions,
    memberships = memberships,
    sizes = sizes,
    overlap = overlap
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
                                   trait_first_module,
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
      m <- trait_first_module[d]
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