#' @title Simulate a Trait Object With Planted SNP Modules
#' @description Generate a synthetic trait object (same shape as
#' `trait(id, include_associations = TRUE)`) with an optional planted latent
#' module structure over SNPs, for testing the univariate clustering pipeline.
#' The returned object feeds directly into `run_univariate_clustering()`, and
#' the accompanying ground truth supports recovery scoring via
#' `evaluate_univariate_simulation()`.
#'
#' Structure: `n_coloc_groups` SNPs (one coloc group / locus each). `K`
#' contiguous module regions are laid out left-to-right; adjacent regions share
#' `overlap` SNPs (shared SNPs have driver traits from both modules). SNPs not
#' covered by any module are unstructured background: sparse random pleiotropic
#' profiles with no shared latent structure.
#'
#' Effects: each module has `drivers_per_module` driver traits with
#' `effect_size`-scale effects on their module's SNPs. Driver support is thinned
#' by `p_structural_zero` (structural zeros inside true support). Cells outside
#' true support receive small noise effects with probability `p_spurious`.
#' Background traits activate on each SNP independently with probability
#' `p_active_background`, giving realistic marginal sparsity without structure.
#'
#' Annotations: `module_annotations` optionally plants genes / tissues /
#' trait categories on each module's driver traits so enrichment machinery can
#' be verified; unplanted components draw randomly from built-in pools, and
#' `annotation_noise` controls contamination of planted labels.
#' @param n_coloc_groups Integer number of SNPs (= coloc groups) for the target trait.
#' @param K Integer number of planted modules. `K = 0` gives the pure null (all
#'   SNPs unstructured).
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
#' @param drivers_per_module Number of driver traits per module.
#' @param n_background_traits Number of unstructured background traits.
#' @param effect_size Mean absolute z-score of driver effects inside their module.
#' @param noise_sd Standard deviation of noise effects added everywhere.
#' @param sign_pattern Driver-sign regime: `"coherent"` (all drivers positive),
#'   `"flipped"` (signs alternate across drivers within a module, so
#'   anti-correlated profiles share a module), `"random"` (independent signs).
#' @param p_structural_zero Probability that a cell inside a module's true
#'   support is absent entirely (structural zero).
#' @param p_spurious Probability that a cell outside true support carries a
#'   small noise effect rather than being absent.
#' @param p_active_background Per-SNP activation probability for background
#'   trait profiles.
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
#' @param seed Optional RNG seed; a random one is drawn and recorded when `NULL`.
#' @return A list with:
#'   \itemize{
#'     \item trait_object: simulated trait result usable by
#'       `run_univariate_clustering()`
#'     \item ground_truth: list with `module_of_snp` (named vector, `0` =
#'       background; shared SNPs take their earlier module), `multi_module_snps`
#'       (shared-SNP names), `module_memberships` (SNP index list, shared SNPs
#'       appear in both), `driver_traits` (per-module trait names),
#'       `planted_annotations`, `seed`, and `parameters`
#'   }
#' @export
simulate_trait <- function(n_coloc_groups = 100,
                           K = 0,
                           module_sizes = NULL,
                           n_background_snps = 0,
                           overlap = "disjoint",
                           drivers_per_module = 5L,
                           n_background_traits = 40L,
                           effect_size = 3,
                           noise_sd = 0.5,
                           sign_pattern = c("coherent", "flipped", "random"),
                           p_structural_zero = 0.4,
                           p_spurious = 0.05,
                           p_active_background = 0.08,
                           snps_per_trait = NULL,
                           traits_per_snp = NULL,
                           module_annotations = NULL,
                           annotation_noise = 0.1,
                           seed = NULL) {
  sign_pattern <- match.arg(sign_pattern)

  if (!is.numeric(n_coloc_groups) || length(n_coloc_groups) != 1 ||
        n_coloc_groups < 2 || n_coloc_groups != floor(n_coloc_groups)) {
    stop("n_coloc_groups must be an integer >= 2")
  }
  if (!is.numeric(K) || length(K) != 1 || K < 0 || K != floor(K)) {
    stop("K must be a non-negative integer")
  }
  n_snps <- as.integer(n_coloc_groups)
  K <- as.integer(K)

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

  n_drivers <- K * drivers_per_module
  n_bg_traits <- as.integer(n_background_traits)
  trait_ids <- 1:(1 + n_drivers + n_bg_traits)
  target_tid <- 1L
  driver_tids <- if (n_drivers > 0) 2:(1 + n_drivers) else integer(0)
  bg_tids <- if (n_bg_traits > 0) seq_len(n_bg_traits) + 1L + n_drivers else integer(0)
  driver_module <- if (n_drivers > 0) rep(seq_len(K), each = drivers_per_module) else integer(0)

  M <- matrix(NA_real_, nrow = length(trait_ids), ncol = n_snps,
              dimnames = list(as.character(trait_ids), .sim_variant_ids(n_snps)))
  snp_module <- module_of

  M[as.character(target_tid), ] <- ifelse(
    snp_module > 0,
    effect_size * ifelse(snp_module %% 2 == 1, 1, -1),
    stats::rnorm(n_snps, 0, noise_sd)
  )

  driver_signs <- switch(sign_pattern,
    coherent = rep(1, n_drivers),
    flipped = ifelse(seq_len(max(n_drivers, 1)) %% 2 == 1, 1, -1)[seq_len(n_drivers)],
    random = sample(c(-1, 1), n_drivers, replace = TRUE)
  )
  if (n_drivers > 0) {
    for (d in seq_len(n_drivers)) {
      m <- driver_module[d]
      row <- as.character(driver_tids[d])
      in_module <- geometry$memberships[[m]]
      M[row, in_module] <- ifelse(
        stats::runif(length(in_module)) < (1 - p_structural_zero),
        effect_size * driver_signs[d] + stats::rnorm(length(in_module), 0, noise_sd),
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
      act <- stats::runif(n_snps) < p_active_background
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
    K = K
  )

  cg <- .sim_coloc_groups_from_matrix(M, annotations)
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

  return(list(
    trait_object = trait_object,
    ground_truth = list(
      module_of_snp = stats::setNames(module_of, colnames(M)),
      multi_module_snps = multi_snps,
      module_memberships = memberships_named,
      driver_traits = stats::setNames(driver_names, sprintf("module_%d", seq_len(K))),
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
        n_background_traits = n_bg_traits,
        effect_size = effect_size,
        noise_sd = noise_sd,
        sign_pattern = sign_pattern,
        p_structural_zero = p_structural_zero,
        p_spurious = p_spurious,
        p_active_background = p_active_background,
        annotation_noise = annotation_noise
      )
    )
  ))
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
    if (length(traits_per_snp) != 2 || any(traits_per_snp < 0)) {
      stop("traits_per_snp must be c(min, max)")
    }
    counts <- colSums(!is.na(M))
    over <- counts > traits_per_snp[2]
    for (c in colnames(M)[over]) {
      drop <- sample(which(!is.na(M[, c])), counts[c] - traits_per_snp[2])
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
      if (!is.null(spec) && !is.null(spec$genes) && length(spec$genes) > 0) {
        sym <- spec$genes[(d - 1) %% length(spec$genes) + 1]
        out$gene[tid_row] <- sym
        out$gene_id[tid_row] <- register_gene(sym)
      } else {
        sym <- draw(pools$genes, 0.3)
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
      sym <- draw(pools$genes, 0.2)
      out$gene[tid_row] <- sym
      out$gene_id[tid_row] <- register_gene(sym)
      out$tissue[tid_row] <- draw(pools$tissues, 0.5)
      out$trait_category[tid_row] <- draw(pools$trait_categories, 0.5)
    }
  }

  return(out)
}


.sim_coloc_groups_from_matrix <- function(M, annotations) {
  observed <- which(!is.na(M), arr.ind = TRUE)
  if (nrow(observed) == 0) {
    stop("simulated matrix has no observations")
  }
  trait_idx <- as.integer(rownames(M))[observed[, "row"]]
  snp_idx <- observed[, "col"]
  z <- M[observed]

  anno <- annotations[match(trait_idx, annotations$trait_id), , drop = FALSE]
  beta <- z
  se <- 1
  p <- pmax(2 * stats::pnorm(-abs(z)), 1e-300)

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
