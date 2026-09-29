#' @title Simulate Two Target Traits Sharing Planted Programs
#' @description Generate a pair of trait objects that share planted programs,
#' with **locus overlap** (what fraction of a shared program's SNPs are the same
#' coloc group in both traits) and **profile overlap** (what fraction of its
#' driver studies are the same) controlled separately. This is the generative
#' model for the multi-trait correspondence study.
#'
#' Effects are generated in a **common allele frame**, as summary statistics
#' are reported. Each program locus `s` carries an effect `a_s` on its program,
#' drawn once and reused in both traits at a shared locus. A driver study's
#' effect at `s` is proportional to `a_s` and has the same sign in both traits'
#' matrices, because it is the same GWAS. The target trait's own effect at `s`
#' is `sign * a_s`, where `sign` is `+1` in the first trait and the program's
#' planted direction in the second, so an antagonistic program has opposite-sign
#' target effects at its loci. `run_univariate_clustering()` then orients each
#' matrix to its own target's risk allele, exactly as on real data, which is
#' where antagonism reaches the fitted programs. Because the target effects
#' vary with `a_s`, [module_rg()] can recover the planted direction from the
#' two targets' own effects.
#'
#' Set `K_shared = 0` for the multi-trait null: two traits with programs of their
#' own and no correspondence between them, so any link reported is a false
#' discovery.
#' @param n_loci_per_trait Loci each target trait carries.
#' @param n_shared_loci Size of the pool of loci both traits carry. Caps how much
#'   locus overlap is achievable; must be at least
#'   `K_shared * module_size * locus_overlap`.
#' @param K_shared Number of programs planted in **both** traits.
#' @param K_specific Number of programs planted in only one trait, per trait.
#'   These have no counterpart and must not be linked.
#' @param module_size SNPs per planted program.
#' @param n_drivers_per_program Background studies loading on each program.
#' @param locus_overlap Fraction (0-1) of a shared program's SNPs that are the
#'   same locus in both traits. `0` gives disjoint loci, so the locus axis has
#'   nothing to score.
#' @param profile_overlap Fraction (0-1) of a shared program's driver studies
#'   that are the same in both traits.
#' @param directions Optional `+1`/`-1` per shared program: `-1` gives the
#'   program an opposite-sign effect on the second trait (target effects of
#'   opposite sign at its loci, in the common allele frame), making it
#'   antagonistic. Defaults to alternating, so sign accuracy is scoreable in
#'   both directions.
#' @param n_background_traits Unstructured background studies, shared between
#'   the two traits so the feature axis intersects.
#' @param effect_size Mean absolute z of driver and target effects inside a
#'   program.
#' @param locus_effect_sd Spread of the per-locus program effect `a_s`, drawn as
#'   `|N(1, locus_effect_sd)|`. It is what makes the two targets' effects
#'   correlated across a shared program's loci, so it sets how well
#'   [module_rg()] can determine direction. Defaults to `0.5`.
#' @param noise_sd Noise added to every observed cell.
#' @param p_structural_zero Probability a cell inside a program's true support is
#'   absent entirely.
#' @param p_active_background,background_sparsity_sd,background_effect_scale
#'   Background realism, as in `simulate_trait()`.
#' @param background_corr,n_bg_factors,snp_pleiotropy_sd Similarity-graph
#'   realism, as in `simulate_trait()`. Defaults reproduce the historical
#'   independent, random-signed background.
#' @param trait_ids Length-2 integer ids for the two target traits.
#' @param seed RNG seed.
#' @return A list with:
#'   \itemize{
#'     \item traits: named list of two trait objects, each consumable by
#'       `run_univariate_clustering()`
#'     \item ground_truth: `correspondence` (the planted cross-trait links and
#'       their direction), `programs` (per-trait planted SNP and driver sets),
#'       `program_snps`, `program_drivers`, `program_effects` (each program's
#'       per-locus effect `a_s`, named by SNP), `shared_loci`, and `parameters`
#'   }
#' @export
simulate_trait_pair <- function(n_loci_per_trait = 200L,
                                n_shared_loci = 80L,
                                K_shared = 3L,
                                K_specific = 2L,
                                module_size = 20L,
                                n_drivers_per_program = 8L,
                                locus_overlap = 0.5,
                                profile_overlap = 0.8,
                                directions = NULL,
                                n_background_traits = 300L,
                                effect_size = 6,
                                locus_effect_sd = 0.5,
                                noise_sd = 1,
                                p_structural_zero = 0,
                                p_active_background = 0.02,
                                background_sparsity_sd = 1.2,
                                background_effect_scale = 0.5,
                                background_corr = 0,
                                n_bg_factors = 4L,
                                snp_pleiotropy_sd = 0,
                                trait_ids = c(9001L, 9002L),
                                seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  K_shared <- as.integer(K_shared)
  K_specific <- as.integer(K_specific)
  module_size <- as.integer(module_size)
  n_shared_loci <- as.integer(n_shared_loci)
  n_loci_per_trait <- as.integer(n_loci_per_trait)
  if (length(trait_ids) != 2L || anyDuplicated(trait_ids) > 0) {
    stop("trait_ids must be two distinct ids")
  }
  if (locus_overlap < 0 || locus_overlap > 1) {
    stop("locus_overlap must be between 0 and 1")
  }
  if (profile_overlap < 0 || profile_overlap > 1) {
    stop("profile_overlap must be between 0 and 1")
  }
  if (n_shared_loci > n_loci_per_trait) {
    stop("n_shared_loci cannot exceed n_loci_per_trait")
  }

  n_shared_snps <- as.integer(round(module_size * locus_overlap))
  if (K_shared * n_shared_snps > n_shared_loci) {
    stop(
      "shared locus pool too small: need at least ",
      K_shared * n_shared_snps, " shared loci for K_shared = ", K_shared,
      " at locus_overlap = ", locus_overlap
    )
  }
  n_private_needed <- (K_shared + K_specific) * module_size
  n_private_per_trait <- n_loci_per_trait - n_shared_loci
  if (n_private_needed > n_private_per_trait + K_shared * n_shared_snps) {
    stop(
      "not enough private loci: raise n_loci_per_trait or lower ",
      "K_shared / K_specific / module_size"
    )
  }

  # --- locus pool -----------------------------------------------------------
  # coloc_group_id is what makes a locus "the same" locus across traits, so the
  # shared pool takes the lowest ids and each trait's private loci take their own
  # disjoint block.
  shared_loci <- seq_len(n_shared_loci)
  private_a <- n_shared_loci + seq_len(n_private_per_trait)
  private_b <- n_shared_loci + n_private_per_trait + seq_len(n_private_per_trait)
  loci_a <- c(shared_loci, private_a)
  loci_b <- c(shared_loci, private_b)

  # --- trait id pools -------------------------------------------------------
  # Driver studies are drawn from one pool shared by both traits, so a shared
  # driver is literally the same feature id in both fits -- which is what the
  # profile axis needs.
  n_driver_pool <- (K_shared * 2L + K_specific * 2L) * n_drivers_per_program
  driver_pool <- 20000L + seq_len(n_driver_pool)
  background_pool <- 30000L + seq_len(as.integer(n_background_traits))

  # --- plant programs -------------------------------------------------------
  take_private <- function(pool, n, used) {
    available <- setdiff(pool, used)
    if (length(available) < n) {
      stop("private locus pool exhausted")
    }
    return(sample(available, n))
  }
  used_a <- integer(0)
  used_b <- integer(0)
  used_shared <- integer(0)
  driver_cursor <- 0L
  next_drivers <- function(n) {
    idx <- driver_cursor + seq_len(n)
    driver_cursor <<- driver_cursor + n
    return(driver_pool[idx])
  }

  if (is.null(directions)) {
    directions <- rep(c(1L, -1L), length.out = max(K_shared, 1L))
  }
  directions <- as.integer(rep(directions, length.out = max(K_shared, 1L)))

  # Per-locus program effect a_s. A shared locus draws it once, so both traits
  # see the same locus effect there, as they would from one causal variant.
  draw_locus_effect <- function(n) {
    return(abs(stats::rnorm(n, 1, locus_effect_sd)))
  }

  programs <- list()
  correspondence <- list()

  for (k in seq_len(K_shared)) {
    snps_shared <- if (n_shared_snps > 0) {
      s <- sample(setdiff(shared_loci, used_shared), n_shared_snps)
      used_shared <- c(used_shared, s)
      s
    } else {
      integer(0)
    }
    n_own <- module_size - n_shared_snps
    own_a <- take_private(private_a, n_own, used_a)
    own_b <- take_private(private_b, n_own, used_b)
    used_a <- c(used_a, own_a)
    used_b <- c(used_b, own_b)

    n_shared_drivers <- as.integer(round(n_drivers_per_program * profile_overlap))
    drivers_shared <- next_drivers(n_shared_drivers)
    n_own_drivers <- n_drivers_per_program - n_shared_drivers
    drivers_a <- c(drivers_shared, next_drivers(n_own_drivers))
    drivers_b <- c(drivers_shared, next_drivers(n_own_drivers))

    effect_shared <- draw_locus_effect(length(snps_shared))

    label <- sprintf("shared%d", k)
    programs[[length(programs) + 1L]] <- list(
      trait_id = trait_ids[[1]], label = label,
      snps = c(snps_shared, own_a), drivers = drivers_a, sign = 1L,
      effect = c(effect_shared, draw_locus_effect(n_own))
    )
    programs[[length(programs) + 1L]] <- list(
      trait_id = trait_ids[[2]], label = label,
      snps = c(snps_shared, own_b), drivers = drivers_b,
      sign = directions[[k]],
      effect = c(effect_shared, draw_locus_effect(n_own))
    )
    correspondence[[length(correspondence) + 1L]] <- data.frame(
      program_label = label,
      trait_id_a = trait_ids[[1]],
      trait_id_b = trait_ids[[2]],
      direction = if (directions[[k]] > 0) "concordant" else "antagonistic",
      n_loci_shared = length(snps_shared),
      n_drivers_shared = n_shared_drivers,
      stringsAsFactors = FALSE
    )
  }

  for (side in 1:2) {
    for (k in seq_len(K_specific)) {
      pool <- if (side == 1L) private_a else private_b
      used <- if (side == 1L) used_a else used_b
      own <- take_private(pool, module_size, used)
      if (side == 1L) used_a <- c(used_a, own) else used_b <- c(used_b, own)
      programs[[length(programs) + 1L]] <- list(
        trait_id = trait_ids[[side]],
        label = sprintf("specific%d_t%d", k, side),
        snps = own, drivers = next_drivers(n_drivers_per_program), sign = 1L,
        effect = draw_locus_effect(module_size)
      )
    }
  }

  # --- build each trait's matrix -------------------------------------------
  build_one <- function(target_id, loci) {
    progs <- Filter(function(p) p$trait_id == target_id, programs)
    driver_ids <- sort(unique(unlist(lapply(progs, function(p) p$drivers))))
    row_ids <- c(target_id, driver_ids, background_pool)
    M <- matrix(
      NA_real_,
      nrow = length(row_ids), ncol = length(loci),
      dimnames = list(as.character(row_ids), paste0("snp", loci))
    )
    col_of <- stats::setNames(seq_along(loci), as.character(loci))

    # Target row: dense and significant at every locus. Outside any program it
    # is positive; inside a program it follows that program's locus effects and
    # its sign on this target (below).
    M[as.character(target_id), ] <- abs(
      stats::rnorm(length(loci), effect_size, noise_sd)
    )

    # Everything is in the common allele frame: drivers carry the same sign in
    # both traits, and an antagonistic program instead has a negative target
    # effect. run_univariate_clustering() orients each column to the target's
    # risk allele, which is what turns that into negative driver cells.
    for (p in progs) {
      cols <- col_of[as.character(p$snps)]
      keep <- !is.na(cols)
      cols <- cols[keep]
      if (length(cols) == 0) {
        next
      }
      a <- p$effect[keep]
      M[as.character(target_id), cols] <- p$sign * effect_size * a +
        stats::rnorm(length(cols), 0, noise_sd)
      for (d in p$drivers) {
        z <- abs(stats::rnorm(1, 1, 0.2))
        vals <- effect_size * z * a +
          stats::rnorm(length(cols), 0, noise_sd)
        if (p_structural_zero > 0) {
          drop <- stats::runif(length(cols)) < p_structural_zero
          vals[drop] <- NA_real_
        }
        M[as.character(d), cols] <- vals
      }
    }

    # Unstructured background, with the same realism controls as simulate_trait().
    bg_rows <- as.character(background_pool)
    n_bg <- length(bg_rows)
    rate_mult <- if (background_sparsity_sd > 0) {
      stats::rlnorm(n_bg, 0, background_sparsity_sd) /
        exp(background_sparsity_sd^2 / 2)
    } else {
      rep(1, n_bg)
    }
    snp_mult <- if (snp_pleiotropy_sd > 0) {
      stats::rlnorm(ncol(M), 0, snp_pleiotropy_sd) /
        exp(snp_pleiotropy_sd^2 / 2)
    } else {
      rep(1, ncol(M))
    }
    components <- .sim_background_components(
      n_snps = ncol(M), n_bg_traits = n_bg,
      n_bg_factors = n_bg_factors, background_corr = background_corr
    )
    for (i in seq_len(n_bg)) {
      active <- stats::runif(ncol(M)) <
        pmin(1, p_active_background * rate_mult[i] * snp_mult)
      if (!any(active)) {
        next
      }
      M[bg_rows[i], active] <- .sim_background_effects(
        sum(active), background_effect_scale, effect_size,
        shared = if (is.null(components)) NULL else components[i, active]
      )
    }

    annotations <- data.frame(
      trait_id = row_ids,
      trait_name = ifelse(
        row_ids == target_id, sprintf("Simulated trait %d", target_id),
        sprintf("Study %d", row_ids)
      ),
      tissue = NA_character_,
      gene = NA_character_,
      gene_id = NA_character_,
      trait_category = ifelse(
        row_ids %in% driver_pool, "Driver", "Background"
      ),
      stringsAsFactors = FALSE
    )
    annotations$trait_category[annotations$trait_id == target_id] <- "Target"

    cg <- .sim_coloc_groups_from_matrix(M, annotations)
    # coloc_group_id must be the GLOBAL locus id, not the column position, or the
    # two traits cannot be recognised as sharing a locus.
    cg$coloc_group_id <- loci[cg$coloc_group_id]
    cg$bp <- as.integer(cg$coloc_group_id) * 1000L

    return(list(
      trait_object = list(
        trait = list(
          id = target_id,
          trait = "SIMULATED",
          trait_name = sprintf("Simulated trait %d", target_id),
          source_url = NA_character_
        ),
        coloc_groups = cg
      ),
      x_matrix = M
    ))
  }

  built <- lapply(trait_ids, function(tid) {
    build_one(tid, if (tid == trait_ids[[1]]) loci_a else loci_b)
  })
  names(built) <- as.character(trait_ids)

  program_table <- dplyr::bind_rows(lapply(programs, function(p) {
    data.frame(
      trait_id = p$trait_id, program_label = p$label,
      n_snps = length(p$snps), n_drivers = length(p$drivers),
      sign = p$sign, stringsAsFactors = FALSE
    )
  }))

  return(list
  (
    traits = lapply(built, function(x) x$trait_object),
    x_matrices = lapply(built, function(x) x$x_matrix),
    ground_truth = list(
      correspondence = if (length(correspondence) > 0) {
        dplyr::bind_rows(correspondence)
      } else {
        data.frame(
          program_label = character(0), trait_id_a = integer(0),
          trait_id_b = integer(0), direction = character(0),
          n_loci_shared = integer(0), n_drivers_shared = integer(0),
          stringsAsFactors = FALSE
        )
      },
      programs = program_table,
      program_snps = stats::setNames(
        lapply(programs, function(p) paste0("snp", p$snps)),
        vapply(programs, function(p) {
          paste0(p$trait_id, ":", p$label)
        }, character(1))
      ),
      program_drivers = stats::setNames(
        lapply(programs, function(p) as.character(p$drivers)),
        vapply(programs, function(p) {
          paste0(p$trait_id, ":", p$label)
        }, character(1))
      ),
      program_effects = stats::setNames(
        lapply(programs, function(p) {
          stats::setNames(p$effect, paste0("snp", p$snps))
        }),
        vapply(programs, function(p) {
          paste0(p$trait_id, ":", p$label)
        }, character(1))
      ),
      shared_loci = paste0("snp", shared_loci),
      seed = seed,
      parameters = list(
        n_loci_per_trait = n_loci_per_trait,
        n_shared_loci = n_shared_loci,
        K_shared = K_shared, K_specific = K_specific,
        module_size = module_size,
        n_drivers_per_program = as.integer(n_drivers_per_program),
        locus_overlap = locus_overlap, profile_overlap = profile_overlap,
        n_background_traits = as.integer(n_background_traits),
        effect_size = effect_size, locus_effect_sd = locus_effect_sd,
        noise_sd = noise_sd,
        p_structural_zero = p_structural_zero,
        p_active_background = p_active_background,
        background_sparsity_sd = background_sparsity_sd,
        background_effect_scale = background_effect_scale,
        background_corr = background_corr,
        snp_pleiotropy_sd = snp_pleiotropy_sd,
        trait_ids = trait_ids
      )
    )
  ))
}
