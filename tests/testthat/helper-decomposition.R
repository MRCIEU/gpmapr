# Small simulated decompositions and correspondence shared by the
# decompose_trait(), trait_correspondence() and create_report() tests. Built
# once per test run: the fits take a few seconds.
.decomposition_cache <- new.env(parent = emptyenv())

fast_settings <- function(...) {
  return(gpmap_settings(similarity_n_perm = 200, n_rep = 0, n_perm = 200, ...))
}

# Pathway mappings covering a trait object's genes (no API call).
fake_pathway_mappings <- function(coloc_groups) {
  genes <- unique(coloc_groups[!is.na(coloc_groups$gene_id), c("gene_id", "gene")])
  if (nrow(genes) == 0) {
    genes <- data.frame(gene_id = 1L, gene = "G1")
  }
  return(list(mappings = data.frame(
    gene_id = genes$gene_id,
    term_id = paste0("path", (seq_len(nrow(genes)) %% 4) + 1),
    source = rep(c("KEGG", "Reactome"), length.out = nrow(genes)),
    description = paste("Pathway", (seq_len(nrow(genes)) %% 4) + 1),
    stringsAsFactors = FALSE
  )))
}

# A duplicates table in the API 1.1.0 shape, with a duplicate group and a proxy
# group around trait "100" and an unrelated group.
fake_duplicates <- function() {
  return(data.frame(
    trait_id = c(100L, 101L, 102L, 103L, 104L, 200L, 201L),
    trait_name = c("Anchor", "Anchor (UKB)", "Anchor (SPA)", "Anchor diagnosis",
                   "Anchor medication", "Other", "Other copy"),
    duplicate_group = c(1L, 1L, 1L, 2L, 3L, 4L, 4L),
    is_representative = c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
    match_rules = c(NA, "ukb_field", "firth_spa", NA, NA, NA, "synonym"),
    proxy_group_id = c(10L, 10L, 10L, 10L, 10L, NA, NA),
    proxy_group_name = c(rep("Anchor", 5), NA, NA),
    proxy_rules = c(NA, NA, NA, "icd10", "wrapper;concept_map", NA, NA),
    stringsAsFactors = FALSE
  ))
}

simulated_pair <- function() {
  if (is.null(.decomposition_cache$pair)) {
    .decomposition_cache$pair <- simulate_trait_pair(seed = 1)
  }
  return(.decomposition_cache$pair)
}

simulated_decompositions <- function() {
  if (is.null(.decomposition_cache$decomps)) {
    sim <- simulated_pair()
    .decomposition_cache$decomps <- lapply(sim$traits, function(obj) {
      return(decompose_trait(
        obj, settings = fast_settings(), proxies = FALSE, remove_duplicates = FALSE,
        annotate = TRUE, mappings = fake_pathway_mappings(obj$coloc_groups),
        cores = 1, verbose = FALSE
      ))
    })
  }
  return(.decomposition_cache$decomps)
}

simulated_correspondence <- function() {
  if (is.null(.decomposition_cache$corr)) {
    .decomposition_cache$corr <- trait_correspondence(
      simulated_decompositions(), settings = fast_settings(), annotate = TRUE, verbose = FALSE
    )
  }
  return(.decomposition_cache$corr)
}
