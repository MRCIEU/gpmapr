# gpmapr 0.1.1

First CRAN release. Works with versions 1.0.1 and 1.1.0 of the GPMap API.

* An R client for the 'OpenGWAS' Genotype-Phenotype Map API
  (<https://gpmap.opengwas.io>): look up traits, genes, variants and regions
  with their colocalisation groups, rare-variant results and associations;
  search the map with `search_gpmap()`; fetch LD proxies, pleiotropy summaries
  and duplicate traits; upload, fetch and delete your own GWAS.
* Requests that hit the API rate limit (HTTP 429) are retried with exponential
  backoff.
* The API version is checked on the first request of a session rather than
  when the package loads, so `library(gpmapr)` needs no network access. If the
  installed gpmapr does not support the API version, the first request stops
  with instructions to update the package.
* API errors and connection failures now give an informative error message.
* New `gpmap_available()` reports whether the API can be reached.
* `search_gpmap()` for a variant now reports a per-row number of colocalising
  studies and rare results (previously every row showed the first row's
  counts), and gains an `rsq` column with each proxy's r-squared with the
  searched variant. Variant searches no longer download the full list of
  traits and genes.
* Coloc pairs from GWAS uploads use consistent column names:
  `upload_study_extraction_a_id` / `_b_id` for the upload's study extractions
  and `existing_study_extraction_a_id` / `_b_id` for those already in the map
  (previously `study_extraction_id_a` and `existing_study_extraction_id_a`).
  The old names sent by earlier versions of the API are renamed automatically.
* Coloc pairs from `gene()`, `genes()`, `variant()` and `get_gwas()` include the
  trait name, data type, tissue and gene of both sides (`trait_name_a`,
  `trait_name_b`, ...), so they no longer need joining to `study_extractions`.
  Coloc pairs from GWAS uploads also say whether each side is the upload or a
  study already in the map (`source_a`, `source_b`).
* Breaking: column names follow version 1.1.0 of the GPMap API, whichever
  version the API is. Names sent by API 1.0.1 are renamed automatically.
  * `trait` is now `pipeline_id` (the GPMap pipeline's id for the trait, such
    as its OpenGWAS dataset id) and `category` is now `variable_type`
    (continuous or categorical), in `all_traits()`, `trait()`, `traits()` and
    GWAS uploads.
  * In GWAS upload coloc groups, `study_extraction_id` is now
    `upload_study_extraction_id`.
* Breaking: `$associations` is no longer returned. With
  `include_associations = TRUE`, each study's association at its variant
  (`beta`, `se`, `p`, `eaf`, `imputed`) is added to `coloc_groups`,
  `rare_results` and `study_extractions` (and, for GWAS uploads, to
  `upload_study_extractions`). `get_gwas()` now adds them too, as
  `trait(<guid>)` already did.
* `upload_gwas()` takes `variable_type` in place of `category`, which still
  works but is deprecated.
* GWAS upload coloc groups gain a `source` column saying whether each row is
  one of the upload's study extractions ("upload") or a study already in the
  map ("existing").
* `ld_proxies()` and `ld_matrix()` gain an `rsq` column. `ld_matrix()` returns
  an empty table, rather than an error, when none of the variants are a lead
  variant and its proxy.
* `genes()` and `traits()` drop empty elements, as `gene()` and `trait()` do,
  and `genes()` no longer returns `tissues`.
* `region()`, `gene()`, `trait()`, `variant()` and their plural forms no longer
  error with `include_associations = TRUE` when the API returns no
  associations.
