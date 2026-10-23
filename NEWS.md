# gpmapr 1.0.0

First CRAN release. Requires version 1.1.0 of the GPMap API.

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
* Breaking: GWAS upload results (`get_gwas()`, and `trait()` / `traits()` with
  a GUID) have the same tables and column names as other results.
  * One `study_extractions` table holds the upload's study extractions and the
    existing ones it colocalises with; `upload_study_extractions` is gone.
  * Rows of `study_extractions` and `coloc_groups` have a `source` column
    ("upload" or "existing") and a `gwas_upload_id` (the upload the row comes
    from; `NA` for existing rows). Upload and existing study extractions have
    separate ids, so identify one by `source` and its id together. The
    upload's own rows are `source == "upload" & gwas_upload_id == result$trait$id`.
  * Coloc groups use `study_extraction_id` for both sides, and coloc pairs use
    `study_extraction_a_id` / `study_extraction_b_id` with `source_a` /
    `source_b`, replacing the `upload_*` / `existing_*` id columns.
  * Upload rows are named after their upload (`trait_name`), have `data_type`
    "Phenotype" and no `trait_id` or `study_id`.
* New `compare_upload()` lines up a GWAS upload with another trait (or a
  compared upload): each one's effect in every coloc group they share, and the
  coloc pairs between them.
* Coloc pairs from `gene()`, `genes()`, `variant()` and `get_gwas()` include the
  trait name, data type, tissue and gene of both sides (`trait_name_a`,
  `trait_name_b`, ...), so they no longer need joining to `study_extractions`.
* Breaking: `trait` is now `pipeline_id` (the GPMap pipeline's id for the
  trait, such as its OpenGWAS dataset id) and `category` is now
  `variable_type` (continuous or categorical), in `all_traits()`, `trait()`,
  `traits()` and GWAS uploads.
* Breaking: `$associations` is no longer returned. With
  `include_associations = TRUE`, each study's association at its variant
  (`beta`, `se`, `p`, `eaf`, `imputed`) is added to `coloc_groups`,
  `rare_results` and `study_extractions`. `get_gwas()` now adds them too, as
  `trait(<guid>)` already did, and so does `traits()` for GUIDs.
* Associations of variants a study did not report are placeholders in the
  GPMap (beta 0, se 1, p 1); they are now returned with `beta`, `se`, `p` and
  `eaf` set to `NA`.
* `upload_gwas()` takes `variable_type` in place of `category`, which still
  works but is deprecated.
* `traits()` no longer errors when given a GWAS upload GUID alongside trait ids.
* `ld_proxies()` and `ld_matrix()` gain an `rsq` column. `ld_matrix()` returns
  an empty table when none of the variants are a lead variant and its proxy.
* `genes()` and `traits()` drop empty elements, as `gene()` and `trait()` do,
  and `genes()` no longer returns `tissues`.
* `region()`, `gene()`, `trait()`, `variant()` and their plural forms no longer
  error with `include_associations = TRUE` when the API returns no
  associations.
