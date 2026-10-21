# gpmapr 0.1.0

First CRAN release.

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
* `region()`, `gene()`, `trait()`, `variant()` and their plural forms no longer
  error with `include_associations = TRUE` when the API returns no
  associations.
