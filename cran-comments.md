## Update

This is a patch release (1.0.0 -> 1.0.1). gpmapr now requests coloc pairs and full associations in
the API's columnar format, which is smaller and faster for large traits. See NEWS.md.

## R CMD check results

0 errors | 0 warnings | 1 note

* Possibly misspelled words in DESCRIPTION: "colocalization" and "colocalizing" are spelled
  correctly (standard genetics terms, as in the associated paper).

## Test environments

* CRAN incoming checks (R-devel, Windows and Debian Linux), on the previous upload of this
  version: 0 errors | 0 warnings | the misspelling note above. The Windows check also noted that the
  `all_genes()` and `get_all_variant_pleiotropies()` examples took over 10 s, because each
  downloads a whole table from the API. Both examples are now in `\donttest{}`.
* local macOS, R release: 0 errors | 0 warnings | no notes from the package (the other notes were
  an outdated local HTML Tidy, and a 403 from mrcieu.r-universe.dev, which was refusing all requests
  at the time).
* R-hub: linux, macos, macos-arm64, windows (R-devel): 0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

## Notes for the reviewer

* gpmapr is a client for a web API (<https://gpmap.opengwas.io>).
  - Loading the package makes no network request.
  - Functions that reach the API fail with an informative error when it is unavailable.
* Tests that need the API are skipped on CRAN with `skip_on_cran()`. The tests that run on CRAN use
  recorded API responses (tests/testthat/fixtures) and need no network access.
* Examples that call the API are wrapped in `@examplesIf gpmap_available()`, so they are skipped
  when the API can't be reached.
  - `upload_gwas()` and `delete_gwas()` are in `\dontrun{}` because they change data on the server
    and need the user's email address.
* The vignettes are precomputed, so building them needs no network access.
