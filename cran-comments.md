## Update

This is a patch release (1.0.0 -> 1.0.1). gpmapr now requests coloc pairs and full associations in
the API's columnar format, which is smaller and faster for large traits. See NEWS.md.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local macOS, R release: 0 errors | 0 warnings | 0 notes (with the network blocked, 1 note:
  "unable to verify current time")
* win-builder, R-devel: TODO
* R-hub: linux, macos, macos-arm64, windows (R-devel): TODO

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
