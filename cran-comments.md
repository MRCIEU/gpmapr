## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Test environments

* local macOS, R (release)
* win-builder (devel)
* R-hub (linux, macos, windows)

## Notes for the reviewer

* gpmapr is a client for a web API (<https://gpmap.opengwas.io>).
  - Loading the package makes no network request.
  - Functions that reach the API fail with an informative error when it is unavailable.
* Tests that need the API are skipped on CRAN with `skip_on_cran()`. The tests that run on CRAN use recorded API responses (tests/testthat/fixtures) and need no network access.
* Examples that call the API are wrapped in `@examplesIf gpmap_available()`.
  - `upload_gwas()` and `delete_gwas()` are in `\dontrun{}` because they change data on the server and need the user's email address.
  - `ld_matrix()` is also in `\dontrun{}`.
* The vignettes are precomputed, so building them needs no network access.
