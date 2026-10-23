## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* Possibly misspelled words in DESCRIPTION: "colocalisation", "colocalising" and "pleiotropy" are
  spelled correctly (British spelling of standard genetics terms).

## Test environments

* local macOS, R release: 0 errors | 0 warnings | 1 note ("unable to verify current time", from the
  local network)
* win-builder, R-devel (2026-10-22 r90748 ucrt), x86_64-w64-mingw32: 0 errors | 0 warnings | 1 note
  (the note above)
* R-hub: linux, macos, macos-arm64, windows (R-devel)

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
