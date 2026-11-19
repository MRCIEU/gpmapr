# Releasing gpmapr

Maintainer notes; not shipped in the package (see `.Rbuildignore`).

## How the API version gate works

`api_to_package_version` in `R/common.R` maps each GPMap API version to the
gpmapr versions that support it. On the first request of a session,
`ensure_api_compatible()` (in `R/api.r`) fetches `/v1/info/version`. If the
installed gpmapr isn't listed for that version, it stops and tells the user to
update. So **every API version change blocks every released gpmapr that isn't
listed for it**.

## Before deploying an API version change

1. Add the new API version to `api_to_package_version`, mapped to every
   gpmapr version that works with it. For a non-breaking change, that is all
   the versions listed for the previous API version.
2. Release that gpmapr to CRAN (below), and wait for acceptance.
3. Only then deploy the API. CRAN versions released before step 2 stop working
   if the new API version isn't listed for them.

Prefer a new path (e.g. `/v2/`) for breaking changes, so older clients keep
working on `/v1/`.

## Breaking API changes

Once gpmapr is on CRAN, a breaking API change (renamed or restructured fields)
should not stop released versions before a fixed one is accepted. Either:

- make the new gpmapr read both the old and the new responses (map old field
  names to new ones as the responses are parsed; for fields gpmapr sends, send
  both names, as the API ignores fields it doesn't know), list it under both API
  versions, release it, then deploy; or
- serve the change under a new path (e.g. `/v2/`) so older clients keep working.

1.0.0, the first CRAN release, did neither: it requires API 1.1.0, and the
version gate stops every earlier gpmapr. With no CRAN users yet, the order was
API 1.1.0 deployed first, then 1.0.0 released.

## Testing against a local API

The live tests and the fixture recorder call production unless `GPMAP_URL` is
set, e.g. `GPMAP_URL=http://localhost:8000`. Use it to test against an API
version that isn't deployed yet (gpmapr refuses an API version it isn't listed
for, so live tests against production fail until the matching API is deployed).

## Each release

1. `usethis::use_release_issue()` for the checklist.
2. `usethis::use_version()` (major / minor / patch).
3. Update `NEWS.md`.
4. Add the new package version to `api_to_package_version` for every API
   version it supports.
5. Re-knit the precomputed vignettes and review the diffs:
   `cd vignettes && Rscript precompile.R`.
   Edit the `.Rmd.orig` files, never the generated `.Rmd` files.
6. Re-record the mocked-test fixtures, then run both test groups against
   production:
   ```sh
   Rscript data-raw/record_api_fixtures.R
   NOT_CRAN=true Rscript -e 'devtools::test()'
   ```
7. `devtools::document()`.
8. `devtools::check(remote = TRUE, manual = TRUE)`, and the same check with
   the network blocked (e.g. `https_proxy=http://127.0.0.1:9`). Expect only
   notes that are about the machine ("unable to verify current time", an old
   HTML Tidy) and the "possibly misspelled words" note for "colocalization" and
   "colocalizing": no dictionary knows them, so explain it in `cran-comments.md`.
   CRAN's own check also flags examples that take over 10 s (Windows is
   slower), so put any example that downloads a whole table in `\donttest{}`.
   Explain every acronym in the DESCRIPTION text (API, QTL, GWAS, ...).
9. `urlchecker::url_check()`.
10. `devtools::check_win_devel()` and `rhub::rhub_check()`.
11. Update `cran-comments.md` with the environments and results.
12. Merge the pull request into `main` (these repos squash-merge), then
    `git checkout main && git pull`.
13. From `main`, `devtools::submit_cran()`, then confirm the submission from the
    maintainer's email. Submitting from `main` means the commit written to
    `CRAN-SUBMISSION` is on `main`, so the release tag lands on it. (`CRAN-SUBMISSION`
    is not tracked, so it doesn't exist in other checkouts or worktrees.)
14. Answer reviewer comments. Fix on a branch, merge to `main` and resubmit from
    `main`, with an incremented version if asked.

## After acceptance

Don't tag before CRAN accepts, so the tag always matches what is on CRAN;
nothing tags automatically.

1. `git checkout main && git pull`, so the release is made from `main`.
   Don't tag from a feature branch: after a squash merge, the branch commits
   aren't on `main`, so the tag would point at a commit outside its history.
2. `usethis::use_github_release()`. This tags `v<version>` on the commit recorded
   in `CRAN-SUBMISSION` (or the current HEAD, if there is no such file) and creates
   the GitHub release from `NEWS.md`. Answer "Definitely" only when the branch is
   `main`; it asks otherwise.
3. Optional: `usethis::use_dev_version()` (e.g. 1.0.1 -> 1.0.1.9000), which lets
   a development install be told apart from the release. If you do, add the dev
   version to `api_to_package_version`, or development installs refuse every
   request. Skip both to stay on the released version number.
4. After the first release, add the CRAN badge to the README:
   `usethis::use_cran_badge()`.

## Tests

- **Mocked tests** (`test_mocked_api.R` and others) answer every request
  from `tests/testthat/fixtures`, so they run on CRAN.
- **Live tests** (`test_live_api.R`) call production. They're skipped on CRAN
  and when the API is unreachable.

When a live test fails, the API has changed: update the code, then re-record
the fixtures so the mocked tests match.

## How users upgrade

```r
install.packages("gpmapr")   # or update.packages()
pak::pak("MRCIEU/gpmapr")    # development version
```
