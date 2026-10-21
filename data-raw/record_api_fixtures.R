# Record the API responses used by the mocked tests (tests/testthat/test_mocked_api.R).
#
# Each exported function below is called once against the production API while
# api_get() / api_post() are wrapped to save every raw response. Responses are
# trimmed (every array of records cut to its first few, all fields kept) so
# the fixtures stay small, then written to tests/testthat/fixtures/ with an
# index mapping each request URL (relative to the API root) to its file.
#
# Re-run from the package root whenever the live tests (test_live_api.R) show
# that the API has changed:
#   Rscript data-raw/record_api_fixtures.R

pkgload::load_all(".", quiet = TRUE)
options(gpmap_url = "https://gpmap.opengwas.io/api")
ns <- asNamespace("gpmapr")

fixture_dir <- file.path("tests", "testthat", "fixtures")
unlink(fixture_dir, recursive = TRUE)
dir.create(fixture_dir, recursive = TRUE)

max_elements <- 3

# Cut every JSON array of records (objects or arrays) to its first
# `max_elements` elements, recursively. Arrays of plain values, such as column
# names or the values of one row, are kept whole so every field survives.
trim <- function(x, keep = max_elements) {
  if (!is.list(x)) {
    return(x)
  }
  is_record_array <- is.null(names(x)) && length(x) > 0 && all(vapply(x, is.list, logical(1)))
  if (is_record_array && length(x) > keep) {
    x <- x[seq_len(keep)]
  }
  return(lapply(x, trim, keep = keep))
}

recorded <- list()
recording <- function(verb) {
  return(function(url, ..., check_version = TRUE) {
    response <- ns$api_request(verb, url, ..., check_version = check_version)
    relative_url <- sub(getOption("gpmap_url"), "", url, fixed = TRUE)
    recorded[[relative_url]] <<- httr::content(response, "text", encoding = "UTF-8")
    return(response)
  })
}

calls <- list(
  quote(health_api()),
  quote(all_traits()),
  quote(trait_duplicates()),
  quote(trait(5020)),
  quote(trait(5020, include_associations = TRUE)),
  quote(trait(5020, include_full_associations = TRUE)),
  quote(traits(c(4405, 4872), include_associations = TRUE)),
  quote(all_genes()),
  quote(gene("WNT7B")),
  quote(genes(c("WNT7B", "WNT7A"), include_associations = TRUE)),
  quote(variant(5553693)),
  quote(variants(c(5553693, 5553694), include_associations = TRUE, expand = TRUE)),
  quote(region(459)),
  quote(region(459, include_associations = TRUE)),
  quote(associations(80717, c(926, 927))),
  quote(ld_proxies(c(5553693, 5553694))),
  quote(get_all_gene_pleiotropies()),
  quote(get_all_variant_pleiotropies()),
  quote(get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155")),
  quote(search_gpmap("rs9939609")),
  quote(search_gpmap("haemoglobin"))
)

assignInNamespace("api_get", recording(httr::GET), ns = "gpmapr")
assignInNamespace("api_post", recording(httr::POST), ns = "gpmapr")
for (call in calls) {
  message("Recording ", deparse(call))
  eval(call, envir = ns)
}
pkgload::load_all(".", quiet = TRUE)

index <- list()
for (i in seq_along(recorded)) {
  relative_url <- names(recorded)[i]
  body <- jsonlite::fromJSON(recorded[[i]], simplifyVector = FALSE)
  if (relative_url == "/v1/search/options") {
    # Thousands of terms: keep those the mocked text search matches, rather
    # than the first few.
    matches <- vapply(body$search_terms, function(term) {
      return(grepl("haemoglobin", tolower(term$name)))
    }, logical(1))
    body$search_terms <- body$search_terms[matches]
  } else {
    body <- trim(body)
  }
  file <- sprintf("response-%02d.json", i)
  writeLines(jsonlite::toJSON(body, auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = TRUE),
             file.path(fixture_dir, file))
  index[[relative_url]] <- file
}
writeLines(jsonlite::toJSON(index, auto_unbox = TRUE, pretty = TRUE), file.path(fixture_dir, "index.json"))
message("Wrote ", length(index), " fixtures to ", fixture_dir)
