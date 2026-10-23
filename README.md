
# gpmapr

<!-- badges: start -->
[![R-CMD-check](https://github.com/MRCIEU/gpmapr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MRCIEU/gpmapr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

gpmapr is an R client for the [OpenGWAS Genotype-Phenotype Map](https://gpmap.opengwas.io).
You can:
- search the map for traits, genes and variants;
- fetch their colocalisation groups, rare-variant results and associations;
- upload your own GWAS to compare against the map.

## Installation


Install the package from GitHub:

``` r
# install.packages("pak")
pak::pak("MRCIEU/gpmapr")
```

## Usage

``` r
library(gpmapr)

# Find a trait, gene or variant
search_gpmap("haemoglobin")

# Then fetch it with the call shown in the `call` column
result <- trait(5020)
head(result$coloc_groups)
```

See `vignette("gpmapr")` to get started.

## Updating

gpmapr checks on the first request of each session that it supports the
version of the GPMap API. If the API has changed, you'll see an error asking
you to update. Update with:

``` r
pak::pak("MRCIEU/gpmapr")
```

then restart R.
