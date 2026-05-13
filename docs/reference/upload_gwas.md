# Upload a GWAS to the API

Upload a GWAS to the API

## Usage

``` r
upload_gwas(
  file,
  name,
  p_value_threshold = 5e-08,
  column_names = list(),
  email = NA,
  category = "continuous",
  is_published = FALSE,
  doi = NA,
  should_be_added = FALSE,
  ancestry = "EUR",
  sample_size = NA,
  reference_build = "GRCh38",
  compare_with_upload_guids = NA
)
```

## Arguments

- file:

  The path to the GWAS file, maximum size is 1GB

- name:

  The name of the GWAS

- p_value_threshold:

  The p-value threshold for the GWAS

- column_names:

  A list of column names in the format of: list(CHR = "chr", BP =
  "pos"...)

  - CHR: chromosome

  - BP: base pair position

  - P: p-value

  - EA: allele 1

  - OA: allele 2

  - EAF: allele frequency And either BETA and SE, or OR, LB, and UB

  - BETA: beta

  - SE: standard error

  - OR: odds ratio

  - LB: lower bound of the confidence interval

  - UB: upper bound of the confidence interval

- email:

  The email of the user

- category:

  The category of the GWAS. Only "continuous" and "categorical" are
  accepted.

- is_published:

  Whether the GWAS is published

- doi:

  The DOI of the GWAS

- should_be_added:

  Whether the GWAS should be added to the API

- ancestry:

  The ancestry of the GWAS. Currently only "EUR" is accepted.

- sample_size:

  The sample size of the GWAS

- reference_build:

  The reference build of the GWAS. Only "GRCh37" and "GRCh38" are
  accepted.

- compare_with_upload_guids:

  A vector of GUIDs of uploads to compare with

## Value

A list containing the GWAS information
