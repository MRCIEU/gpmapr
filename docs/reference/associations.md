# Get Associations by SNP ID and Study ID

Get associations from the API by SNP id and study id

## Usage

``` r
associations(variant_ids, study_ids)
```

## Arguments

- variant_ids:

  A vector of numeric values specifying the SNP IDs

- study_ids:

  A vector of numeric values specifying the Study IDs

## Value

A dataframe containing the associations

## associations_dataframe

The associations dataframe contains information about which studies have
association results. It has the following columns:

- variant_id: the id of the SNP associated with this association

- study_id: the id of the study associated with this association

- beta: the beta value of the association

- se: the standard error of the association

- p: the p-value of the association

- eaf: the estimated allele frequency of the association

- imputed: whether the association is imputed
