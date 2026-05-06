# Associations Dataframe Documentation

Shared documentation for associations dataframe columns

## Usage

``` r
associations_doc()
```

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
