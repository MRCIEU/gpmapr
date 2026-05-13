# All traits

Get all traits from the API

## Usage

``` r
all_traits()
```

## Value

A dataframe containing all traits with the following columns:

- id: the id of the trait

- data_type: the data type of the trait

- trait: the internal string id of the trait

- trait_name: the name of the trait

- trait_category: the trait category of the trait

- variant_type: the type of variant

- sample_size: the sample size of the trait

- category: the category of the trait (continuous, categorical)

- ancestry: the ancestry of the trait

- heritability: the LDSC heritability score of the trait

- heritability_se: the standard error of the LDSC heritability score of
  the trait

- num_study_extractions: the number of study extractions for this trait

- num_coloc_groups: the number of coloc groups for this trait

- num_coloc_studies: the number of studies that have coloc results for
  this trait

- num_rare_results: the number of rare results for this trait
