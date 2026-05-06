# SNP Pleiotropies Dataframe Documentation

Shared documentation for variant_pleiotropies dataframe columns

## Usage

``` r
variant_pleiotropies_doc()
```

## variant_pleiotropies_dataframe

The variant_pleiotropies dataframe contains information about which SNPs
are pleiotropic. It has the following columns:

- variant_id: the id of the SNP

- snp: the name of the SNP

- : distinct_trait_categories the number of trait categories that the
  SNP is associated with via coloc groups

- : distinct_protein_coding_genes the number of genes that the SNP is
  associated with via coloc groups
