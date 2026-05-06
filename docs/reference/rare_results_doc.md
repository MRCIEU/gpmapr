# Rare Results Dataframe Documentation

Shared documentation for rare_results dataframe columns

## Usage

``` r
rare_results_doc()
```

## rare_results_dataframe

The rare_results dataframe contains information about which studies have
coloc results. It has the following columns:

- rare_result_group_id: the unique id for this rare result group

- study_id: the id of the study associated with this rare result

- study_extraction_id: the id of the study extraction associated with
  this rare result

- variant_id: the id of the SNP

- ld_block_id: the id of the LD block

- chr: the chromosome of the SNP

- bp: the base pair position of the SNP

- min_p: the minimum p-value related to the study_extraction_id

- display_snp: the display SNP name

- gene: the gene associated with the SNP

- gene_id: the id of the gene

- trait_id: the id of the trait

- trait_name: the name of the trait

- trait_category: the category of the trait

- data_type: the data type of the trait

- tissue: the tissue of the trait

- ld_block: the LD block of the SNP
