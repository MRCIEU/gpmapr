# Coloc Groups Dataframe Documentation

Shared documentation for coloc_groups dataframe columns

## Usage

``` r
coloc_groups_doc()
```

## coloc_groups_dataframe

The coloc_groups dataframe contains information about which studies have
coloc results. It has the following columns:

- coloc_group_id: the unique id for this group of colocalised results

- study_id: the id of the study

- study_extraction_id: the id of the study extraction

- variant_id: the id of the SNP

- ld_block_id: the id of the LD block

- chr: the chromosome of the SNP

- bp: the base pair position of the SNP

- min_p: the minimum p-value related to the study_extraction_id

- cis_trans: the cis/trans status of the SNP

- ld_block: the LD block of the SNP

- display_snp: the display SNP name

- gene: the gene associated with the SNP

- gene_id: the id of the gene

- trait_id: the id of the trait

- trait_name: the name of the trait

- trait_category: the category of the trait

- data_type: the data type of the trait

- tissue: the tissue of the trait
