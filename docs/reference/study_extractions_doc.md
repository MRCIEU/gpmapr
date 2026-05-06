# Study Extractions Dataframe Documentation

Shared documentation for study_extractions dataframe columns

## Usage

``` r
study_extractions_doc()
```

## study_extractions_dataframe

The study_extractions dataframe contains information about which studies
have coloc results. It has the following columns:

- id: the unique id for this study extraction

- study_id: the id of the study associated with this study extraction

- variant_id: the id of the SNP

- snp: the SNP name

- ld_block_id: the id of the LD block

- unique_study_id: the unique id for this study

- study: the study name

- file: the file name

- svg_file: the SVG file name

- file_with_lbfs: the file name with lbfs

- chr: the chromosome of the SNP

- bp: the base pair position of the SNP

- min_p: the minimum p-value related to the study_extraction_id

- cis_trans: the cis/trans status of the SNP

- ld_block: the LD block of the SNP

- gene: the gene associated with the SNP

- gene_id: the id of the gene

- trait_id: the id of the trait

- trait_name: the name of the trait

- trait_category: the category of the trait

- data_type: the data type of the trait

- tissue: the tissue of the trait
