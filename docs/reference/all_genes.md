# All genes

Get all genes from the API

## Usage

``` r
all_genes()
```

## Value

A dataframe containing all genes with the following columns:

- id: the id of the gene

- gene: the name of the gene

- description: the description of the gene

- gene_biotype: the gene biotype

- chr: the chromosome of the gene

- start: the start position of the gene

- stop: the end position of the gene

- strand: the strand of the gene

- source: the source of the gene

- distinct_trait_categories: the number of trait categories that the
  gene is associated with via coloc groups

- distinct_protein_coding_genes: the number of genes that the gene is
  associated with via coloc groups

- num_study_extractions: the number of study extractions for this gene

- num_coloc_groups: the number of coloc groups for this gene

- num_coloc_studies: the number of studies that have coloc results for
  this gene

- num_rare_groups: the number of rare groups for this gene
