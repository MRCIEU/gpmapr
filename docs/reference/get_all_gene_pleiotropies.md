# Get All Gene Pleiotropies

Get gene pleiotropy from the API by gene id

## Usage

``` r
get_all_gene_pleiotropies()
```

## Value

A list containing the gene pleiotropy

- gene_id: the id of the gene

- gene: the name of the gene

- distinct_trait_categories: the number of trait categories that the
  gene is associated with via coloc groups

- distinct_protein_coding_genes: the number of genes that the gene is
  associated with via coloc groups
