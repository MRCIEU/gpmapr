# Gene

A collection of studies that are associated with a particular gene.

## Usage

``` r
gene(
  gene_id,
  include_associations = FALSE,
  include_coloc_pairs = FALSE,
  include_trans = TRUE,
  h4_threshold = 0.8
)
```

## Arguments

- gene_id:

  A numeric value specifying the gene id

- include_associations:

  A logical value specifying whether to include associations (BETA, SE,
  P), defaults to FALSE

- include_coloc_pairs:

  A logical value specifying whether to include coloc pairs, defaults to
  FALSE

- include_trans:

  A logical value specifying whether to include trans genetic effects,
  defaults to TRUE

- h4_threshold:

  A numeric value specifying the h4 threshold for coloc pairs, defaults
  to 0.8

## Value

A list which contains the following elements:

- gene: A list containing metadata about the gene, including region, and
  neighboring genes.

- coloc_groups: a dataframe containing information about which studies
  have coloc results for this gene. See below for details.

- study_extractions: a list of dataframes containing the study
  extractions for this trait. See below for details.

- rare_results: (optional) a list of dataframes containing the rare
  results for this trait

- coloc_pairs: (optional) a dataframe containing all pairwise coloc
  results for this trait.

- variants: a dataframe containing the variants for each associated
  coloc group or rare group.

See below for details.

## Details

The dataframes returned by this function are as follows:

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

## coloc_pairs_dataframe

The coloc_pairs dataframe contains information about which studies have
coloc pairs. It has the following columns:

- study_extraction_a_id: the id of the study extraction associated with
  this coloc pair

- study_extraction_b_id: the id of the study extraction associated with
  this coloc pair

- ld_block_id: the id of the LD block

- h3: the h3 value for this coloc pair

- h4: the h4 value for this coloc pair

- spurious: whether this coloc pair is spurious

## variants_dataframe

The variants dataframe contains variant information that is pulled from
the Variant Effect Predictor (VEP) database. It has the following
columns, along side many more columns from VEP:

- id: the id of the SNP

- gene_id: the id of the gene as predicted by VEP

- gene: the gene name as predicted by VEP
