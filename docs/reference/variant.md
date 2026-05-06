# Variant

A collection of studies that are associated with a particular variant.

## Usage

``` r
variant(
  variant_id,
  include_coloc_pairs = FALSE,
  h4_threshold = 0.8,
  include_summary_stats = FALSE
)
```

## Arguments

- variant_id:

  A character string specifying the SNP ID

- include_coloc_pairs:

  A logical value specifying whether to include coloc pairs

- h4_threshold:

  A numeric value specifying the cutoff for included coloc pairs,
  defaults to 0.8. Only used if include_coloc_pairs is TRUE.

- include_summary_stats:

  A logical value specifying whether to include summary stats

## Value

A list which contains the following elements:

- variant: named list containing the variant information

- coloc_groups: a dataframe containing information about which studies
  have coloc results for this variant

- rare_results: a list of dataframes containing the rare variants

- study_extractions: a list of dataframes containing the study
  extractions

- summary_stats (optional): a list of dataframes containing the summary
  stats for each study, where the name of each element is the study_id.
  Column names are uppercase (e.g. SNP, BP, BETA, SE, LBF_1).

- coloc_pairs (optional): a dataframe containing information about which
  studies have coloc pairs for this variant where the
  study_extraction_a_id and study_extraction_b_id are the
  study_extraction_ids of the two studies. h4_threshold is the cutoff
  for included coloc pairs, defaults to 0.8

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

## summary_statistics_dataframe

The summary_statistics dataframe contains information about which
studies have summary statistics. From the API, column names are
typically uppercase (SNP, CHR, BP, EA, OA, EAF, Z, BETA, SE, P, LBF_1,
etc.). It has the following columns (names may be upper or lower case
depending on source):

- SNP / variant_id: the id of the SNP

- CHR / chr: the chromosome of the SNP

- BP / bp: the base pair position of the SNP

- EA / ea: the effect allele

- OA / oa: the other allele

- EAF / eaf: the estimated allele frequency

- Z / z: the z-score

- BETA / beta: the beta value

- SE / se: the standard error

- P / p: the p-value

- imputed: whether the summary statistics are imputed

- LBF\_\* / lbf\_\*: all different finemapped log-bayes factors for each
  credible set. Each credible set is numbered from 1 to 10. If
  finemapped failed or only returned 1 credible set, the LBF_1 column is
  just converted directly from the z-score.

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
