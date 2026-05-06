# Summary Statistics Dataframe Documentation

Shared documentation for summary_statistics dataframe columns

## Usage

``` r
summary_statistics_doc()
```

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
