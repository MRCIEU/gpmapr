# Suzuki et al. (2024) T2D mechanistic clusters

Reference tables used by `vignettes/investigation-suzuki-t2d.Rmd`, extracted from:

> Suzuki K, Hatzikotoulas K, Southam L, Taylor HJ, Yin X, Lorenz KM, Mandla R, et al.
> Genetic drivers of heterogeneity in type 2 diabetes pathophysiology.
> *Nature* 627, 347–357 (2024). <https://doi.org/10.1038/s41586-024-07019-6> (PMC10937372)

The article and its supplementary information are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The files below are
reshaped extracts. Values are unchanged, apart from the GPMap columns added to
`trait_map.csv` and the GRCh38 columns added to `index_snv_clusters.csv`.

| File | Source | Contents |
|---|---|---|
| `clusters.csv` | Supplementary Table 9; main-text Table 1, Fig. 2, Extended Data Fig. 4 | One row per cluster: number of signals, cluster T2D log-OR, defining profile, example loci, insulin secretion/sensitivity direction, open-chromatin enrichment |
| `cluster_trait_profiles.csv` | Supplementary Table 7 | Cluster × trait association of T2D risk alleles (beta, SE, P) for the 37 clustering traits. The beta is the cluster mean of the sample-size-scaled z, `beta / (sqrt(N) * se)` |
| `trait_map.csv` | Supplementary Table 5, plus curated GPMap mapping | The 37 clustering traits, with the GPMap trait row each corresponds to and a `match_type` of `same_gwas`, `proxy`, `unadjusted_proxy` or `none` |
| `index_snv_clusters.csv` | Supplementary Table 6, plus Ensembl | Cluster assignment of the 1,289 index SNVs (locus, chr, rsid, GRCh37 position, cluster, distance to centroid), with each SNV's GRCh38 position (`pos_b38`) and current rsid (`rsid_current`) from Ensembl |

Regenerate with `Rscript data-raw/suzuki2024.R` from the package root. The
script downloads the supplementary workbook from the Europe PMC API and checks
the tables against each other. It also checks the GPMap ids against the
production API.

`pos_b38` is looked up by rsid in the Ensembl REST API (dbSNP mappings, not a
coordinate liftover), because GPMap is on GRCh38. The script also looks each
rsid up on GRCh37 and stops unless that reproduces `pos_b37`. Four rsids have
since been merged into another by dbSNP; `rsid_current` gives the rsid they
were merged into, and equals `rsid` otherwise.
