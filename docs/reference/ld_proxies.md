# LD Proxies

Get LD proxies from the API by Variant ID

## Usage

``` r
ld_proxies(variant_ids = c())
```

## Arguments

- variant_ids:

  A character string specifying the Variant ID. Variant IDs can be SNP
  IDs or variant IDs.

## Value

A list containing the LD proxies

## ld_dataframe

The ld dataframe contains information about the LD matrix. It has the
following columns:

- lead_variant_id: the id of the lead SNP

- proxy_variant_id: the id of the variant SNP

- ld_block_id: the id of the LD block

- r: the r value between the lead and variant SNPs
