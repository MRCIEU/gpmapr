# Search the Genotype-Phenotype Map

Search the GP Map for Traits, Genes or Variants

## Usage

``` r
search_gpmap(search_text, rsquared_threshold = 0.8)
```

## Arguments

- search_text:

  A character string specifying the search text

- rsquared_threshold:

  A numeric value specifying the rsquared threshold for proxy variants,
  defaults to 0.8

## Value

A dataframe containing the search results with the following columns:

- type: the type of the search result: "original_variant",
  "proxy_variant", "trait", "gene"

- name: the name of the search result

- type_id: the type_id of the search result. This is the internal id in
  which the data can be accessed.

- call: the call to get the search result: "variant(type_id)",
  "trait(type_id)", "gene(type_id)"

- info: a string containing informaiton about the search result, which
  may include:

  - Extractions: the number of extractions

  - Colocalisation Groups: the number of colocalisation groups

  - Colocalisation Studies: the number of colocalisation studies

  - Rare Results: the number of rare results

  - Rsquared: the rsquared of the proxy variant compared to the original
    variant

## Details

After calling search, you can use call the subsequent data as described
in the `call` column of the search results.
