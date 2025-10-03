.search_cache <- new.env()

#' @title Search the Genotype-Phenotype Map
#' @description Search the GP Map for Traits, Genes or Variants
#' @param search_text A character string specifying the search text
#' @details After calling search, you can use call the subsequent data as described in the `call` column of the
#' search results.
#' @return A dataframe containing the search results with the following columns:
#' \itemize{
#'   \item type: the type of the search result: "original_variant", "proxy_variant", "trait", "gene"
#'   \item name: the name of the search result
#'   \item type_id: the type_id of the search result.  This is the internal id in which the data can be accessed.
#'   \item call: the call to get the search result: "variant(type_id)", "trait(type_id)", "gene(type_id)"
#'   \item info: a string containing informaiton about the search result, which may include:
#'     \itemize{
#'       \item Extractions: the number of extractions
#'       \item Colocalisation Groups: the number of colocalisation groups
#'       \item Colocalisation Studies: the number of colocalisation studies
#'       \item Rare Results: the number of rare results
#'       \item Rsquared: the rsquared of the proxy variant compared to the original variant
#'     }
#' }
#' @export
search_gpmap <- function(search_text) {
  if (is.null(search_text) || is.na(search_text) || search_text == "") {
    stop("search_text is required")
  }

  #get search results from cache
  cache_key <- "search_results"
  if (exists(cache_key, envir = .search_cache)) {
    search_results <- get(cache_key, envir = .search_cache)
  } else {
    search_results <- search_options_api()
    search_results <- search_results$search_terms
    search_results$name_lower <- tolower(search_results$name)
    search_results$alt_name_lower <- tolower(search_results$alt_name)
    assign(cache_key, search_results, envir = .search_cache)
  }

  search_text <- tolower(search_text)

  # if the search text is a variant, get the variant results
  if (grepl("^rs", search_text) || grepl("^\\d+:\\d+", search_text)) {
    variant_results <- search_variants_api(search_text)

    if (length(variant_results$original_variants) > 0) {
      original_variant <- variant_results$original_variants |>
        dplyr::rename(type_id = id, name = display_snp) |>
        dplyr::mutate(
          type = "original_variant",
          rsq = NA,
          num_coloc_groups = num_colocs,
          num_coloc_studies = ifelse(is.data.frame(coloc_groups[[1]]), nrow(coloc_groups[[1]]), 0),
          num_rare_results = ifelse(is.data.frame(rare_results[[1]]), nrow(rare_results[[1]]), 0),
          num_study_extractions = NA) |>
        dplyr::select(type,
          name,
          type_id,
          rsq,
          num_coloc_groups,
          num_coloc_studies,
          num_rare_results,
          num_study_extractions
        )
    } else {
      original_variant <- data.frame()
    }

    if (length(variant_results$proxy_variants) > 0) {
      proxy_variants <- variant_results$proxy_variants |>
        dplyr::rename(type_id = id, name = display_snp) |>
        dplyr::mutate(
          type = "proxy_variant",
          rsq = round(ld_proxies[[1]]$r^2, 2),
          num_coloc_groups = num_colocs,
          num_coloc_studies = ifelse(is.data.frame(coloc_groups[[1]]), nrow(coloc_groups[[1]]), 0),
          num_rare_results = ifelse(is.data.frame(rare_results[[1]]), nrow(rare_results[[1]]), 0),
          num_study_extractions = NA
        ) |>
        dplyr::select(type,
          name,
          type_id,
          rsq,
          num_coloc_groups,
          num_coloc_studies,
          num_rare_results,
          num_study_extractions
        )
    } else {
      proxy_variants <- data.frame()
    }

    search_results <- dplyr::bind_rows(original_variant, proxy_variants) |>
      dplyr::distinct()
    if (nrow(search_results) == 0) return(search_results)
  } else {
    # if the search text is a trait, get the trait results
    search_results <- search_results[
      grepl(search_text, search_results$name_lower) | grepl(search_text, search_results$alt_name_lower),
    ]
    search_results$name <- ifelse(
      is.na(search_results$alt_name), search_results$name,
      paste0(search_results$name, " (", search_results$alt_name, ")")
    )
  }

  if (nrow(search_results) > 0) {
    search_results$call <- dplyr::case_when(
      search_results$type == "original_variant" ~ paste0("variant(", search_results$type_id, ")"),
      search_results$type == "proxy_variant" ~ paste0("variant(", search_results$type_id, ")"),
      search_results$type == "trait" ~ paste0("trait(", search_results$type_id, ")"),
      search_results$type == "gene" ~ paste0("gene('", search_results$type_id, "')"),
      TRUE ~ NA_character_
    )

    search_results$info <- paste0(
      ifelse(search_results$num_extractions > 0, paste0("Extractions: ", search_results$num_extractions, ", "), ""),
      "Colocalisation Groups: ", search_results$num_coloc_groups, ", ",
      "Colocalising Traits: ", search_results$num_coloc_studies, ", ",
      "Rare Results: ", search_results$num_rare_results,
      ifelse(search_results$type == "proxy_variant", paste0(", Rsquared: ", search_results$rsq), "")
    )
    search_results$rsq <- NULL

    search_results <- dplyr::select(search_results, call, name, type, type_id, info)
  }
  return(search_results)
}