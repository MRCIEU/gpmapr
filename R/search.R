.search_cache <- new.env()

variant_result_columns <- c(
  "type", "name", "type_id", "rsq", "num_coloc_groups", "num_coloc_studies",
  "num_rare_results", "num_study_extractions"
)

#' @title Search the Genotype-Phenotype Map
#' @description Search the GP Map for Traits, Genes or Variants
#' @param search_text A character string specifying the search text
#' @param rsquared_threshold A numeric value specifying the rsquared threshold for proxy variants, defaults to 0.8
#' @details After calling search, you can use call the subsequent data as described in the `call` column of the
#' search results.
#' @return A dataframe containing the search results, ordered by the number of colocalisation groups plus rare
#' results (most first).  For variant searches the searched variant is always first.  Columns include:
#' \itemize{
#'   \item type: the type of the search result: "original_variant", "proxy_variant", "trait", "gene"
#'   \item name: the name of the search result
#'   \item type_id: the type_id of the search result.  This is the internal id in which the data can be accessed.
#'   \item rsq: (variant searches only) the r-squared between a proxy variant and the searched variant; `NA` for
#'     the searched variant itself
#'   \item num_study_extractions: the number of study extractions (`NA` for variants)
#'   \item num_coloc_groups: the number of colocalisation groups
#'   \item num_coloc_studies: the number of colocalising studies
#'   \item num_rare_results: the number of rare results
#'   \item call: the call to get the search result: "variant(type_id)", "trait(type_id)", "gene(type_id)"
#' }
#' Trait and gene results also carry the descriptive columns of the search options (e.g. `aliases`,
#' `sample_size`, `ancestry`).
#' @examplesIf gpmap_available()
#' search_gpmap("haemoglobin")
#' search_gpmap("rs9939609")
#' @export
search_gpmap <- function(search_text, rsquared_threshold = 0.8) {
  if (is.null(search_text) || is.na(search_text) || search_text == "") {
    stop("search_text is required")
  }

  search_text <- tolower(search_text)

  # if the search text is a variant, get the variant results
  if (grepl("^rs", search_text) || grepl("^\\d+:\\d+", search_text)) {
    variant_results <- search_variants_api(search_text, rsquared_threshold)

    if (length(variant_results$original_variants) > 0) {
      original_variant <- variant_results$original_variants |>
        dplyr::rename(type_id = "id", name = "display_snp") |>
        dplyr::mutate(
          type = "original_variant",
          rsq = NA_real_,
          num_coloc_groups = .data$num_colocs,
          num_coloc_studies = .count_rows(.data$coloc_groups),
          num_rare_results = .count_rows(.data$rare_results),
          num_study_extractions = NA) |>
        dplyr::select(dplyr::all_of(variant_result_columns))
    } else {
      original_variant <- data.frame()
    }

    if (length(variant_results$proxy_variants) > 0) {
      proxy_variants <- variant_results$proxy_variants |>
        dplyr::rename(type_id = "id", name = "display_snp") |>
        dplyr::mutate(
          type = "proxy_variant",
          rsq = .proxy_rsq(.data$ld_proxies),
          num_coloc_groups = .data$num_colocs,
          num_coloc_studies = .count_rows(.data$coloc_groups),
          num_rare_results = .count_rows(.data$rare_results),
          num_study_extractions = NA
        ) |>
        dplyr::select(dplyr::all_of(variant_result_columns))
    } else {
      proxy_variants <- data.frame()
    }

    search_results <- dplyr::bind_rows(original_variant, proxy_variants) |>
      dplyr::distinct()
    if (nrow(search_results) == 0) return(search_results)
  } else {
    # if the search text is a trait, get the trait results
    search_results <- .search_options()
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

    columns_to_remove <- c("name_lower", "alt_name", "alt_name_lower")
    search_results <- dplyr::select(search_results, -dplyr::any_of(columns_to_remove))

    # most important results first, keeping an exact variant match above its proxies
    num_coloc_groups <- search_results$num_coloc_groups
    num_rare_results <- search_results$num_rare_results
    importance <- ifelse(is.na(num_coloc_groups), 0, num_coloc_groups) +
      ifelse(is.na(num_rare_results), 0, num_rare_results)
    is_original_variant <- search_results$type == "original_variant"
    search_results <- search_results[order(!is_original_variant, -importance, method = "radix"), ]
    rownames(search_results) <- NULL
  }
  return(search_results)
}


# The full list of searchable traits and genes, downloaded once per session.
.search_options <- function() {
  cache_key <- "search_results"
  if (exists(cache_key, envir = .search_cache)) {
    return(get(cache_key, envir = .search_cache))
  }
  search_results <- search_options_api()$search_terms
  search_results$name_lower <- tolower(search_results$name)
  search_results$alt_name_lower <- tolower(search_results$alt_name)
  assign(cache_key, search_results, envir = .search_cache)
  return(search_results)
}


# Rows in each element of a list column of data frames (0 for a missing one).
.count_rows <- function(tables) {
  return(vapply(tables, function(x) if (is.data.frame(x)) nrow(x) else 0L, integer(1)))
}


# r-squared of each proxy with the searched variant, from its ld_proxies table
# (the strongest pair when there are several; NA when there is none).
.proxy_rsq <- function(ld_proxies) {
  return(vapply(ld_proxies, function(x) {
    if (!is.data.frame(x) || nrow(x) == 0 || !"r" %in% names(x)) {
      return(NA_real_)
    }
    return(round(max(x$r^2, na.rm = TRUE), 2))
  }, numeric(1)))
}
