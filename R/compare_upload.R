#' @title Compare a GWAS upload with a trait
#' @description Line up a GWAS upload with another trait (or another uploaded GWAS it was compared with) at
#'   the places they colocalise.
#'
#' For each coloc group containing both the upload and the trait, it gives the effect of each at the group's
#' variant.  Both effects are read from the upload's own result, where they are on the same alleles, so
#' their signs can be compared directly.
#' @param upload A GWAS upload result from [get_gwas()] or [trait()] with a GUID, fetched with
#'   `include_associations = TRUE`
#' @param trait_id The id of a trait in the GPMap, or the GUID of another GWAS upload that `upload` was
#'   compared with (`compare_with_upload_guids` in [upload_gwas()])
#' @return A list with two dataframes:
#' \itemize{
#'   \item coloc_groups: one row per coloc group containing both the upload and the trait, with
#'     `coloc_group_id`, `ld_block_id`, `ld_block`, `variant_id`, `display_snp`, `h4_connectedness`, the
#'     upload's lead study extraction in the group (`upload_study_extraction_id`) and its `upload_beta`,
#'     `upload_se` and `upload_p`, and the same for the trait (`trait_name`, `trait_study_extraction_id`,
#'     `trait_beta`, `trait_se`, `trait_p`).  The lead is the study extraction with the smallest p-value.
#'   \item coloc_pairs: the coloc pairs between one of the upload's study extractions and one of the
#'     trait's, with `ld_block_id`, `upload_study_extraction_id`, `upload_display_snp`, `trait_name`,
#'     `trait_study_extraction_id`, `trait_display_snp`, `h3`, `h4` and `in_shared_coloc_group` (whether
#'     the pair's LD block has a coloc group in `coloc_groups`).  Pairs have no effects, as the two study
#'     extractions are at different variants.
#' }
#' Both are empty when the upload and the trait don't colocalise.
#' @examplesIf gpmap_available()
#' upload <- get_gwas("7a289615-c1b4-91f3-3d97-887f60de9155", include_associations = TRUE)
#' comparison <- compare_upload(upload, 2527)
#' head(comparison$coloc_groups)
#' @export
compare_upload <- function(upload, trait_id) {
  groups <- upload$coloc_groups
  if (!is.list(upload$trait) || is.null(upload$trait$id) || !is.data.frame(groups) ||
        !"source" %in% names(groups)) {
    stop("upload must be a GWAS upload result from get_gwas() or trait() with a GUID", call. = FALSE)
  }
  if (!"beta" %in% names(groups)) {
    stop(
      "upload has no associations: fetch it with include_associations = TRUE, ",
      "e.g. get_gwas(guid, include_associations = TRUE)",
      call. = FALSE
    )
  }
  if (length(trait_id) != 1 || is.na(trait_id)) {
    stop("trait_id must be a single trait id or GWAS upload GUID", call. = FALSE)
  }

  extractions <- as_table(upload$study_extractions)
  is_own_group <- groups$source %in% "upload" & groups$gwas_upload_id %in% upload$trait$id
  is_own_extraction <- row_sources(extractions) %in% "upload" & extractions$gwas_upload_id %in% upload$trait$id

  if (is_guid(trait_id)) {
    compared_upload_id <- unique(c(
      groups$gwas_upload_id[groups$source %in% "upload" & groups$study %in% trait_id],
      extractions$gwas_upload_id[row_sources(extractions) %in% "upload" & extractions$study %in% trait_id]
    ))
    if (length(compared_upload_id) == 0) {
      stop("trait_id ", trait_id, " is not a GWAS upload this upload was compared with", call. = FALSE)
    }
    is_trait_group <- groups$source %in% "upload" & groups$gwas_upload_id %in% compared_upload_id
    is_trait_extraction <- row_sources(extractions) %in% "upload" &
      extractions$gwas_upload_id %in% compared_upload_id
  } else {
    is_trait_group <- groups$source %in% "existing" & groups$trait_id %in% trait_id
    is_trait_extraction <- row_sources(extractions) %in% "existing" & extractions$trait_id %in% trait_id
  }

  upload_leads <- lead_per_coloc_group(groups[is_own_group, , drop = FALSE])
  trait_leads <- lead_per_coloc_group(groups[is_trait_group, , drop = FALSE])
  shared <- upload_leads[upload_leads$coloc_group_id %in% trait_leads$coloc_group_id, , drop = FALSE]
  trait_leads <- trait_leads[match(shared$coloc_group_id, trait_leads$coloc_group_id), , drop = FALSE]
  compared_groups <- data.frame(
    coloc_group_id = shared$coloc_group_id,
    ld_block_id = shared$ld_block_id,
    ld_block = column_or_na(shared, "ld_block"),
    variant_id = shared$variant_id,
    display_snp = column_or_na(shared, "display_snp"),
    h4_connectedness = column_or_na(shared, "h4_connectedness"),
    upload_study_extraction_id = shared$study_extraction_id,
    upload_beta = shared$beta,
    upload_se = shared$se,
    upload_p = column_or_na(shared, "p"),
    trait_name = column_or_na(trait_leads, "trait_name"),
    trait_study_extraction_id = trait_leads$study_extraction_id,
    trait_beta = trait_leads$beta,
    trait_se = trait_leads$se,
    trait_p = column_or_na(trait_leads, "p"),
    stringsAsFactors = FALSE
  )

  compared_pairs <- compare_upload_pairs(
    upload$coloc_pairs, extractions, is_own_extraction, is_trait_extraction, compared_groups$ld_block_id
  )
  return(list(coloc_groups = compared_groups, coloc_pairs = compared_pairs))
}

# The row of each coloc group with the smallest p-value (min_p when p is missing), among rows with an effect.
lead_per_coloc_group <- function(rows) {
  rows <- rows[!is.na(rows$beta) & !is.na(rows$se), , drop = FALSE]
  p <- column_or_na(rows, "p")
  p <- ifelse(is.na(p), column_or_na(rows, "min_p"), p)
  rows <- rows[order(rows$coloc_group_id, p, na.last = TRUE), , drop = FALSE]
  rows <- rows[!duplicated(rows$coloc_group_id), , drop = FALSE]
  rownames(rows) <- NULL
  return(rows)
}

# The coloc pairs between the upload's study extractions and the trait's, with the upload on the first side.
compare_upload_pairs <- function(pairs, extractions, is_own_extraction, is_trait_extraction, shared_ld_blocks) {
  empty <- data.frame(
    ld_block_id = integer(0), upload_study_extraction_id = integer(0), upload_display_snp = character(0),
    trait_name = character(0), trait_study_extraction_id = integer(0), trait_display_snp = character(0),
    h3 = numeric(0), h4 = numeric(0), in_shared_coloc_group = logical(0), stringsAsFactors = FALSE
  )
  if (!is.data.frame(pairs) || nrow(pairs) == 0 || nrow(extractions) == 0) {
    return(empty)
  }
  extraction_key <- paste(row_sources(extractions), extractions$id)
  own_keys <- extraction_key[is_own_extraction]
  trait_keys <- extraction_key[is_trait_extraction]
  key_a <- paste(pairs$source_a, pairs$study_extraction_a_id)
  key_b <- paste(pairs$source_b, pairs$study_extraction_b_id)
  upload_first <- key_a %in% own_keys & key_b %in% trait_keys
  upload_second <- key_b %in% own_keys & key_a %in% trait_keys
  keep <- upload_first | upload_second
  if (!any(keep)) {
    return(empty)
  }
  upload_key <- ifelse(upload_first, key_a, key_b)[keep]
  trait_key <- ifelse(upload_first, key_b, key_a)[keep]
  upload_row <- match(upload_key, extraction_key)
  trait_row <- match(trait_key, extraction_key)
  compared <- data.frame(
    ld_block_id = pairs$ld_block_id[keep],
    upload_study_extraction_id = extractions$id[upload_row],
    upload_display_snp = column_or_na(extractions, "display_snp")[upload_row],
    trait_name = column_or_na(extractions, "trait_name")[trait_row],
    trait_study_extraction_id = extractions$id[trait_row],
    trait_display_snp = column_or_na(extractions, "display_snp")[trait_row],
    h3 = pairs$h3[keep],
    h4 = pairs$h4[keep],
    stringsAsFactors = FALSE
  )
  compared$in_shared_coloc_group <- compared$ld_block_id %in% shared_ld_blocks
  rownames(compared) <- NULL
  return(compared)
}

# A column of a dataframe, or NA for each row when it has no such column.
column_or_na <- function(table, column) {
  if (column %in% names(table)) {
    return(table[[column]])
  }
  return(rep(NA, nrow(table)))
}
