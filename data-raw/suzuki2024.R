# Build inst/extdata/suzuki2024/ from the supplementary material of
# Suzuki et al. (2024) Genetic drivers of heterogeneity in type 2 diabetes
# pathophysiology. Nature 627:347-357. doi:10.1038/s41586-024-07019-6
#
# The article and its supplementary tables are CC BY 4.0. The supplementary
# workbook is downloaded from the Europe PMC API (PMC10937372) rather than
# nature.com, which does not serve it to scripted clients.
#
# Run from the package root:  Rscript data-raw/suzuki2024.R

out_dir <- file.path("inst", "extdata", "suzuki2024")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

supp_zip <- tempfile(fileext = ".zip")
utils::download.file(
  "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10937372/supplementaryFiles",
  supp_zip,
  mode = "wb",
  quiet = TRUE
)
supp_dir <- tempfile("suzuki2024_")
utils::unzip(supp_zip, exdir = supp_dir)
xlsx <- file.path(supp_dir, "41586_2024_7019_MOESM3_ESM.xlsx")
stopifnot(file.exists(xlsx))

read_sheet <- function(sheet) {
  return(as.data.frame(suppressMessages(readxl::read_excel(
    xlsx,
    sheet = sheet,
    col_names = FALSE,
    col_types = "text"
  ))))
}

fill_down <- function(x) {
  for (i in seq_along(x)[-1]) {
    if (is.na(x[i])) x[i] <- x[i - 1]
  }
  return(x)
}

cluster_levels <- c(
  "Beta cell +PI", "Beta cell -PI", "Residual glycaemic", "Body fat",
  "Metabolic syndrome", "Obesity", "Lipodystrophy", "Liver/lipid metabolism"
)

# --- Supplementary Table 9: cluster sizes and T2D effects -------------------
st9 <- read_sheet("ST9")[4:11, 1:4]
names(st9) <- c("cluster", "n_signals", "t2d_log_or", "t2d_log_or_se")
st9$n_signals <- as.integer(st9$n_signals)
st9$t2d_log_or <- as.numeric(st9$t2d_log_or)
st9$t2d_log_or_se <- as.numeric(st9$t2d_log_or_se)
stopifnot(identical(st9$cluster, cluster_levels))

# Main text (Table 1, Extended Data Fig. 4, Fig. 2): defining profile, example
# loci, physiology and single-cell chromatin enrichment of each cluster.
cluster_annotation <- data.frame(
  cluster = cluster_levels,
  previously_reported = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
  cardiometabolic_profile = c(
    "+FG*, +2hG*, +HbA1c, +PI*",
    "+FG*, +2hG*, +HbA1c, -PI*",
    "+FG*, +HbA1c",
    "+Body fat, +ASAT*",
    "+FG*, +FI*, +WHR, +VAT*, -GFAT*, +TG, -HDL, +BP",
    "+BMI, +WHR, +body fat, +BMR, +TG, -HDL",
    "+FI*, +WHR, -body fat, -GFAT*, +TG, -HDL, +BP",
    "-LDL, -TC, +liver fat, +liver biomarkers"
  ),
  example_loci = c(
    "TCF7L2; KCNQ1; CDKAL1; CDKN2A-CDKN2B; SLC30A8",
    "CDC123-CAMK1D; HNF1B; KCNJ11-ABCC8; HNF4A; HNF1A",
    "GCC1-PAX4-LEP; ANKRD55; GCKR; UBE2E2",
    "ZMIZ1; HMGA2; CTBP1",
    "IGF2BP2; CCND2; HHEX-IDE; JAZF1; GPSM1",
    "FTO; MC4R; MACF1; TMEM18",
    "IRS1; GRB14-COBLL1; PPARG",
    "TOMM40-APOE-GIPR; TM6SF2; PNPLA3"
  ),
  insulin_secretion = c("-", "-", "-", "+", "+", "+", "+", "-"),
  insulin_sensitivity = c("+", "+", "-", "-", "-", "-", "-", "-"),
  open_chromatin_enrichment = c(
    "Fetal islets; adult islet alpha/beta/gamma/delta cells",
    "Fetal islets; adult islet endocrine cells; adult enterochromaffin cells",
    "Fetal islets; adult islet endocrine cells; fetal and adult pancreatic ductal cells",
    "None",
    "Adult pericytes; fetal endothelial cells; fetal mesangial cells; fetal fibroblasts; brain IT and SST+ neurons",
    paste(
      "Adult islet alpha/gamma/delta cells; fetal adrenal, heart and kidney cells;",
      "brain IT, SST+ and D1 medium spiny neurons"
    ),
    "Adult adipocytes",
    "Not estimated (3 signals)"
  ),
  stringsAsFactors = FALSE
)
clusters <- merge(st9, cluster_annotation, by = "cluster", sort = FALSE)
clusters <- clusters[match(cluster_levels, clusters$cluster), ]

# --- Supplementary Table 5: the 37 clustering traits ------------------------
st5 <- read_sheet("ST5")[4:40, 1:6]
names(st5) <- c(
  "trait_category", "trait_acronym", "trait", "sample_size", "ancestry",
  "reference"
)
st5$trait_category <- fill_down(st5$trait_category)
st5$sample_size <- as.integer(st5$sample_size)
stopifnot(nrow(st5) == 37, !anyNA(st5$trait_acronym))

# Curated lookup from each Suzuki trait to the GPMap trait row it corresponds
# to. This is used only to read existing rows of the pleiotropy matrix on the
# same traits Suzuki clustered on; no extra analyses are run on these traits.
# match_type:
#   same_gwas         the GWAS Suzuki used (same source and sample size)
#   proxy             a different GWAS of the same trait
#   unadjusted_proxy  Suzuki used a BMI-adjusted trait; GPMap has it unadjusted
#   none              no GPMap equivalent
# Where GPMap holds several GWAS of a trait, the one with the most coloc
# groups is used, because sparse rows are dropped by the univariate filters.
gpmap_map <- data.frame(
  trait_acronym = c(
    "FGadjBMI", "2hGadjBMI", "HbA1c", "PIadjBMI", "FIadjBMI", "ISIadjBMI",
    "IFCadjBMI", "SBPadjBMI", "DBPadjBMI", "PPadjBMI", "BW", "BMR", "BMI",
    "WHR", "WC", "HC", "BFP", "TFP", "GFATadjBMI", "VATadjBMI", "ASATadjBMI",
    "ALP", "GGT", "Tbil", "Dbil", "AST", "ALT", "LiverFat", "TG", "LDL",
    "HDL", "TC", "nonHDL", "APOA", "APOB", "CRP", "IGF1"
  ),
  gpmap_trait_id = c(
    913L, 912L, 1640L, 3L, 914L, 2339L,
    NA, 1996L, 1995L, 1488L, 134L, 2006L, 1992L,
    1994L, 4503L, 3342L, 2581L, 3424L, NA, 1324L, 1325L,
    1624L, 1634L, 1639L, NA, 1643L, 1642L, 1326L, 1631L, 939L,
    1630L, 1627L, NA, 1629L, 1626L, 2644L, 1647L
  ),
  match_type = c(
    "same_gwas", "same_gwas", "proxy", "proxy", "same_gwas", "proxy",
    "none", "unadjusted_proxy", "unadjusted_proxy", "unadjusted_proxy",
    "proxy", "proxy", "proxy",
    "proxy", "proxy", "proxy", "proxy", "proxy", "none", "unadjusted_proxy",
    "unadjusted_proxy",
    "proxy", "proxy", "proxy", "none", "proxy", "proxy", "same_gwas", "proxy",
    "proxy", "proxy", "proxy", "none", "proxy", "proxy", "proxy", "proxy"
  ),
  stringsAsFactors = FALSE
)
stopifnot(setequal(gpmap_map$trait_acronym, st5$trait_acronym))

pkgload::load_all(quiet = TRUE)
select_api("production")
gpmap_traits <- all_traits()
gpmap_traits <- gpmap_traits[!duplicated(gpmap_traits$id), ]
mapped_ids <- stats::na.omit(gpmap_map$gpmap_trait_id)
stopifnot(all(mapped_ids %in% gpmap_traits$id))

trait_map <- merge(st5, gpmap_map, by = "trait_acronym", sort = FALSE)
trait_map <- trait_map[match(st5$trait_acronym, trait_map$trait_acronym), ]
gpmap_idx <- match(trait_map$gpmap_trait_id, gpmap_traits$id)
trait_map$gpmap_study <- gpmap_traits$trait[gpmap_idx]
trait_map$gpmap_trait_name <- gpmap_traits$trait_name[gpmap_idx]
trait_map$gpmap_sample_size <- gpmap_traits$sample_size[gpmap_idx]
trait_map$gpmap_num_coloc_groups <- gpmap_traits$num_coloc_groups[gpmap_idx]

# --- Supplementary Table 7: cluster x trait profiles -------------------------
st7 <- read_sheet("ST7")
st7_header <- st7[3, ]
st7_body <- st7[5:41, ]
stopifnot(nrow(st7_body) == 37)
cluster_cols <- which(!is.na(st7_header[1, ]) & seq_along(st7_header) > 3)
profiles <- do.call(rbind, lapply(cluster_cols, function(j) {
  return(data.frame(
    cluster = as.character(st7_header[1, j]),
    trait_acronym = st7_body[[2]],
    beta = as.numeric(st7_body[[j]]),
    se = as.numeric(st7_body[[j + 1]]),
    p = as.numeric(st7_body[[j + 2]]),
    stringsAsFactors = FALSE
  ))
}))
stopifnot(
  setequal(profiles$cluster, cluster_levels),
  nrow(profiles) == 37 * 8,
  setequal(profiles$trait_acronym, st5$trait_acronym)
)

# --- Supplementary Table 6: index SNV cluster assignments --------------------
st6 <- read_sheet("ST6")[-(1:4), 1:6]
names(st6) <- c("locus", "chr", "rsid", "pos_b37", "cluster", "distance_to_centroid")
st6 <- st6[!is.na(st6$rsid), ]
st6$locus <- fill_down(st6$locus)
st6$chr <- as.integer(fill_down(st6$chr))
st6$pos_b37 <- as.integer(st6$pos_b37)
st6$distance_to_centroid <- as.numeric(st6$distance_to_centroid)
stopifnot(
  nrow(st6) == 1289,
  !anyNA(st6$locus),
  all(st6$cluster %in% cluster_levels)
)
cluster_counts <- table(factor(st6$cluster, levels = cluster_levels))
stopifnot(all(as.integer(cluster_counts) == clusters$n_signals))

# --- GRCh38 positions of the index SNVs --------------------------------------
# GPMap is on GRCh38 and Supplementary Table 6 gives GRCh37 positions. Each
# rsid is looked up in Ensembl (dbSNP mappings, no coordinate liftover) on both
# builds; the GRCh37 lookup must reproduce pos_b37, so a wrong or re-used rsid
# fails here. Retired rsids that dbSNP has merged into another are returned
# under the current rsid, which is kept in rsid_current.
ensembl_positions <- function(rsids, host) {
  post_batch <- function(ids) {
    resp <- httr::RETRY(
      "POST",
      paste0(host, "/variation/homo_sapiens"),
      httr::add_headers("Content-Type" = "application/json", Accept = "application/json"),
      body = jsonlite::toJSON(list(ids = ids)),
      encode = "raw",
      httr::timeout(120),
      times = 8,
      pause_cap = 120,
      quiet = TRUE
    )
    httr::stop_for_status(resp)
    return(jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE))
  }
  get_one <- function(id) {
    resp <- httr::RETRY(
      "GET",
      paste0(host, "/variation/human/", id),
      httr::add_headers(Accept = "application/json"),
      httr::timeout(120),
      times = 8,
      pause_cap = 120,
      quiet = TRUE
    )
    httr::stop_for_status(resp)
    return(jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE))
  }
  to_row <- function(query, record) {
    mappings <- Filter(function(m) m$seq_region_name %in% c(1:22, "X"), record$mappings)
    stopifnot(length(mappings) == 1)
    return(data.frame(
      rsid = query,
      rsid_current = record$name,
      chr = as.integer(mappings[[1]]$seq_region_name),
      pos = as.integer(mappings[[1]]$start),
      stringsAsFactors = FALSE
    ))
  }
  batches <- lapply(split(rsids, ceiling(seq_along(rsids) / 200)), post_batch)
  records <- do.call(c, batches)
  found <- intersect(rsids, names(records))
  out <- do.call(rbind, lapply(found, function(id) return(to_row(id, records[[id]]))))
  # Merged rsids come back keyed by their current rsid; look them up one by one.
  merged <- setdiff(rsids, found)
  if (length(merged) > 0) {
    out <- rbind(out, do.call(rbind, lapply(merged, function(id) return(to_row(id, get_one(id))))))
  }
  return(out[match(rsids, out$rsid), ])
}
b37 <- ensembl_positions(st6$rsid, "https://grch37.rest.ensembl.org")
b38 <- ensembl_positions(st6$rsid, "https://rest.ensembl.org")
stopifnot(
  identical(b37$rsid, st6$rsid),
  identical(b38$rsid, st6$rsid),
  all(b37$chr == st6$chr),
  all(b37$pos == st6$pos_b37),
  all(b38$chr == st6$chr),
  identical(b37$rsid_current, b38$rsid_current)
)
st6$pos_b38 <- b38$pos
st6$rsid_current <- b38$rsid_current
st6 <- st6[, c(
  "locus", "chr", "rsid", "rsid_current", "pos_b37", "pos_b38", "cluster",
  "distance_to_centroid"
)]

utils::write.csv(clusters, file.path(out_dir, "clusters.csv"), row.names = FALSE)
utils::write.csv(trait_map, file.path(out_dir, "trait_map.csv"), row.names = FALSE)
utils::write.csv(
  profiles,
  file.path(out_dir, "cluster_trait_profiles.csv"),
  row.names = FALSE
)
utils::write.csv(
  st6,
  file.path(out_dir, "index_snv_clusters.csv"),
  row.names = FALSE
)
