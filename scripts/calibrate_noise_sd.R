#!/usr/bin/env Rscript
# Empirical noise calibration for simulate_trait() from the real BMI matrix.
#
# Rationale: simulate_trait()'s noise_sd is the SD of the N(0, noise_sd) noise
# added to trait z-score observations (default log_se_sd = 0, so beta == z).
# We estimate a plausible value from the real trait x SNP z-score matrix using
# the per-SNP scaled MAD (stats::mad, i.e. a robust SD analogue) of the
# across-trait z-scores, computed separately for:
#   - "program"    : SNPs assigned to any EBMF program
#   - "background" : SNPs not assigned to any program
#   - "all"        : every SNP in the matrix
# Comparing the three groups shows whether the noise scale differs between
# structured and unstructured SNPs.
#
# Usage: Rscript scripts/calibrate_noise_sd.R
# Set `output_dir` below to keep the plot + tables somewhere permanent.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
pkg_root <- normalizePath(file.path(dirname(script_path), ".."))

# ---- Settings ---------------------------------------------------------------
trait_id <- 1992L                 # BMI
output_dir <- file.path(pkg_root, "analysis", "noise_calibration")
min_obs <- 10L                    # min non-missing traits per SNP
min_snp_signals <- 5L
compress_method <- "asinh"
compress_scale <- 2
similarity_threshold <- 0.2
ebmf_lfsr_threshold <- 0.05
ebmf_magnitude_threshold <- 0.5
ebmf_drop_global <- TRUE
ebmf_prior <- "point_normal"
ebmf_backfit <- TRUE
ebmf_hard_assignment <- FALSE
min_module_size <- 3L
n_null <- 5L
cores <- 5L
multipliers <- c(clean = 0.5, moderate = 1.0, noisy = 1.5, very_noisy = 2.0)
group_order <- c("background", "program", "all")
group_cols <- c(
  background = "#7570b3", program = "#1b9e77", all = "#d95f02"
)

# ---- Load package from source ----------------------------------------------
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(gpmapr)
}

# ---- Fetch real BMI data and run the EBMF pipeline --------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
select_api("production")
obj <- trait(trait_id, include_associations = TRUE)

res <- run_univariate_clustering(
  obj,
  cluster_type = "ebmf",
  min_snp_signals = min_snp_signals,
  compress_method = compress_method,
  compress_scale = compress_scale,
  similarity_threshold = similarity_threshold,
  ebmf_lfsr_threshold = ebmf_lfsr_threshold,
  ebmf_magnitude_threshold = ebmf_magnitude_threshold,
  ebmf_drop_global = ebmf_drop_global,
  ebmf_prior = ebmf_prior,
  ebmf_backfit = ebmf_backfit,
  ebmf_hard_assignment = ebmf_hard_assignment,
  min_module_size = min_module_size
)

# Program membership: SNPs gated by lFSR + magnitude (stability skipped, n_rep = 0).
smry <- summarise_ebmf_programs(
  res,
  s_matrix = res$s_matrix,
  edge_threshold = similarity_threshold,
  min_module_size = min_module_size,
  n_null = n_null,
  n_rep = 0,
  cores = cores,
  verbose = FALSE
)
assigned_snps <- unique(smry$assigned$snp_id)

# Raw traits x SNPs z-score matrix (after trait filtering, before compression).
X <- res$x_matrix
target_row <- as.character(trait_id)
bg_rows <- setdiff(rownames(X), target_row)   # exclude the target trait itself
bg_snps <- setdiff(colnames(X), assigned_snps)

# ---- Per-SNP robust across-trait scale, per SNP group ------------------------
mad_by_snp <- function(snp) {
  vals <- X[bg_rows, snp]
  vals <- vals[!is.na(vals)]
  if (length(vals) < min_obs) {
    return(NA_real_)
  }
  stats::mad(vals)   # scaled MAD (constant 1.4826) = robust SD analogue
}

compute_group_mad <- function(snps, label) {
  m <- vapply(snps, mad_by_snp, numeric(1))
  data.frame(
    group = label,
    snp_id = names(m),
    mad = unname(m),
    stringsAsFactors = FALSE
  )
}

mad_df <- dplyr::bind_rows(
  compute_group_mad(bg_snps, "background"),
  compute_group_mad(assigned_snps, "program"),
  compute_group_mad(colnames(X), "all")
) |>
  dplyr::filter(!is.na(mad)) |>
  dplyr::mutate(group = factor(group, levels = group_order))

group_totals <- data.frame(
  group = c("background", "program", "all"),
  total = c(length(bg_snps), length(assigned_snps), length(colnames(X))),
  stringsAsFactors = FALSE
)

group_counts <- mad_df |>
  dplyr::group_by(group) |>
  dplyr::summarise(used = dplyr::n(), .groups = "drop")

summary_df <- mad_df |>
  dplyr::group_by(group) |>
  dplyr::summarise(
    P10 = stats::quantile(mad, 0.10),
    P25 = stats::quantile(mad, 0.25),
    P50 = stats::quantile(mad, 0.50),
    P75 = stats::quantile(mad, 0.75),
    P90 = stats::quantile(mad, 0.90),
    .groups = "drop"
  ) |>
  dplyr::left_join(group_totals, by = "group") |>
  dplyr::left_join(group_counts, by = "group") |>
  dplyr::mutate(excluded = total - used) |>
  dplyr::relocate(group, total, used, excluded)

noise_sd_real <- stats::setNames(summary_df$P50, summary_df$group)

calibration_df <- do.call(rbind, lapply(group_order, function(g) {
  data.frame(
    group = g,
    scenario = names(multipliers),
    multiplier = unname(multipliers),
    noise_sd = round(unname(multipliers) * noise_sd_real[[g]], 4),
    stringsAsFactors = FALSE
  )
})) |>
  dplyr::mutate(group = factor(group, levels = group_order))

# ---- Report -----------------------------------------------------------------
cat("==============================================================\n")
cat("noise_sd calibration from real BMI (trait ", trait_id, ") matrix\n", sep = "")
cat("==============================================================\n")
cat("trait x SNP matrix: ", nrow(X), " traits x ", ncol(X), " SNPs\n", sep = "")
cat("MAD is the scaled MAD (stats::mad default), i.e. a robust SD analogue\n",
    "directly comparable to simulate_trait()'s noise_sd.\n", sep = "")
cat("SNPs are split by EBMF program membership (min_obs = ", min_obs, ").\n\n", sep = "")

cat("Per-SNP MAD distribution by SNP group:\n")
print(summary_df |> dplyr::mutate(dplyr::across(
  dplyr::where(is.numeric), ~ round(.x, 4)
)), row.names = FALSE)

cat("\nEmpirical noise scale (P50 per-SNP MAD) by group:\n")
print(data.frame(
  group = group_order,
  noise_sd_real = round(unname(noise_sd_real), 4),
  P10 = round(summary_df$P10[match(group_order, summary_df$group)], 4),
  P90 = round(summary_df$P90[match(group_order, summary_df$group)], 4)
), row.names = FALSE)

cat("\nSuggested simulation noise_sd values by group:\n")
print(calibration_df, row.names = FALSE)

# ---- Diagnostic plot ---------------------------------------------------------
p <- ggplot(mad_df, aes(mad, colour = group, fill = group)) +
  geom_density(alpha = 0.35) +
  geom_vline(
    aes(xintercept = P50, colour = group),
    data = summary_df, linewidth = 1, linetype = "solid"
  ) +
  scale_colour_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  coord_cartesian(xlim = stats::quantile(mad_df$mad, c(0, 0.99))) +
  labs(
    x = "per-SNP MAD of trait z-scores",
    y = "density",
    title = "Empirical noise scale from BMI background vs program SNPs",
    subtitle = paste0(
      "solid lines = group medians (background noise_sd_real = ",
      round(noise_sd_real[["background"]], 3), ", program = ",
      round(noise_sd_real[["program"]], 3), ", all = ",
      round(noise_sd_real[["all"]], 3), ")"
    ),
    colour = "SNP group",
    fill = "SNP group"
  ) +
  theme_minimal()

plot_path <- file.path(output_dir, "noise_calibration.png")
ggsave(plot_path, p, width = 8, height = 5.5, dpi = 150)

# ---- Write outputs ------------------------------------------------------------
write.csv(summary_df, file.path(output_dir, "noise_calibration_summary.csv"),
          row.names = FALSE)
write.csv(calibration_df, file.path(output_dir, "noise_calibration.csv"),
          row.names = FALSE)
write.csv(
  mad_df |> dplyr::arrange(group, snp_id),
  file.path(output_dir, "noise_calibration_snp_mad.csv"),
  row.names = FALSE
)

cat("\nOutputs written to: ", output_dir, "\n", sep = "")
cat("  plot:   noise_calibration.png\n")
cat("  tables: noise_calibration_summary.csv, noise_calibration.csv,\n")
cat("          noise_calibration_snp_mad.csv\n")