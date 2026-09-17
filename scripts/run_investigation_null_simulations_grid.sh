# Render the null simulation vignette across the background-realism parameter
# grids. Every render is a pure null (K = 0): EBMF tries to discover programs
# from unstructured background, and the output is how often each validation
# check is fooled by randomness.
#
# The vignette's internal null versions all shift with the baseline parameters
# overridden below.
#
# Production settings for the HPC run. quick = FALSE is passed below so the
# stability gate is actually exercised -- it was never evaluated in any
# simulation, while being the most influential gate on real data. The stability
# threshold should be read off this study's replication distribution rather than
# assumed.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIGNETTES_DIR="$(cd "${SCRIPT_DIR}/../vignettes" && pwd)"
RMD="investigation-univariate-null-simulations.Rmd"

# ---------------------------------------------------------------------------
# Replication settings
# ---------------------------------------------------------------------------

N_REPS=5
N_NULL=5
N_REP=5

# ---------------------------------------------------------------------------
# Background density (p_active_background)
# ---------------------------------------------------------------------------

background_rates=(0 0.02 0.06 0.1 0.2 0.4)

# ---------------------------------------------------------------------------
# Per-trait sparsity heterogeneity (background_sparsity_sd)
# ---------------------------------------------------------------------------

background_sds=(0 0.6 1.2 1.8)

# ---------------------------------------------------------------------------
# Background magnitude relative to the module effect size
# ---------------------------------------------------------------------------

background_scales=(0.25 0.5 1 1.5)

# ---------------------------------------------------------------------------
# Generative realism of the similarity graph
#
# background_corr = 0 (the old default) gives background traits independent,
# random-signed effects, which contribute nothing to SNP-SNP similarity: the
# simulated graph sits near zero while a real one sits around 0.35, and a random
# SNP set clears the old 0.3 gate ~0% of the time in simulation against ~94% on
# real BMI data. Under that regime a null study cannot fail, so the reported
# FDR = 0 was close to un-failable by construction.
#
# snp_pleiotropy_sd is the per-SNP analogue of background_sparsity_sd and is
# what brings the median pairwise trait overlap down from ~52 to the real ~3.
#
# Each render reports simulated_graph_diagnostics() against
# real_bmi_graph_targets(). The null FDR is only meaningful for the settings
# that land near those targets.
# ---------------------------------------------------------------------------

# CROSSED, not one-at-a-time. Varying either axis alone from the baseline
# p_active_background = 0.1 never reaches the target regime: every one-at-a-time
# arm gives random_set_pass_rate between 0.000 and 0.245 against the real 0.935,
# because the baseline is far too dense (median pairwise trait overlap ~15
# against the real ~3). The three axes have to move together. Each entry below is
# "p_active_background:snp_pleiotropy_sd:background_corr".
realism_grid=(
  0.10:0:0
  0.02:1.0:0.6
  0.02:1.5:0.9
  0.02:2.0:0.9
  0.01:1.5:0.9
  0.01:2.0:0.9
)

# ---------------------------------------------------------------------------
# Noise (noise_sd)
# ---------------------------------------------------------------------------

noise_sds=(0.25 0.5 1 2 4)

# ---------------------------------------------------------------------------
# Rendering helper
# ---------------------------------------------------------------------------

render_one() {
  local label=$1
  shift

  local out="investigation-univariate-null-simulations_${label}.html"

  echo "$(date +%Y-%m-%d\ %H:%M:%S) >>> Rendering ${label} -> ${out}"

  rm -f "${RMD%.Rmd}.html"

  Rscript - "$RMD" "$N_REPS" "$N_NULL" "$N_REP" "$@" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

rmd <- args[[1]]
n_reps <- as.integer(args[[2]])
n_null <- as.integer(args[[3]])
n_rep <- as.integer(args[[4]])

# Remaining arguments are name=value pairs.
overrides <- args[-(1:4)]

params <- list(
  # The vignettes declare n_sim / n_stability_rep; n_reps / n_rep are the old
  # names and rmarkdown rejects params it has not declared.
  n_sim = n_reps,
  n_null = n_null,
  n_stability_rep = n_rep,
  quick = FALSE
)

for (x in overrides) {
  parts <- strsplit(x, "=", fixed = TRUE)[[1]]
  name <- parts[[1]]
  value <- parts[[2]]

  # Numeric vectors are supplied as comma-separated values.
  if (grepl(",", value, fixed = TRUE)) {
    params[[name]] <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
  } else if (value %in% c("TRUE", "FALSE")) {
    params[[name]] <- as.logical(value)
  } else if (grepl("^-?[0-9.]+$", value)) {
    params[[name]] <- as.numeric(value)
  } else {
    params[[name]] <- value
  }
}

rmarkdown::render(
  rmd,
  params = params,
  envir = new.env(),
  quiet = TRUE
)
RSCRIPT

  mv -f "${RMD%.Rmd}.html" "${out}"

  echo "    -> ${out}"
}

cd "${VIGNETTES_DIR}"

# ===========================================================================
# Background density
# ===========================================================================

for value in "${background_rates[@]}"; do

  render_one \
    "bgrate${value}" \
    "p_active_background=${value}"

done

# ===========================================================================
# Per-trait sparsity heterogeneity
# ===========================================================================

for value in "${background_sds[@]}"; do

  render_one \
    "bgsd${value}" \
    "background_sparsity_sd=${value}"

done

# ===========================================================================
# Background magnitude
# ===========================================================================

for value in "${background_scales[@]}"; do

  render_one \
    "bgscale${value}" \
    "background_effect_scale=${value}"

done

# ===========================================================================
# Noise
# ===========================================================================

for value in "${noise_sds[@]}"; do

  render_one \
    "noise${value}" \
    "noise_sd=${value}"

done

for combo in "${realism_grid[@]}"; do

  IFS=: read -r p_act snp_sd bg_corr <<< "${combo}"

  render_one \
    "realism_p${p_act}_snpsd${snp_sd}_corr${bg_corr}" \
    "p_active_background=${p_act}" \
    "snp_pleiotropy_sd=${snp_sd}" \
    "background_corr=${bg_corr}"

done

echo ""
echo "Done."
echo ""
echo "Null background density:"
ls -1 investigation-univariate-null-simulations_bgrate*.html

echo ""
echo "Null background sparsity:"
ls -1 investigation-univariate-null-simulations_bgsd*.html

echo ""
echo "Null background magnitude:"
ls -1 investigation-univariate-null-simulations_bgscale*.html

echo ""
echo "Null noise:"
ls -1 investigation-univariate-null-simulations_noise*.html

echo ""
echo "Similarity-graph realism (crossed density x per-SNP heterogeneity x background correlation):"
echo "Read random_set_pass_rate in each render's graph-diagnostics table; the"
echo "null FDR is only meaningful for arms that land near the real targets."
ls -1 investigation-univariate-null-simulations_realism*.html