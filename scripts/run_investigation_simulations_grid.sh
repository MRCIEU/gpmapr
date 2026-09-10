# Render the univariate simulation vignette across the EBMF clustering and
# program validation parameter grids.
#
# Investigation E: EBMF clustering parameters
#   compress_method: none vs asinh (compress_scale = 5)
#   ebmf_prior: point_normal vs point_laplace
#   ebmf_magnitude_threshold: 0.25 / 0.5 / 0.75
#
# Investigation V: program validation parameters
#   similarity_threshold: 0.2 vs 0.5
#   min_module_size: 2 vs 5
#   min_connectedness: 0.25 vs 0.5
#   stability_threshold: 0.3 vs 0.7
#
# Each parameter is varied independently while the others remain at their
# values in the Rmd. The vignette itself retains its generative input-version
# tiers (background density x trait overlap, driver architecture, spurious
# hits) for every run, so each analysis parameter is evaluated across the full
# battery of input versions.
#
# Development settings below are intentionally small. For final results,
# increase n_reps / n_null / n_rep as appropriate.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIGNETTES_DIR="$(cd "${SCRIPT_DIR}/../vignettes" && pwd)"
RMD="investigation-univariate-simulations.Rmd"

# ---------------------------------------------------------------------------
# Replication settings
# ---------------------------------------------------------------------------

N_REPS=5
N_NULL=5
N_REP=5

# ---------------------------------------------------------------------------
# Investigation E: EBMF clustering parameters
#
# Rmd baselines for reference: compress_method = asinh (compress_scale = 5),
# ebmf_prior = point_normal, ebmf_magnitude_threshold = 0.3.
# ---------------------------------------------------------------------------

# Paired arrays: compress_scale only applies to the asinh arm (ignored by
# "none", which returns the matrix unchanged).
compress_methods=(none asinh)
compress_scales=(5 5)

ebmf_priors=(point_normal point_laplace)

magnitude_thresholds=(0.25 0.5 0.75)

# ---------------------------------------------------------------------------
# Investigation V: program validation parameters
#
# Rmd baselines for reference: similarity_threshold = 0.2,
# min_module_size = 3, min_connectedness = 0.25, stability_threshold = 0.5.
# ---------------------------------------------------------------------------

similarity_thresholds=(0.2 0.5)

module_size_thresholds=(2 5)

connectedness_thresholds=(0.25 0.5)

stability_thresholds=(0.3 0.7)

# ---------------------------------------------------------------------------
# Rendering helper
# ---------------------------------------------------------------------------

render_one() {
  local label=$1
  shift

  local out="investigation-univariate-simulations_${label}.html"

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
  n_reps = n_reps,
  n_null = n_null,
  n_rep = n_rep
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
# Investigation E: EBMF clustering parameters
# ===========================================================================

for i in "${!compress_methods[@]}"; do

  render_one \
    "E_compress_${compress_methods[$i]}" \
    "compress_method=${compress_methods[$i]}" \
    "compress_scale=${compress_scales[$i]}"

done

for prior in "${ebmf_priors[@]}"; do

  render_one \
    "E_prior_${prior}" \
    "ebmf_prior=${prior}"

done

for value in "${magnitude_thresholds[@]}"; do

  render_one \
    "E_magnitude_${value}" \
    "ebmf_magnitude_threshold=${value}"

done

# ===========================================================================
# Investigation V: program validation parameters
# ===========================================================================

for value in "${similarity_thresholds[@]}"; do

  render_one \
    "V_simthr_${value}" \
    "similarity_threshold=${value}"

done

for value in "${module_size_thresholds[@]}"; do

  render_one \
    "V_modulesize_${value}" \
    "min_module_size=${value}"

done

for value in "${connectedness_thresholds[@]}"; do

  render_one \
    "V_connectedness_${value}" \
    "min_connectedness=${value}"

done

for value in "${stability_thresholds[@]}"; do

  render_one \
    "V_stability_${value}" \
    "stability_threshold=${value}"

done

echo ""
echo "Done."
echo ""
echo "Investigation E: EBMF clustering parameters:"
ls -1 investigation-univariate-simulations_E_*.html

echo ""
echo "Investigation V: program validation parameters:"
ls -1 investigation-univariate-simulations_V_*.html
