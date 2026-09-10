# Render the null simulation vignette across the background-realism parameter
# grids. Every render is a pure null (K = 0): EBMF tries to discover programs
# from unstructured background, and the output is how often each validation
# check is fooled by randomness.
#
# The vignette's internal null versions all shift with the baseline parameters
# overridden below.
#
# Development settings below are intentionally small. For final results,
# increase n_reps / n_null / n_rep as appropriate.

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
# Noise (noise_sd)
# ---------------------------------------------------------------------------

noise_sds=(0.25 0.5 1 2 4)

# ---------------------------------------------------------------------------
# Pleiotropic hub traits (n_hub_traits)
# ---------------------------------------------------------------------------

hub_traits=(0 5 15 30)

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

# ===========================================================================
# Pleiotropic hub traits
# ===========================================================================

for value in "${hub_traits[@]}"; do

  render_one \
    "hubs${value}" \
    "n_hub_traits=${value}"

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
echo "Null hub traits:"
ls -1 investigation-univariate-null-simulations_hubs*.html