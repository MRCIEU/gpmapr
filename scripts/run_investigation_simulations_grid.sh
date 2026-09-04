# Render the univariate simulation vignette across three simulation
# investigations.
#
# Investigation A: overall signal-to-noise
#   noise_sd = 0.25 / 0.5 / 1 / 2 / 4 / 6
#
# Investigation A2: individual sources of realism/noise
#   p_structural_zero
#   p_spurious
#   p_active_background
#   background_sparsity_sd
#   n_hub_traits
#   effect_tail
#
# Investigation B: program size and number of defining traits
#   module_sizes x n_traits_per_module
#
# The vignette itself retains its existing trait-overlap tiers for each run.
# All other simulation parameters remain at their values in the Rmd unless
# explicitly overridden below.
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
# Investigation A: how noisy?
#
# effect_size stays fixed at the Rmd default (6).
# noise_sd is the primary signal-to-noise parameter.
# ---------------------------------------------------------------------------

noise_sds=(0.25 0.5 1 2 4 6)

# ---------------------------------------------------------------------------
# Investigation A2: what type of noise/realism matters?
#
# Each parameter is varied independently while the others remain at their
# realistic baseline values.
# ---------------------------------------------------------------------------

structural_zeros=(0 0.1 0.2 0.4 0.6)
spurious_rates=(0 0.005 0.02 0.05)
background_rates=(0.002 0.006 0.02 0.05)
background_sds=(0 0.6 1.2 1.8)
hub_traits=(0 5 15 30)
effect_tails=(0 0.2 0.4 0.8)

# ---------------------------------------------------------------------------
# Investigation B: module size x number of defining traits
#
# Keep the three modules equal-sized for this experiment so that the
# detection boundary is easier to interpret.
# ---------------------------------------------------------------------------

module_sizes=(5 10 20 40 80)
traits_per_module=(2 3 5 10 20)

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
# Investigation A: How noisy?
# ===========================================================================

for noise in "${noise_sds[@]}"; do

  render_one \
    "A_noise${noise}" \
    "noise_sd=${noise}"

done

# ===========================================================================
# Investigation A2: Individual noise / realism parameters
# ===========================================================================

for value in "${structural_zeros[@]}"; do

  render_one \
    "A2_structuralzero${value}" \
    "p_structural_zero=${value}"

done

for value in "${spurious_rates[@]}"; do

  render_one \
    "A2_spurious${value}" \
    "p_spurious=${value}"

done

for value in "${background_rates[@]}"; do

  render_one \
    "A2_background${value}" \
    "p_active_background=${value}"

done

for value in "${background_sds[@]}"; do

  render_one \
    "A2_backgroundsd${value}" \
    "background_sparsity_sd=${value}"

done

for value in "${hub_traits[@]}"; do

  render_one \
    "A2_hubs${value}" \
    "n_hub_traits=${value}"

done

for value in "${effect_tails[@]}"; do

  render_one \
    "A2_effecttail${value}" \
    "effect_tail=${value}"

done

# ===========================================================================
# Investigation B: Module size x number of defining traits
# ===========================================================================

for module_size in "${module_sizes[@]}"; do
  for n_traits in "${traits_per_module[@]}"; do

    render_one \
      "B_snps${module_size}_traits${n_traits}" \
      "module_sizes=${module_size},${module_size},${module_size}" \
      "n_traits_per_module=${n_traits},${n_traits},${n_traits}"

  done
done

echo ""
echo "Done."
echo ""
echo "Investigation A:"
ls -1 investigation-univariate-simulations_A_noise*.html

echo ""
echo "Investigation A2:"
ls -1 investigation-univariate-simulations_A2_*.html

echo ""
echo "Investigation B:"
ls -1 investigation-univariate-simulations_B_*.html
