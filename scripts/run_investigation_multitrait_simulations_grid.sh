# Render the multi-trait correspondence simulation vignette across its
# parameter grids.
#
# The vignette already crosses locus_overlap x profile_overlap internally -- that
# cross IS the study, and it is what decides whether the trait-profile axis
# recovers shared programs the locus axis cannot. The grids below therefore vary
# things ORTHOGONAL to that cross, so each render evaluates the full
# locus x profile design under a different condition.
#
# Investigation R: similarity-graph realism
#   The headline contrast between the two axes is structural and largely
#   insensitive to the background regime, but the NULL arm is not: a
#   false-discovery rate measured in an unrealistic noise regime is as weak here
#   as in the univariate null study. Do not quote an FDR from the
#   snpsd0_corr0 arm.
#
# Investigation S: how much signal each axis is given
#   n_drivers_per_program is the size of the profile signal; module_size the size
#   of the locus signal; n_shared_loci caps how much locus overlap is achievable
#   at all. Shrinking each in turn shows which axis degrades first.
#
# Investigation T: thresholds and contamination
#
# These are production settings for the HPC run. quick = FALSE is passed on every
# render: quick caps n_sim at 2, at which recall saturates at 1.00 and the two
# profile_overlap levels cannot be told apart.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIGNETTES_DIR="$(cd "${SCRIPT_DIR}/../vignettes" && pwd)"
RMD="investigation-multitrait-simulations.Rmd"

# ---------------------------------------------------------------------------
# Replication settings
# ---------------------------------------------------------------------------

N_SIM=50
N_PERM=1000

# ---------------------------------------------------------------------------
# Investigation R: similarity-graph realism
#
# CROSSED, not one-at-a-time: varying either axis alone never reaches the regime
# real data is in (see simulated_graph_diagnostics() and
# real_bmi_graph_targets()). Entries are
# "p_active_background:snp_pleiotropy_sd:background_corr".
# ---------------------------------------------------------------------------

realism_grid=(
  0.02:0:0
  0.02:1.0:0.6
  0.02:1.5:0.9
  0.02:2.0:0.9
)

# ---------------------------------------------------------------------------
# Investigation S: how much signal each axis is given
# ---------------------------------------------------------------------------

drivers_per_program=(3 5 8)

module_sizes=(8 15 25)

shared_locus_pools=(20 80 160)

# ---------------------------------------------------------------------------
# Investigation T: thresholds and contamination
# ---------------------------------------------------------------------------

fdr_thresholds=(0.01 0.05 0.10)

structural_zeros=(0 0.2 0.4)

# ---------------------------------------------------------------------------
# Rendering helper
# ---------------------------------------------------------------------------

render_one() {
  local label=$1
  shift

  local out="investigation-multitrait-simulations_${label}.html"

  echo "$(date +%Y-%m-%d\ %H:%M:%S) >>> Rendering ${label} -> ${out}"

  rm -f "${RMD%.Rmd}.html"

  Rscript - "$RMD" "$N_SIM" "$N_PERM" "$@" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

rmd <- args[[1]]
n_sim <- as.integer(args[[2]])
n_perm <- as.integer(args[[3]])

overrides <- args[-(1:3)]

params <- list(
  n_sim = n_sim,
  n_perm = n_perm,
  # Full-depth run. Under quick = TRUE recall saturates at 1.00 on 2 replicates
  # and the input versions cannot be separated.
  quick = FALSE
)

for (x in overrides) {
  parts <- strsplit(x, "=", fixed = TRUE)[[1]]
  name <- parts[[1]]
  value <- parts[[2]]

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

rmarkdown::render(rmd, params = params, envir = new.env())
RSCRIPT

  mv "${RMD%.Rmd}.html" "$out"
}

cd "$VIGNETTES_DIR"

# ---------------------------------------------------------------------------
# Investigation R
# ---------------------------------------------------------------------------

for combo in "${realism_grid[@]}"; do

  IFS=: read -r p_act snp_sd bg_corr <<< "${combo}"

  render_one \
    "R_realism_p${p_act}_snpsd${snp_sd}_corr${bg_corr}" \
    "p_active_background=${p_act}" \
    "snp_pleiotropy_sd=${snp_sd}" \
    "background_corr=${bg_corr}"

done

# ---------------------------------------------------------------------------
# Investigation S
# ---------------------------------------------------------------------------

for value in "${drivers_per_program[@]}"; do

  render_one \
    "S_drivers_${value}" \
    "n_drivers_per_program=${value}"

done

for value in "${module_sizes[@]}"; do

  render_one \
    "S_modulesize_${value}" \
    "module_size=${value}"

done

for value in "${shared_locus_pools[@]}"; do

  render_one \
    "S_sharedloci_${value}" \
    "n_shared_loci=${value}"

done

# ---------------------------------------------------------------------------
# Investigation T
# ---------------------------------------------------------------------------

for value in "${fdr_thresholds[@]}"; do

  render_one \
    "T_fdr_${value}" \
    "fdr_threshold=${value}"

done

for value in "${structural_zeros[@]}"; do

  render_one \
    "T_structzero_${value}" \
    "p_structural_zero=${value}"

done

echo ""
echo "Done."
echo ""
echo "Investigation R: similarity-graph realism (FDR only meaningful here):"
ls -1 investigation-multitrait-simulations_R_*.html

echo ""
echo "Investigation S: signal available to each axis:"
ls -1 investigation-multitrait-simulations_S_*.html

echo ""
echo "Investigation T: thresholds and contamination:"
ls -1 investigation-multitrait-simulations_T_*.html
