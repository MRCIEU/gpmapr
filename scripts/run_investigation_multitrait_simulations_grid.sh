# Render the multi-trait correspondence simulation vignette across its
# parameter grids.
#
# The vignette already crosses SNP overlap x planted direction (concordant /
# antagonistic) internally, plus a null -- that cross IS the study. The grids
# below vary things ORTHOGONAL to it, so each render evaluates the full design
# under a different condition.
#
# Investigation R: similarity-graph realism
#   Discovery and linking of the shared program are largely insensitive to the
#   background regime, but the NULL arm is not: a false-link rate measured in an
#   unrealistic noise regime is as weak here as in the univariate null study.
#   Do not quote a false-link rate from the snpsd0_corr0 arm.
#
# Investigation S: how many background studies a program needs
#   n_drivers_per_program is the number of driver studies behind each planted
#   program. Shrinking it shows when the shared program stops being found (and
#   by which validation gate it is filtered), and so stops being linked.
#
# These are production settings for the HPC run. quick = FALSE is passed on every
# render: quick caps n_sim at 2.

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
# Investigation S: driver studies per program
# ---------------------------------------------------------------------------

drivers_per_program=(3 5 8 15)

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
  # Full-depth run. Under quick = TRUE only 2 replicates are run, too few to
  # separate the versions.
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

echo ""
echo "Done."
echo ""
echo "Investigation R: similarity-graph realism (false-link rate only meaningful here):"
ls -1 investigation-multitrait-simulations_R_*.html

echo ""
echo "Investigation S: driver studies per program:"
ls -1 investigation-multitrait-simulations_S_*.html
