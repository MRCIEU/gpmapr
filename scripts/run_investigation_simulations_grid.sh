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
# These are production settings, intended for the HPC run. The vignette must
# also be rendered with quick = FALSE (passed below) so stability is actually
# exercised: it was disabled in every simulation run so far, while being the
# most influential gate on real data.

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

module_size_thresholds=(2 3 5)

# The coherence gate replaces min_mean_internal / min_connectedness, which were
# absolute thresholds sitting below the real similarity graph's own baseline.
# What needs calibrating now is the BH level and the permutation count.
coherence_q_thresholds=(0.01 0.05 0.10)

stability_thresholds=(0.3 0.5 0.7)

# greedy_Kmax was binding at 50 (null fits returned 49-50 factors, real BMI
# returned exactly 50). Confirm where it stops binding.
greedy_kmax_values=(50 100 200)

# Background realism. background_corr = 0 is the old independent, random-signed
# background, which contributes nothing to SNP-SNP similarity; snp_pleiotropy_sd
# is the per-SNP analogue of background_sparsity_sd. Both are needed to move the
# simulated similarity graph towards the real one -- see
# simulated_graph_diagnostics() and real_bmi_graph_targets().
# CROSSED, not one-at-a-time: varying either axis alone never reaches the target
# regime (every one-at-a-time arm gives random_set_pass_rate <= 0.245 against the
# real 0.935). Entries are "snp_pleiotropy_sd:background_corr".
#
# Unlike the null study, p_active_background is NOT overridden here: this
# vignette already iterates it as a generative axis crossed with trait_overlap
# ([0.00, 0.02, 0.05, 0.07]), so each render below evaluates the realism setting
# at every background density. Overriding it would collapse that grid. The
# p_active_background = 0.02 versions are the ones that land nearest the real
# similarity graph -- check each render's graph-diagnostics table.
realism_grid=(
  0:0
  1.0:0.6
  1.5:0.9
  2.0:0.9
)

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
  # The vignettes declare n_sim / n_stability_rep; n_reps / n_rep are the old
  # names and rmarkdown rejects params it has not declared.
  n_sim = n_reps,
  n_null = n_null,
  n_stability_rep = n_rep,
  # Full-depth run: quick = TRUE caps n_sim at 3 and disables the stability
  # subsampling entirely, which is how every simulation so far was rendered.
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

for value in "${coherence_q_thresholds[@]}"; do

  render_one \
    "V_coherenceq_${value}" \
    "coherence_q=${value}"

done

for value in "${stability_thresholds[@]}"; do

  render_one \
    "V_stability_${value}" \
    "stability_threshold=${value}"

done

for value in "${greedy_kmax_values[@]}"; do

  render_one \
    "V_kmax_${value}" \
    "ebmf_greedy_Kmax=${value}"

done

# ---------------------------------------------------------------------------
# Investigation G: generative realism of the similarity graph
#
# The simulated similarity graph does not resemble a real one: a random SNP set
# clears the old 0.3 internal-similarity gate essentially never in simulation
# and about 94% of the time on real BMI data, so the null study could not
# detect that the gate was uninformative. These two axes are what move it.
# Every render reports simulated_graph_diagnostics() against
# real_bmi_graph_targets(); pick the pairing that lands nearest the targets.
# ---------------------------------------------------------------------------

for combo in "${realism_grid[@]}"; do

  IFS=: read -r snp_sd bg_corr <<< "${combo}"

  render_one \
    "G_realism_snpsd${snp_sd}_corr${bg_corr}" \
    "snp_pleiotropy_sd=${snp_sd}" \
    "background_corr=${bg_corr}"

done

echo ""
echo "Done."
echo ""
echo "Investigation E: EBMF clustering parameters:"
ls -1 investigation-univariate-simulations_E_*.html

echo ""
echo "Investigation V: program validation parameters:"
ls -1 investigation-univariate-simulations_V_*.html

echo ""
echo "Investigation G: generative realism:"
ls -1 investigation-univariate-simulations_G_*.html
