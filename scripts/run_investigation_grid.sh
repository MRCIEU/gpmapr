#!/usr/bin/env bash
# Render vignettes/investigation-univariate-analysis.Rmd across a grid of
# parameter settings, saving each HTML under a filename that encodes the run.
#
# Grid (edit the arrays below to change what gets tested):
#   louvain_gamma           2 vs 3
#   min_snp_signals         3 vs 5
#   similarity_threshold    0 vs 0.5
#   associations            coloc vs full  (coloc = colocalisation rows only,
#                            full = all /associations-full rows)
#
# Output files land in vignettes/ as
#   investigation-univariate_g<gamma>_mss<min_snp_signals>_thr<threshold>_<associations>.html
#
# Requirements: Rscript on PATH with the gpmapr deps (dplyr, tidyr, knitr,
# rmarkdown) installed. Each run makes live API calls, so the full 16-combo
# grid can take a while.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIGNETTES_DIR="$(cd "${SCRIPT_DIR}/../vignettes" && pwd)"
RMD="investigation-univariate-analysis.Rmd"

gammas=(2 3)
min_snp_signals=(3 5)
similarity_thresholds=(0 0.2 0.5)
associations=(coloc)

render_one() {
  local gamma=$1 mss=$2 thr=$3 assoc=$4
  local out="investigation-univariate_g${gamma}_mss${mss}_thr${thr}_${assoc}.html"
  echo ">>> Rendering gamma=${gamma} min_snp_signals=${mss} similarity_threshold=${thr} associations=${assoc}"

  rm -f "${RMD%.Rmd}.html"

  Rscript -e "rmarkdown::render(
    '${RMD}',
    params = list(
      louvain_gamma = ${gamma},
      min_snp_signals = ${mss},
      similarity_threshold = ${thr},
      associations = '${assoc}'
    ),
    quiet = TRUE
  )"

  mv -f "${RMD%.Rmd}.html" "${out}"
  echo "    -> ${out}"
}

cd "${VIGNETTES_DIR}"

for gamma in "${gammas[@]}"; do
  for mss in "${min_snp_signals[@]}"; do
    for thr in "${similarity_thresholds[@]}"; do
      for assoc in "${associations[@]}"; do
        render_one "${gamma}" "${mss}" "${thr}" "${assoc}"
      done
    done
  done
done

echo ""
echo "Done. ${#gammas[@]}x${#min_snp_signals[@]}x${#similarity_thresholds[@]}x${#associations[@]} runs:"
ls -1 investigation-univariate_g*.html
