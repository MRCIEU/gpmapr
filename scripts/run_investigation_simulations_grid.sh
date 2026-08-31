#!/usr/bin/env bash
# Render vignettes/investigation-univariate-simulations.Rmd across a grid of
# parameter settings, saving each HTML under a filename that encodes the run.
#
# Grid (edit the arrays below to change what gets tested):
#   ebmf_lfsr_threshold      0.02 / 0.05 / 0.10 / 0.20
#   ebmf_magnitude_threshold 0.10 / 0.25 / 0.50 / 0.75
#   compress_scale           none (no compression) / 2 / 5  (asinh scale;
#                            "none" maps to compress_method = "none")
#   min_module_size          1 / 2 / 5
#
# Output files land in vignettes/ as
#   investigation-univariate-simulations_lfsr<lfsr>_mag<magnitude>_cs<scale>_mms<min_module_size>.html
#
# Requirements: Rscript on PATH with the gpmapr deps (dplyr, ggplot2, knitr,
# rmarkdown) installed. Each run re-simulates the data and refits EBMF
# (including stability), so the full 4x4x3x3 grid can take a while.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIGNETTES_DIR="$(cd "${SCRIPT_DIR}/../vignettes" && pwd)"
RMD="investigation-univariate-simulations.Rmd"

lfsrs=(0.02 0.05 0.10 0.20)
magnitudes=(0.10 0.25 0.50 0.75)
compress_scales=(none 2 5)
min_module_sizes=(1 2 5)

render_one() {
  local lfsr=$1 mag=$2 cs=$3 mms=$4
  local method="asinh"
  local scale="$cs"
  local cslabel="$cs"
  if [ "$cs" = "none" ]; then
    method="none"
    scale=1
    cslabel="none"
  fi
  local out="investigation-univariate-simulations_lfsr${lfsr}_mag${mag}_cs${cslabel}_mms${mms}.html"
  echo ">>> Rendering ebmf_lfsr_threshold=${lfsr} ebmf_magnitude_threshold=${mag} compress_method=${method} compress_scale=${scale} min_module_size=${mms}"

  rm -f "${RMD%.Rmd}.html"

  Rscript -e "rmarkdown::render(
    '${RMD}',
    params = list(
      ebmf_lfsr_threshold = ${lfsr},
      ebmf_magnitude_threshold = ${mag},
      compress_method = '${method}',
      compress_scale = ${scale},
      min_module_size = ${mms}
    ),
    quiet = TRUE
  )"

  mv -f "${RMD%.Rmd}.html" "${out}"
  echo "    -> ${out}"
}

cd "${VIGNETTES_DIR}"

for lfsr in "${lfsrs[@]}"; do
  for mag in "${magnitudes[@]}"; do
    for cs in "${compress_scales[@]}"; do
      for mms in "${min_module_sizes[@]}"; do
        render_one "${lfsr}" "${mag}" "${cs}" "${mms}"
      done
    done
  done
done

echo ""
echo "Done. ${#lfsrs[@]}x${#magnitudes[@]}x${#compress_scales[@]}x${#min_module_sizes[@]} runs:"
ls -1 investigation-univariate-simulations_lfsr*.html