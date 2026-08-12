#!/usr/bin/env bash
# Knit investigation-bivariate-analysis.Rmd across the parameter grid and
# save HTML outputs with descriptive filenames for side-by-side comparison.
#
# Grid factors:
#   api:                   production | local
#   association_source:    coloc | full
#   min_snp_signals:       1 | 3 | 5
#   louvain_gamma:         2 | 3
#   compress_method:       none | asinh
#   similarity_threshold:  0 | 0.5
#   collapse_gene_tissue:  FALSE | TRUE
#
# Usage:
#   ./scripts/run_bivariate_param_grid.sh
#   ./scripts/run_bivariate_param_grid.sh --apis production --assocs coloc
#   ./scripts/run_bivariate_param_grid.sh --dry-run
#
# Optional filters (space-separated lists):
#   --apis --assocs --min-snps --gammas --compress --sim-thresholds --collapse-gene-tissue
#   --out-dir DIR   (default: vignettes/bivariate-param-grid)
#   --continue      skip runs whose HTML already exists
#   --dry-run       print planned renders only

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RMD="${ROOT}/vignettes/investigation-bivariate-analysis.Rmd"
OUT_DIR="${ROOT}/vignettes/bivariate-param-grid"
TMP_DIR="${OUT_DIR}/.knit-tmp"

APIS=(production local)
APIS=(local)
# ASSOCS=(coloc full)
ASSOCS=(coloc)
MIN_SNPS=(5)
GAMMAS=(2)
COMPRESS=(none)
SIM_THRESHOLDS=(0.5)
# COLLAPSE_GENE_TISSUE=(FALSE TRUE)
COLLAPSE_GENE_TISSUE=(FALSE TRUE)
CONTINUE=0
DRY_RUN=0

usage() {
  sed -n '2,26p' "$0" | sed 's/^# \?//'
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --apis)
      shift
      APIS=()
      while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do APIS+=("$1"); shift; done
      ;;
    --assocs)
      shift
      ASSOCS=()
      while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do ASSOCS+=("$1"); shift; done
      ;;
    --min-snps)
      shift
      MIN_SNPS=()
      while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do MIN_SNPS+=("$1"); shift; done
      ;;
    --gammas)
      shift
      GAMMAS=()
      while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do GAMMAS+=("$1"); shift; done
      ;;
    --compress)
      shift
      COMPRESS=()
      while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do COMPRESS+=("$1"); shift; done
      ;;
    --sim-thresholds)
      shift
      SIM_THRESHOLDS=()
      while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do SIM_THRESHOLDS+=("$1"); shift; done
      ;;
    --collapse-gene-tissue)
      shift
      COLLAPSE_GENE_TISSUE=()
      while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do COLLAPSE_GENE_TISSUE+=("$1"); shift; done
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --continue)
      CONTINUE=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

n_total=$(( ${#APIS[@]} * ${#ASSOCS[@]} * ${#MIN_SNPS[@]} * ${#GAMMAS[@]} * ${#COMPRESS[@]} * ${#SIM_THRESHOLDS[@]} * ${#COLLAPSE_GENE_TISSUE[@]} ))
n_done=0
n_skip=0
n_fail=0

echo "Root: ${ROOT}"
echo "Output: ${OUT_DIR}"
echo "Planned runs: ${n_total}"
echo

for api in "${APIS[@]}"; do
  for assoc in "${ASSOCS[@]}"; do
    for min_snp in "${MIN_SNPS[@]}"; do
      for gamma in "${GAMMAS[@]}"; do
        for compress in "${COMPRESS[@]}"; do
          for sim in "${SIM_THRESHOLDS[@]}"; do
            for collapse in "${COLLAPSE_GENE_TISSUE[@]}"; do
              tag="api-${api}_assoc-${assoc}_minsnp-${min_snp}_gamma-${gamma}_compress-${compress}_sim-${sim}_collapse-${collapse}"
              out_html="${OUT_DIR}/${tag}.html"
              knit_html="${TMP_DIR}/${tag}.html"

              if [[ "${CONTINUE}" -eq 1 && -f "${out_html}" ]]; then
                echo "[skip] ${tag}"
                n_skip=$((n_skip + 1))
                continue
              fi

              echo "[run ] ${tag}"
              if [[ "${DRY_RUN}" -eq 1 ]]; then
                continue
              fi

              if Rscript --vanilla -e "
                setwd('${ROOT}')
                rmarkdown::render(
                  input = '${RMD}',
                  output_file = '${tag}.html',
                  output_dir = '${TMP_DIR}',
                  intermediates_dir = '${TMP_DIR}/${tag}-int',
                  quiet = TRUE,
                  params = list(
                    api = '${api}',
                    association_source = '${assoc}',
                    min_snp_signals = as.integer(${min_snp}),
                    louvain_gamma = as.numeric(${gamma}),
                    compress_method = '${compress}',
                    similarity_threshold = as.numeric(${sim}),
                    collapse_gene_tissue = as.logical('${collapse}')
                  )
                )
              "; then
                mv -f "${knit_html}" "${out_html}"
                n_done=$((n_done + 1))
                echo "[ok  ] ${out_html}"
              else
                n_fail=$((n_fail + 1))
                echo "[fail] ${tag}" >&2
              fi
            done
          done
        done
      done
    done
  done
done

echo
echo "Finished. ok=${n_done} skip=${n_skip} fail=${n_fail} planned=${n_total}"
if [[ "${n_fail}" -gt 0 ]]; then
  exit 1
fi
