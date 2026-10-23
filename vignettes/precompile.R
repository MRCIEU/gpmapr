# Precompile the vignettes, which need the live GPMap API.
#
# CRAN builds vignettes without network access, so each vignette is written as
# <name>.Rmd.orig and knitted here, against the production API, into
# <name>.Rmd with its results (and figures in figures/) baked in. CRAN then
# rebuilds the static .Rmd without running any code.
#
# Run from the package root before every release, and commit the regenerated
# .Rmd files and figures:
#   Rscript vignettes/precompile.R

if (basename(getwd()) != "vignettes") {
  setwd("vignettes")
}
pkgload::load_all("..", quiet = TRUE)
options(gpmap_url = "https://gpmap.opengwas.io/api")

for (name in c("gpmapr", "tutorial-case-study", "tutorial-gwas-upload")) {
  message("Knitting ", name)
  knitr::opts_chunk$set(fig.path = paste0("figures/", name, "-"), dpi = 96)
  # knitr::knit() doesn't set `params` as rmarkdown::render() does: use the vignette's defaults.
  params <- lapply(knitr::knit_params(readLines(paste0(name, ".Rmd.orig"))), function(param) param$value)
  knitr::knit(paste0(name, ".Rmd.orig"), paste0(name, ".Rmd"), quiet = TRUE)
  knitr::opts_chunk$restore()
}
