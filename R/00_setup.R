# =============================================================================
# 00_setup.R  —  Packages and global options
# -----------------------------------------------------------------------------
# Run this ONCE per machine to install dependencies, then `source()` it at the
# top of every analysis script. It is safe to source repeatedly: it installs
# only what is missing.
# =============================================================================

# ---- Required packages ------------------------------------------------------
.required <- c(
  "readxl",     # read the .xlsx extraction sheets
  "metafor",    # SMD meta-analysis (escalc + rma)
  "mada",       # diagnostic test accuracy meta-analysis (Reitsma / SROC)
  "ggplot2",    # figures
  "dplyr",      # data wrangling
  "tidyr",      # data wrangling
  "stringr",    # string cleaning of sheet names
  "purrr",      # functional iteration over biomarkers
  "scales",     # axis / weight rescaling
  "patchwork",  # compose multi-panel figures
  "officer"     # build Word (.docx) tables (scripts/07)
)

.install_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(TRUE)
}
.install_missing(.required)

# Load for interactive use
suppressPackageStartupMessages({
  library(readxl);  library(metafor); library(mada)
  library(ggplot2); library(dplyr);   library(tidyr)
  library(stringr); library(purrr);   library(scales); library(patchwork)
})

# ---- Global options ---------------------------------------------------------
options(stringsAsFactors = FALSE)

# Project paths (relative to the .Rproj root)
PATHS <- list(
  data    = "data",
  figures = "output/figures",
  tables  = "output/tables"
)
for (p in PATHS) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

# ---- Analysis configuration -------------------------------------------------
# Random-effects estimator. Stata's older `metan` defaults to "DL"; modern
# practice (and metafor's default) is "REML". Set this to match whatever your
# Stata run used so the reproduction comparison is apples-to-apples.
RE_METHOD <- "REML"          # alternatives: "DL", "PM", "SJ"
MIN_K <- 3
# Prediction intervals (t distribution, k - 2 df) are computed for every pooled
# estimate but shown in the Word tables only when k >= PI_MIN_K; with fewer
# studies they are too wide to be informative.
PI_MIN_K <- 5
# Knapp-Hartung small-sample adjustment for the pooled CI. This analysis uses
# z-based CIs (USE_KNHA = FALSE) to match the original (Stata) analysis being
# reproduced. Do not change this without re-running scripts 02-07.
USE_KNHA  <- FALSE

# Reproducibility
set.seed(20260519)

message("Setup complete. RE_METHOD = ", RE_METHOD,
        " | Knapp-Hartung = ", USE_KNHA)
