# =============================================================================
# Step 08 - Funnel plots and Egger tests for small-study effects
# Script: 08_publication_bias.R
# =============================================================================
#
# Project
# -------
# Cerebrospinal fluid biomarkers distinguishing idiopathic normal pressure
# hydrocephalus (iNPH) from healthy controls and Alzheimer disease:
# a systematic review and meta-analysis (41 studies).
#
# Purpose
# -------
# Assess small-study effects for biomarker x comparison pairs with at least
# 10 studies (PB_MIN_K = 10): AB42, p-tau181, and total tau, both comparisons.
#
# Methods
# -------
# - Funnel plot of SMD against standard error.
# - PRIMARY: modified test of Pustejovsky & Rodgers (2019), predictor
#   sqrt(1/n1 + 1/n2), which avoids the artefactual SMD-SE correlation.
# - SENSITIVITY: conventional Egger test, metafor::regtest() (SE predictor).
# - Interpretation: with I^2 of 83-96%, asymmetry may reflect publication
#   bias, small-study effects, or heterogeneity-induced distortion; these are
#   not distinguishable and asymmetry is not attributed to publication bias
#   alone.
#
# Inputs
# ------
# output/tables/smd_results.rds          (from Step 02)
# R/funnel_pb.R, R/theme_bist.R
#
# Outputs
# -------
# output/figures/funnel_<biomarker>_<comparison>.{pdf,tiff}   (supplementary funnel plots)
# output/tables/pb_egger.csv
#
# Run order
# ---------
# Requires Step 02.
#
# publication_bias—  Funnel plots + Egger's test (small-study effects)
# --------------------------------------------------------------------
# Assesses publication bias for the biomarker x comparison analyses with
# k >= PB_MIN_K (10). For this dataset that is Abeta42, p-tau181 and total tau,
# in BOTH comparisons - up to six analyses.
#
# Outputs (output/figures/ and output/tables/):
#   funnel_<biomarker>_<comparison>.{pdf,tiff}  - one funnel per analysis
#   pb_egger.csv                                - Egger's test results, all
#                                                 biomarkers (k>=10 testable)
#
# PREREQUISITE:
# --------------
# run scripts 01-02 first - this reads smd_results.rds.
# =============================================================================

source("R/00_setup.R")
source("R/load_data.R")
source("R/meta_smd.R")
source("R/theme_bist.R")        # RASTER_DPI, biomarker_label, comparison_style
source("R/funnel_pb.R")

smd_results <- readRDS(file.path(PATHS$tables, "smd_results.rds"))

# Biomarkers to assess (raw names). Only those reaching PB_MIN_K in at least
# one comparison are worth listing; others are reported as "not assessable".
PB_BIOMARKERS <- c("AB42", "pTau181", "tTau")

# ---- Run funnels + Egger across the six biomarker x comparison cells --------
message("Publication-bias assessment (k >= ", PB_MIN_K, ") ...")
egger_rows <- list()

for (bm in PB_BIOMARKERS) {
  for (cmp in c("iNPH_vs_HC", "iNPH_vs_AD")) {
    key <- paste0(cmp, "::", bm)                 # keys are comparison-prefixed
    r   <- smd_results[[key]]
    if (is.null(r)) {
      message("  (", key, " not found - skipped)")
      next
    }
    bm_lab  <- biomarker_label(bm)
    cs      <- comparison_style(cmp)
    cmp_lab <- cs$label

    # Egger's test (records a non-testable row if k < PB_MIN_K).
    egger_rows[[key]] <- egger_test(r, bm_lab, cmp_lab)

    # Funnel plot (skipped internally if k < PB_MIN_K). The header uses the
    # comparison's accent colour, matching the forest and LOO figures.
    funnel_bist(r$fit, bm_lab, cmp_lab,
                file = sprintf("funnel_%s_%s", bm, cmp),
                header_color = cs$color)
  }
}

# ---- Egger's results table --------------------------------------------------
egger <- do.call(rbind, egger_rows)
for (v in c("egger_z", "limit_est", "limit_lb", "limit_ub", "mod_egger_z"))
  egger[[v]] <- round(egger[[v]], 3)
for (v in c("egger_p", "mod_egger_p"))
  egger[[v]] <- round(egger[[v]], 4)
# Primary = modified test (mod_egger_*); conventional Egger = sensitivity.
egger <- egger[, c("biomarker", "comparison", "k",
                   "mod_egger_z", "mod_egger_p",
                   "egger_z", "egger_p",
                   "limit_est", "limit_lb", "limit_ub", "testable")]
write.csv(egger, file.path(PATHS$tables, "pb_egger.csv"), row.names = FALSE)

message("\nFunnel-asymmetry tests (primary: modified test; sensitivity: conventional Egger):")
print(egger[, c("biomarker", "comparison", "k", "mod_egger_p", "egger_p",
                "limit_est", "testable")], row.names = FALSE)

n_test <- sum(egger$testable)
message("\n08_publication_bias.R complete: ", n_test,
        " analyses had k >= ", PB_MIN_K, " and were assessed.")
message("  funnels  -> output/figures/funnel_*.{pdf,tiff}")
message("  Egger's  -> output/tables/pb_egger.csv")
