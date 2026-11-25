# =============================================================================
# Step 10 - Age moderator analyses for CSF AB42 (subgroup + meta-regression)
# Script: 10_age_metaregression.R
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
# Exploratory, study-level investigation of whether iNPH-group mean age
# explains between-study heterogeneity in the CSF AB42 SMD, in both
# comparisons, using two complementary analyses:
#   1. Subgroup analysis by decade band (60-69, 70-79, 80-89), each band
#      pooled separately, with a test for subgroup differences (Q_between).
#   2. Continuous meta-regression: rma(yi, vi, mods = ~ iNPH_age).
#
# Age variable
# ------------
# iNPH-group central age per study: reported mean; median where only the
# median was reported; sample-size-weighted mean where age was reported by
# subgroup. Complete-case analysis; missing age is never imputed.
#
# Reproducibility check
# ---------------------
# The iNPH vs Alzheimer disease subgroup analysis reproduces the original
# Stata analysis (60-69: 0.28; 70-79: 0.53; 80-89: 1.17; overall 0.57).
#
# Interpretation caveats
# ----------------------
# Ecological (study-level) association only; mean age is confounded with
# publication era, assay generation, and diagnostic criteria. Bands with few
# studies are not reliably poolable.
#
# Inputs
# ------
# output/tables/clean_data.rds   (from Step 01; requires iNPH_age, column R)
# R/meta_regression.R, R/theme_bist.R
#
# Outputs
# -------
# output/figures/age_subgroup_forest_AB42_vs_{AD,HC}.{pdf,tiff}
# output/figures/age_bubble_binned_AB42_vs_{AD,HC}.{pdf,tiff}
# output/figures/age_metareg_bubble_AB42_vs_{AD,HC}.{pdf,tiff}
#   (Step 11 composes Supplementary Figure S9 from the same plot builders)
# output/tables/age_subgroup_AB42_vs_{AD,HC}.csv
# output/tables/age_metareg_AB42_vs_{AD,HC}.csv
# output/tables/age_metareg_loo_AB42_vs_{AD,HC}.csv
#
# Run order
# ---------
# Requires Step 01 (re-run after adding iNPH_age to the source workbooks).
# =============================================================================
#
# age_metaregression  —  Age moderator analyses (AB42, both comparisons)
# ----------------------------------------------------------------------
# At the lead author's request, runs BOTH:
#   (1) age subgroup analysis (decade bands 60-69 / 70-79 / 80-89), reproducing
#       Dr. Nila's Stata analysis, and
#   (2) continuous-age meta-regression (the statistically sounder tool),
# for CSF AB42 in BOTH the iNPH vs AD and iNPH vs HC comparisons.
#
# Requires iNPH_age (column R / 18) on the AB42 sheet of:
#   - CSF_iNPH_vs_AD.xlsx   (you have added this)
#   - CSF_iNPH_vs_HC.xlsx   (must ALSO be added for the HC analysis to run)
# The script runs whichever comparison has age data and clearly reports if HC
# age is absent.
#
# Outputs (output/figures/, output/tables/):
#   bubble plots          age_metareg_bubble_AB42_vs_{AD,HC}.{pdf,tiff}
#   meta-reg results      age_metareg_AB42_vs_{AD,HC}.csv
#   meta-reg LOO          age_metareg_loo_AB42_vs_{AD,HC}.csv
#   subgroup results      age_subgroup_AB42_vs_{AD,HC}.csv
#
# PREREQUISITE: re-run script 01 first so clean_data.rds carries iNPH_age.
# =============================================================================

source("R/00_setup.R")
source("R/load_data.R")
source("R/meta_smd.R")
source("R/theme_bist.R")
source("R/meta_regression.R")
library(ggplot2)

clean_smd <- readRDS(file.path(PATHS$tables, "clean_data.rds"))$smd

# Run both analyses for one comparison; returns invisibly, prints + writes.
run_one_comparison <- function(cmp) {
  tag <- comparison_style(cmp)$tag                       # vs_AD / vs_HC
  sub_data <- clean_smd[clean_smd$biomarker == "AB42" &
                        clean_smd$comparison == cmp, ]
  n_age <- sum(!is.na(sub_data$iNPH_age))
  message("\n========================================================")
  message("AB42  ", comparison_style(cmp)$label)
  message("  ", nrow(sub_data), " studies, ", n_age, " with iNPH_age.")

  if (n_age < 3) {
    message("  >>> Too few studies with age (", n_age,
            ") - skipping this comparison.")
    if (cmp == "iNPH_vs_HC")
      message("      If you want the HC analysis, add iNPH_age (column R) to ",
              "the AB42 sheet of CSF_iNPH_vs_HC.xlsx and re-run script 01.")
    return(invisible(NULL))
  }

  # ---- (1) Subgroup analysis (decade bands) --------------------------------
  sg <- run_age_subgroup(clean_smd, "AB42", cmp)
  message("\n  -- Age subgroup analysis (decade bands) --")
  print(sg$bands, row.names = FALSE)
  message(sprintf("  Overall: %.2f (%.2f to %.2f)",
                  as.numeric(sg$overall$b), sg$overall$ci.lb, sg$overall$ci.ub))
  message(sprintf("  Test of group differences: Q=%.2f, df=%s, P=%.3f",
                  sg$qbet$QM, sg$qbet$df, sg$qbet$QMp))
  write.csv(sg$bands,
            file.path(PATHS$tables, sprintf("age_subgroup_AB42_%s.csv", tag)),
            row.names = FALSE)

  # subgroup forest plot (BIST theme), for Dr. Nila's analysis
  sgf <- subgroup_forest(sg)
  save_bist(sgf, sprintf("age_subgroup_forest_AB42_%s", tag),
            width = 6, height = max(5, 2 + 0.34 * (sg$data |> nrow())))
  message("  saved subgroup forest: age_subgroup_forest_AB42_", tag)

  # binned bubble plot (Stata-style, companion to the subgroup analysis)
  bbin <- bubble_plot_binned(sg)
  save_bist(bbin, sprintf("age_bubble_binned_AB42_%s", tag),
            width = 7, height = 5.5)
  message("  saved binned bubble: age_bubble_binned_AB42_", tag)

  # ---- (2) Continuous meta-regression --------------------------------------
  reg <- run_age_metareg(clean_smd, "AB42", cmp)
  message("\n  -- Continuous age meta-regression --")
  message(sprintf("  k=%d (%d dropped); slope=%.4f/yr (95%% CI %.4f to %.4f); P=%.4f; R2=%s",
                  reg$k, reg$n_dropped, reg$slope, reg$ci[1], reg$ci[2], reg$p,
                  ifelse(is.na(reg$R2), "NA", sprintf("%.1f%%", reg$R2))))

  loo <- age_metareg_loo(reg)
  write.csv(loo,
            file.path(PATHS$tables, sprintf("age_metareg_loo_AB42_%s.csv", tag)),
            row.names = FALSE)

  res <- data.frame(
    biomarker = "AB42", comparison = cmp, k = reg$k, n_dropped = reg$n_dropped,
    slope_per_year = round(reg$slope, 4),
    ci_lb = round(reg$ci[1], 4), ci_ub = round(reg$ci[2], 4),
    p_value = signif(reg$p, 4),
    R2_percent = ifelse(is.na(reg$R2), NA, round(reg$R2, 1)))
  write.csv(res,
            file.path(PATHS$tables, sprintf("age_metareg_AB42_%s.csv", tag)),
            row.names = FALSE)

  p <- bubble_plot(reg)
  save_bist(p, sprintf("age_metareg_bubble_AB42_%s", tag),
            width = 7, height = 5.5)
  message("  saved meta-regression bubble + tables for ", tag)
  invisible(NULL)
}

for (cmp in c("iNPH_vs_AD", "iNPH_vs_HC")) run_one_comparison(cmp)

message("\n10_age_metaregression.R complete.")
