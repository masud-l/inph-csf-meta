# =============================================================================
# Step 02 - Random-effects SMD meta-analysis for all biomarkers
# Script: 02_smd_metaanalysis.R
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
# Pool the standardised mean difference (Hedges g) for every biomarker x
# comparison (iNPH vs healthy controls; iNPH vs Alzheimer disease), quantify
# heterogeneity, derive exploratory AUCs, and run leave-one-out analyses for
# the principal biomarkers.
#
# Methods (locked settings in R/00_setup.R)
# -------
# - Effect size: Hedges g via metafor::escalc(measure = "SMD").
#   Negative SMD = lower CSF concentration in iNPH than in the comparator.
# - Model: random effects, REML estimator (RE_METHOD = "REML").
# - Confidence intervals: normal (z) distribution (USE_KNHA = FALSE),
#   matching the original Stata analysis.
# - Pooling threshold: k >= 3 studies (MIN_K); k < 3 reported descriptively.
# - Heterogeneity: I^2 and tau^2.
# - Derived AUC (exploratory): AUC = pnorm(|SMD| / sqrt(2)), assuming normal
#   distributions with equal variance; lower bound set to 0.5 when the SMD
#   CI spans zero.
# - Leave-one-out: AB42, NFL, PTAU181 in both comparisons.
#
# Inputs
# ------
# output/tables/clean_data.rds          (from Step 01)
# R/meta_smd.R
#
# Outputs
# -------
# output/tables/smd_pooled_estimates.csv / .xlsx
# output/tables/smd_results.rds
# output/tables/smd_derived_auc.csv
# output/tables/loo_<biomarker>_<comparison>.csv
# output/tables/smd_reproduction_DL_vs_REML.csv   [VERIFY: still produced]
#
# Run order
# ---------
# Requires Step 01.
# =============================================================================


source("R/00_setup.R")
source("R/load_data.R")
source("R/meta_smd.R")

clean <- readRDS(file.path(PATHS$tables, "clean_data.rds"))$smd

# ---- Run the primary analysis (estimator = RE_METHOD from 00_setup.R) -------
message("Running SMD meta-analysis (", RE_METHOD, ") ...")
res_hc <- run_smd_all(dplyr::filter(clean, comparison == "iNPH_vs_HC"))
res_ad <- run_smd_all(dplyr::filter(clean, comparison == "iNPH_vs_AD"))
# Prefix each element key with its comparison. 14 biomarkers are measured in
# BOTH comparisons (AB42, pTau181, tTau, NFL, ...); without unique keys
# c(res_hc, res_ad) produces duplicate list names, and list[[name]] then
# silently returns only the FIRST match - so the vs-AD plot for every shared
# biomarker would never render. Unique keys make the script-04 loop correct.
names(res_hc) <- paste0("iNPH_vs_HC::", names(res_hc))
names(res_ad) <- paste0("iNPH_vs_AD::", names(res_ad))
results <- c(res_hc, res_ad)

summary_tbl <- dplyr::bind_rows(
  smd_summary_table(res_hc),
  smd_summary_table(res_ad)
)
write.csv(summary_tbl,
          file.path(PATHS$tables, "smd_pooled_estimates.csv"),
          row.names = FALSE)
message("Pooled estimates written to smd_pooled_estimates.csv")
print(summary_tbl, row.names = FALSE, digits = 3)

# Study-level estimates for outcomes with k < MIN_K (not pooled).
desc_tbl <- smd_descriptive_table(c(res_hc, res_ad))
write.csv(desc_tbl,
          file.path(PATHS$tables, "smd_descriptive_studies.csv"),
          row.names = FALSE)
message("Study-level estimates for non-pooled outcomes written to ",
        "smd_descriptive_studies.csv (", nrow(desc_tbl), " rows)")

# ---- Reproducibility: DL vs REML side by side -------------------------------
# Stata's `metan` historically defaults to DerSimonian-Laird. metafor defaults
# to REML. Running both isolates whether any Stata/R mismatch is just the
# estimator (expected on high-I-squared pools) rather than a data problem.
message("\nEstimator comparison (DL vs REML) ...")
repro <- function(method) {
  r <- c(run_smd_all(dplyr::filter(clean, comparison == "iNPH_vs_HC"), method),
         run_smd_all(dplyr::filter(clean, comparison == "iNPH_vs_AD"), method))
  smd_summary_table(r) |>
    dplyr::transmute(comparison, biomarker, k,
                     !!paste0("est_", method)   := round(estimate, 3),
                     !!paste0("cilb_", method)  := round(ci_lb, 3),
                     !!paste0("ciub_", method)  := round(ci_ub, 3))
}
cmp <- dplyr::full_join(repro("DL"), repro("REML"),
                        by = c("comparison", "biomarker", "k")) |>
  dplyr::mutate(abs_diff = abs(est_DL - est_REML)) |>
  dplyr::arrange(dplyr::desc(abs_diff))
write.csv(cmp, file.path(PATHS$tables, "smd_reproduction_DL_vs_REML.csv"),
          row.names = FALSE)
message("Estimator comparison written to smd_reproduction_DL_vs_REML.csv")
message("Largest DL-vs-REML gap: ", round(max(cmp$abs_diff, na.rm = TRUE), 3),
        " (", cmp$biomarker[1], ", ", cmp$comparison[1], ")")

# ---- Leave-one-out for the dual-purpose markers -----------------------------
for (key in c("AB42", "NFL", "PTAU181")) 
  {
  for (r in results) {
    s <- r$summary
    if (!is.null(s) && toupper(s$biomarker) == key && !is.null(r$fit)) {
      loo <- smd_leave1out(r)
      fn  <- sprintf("loo_%s_%s.csv", key, s$comparison)
      write.csv(loo, file.path(PATHS$tables, fn), row.names = FALSE)
      message("Leave-one-out written to ", fn)
    }
  }
}

saveRDS(results, file.path(PATHS$tables, "smd_results.rds"))
message("\n02_smd_metaanalysis.R complete.")
