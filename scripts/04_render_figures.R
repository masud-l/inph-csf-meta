# =============================================================================
# Step 04 - Render individual forest plots and SROC curves
# Script: 04_render_figures.R
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
# Render one BIST-style forest plot per biomarker x comparison, diagnostic
# accuracy forest plots, and SROC curves. The multi-panel submission figures
# (Figure 3, Supplementary Figures S1-S9) are composed by Step 11.
#
# Figure conventions (R/theme_bist.R)
# -------------------
# - Title = biomarker name; subtitle = "iNPH vs Healthy Controls" (blue) or
#   "iNPH vs Alzheimer Disease" (amber), non-bold.
# - k < 3: individual studies shown, no pooled estimate.
# - Height scales with k (minimum 3.5 in); width 9 in.
# - SROC: curve, 95% confidence and prediction regions computed analytically
#   from the Reitsma parameters (coef, vcov, Psi), back-transformed with plogis.
# - PDF + TIFF at 600 dpi.
#
# Inputs
# ------
# output/tables/smd_results.rds          (from Step 02)
# output/tables/diagnostic_results.rds   (from Step 03)
# R/forest_bist.R, R/forest_diagnostic.R, R/theme_bist.R
#
# Outputs
# -------
# output/figures/forest_<biomarker>_<comparison>.{pdf,tiff}  (Suppl. Figures S1-S6)
# output/figures/dta_forest_<biomarker>_<comparison>.{pdf,tiff}
# output/figures/sroc_<biomarker>.{pdf,tiff}                 (Figure 3A)
#
# Run order
# ---------
# Requires Steps 02 and 03.
# =============================================================================

source("R/00_setup.R")
source("R/load_data.R")
source("R/meta_smd.R")
source("R/meta_diagnostic.R")
source("R/theme_bist.R")
source("R/forest_bist.R")
source("R/forest_diagnostic.R")

# Build a file-name-safe biomarker tag, e.g. "AB42/AB40" -> "AB42_AB40".
bm_tag <- function(x) gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", x))

# ---- SMD forest plots -------------------------------------------------------
smd_results <- readRDS(file.path(PATHS$tables, "smd_results.rds"))

message("Rendering SMD forest plots ...")
for (nm in names(smd_results)) {
  r <- smd_results[[nm]]
  s <- r$summary
  if (is.null(s)) next
  # forest_bist() derives the biomarker title and the colour-coded comparison
  # subtitle automatically from r$summary - nothing to pass.
  p <- tryCatch(forest_bist(r, order_by = "effect"),
                error = function(e) {
                  message("  ! ", nm, ": ", conditionMessage(e)); NULL })
  if (is.null(p)) next
  # Height scales with study count; floor at 3.5" so k=2 plots (no pooled row)
  # still have enough breathing room for headers, labels, and axis.
  height <- max(3.5, 1.6 + 0.32 * (s$k + 1))             # scale to study count
  fn <- sprintf("forest_%s_%s",
                bm_tag(s$biomarker),
                comparison_style(s$comparison)$tag)
  save_bist(p, fn, width = 6.5, height = height)
}

# ---- Diagnostic figures -----------------------------------------------------
dr <- readRDS(file.path(PATHS$tables, "diagnostic_results.rds"))

message("\nRendering diagnostic forest plots ...")
# (Revision) Only the audited iNPH-vs-AD comparisons (Supplementary Table S5) are plotted.
# The single-study rows carried over from the original diagnostic file
# (Kapaki Abeta42 vs HC; Miyajima sAPP vs AD+HC) are not part of the
# submission and are no longer rendered.
DTA_COMPARISONS <- "iNPH_vs_AD"
for (bm in unique(dr$per_study$biomarker)) {
  for (cmp in intersect(unique(dr$per_study$comparison), DTA_COMPARISONS)) {
    ps <- dplyr::filter(dr$per_study, biomarker == bm, comparison == cmp)
    if (nrow(ps) == 0) next
    # Attach a pooled diamond only where a bivariate model was fitted (vs AD).
    pooled <- if (cmp == "iNPH_vs_AD") dr$reitsma[[bm]]$summary else NULL
    p <- tryCatch(paired_dta_forest(ps, pooled, comparison = cmp),
                  error = function(e) {
                    message("  ! ", bm, "/", cmp, ": ", conditionMessage(e))
                    NULL })
    if (is.null(p)) next
    height <- 1.8 + 0.32 * (nrow(ps) + 1)
    fn <- sprintf("dta_forest_%s_%s", bm_tag(bm),
                  comparison_style(cmp)$tag)
    save_bist(p, fn, width = 9, height = height)
  }
}

message("\nRendering SROC plots ...")
for (bm in c("ab42", "p-tau181")) {
  r <- dr$reitsma[[bm]]
  if (!is.null(r) && isTRUE(r$boundary))
    message("  ", bm, ": no SROC curve (bivariate correlation at boundary; ",
            "univariate estimates reported).")
  if (is.null(r) || is.null(r$fit)) next
  sroc_bist(r$fit, dta = r$data,
            title = biomarker_label(bm),
            auc   = r$summary$AUC,
            comp  = "iNPH_vs_AD",
            file  = sprintf("sroc_%s", bm_tag(bm)))
}

message("\n04_render_figures.R complete. See output/figures/.")
