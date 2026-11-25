# =============================================================================
# Step 06 - Preview sheets of supplementary forest panels (optional)
# Script: 06_supplementary_preview.R
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
# Tile compact forest panels into 2x2 preview sheets for a quick visual check
# of every pooled and descriptive biomarker.
#
# Status
# ------
# Superseded for the submission: Figure 3 and Supplementary Figures S1-S9 are
# composed by Step 11 (scripts/11_composite_figures.R). The preview sheets are
# unnumbered and their grouping does not follow the manuscript; they are not
# required to regenerate the submission figures.
#
# Inputs
# ------
# output/tables/smd_results.rds          (from Step 02)
# R/forest_panel.R, R/theme_bist.R
#
# Outputs
# -------
# output/figures/preview_main_vs_HC.{pdf,tiff}
# output/figures/preview_main_vs_AD.{pdf,tiff}
# output/figures/preview_other_vs_HC.{pdf,tiff}
# output/figures/preview_other_vs_AD.{pdf,tiff}
# output/figures/preview_descriptive_<a-e>.{pdf,tiff}
#
# Run order
# ---------
# Optional. Requires Step 02.
# =============================================================================

source("R/00_setup.R")
source("R/load_data.R")
source("R/meta_smd.R")
source("R/theme_bist.R")
source("R/forest_bist.R")     # provides .forest_frame() / .study_year()
source("R/forest_panel.R")

smd_results <- readRDS(file.path(PATHS$tables, "smd_results.rds"))

# ---- Index the results by comparison + biomarker ----------------------------
# Keys in smd_results are "iNPH_vs_HC::AB42" etc. Build a small lookup.
idx <- lapply(smd_results, function(r) {
  s <- r$summary
  data.frame(biomarker = s$biomarker, comparison = s$comparison,
             k = s$k, pooled = !is.null(r$fit), stringsAsFactors = FALSE)
})
idx <- do.call(rbind, Map(function(d, key) { d$key <- key; d },
                          idx, names(smd_results)))

# The four main biomarkers, by their raw (pre-label) names in smd_results.
MAIN_BIOMARKERS <- c("AB42", "NFL", "pTau181", "tTau")

# Helper: pull the result objects for a set of keys, in a given order.
get_results <- function(keys) smd_results[keys]

# Helper: render a group of results across as many 2x2 sheets as needed.
# stem = file stem; title_body = sheet title.
render_group <- function(results, stem, title_body, caption) {
  if (length(results) == 0) {
    message("  (", stem, ": no biomarkers in this group - skipped)")
    return(invisible())
  }
  chunks <- split(results, ceiling(seq_along(results) / 4))
  multi  <- length(chunks) > 1
  for (i in seq_along(chunks)) {
    sheet_stem <- if (multi) paste0(stem, "_", letters[i]) else stem
    full_title <- paste0("Preview. ", title_body,
                         if (multi) paste0("  (sheet ", i, " of ",
                                           length(chunks), ")") else "")
    fig <- assemble_supp_figure(chunks[[i]], title = full_title,
                                caption = caption)
    n_panel <- length(chunks[[i]])
    save_bist(fig, sheet_stem,
              width = 11, height = 3.0 + 3.1 * ceiling(n_panel / 2))
  }
}

# ---- Shared caption text ----------------------------------------------------
cap_pooled <- paste0(
  "Each panel is a random-effects (REML) meta-analysis of one CSF biomarker. ",
  "Squares are individual study standardized mean differences (Hedges' g), ",
  "sized by study weight; horizontal lines are 95% CIs; the diamond is the ",
  "pooled estimate. A negative SMD indicates a lower CSF concentration in ",
  "iNPH. I-squared quantifies between-study heterogeneity. ",
  "HC indicates healthy controls; AD, Alzheimer disease."
)
cap_descr <- paste0(
  "Each panel shows the individual study standardized mean differences ",
  "(Hedges' g, 95% CI) for one CSF biomarker reported in fewer than ",
  MIN_K, " studies. These biomarkers were NOT meta-analysed - too few ",
  "studies for a stable random-effects estimate - and are shown ",
  "descriptively, without a pooled diamond. HC indicates healthy controls; ",
  "AD, Alzheimer disease."
)

# ---- Select keys for each preview sheet -------------------------------------
keys_for <- function(comparison, main, pooled) {
  sel <- idx[idx$comparison == comparison &
             (idx$biomarker %in% MAIN_BIOMARKERS) == main &
             idx$pooled == pooled, ]
  # order by descending |effect| via the summary estimates
  ord <- order(vapply(smd_results[sel$key],
                      function(r) -abs(r$summary$estimate), numeric(1)))
  sel$key[ord]
}

message("Rendering supplementary preview sheets ...")

# Preview - main biomarkers, HC
render_group(get_results(keys_for("iNPH_vs_HC", main = TRUE,  pooled = TRUE)),
             "preview_main_vs_HC",
             "Forest plots of the principal CSF biomarkers: iNPH vs healthy controls",
             cap_pooled)

# Preview - main biomarkers, AD
render_group(get_results(keys_for("iNPH_vs_AD", main = TRUE,  pooled = TRUE)),
             "preview_main_vs_AD",
             "Forest plots of the principal CSF biomarkers: iNPH vs Alzheimer disease",
             cap_pooled)

# Preview - other pooled biomarkers, HC
render_group(get_results(keys_for("iNPH_vs_HC", main = FALSE, pooled = TRUE)),
             "preview_other_vs_HC",
             "Forest plots of additional CSF biomarkers: iNPH vs healthy controls",
             cap_pooled)

# Preview - other pooled biomarkers, AD
render_group(get_results(keys_for("iNPH_vs_AD", main = FALSE, pooled = TRUE)),
             "preview_other_vs_AD",
             "Forest plots of additional CSF biomarkers: iNPH vs Alzheimer disease",
             cap_pooled)

# Preview - descriptive (not pooled, k < MIN_K), both comparisons together
descr_keys <- idx$key[idx$pooled == FALSE]
descr_keys <- descr_keys[order(idx$comparison[match(descr_keys, idx$key)],
                               idx$biomarker[match(descr_keys, idx$key)])]
render_group(get_results(descr_keys),
             "preview_descriptive",
             "CSF biomarkers reported in too few studies to meta-analyse",
             cap_descr)

message("\n06_supplementary_preview.R complete. See output/figures/preview_*.")
