# =============================================================================
# Step 05 - Pooled-effect overview across all biomarkers (Figure 2 data)
# Script: 05_figure2_overview.R
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
# Summarise pooled SMDs (95% CI) for every biomarker with k >= 3 in both
# comparisons, as an overview dot plot and table, and export the underlying
# data. The published Figure 2 was drawn from this exported data.
# [VERIFY: Figure 2 final version prepared in GraphPad Prism]
#
# Inputs
# ------
# output/tables/smd_results.rds          (from Step 02)
# R/figure2_overview.R, R/theme_bist.R
#
# Outputs
# -------
# output/figures/Figure2_dotplot.{pdf,tiff}
# output/figures/Figure2_table.{pdf,tiff}
# output/tables/Figure2_overview_data.csv
#
# Run order
# ---------
# Requires Step 02.
# =============================================================================

source("R/00_setup.R")
source("R/theme_bist.R")
source("R/figure2_overview.R")

ov <- load_overview()
message("Loaded ", nrow(ov), " biomarker x comparison rows.")

# ---- Candidate A: dot plot --------------------------------------------------
message("\nBuilding Figure 2 candidate A (dot plot) ...")
n_max <- max(table(ov$comparison))
dot   <- figure2_dotplot(ov)
# No in-figure caption (it is the manuscript legend); height covers the two
# panels plus the thin significance-key strip beneath them.
save_bist(dot, "Figure2_dotplot",
          width = 8, height = 2.4 + 0.30 * n_max)

# ---- Candidate B: table -----------------------------------------------------
message("Building Figure 2 candidate B (table) ...")
tab <- figure2_table(ov)
save_bist(tab, "Figure2_table",
          width = 8.5, height = 2.4 + 0.32 * nrow(ov))

# ---- Underlying data --------------------------------------------------------
write.csv(figure2_overview_data(ov),
          file.path(PATHS$tables, "Figure2_overview_data.csv"),
          row.names = FALSE)

message("\n05_figure2_overview.R complete. Compare:")
message("  output/figures/Figure2_dotplot.{pdf,tiff}")
message("  output/figures/Figure2_table.{pdf,tiff}")
message("  output/tables/Figure2_overview_data.csv")