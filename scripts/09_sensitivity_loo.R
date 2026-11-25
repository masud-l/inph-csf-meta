# =============================================================================
# Step 09 - Leave-one-out sensitivity plots for AB42 and p-tau181
# Script: 09_sensitivity_loo.R
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
# Plot the pooled SMD re-estimated after omitting each study in turn, to show
# whether any single study determines the principal pooled estimates.
#
# Key behaviour
# -------------
# - Rows sorted by absolute shift from the full estimate (most influential
#   first); full-data estimate shown as the reference.
# - Panels feed Supplementary Figures S7-S8 (composed by Step 11).
#
# Inputs
# ------
# output/tables/loo_AB42_iNPH_vs_{HC,AD}.csv      (from Step 02)
# output/tables/loo_PTAU181_iNPH_vs_{HC,AD}.csv   (from Step 02)
#
# Outputs
# -------
# output/figures/loo_forest_AB42_vs_HC.{pdf,tiff}
# output/figures/loo_forest_AB42_vs_AD.{pdf,tiff}
# output/figures/loo_forest_pTau181_vs_HC.{pdf,tiff}
# output/figures/loo_forest_pTau181_vs_AD.{pdf,tiff}
#
# Run order
# ---------
# Requires Step 02.
# =============================================================================
#
# 09_sensitivity_loo.R  —  Leave-one-out sensitivity analysis
# -----------------------------------------------------------
# Generates BIST-style leave-one-out plots with:
#   - study names on the y-axis
#   - pooled SMD and 95% CI as a right-side text column
#   - full-analysis reference estimate shown as a diamond and dotted line
#
# Outputs (Step 11 composes the panels into the submission figures):
#   Supplementary Figure S7A: CSF Aβ42     — iNPH vs healthy controls
#   Supplementary Figure S7B: CSF Aβ42     — iNPH vs Alzheimer disease
#   Supplementary Figure S8A: CSF p-tau181 — iNPH vs healthy controls
#   Supplementary Figure S8B: CSF p-tau181 — iNPH vs Alzheimer disease
#
# PREREQUISITE:
# --------------
# Run script 02 first so that loo_*.csv and smd_results.rds are available.
# =============================================================================

source("R/00_setup.R")
source("R/theme_bist.R")      # save_bist, BIST colours, biomarker_label

library(ggplot2)


# ---- LOO panel builder ------------------------------------------------------

build_loo_panel <- function(loo, biomarker_lbl, comparison_lbl, comp_color,
                            ref_est, ref_lb, ref_ub, panel_tag) {
  
  # Keep only valid estimates
  d <- loo[!is.na(loo$estimate), , drop = FALSE]
  
  # Order by absolute shift from the full pooled estimate
  d$shift <- abs(d$estimate - ref_est)
  d <- d[order(-d$shift), ]
  
  # Clean study labels
  d$omitted <- gsub("Jurjevi´c", "Jurjevic", d$omitted)
  d$omitted <- gsub("Alzheimer's disease", "Alzheimer disease", d$omitted)
  
  # Full-analysis reference row
  ref_row <- data.frame(
    omitted = "Full analysis",
    estimate = ref_est,
    ci_lb = ref_lb,
    ci_ub = ref_ub,
    I2 = NA_real_,
    shift = 0
  )
  
  # Put reference row first so it appears at the bottom after reversing levels
  plot_d <- rbind(ref_row, d)
  
  # Avoid duplicate factor-level errors if study labels repeat
  plot_d$row_id <- make.unique(plot_d$omitted, sep = " ")
  plot_d$display_label <- plot_d$omitted
  
  # Full-analysis row indicator
  plot_d$is_ref <- plot_d$omitted == "Full analysis"
  
  # Text style for SMD column
  plot_d$text_face <- ifelse(plot_d$is_ref, "bold", "plain")
  
  # SMD text column
  plot_d$estimate_label <- sprintf(
    "%.2f (%.2f to %.2f)",
    plot_d$estimate,
    plot_d$ci_lb,
    plot_d$ci_ub
  )
  
  # Row order: full analysis at bottom; largest shift at top
  plot_d$label <- factor(plot_d$row_id, levels = rev(plot_d$row_id))
  
  # X-axis limits based only on effect estimates and CIs
  x_rng <- range(c(plot_d$ci_lb, plot_d$ci_ub), na.rm = TRUE)
  x_pad <- diff(x_rng) * 0.10
  
  if (!is.finite(x_pad) || x_pad == 0) {
    x_pad <- 0.10
  }
  
  x_lim <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)
  
  # Place SMD text outside the plotting panel.
  # This is intentionally NOT included in scale_x_continuous limits.
  x_text <- x_lim[2] + diff(x_lim) * 0.08
  
  # Use breaks based only on the true plotting range
  x_breaks <- scales::pretty_breaks(n = 5)(x_lim)
  x_breaks <- x_breaks[x_breaks >= x_lim[1] & x_breaks <= x_lim[2]]
  
  # Show null line only if SMD = 0 is inside the true plotting range
  show_null <- x_lim[1] <= 0 && x_lim[2] >= 0
  
  # Position for SMD column header
  header_y <- length(levels(plot_d$label)) + 0.65
  
  ggplot(plot_d, aes(y = label)) +
    
    # Full-analysis reference line
    geom_vline(
      xintercept = ref_est,
      colour = comp_color,
      linetype = "dotted",
      linewidth = 0.75
    ) +
    
    # Null-effect line, only when visible within true x-axis range
    {if (show_null)
      geom_vline(
        xintercept = 0,
        colour = "grey70",
        linewidth = 0.40
      )
    } +
    
    # Confidence intervals
    geom_errorbar(
      aes(
        x = estimate,
        xmin = ci_lb,
        xmax = ci_ub,
        colour = is_ref,
        linewidth = is_ref
      ),
      orientation = "y",
      width = 0.22
    ) +
    
    # Point estimates
    geom_point(
      aes(
        x = estimate,
        shape = is_ref,
        colour = is_ref,
        size = is_ref
      )
    ) +
    
    # Right-side SMD column outside the plotting panel
    geom_text(
      aes(
        x = x_text,
        label = estimate_label,
        fontface = text_face
      ),
      hjust = 0,
      size = 3.1,
      colour = BIST_BLACK
    ) +
    
    # Header for SMD column
    annotate(
      "text",
      x = x_text,
      y = header_y,
      label = "SMD (95% CI)",
      hjust = 0,
      size = 3.3,
      fontface = "bold",
      colour = BIST_BLACK
    ) +
    
    # Manual scales
    scale_colour_manual(
      values = c(`FALSE` = BIST_BLACK, `TRUE` = comp_color),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c(`FALSE` = 15, `TRUE` = 18),
      guide = "none"
    ) +
    scale_size_manual(
      values = c(`FALSE` = 2.4, `TRUE` = 4.0),
      guide = "none"
    ) +
    scale_linewidth_manual(
      values = c(`FALSE` = 0.55, `TRUE` = 1.15),
      guide = "none"
    ) +
    
    # X-axis based only on actual plot range
    scale_x_continuous(
      breaks = x_breaks
    ) +
    
    # Study labels
    scale_y_discrete(
      labels = setNames(plot_d$display_label, plot_d$row_id)
    ) +
    
    # Keep x-axis tight while allowing SMD text in the right margin
    coord_cartesian(
      xlim = x_lim,
      clip = "off"
    ) +
    
    # Titles and labels
    labs(
      title = paste0(panel_tag, ". ", biomarker_lbl),
      subtitle = comparison_lbl,
      x = "Pooled SMD with study omitted (Hedges g)",
      y = NULL
    ) +
    
    theme_bist() +
    theme(
      plot.title.position = "plot",
      plot.subtitle = element_text(
        size = 10,
        colour = comp_color,
        margin = margin(b = 10),
        hjust = 0
      ),
      plot.title = element_text(
        size = 12,
        face = "bold",
        colour = BIST_BLACK,
        margin = margin(b = 3),
        hjust = 0
      ),
      
      # Study names: increase this value if you want even larger labels
      axis.text.y = element_text(
        size = 12,
        colour = BIST_BLACK,
        hjust = 1,
        lineheight = 0.95
      ),
      
      axis.text.x = element_text(
        size = 9.0,
        colour = BIST_BLACK
      ),
      axis.title.x = element_text(
        size = 12,
        margin = margin(t = 8)
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        colour = "grey88",
        linewidth = 0.30
      ),
      
      # Right margin leaves space for SMD column
      plot.margin = margin(12, 125, 10, 18)
    )
}


# ---- Save one LOO panel as individual PDF/TIFF ------------------------------

save_loo_panel <- function(panel, loo, file) {
  
  k <- nrow(loo[!is.na(loo$estimate), , drop = FALSE])
  
  save_bist(
    panel,
    file,
    width  = 7.4,
    height = 3.0 + 0.32 * (k + 1)
  )
  
}


# ---- Load LOO CSVs and full-analysis pooled estimates -----------------------

smd_results <- readRDS(file.path(PATHS$tables, "smd_results.rds"))

get_ref <- function(bm_upper, cmp) {
  
  hits <- names(smd_results)[
    grepl(
      paste0(cmp, "::", bm_upper, "$"),
      names(smd_results),
      ignore.case = TRUE
    )
  ]
  
  if (!length(hits)) {
    return(NULL)
  }
  
  r <- smd_results[[hits[1]]]
  
  if (is.null(r$fit)) {
    return(NULL)
  }
  
  list(
    est = as.numeric(r$fit$b),
    lb  = as.numeric(r$fit$ci.lb),
    ub  = as.numeric(r$fit$ci.ub)
  )
}


read_loo <- function(key_upper, cmp) {
  
  fn <- file.path(
    PATHS$tables,
    sprintf("loo_%s_%s.csv", key_upper, cmp)
  )
  
  if (!file.exists(fn)) {
    message("  LOO file not found: ", fn, " - skipping")
    return(NULL)
  }
  
  read.csv(fn, stringsAsFactors = FALSE)
}


message("Building leave-one-out panels with SMD text column ...")


# =============================================================================
# Supplementary Figure S7 — Aβ42
# =============================================================================

bm <- "AB42"
bm_lbl <- biomarker_label(bm)

loo_hc <- read_loo("AB42", "iNPH_vs_HC")
loo_ad <- read_loo("AB42", "iNPH_vs_AD")

ref_hc <- get_ref("AB42", "iNPH_vs_HC")
ref_ad <- get_ref("AB42", "iNPH_vs_AD")

if (!is.null(loo_hc) && !is.null(ref_hc) &&
    !is.null(loo_ad) && !is.null(ref_ad)) {
  
  cs_hc <- comparison_style("iNPH_vs_HC")
  cs_ad <- comparison_style("iNPH_vs_AD")
  
  p_hc <- build_loo_panel(
    loo = loo_hc,
    biomarker_lbl = bm_lbl,
    comparison_lbl = cs_hc$label,
    comp_color = cs_hc$color,
    ref_est = ref_hc$est,
    ref_lb = ref_hc$lb,
    ref_ub = ref_hc$ub,
    panel_tag = "A"
  )
  
  p_ad <- build_loo_panel(
    loo = loo_ad,
    biomarker_lbl = bm_lbl,
    comparison_lbl = cs_ad$label,
    comp_color = cs_ad$color,
    ref_est = ref_ad$est,
    ref_lb = ref_ad$lb,
    ref_ub = ref_ad$ub,
    panel_tag = "B"
  )
  
  save_loo_panel(p_hc, loo_hc, "loo_forest_AB42_vs_HC")
  save_loo_panel(p_ad, loo_ad, "loo_forest_AB42_vs_AD")
  
} else {
  message("  Supplementary Figure S7: missing LOO data for Aβ42 - skipped.")
}


# =============================================================================
# Supplementary Figure S8 — p-tau181
# =============================================================================

bm <- "pTau181"
bm_lbl <- biomarker_label(bm)

loo_hc <- read_loo("PTAU181", "iNPH_vs_HC")
loo_ad <- read_loo("PTAU181", "iNPH_vs_AD")

if (is.null(loo_hc)) {
  loo_hc <- read_loo("pTau181", "iNPH_vs_HC")
}

if (is.null(loo_ad)) {
  loo_ad <- read_loo("pTau181", "iNPH_vs_AD")
}

ref_hc <- get_ref("pTau181", "iNPH_vs_HC")
ref_ad <- get_ref("pTau181", "iNPH_vs_AD")

if (!is.null(loo_hc) && !is.null(ref_hc) &&
    !is.null(loo_ad) && !is.null(ref_ad)) {
  
  cs_hc <- comparison_style("iNPH_vs_HC")
  cs_ad <- comparison_style("iNPH_vs_AD")
  
  p_hc <- build_loo_panel(
    loo = loo_hc,
    biomarker_lbl = bm_lbl,
    comparison_lbl = cs_hc$label,
    comp_color = cs_hc$color,
    ref_est = ref_hc$est,
    ref_lb = ref_hc$lb,
    ref_ub = ref_hc$ub,
    panel_tag = "A"
  )
  
  p_ad <- build_loo_panel(
    loo = loo_ad,
    biomarker_lbl = bm_lbl,
    comparison_lbl = cs_ad$label,
    comp_color = cs_ad$color,
    ref_est = ref_ad$est,
    ref_lb = ref_ad$lb,
    ref_ub = ref_ad$ub,
    panel_tag = "B"
  )
  
  save_loo_panel(p_hc, loo_hc, "loo_forest_pTau181_vs_HC")
  save_loo_panel(p_ad, loo_ad, "loo_forest_pTau181_vs_AD")
  
} else {
  message(
    "  Supplementary Figure S8: missing LOO data for p-tau181 - check that script 02 ",
    "was re-run after adding p-tau181 to the LOO loop."
  )
}


message("\n09_sensitivity_loo.R complete. Individual PDF/TIFF files saved:")
message("  loo_forest_AB42_vs_HC.{pdf,tiff}")
message("  loo_forest_AB42_vs_AD.{pdf,tiff}")
message("  loo_forest_pTau181_vs_HC.{pdf,tiff}")
message("  loo_forest_pTau181_vs_AD.{pdf,tiff}")