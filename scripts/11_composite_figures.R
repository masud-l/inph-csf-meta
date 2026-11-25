# =============================================================================
# Step 11 - Multi-panel manuscript and supplementary figures (patchwork)
# Script: 11_composite_figures.R
# =============================================================================
#
# Purpose
# -------
# Replaces the PowerPoint assembly step (Figures_v4.pptx). Every multi-panel
# figure is composed here from the same plot builders used in Steps 04, 09 and
# 10, so the figures always match the current results.
#
# One layout rule for every figure
# --------------------------------
# - Final width FINAL_WIDTH = 7.1 in (180 mm, double column), two columns.
# - Each panel is drawn at its usual design size (forest plots 6.5 in wide,
#   as in Step 04) and the whole composite is then scaled down uniformly -
#   the same thing the PPT did when 6.5-in plots were placed at 3.4 in.
# - Panel labels "A. Abeta42" in bold above each panel, same size everywhere.
#   Where a row pairs the HC and AD versions of one biomarker, the label sits
#   on the left panel only (as in the PPT).
# - Panels are top-aligned; shorter panels are padded, never stretched.
# - 600-dpi LZW TIFF at final size (for submission) + vector PDF at design
#   size (for proofing; it scales without loss).
# - Figure legends stay in the manuscript text, not in the image.
#
# Figure map (follows Figures_v4.pptx; content follows the revision rules,
# k >= 3 to pool, so some panels moved to Supplementary Figure S6):
#   Figure 3          A Abeta42 SROC | B p-tau181 sensitivity/specificity
#   Suppl. Figure S1  HC principal: A Abeta42, B p-tau181, C NfL, D sAPPa
#   Suppl. Figure S2  HC other pooled: A Total tau, B Abeta40
#   Suppl. Figure S3  Ratios, HC | AD per row: A p-tau181/Abeta42,
#                     B Abeta42/Abeta40, C Total tau/Abeta42
#   Suppl. Figure S4  AD principal: A p-tau181, B Abeta42, C NfL, D sAPPa
#   Suppl. Figure S5  AD other pooled: A Total tau, B Abeta40
#   Suppl. Figure S6  Fewer than 3 studies (descriptive), two files:
#                     S6a = HC, S6b = AD
#   Suppl. Figure S7  Leave-one-out, Abeta42 (HC | AD)
#   Suppl. Figure S8  Leave-one-out, p-tau181 (HC | AD)
#   Suppl. Figure S9  Age: A subgroup forest (HC | AD), B meta-regression (HC | AD)
# Output files are named Figure3_* and Supplementary_Figure_S<n>_*.
# Figure 2 is made in GraphPad from Effect_size_revised.xlsx (R alternative:
# Step 05).
#
# Inputs : output/tables/smd_results.rds, clean_data.rds,
#          diagnostic_results.rds, loo_*.csv   (Steps 01-03)
# Outputs: output/figures/composite/<name>.{tiff,pdf}
# Run    : after Steps 01-10.
# =============================================================================

source("R/00_setup.R")
source("R/load_data.R")
source("R/meta_smd.R")
source("R/meta_diagnostic.R")
source("R/theme_bist.R")
source("R/forest_bist.R")
source("R/forest_diagnostic.R")
source("R/meta_regression.R")
library(ggplot2)
library(patchwork)

# ---- Layout settings (the single place to change the look) -----------------
FINAL_WIDTH <- 7.1     # inches, final printed width of every figure
MAX_FINAL_H <- 9.6     # inches; a message is printed if a figure is taller
GUTTER      <- 0.35    # design inches between the two columns
ROW_GAP     <- 0.30    # design inches between stacked panels
LABEL_SIZE  <- 20      # design pt for "A. Abeta42" (~10-11 pt at final size)
LABEL_EXTRA <- 0.20    # design inches added to each panel for the larger label
MIN_H_SMD   <- 3.0     # design inches, minimum forest-plot height (k = 1-2)

OUT_DIR <- file.path(PATHS$figures, "composite")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

smd_results <- readRDS(file.path(PATHS$tables, "smd_results.rds"))
clean_smd   <- readRDS(file.path(PATHS$tables, "clean_data.rds"))$smd
dr          <- readRDS(file.path(PATHS$tables, "diagnostic_results.rds"))

# The leave-one-out panel builder lives in Step 09: load only its function
# definitions (build_loo_panel, get_ref, read_loo), without running the script.
local({
  ex <- parse("scripts/09_sensitivity_loo.R", encoding = "UTF-8")
  for (e in ex)
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))
      eval(e, envir = globalenv())
})

# ---- Label styling ----------------------------------------------------------
label_el <- element_text(size = LABEL_SIZE, face = "bold", colour = BIST_BLACK,
                         hjust = 0, margin = margin(b = 4))
# For single ggplot panels: the label starts at the left edge of the whole plot
# (axis title included), exactly where plot_annotation() puts the label of the
# patchwork-built panels. ggplot's default ("panel") would indent it to the
# start of the plotting area, past the y-axis title and tick labels. Captions
# follow the same rule so they line up under the label.
label_theme <- theme(plot.title = label_el, plot.title.position = "plot",
                     plot.caption.position = "plot")

# Annotation theme for patchwork-built panels (forest_bist, paired DTA forest,
# subgroup forest): the builder's own theme plus the larger panel label.
annot_theme <- function(cmp, sub_size = 11, small_caption = FALSE) {
  th <- theme_bist() +
    theme(plot.title    = label_el,
          plot.subtitle = element_text(colour = comparison_style(cmp)$color,
                                       size = sub_size))
  if (small_caption)
    th <- th + theme(plot.caption = element_text(hjust = 0, size = 8,
                                                 colour = BIST_GREY))
  th
}

# A panel = the plot embedded as one element + its design size (inches).
panel <- function(p, w, h) list(p = wrap_elements(full = p), w = w,
                                h = h + LABEL_EXTRA)

# ---- Panel builders ---------------------------------------------------------
find_res <- function(bm, cmp) {
  hit <- Filter(function(r) !is.null(r$summary) &&
                  identical(r$summary$biomarker, bm) &&
                  identical(r$summary$comparison, cmp), smd_results)
  if (!length(hit)) stop("No SMD result for ", bm, " / ", cmp)
  hit[[1]]
}

smd_panel <- function(bm, cmp, label) {
  r <- find_res(bm, cmp)
  p <- forest_bist(r, order_by = "effect") +
    plot_annotation(title = label, theme = annot_theme(cmp))
  k <- r$summary$k
  panel(p, 6.5, max(MIN_H_SMD, 1.6 + 0.32 * (k + 1)))
}

loo_panel <- function(bm_file, bm_ref, cmp, label) {
  loo <- read_loo(bm_file, cmp)
  if (is.null(loo)) loo <- read_loo(bm_ref, cmp)
  ref <- get_ref(bm_ref, cmp)
  if (is.null(loo) || is.null(ref)) stop("LOO inputs missing for ", bm_ref)
  cs <- comparison_style(cmp)
  p <- build_loo_panel(loo = loo, biomarker_lbl = biomarker_label(bm_ref),
                       comparison_lbl = cs$label, comp_color = cs$color,
                       ref_est = ref$est, ref_lb = ref$lb, ref_ub = ref$ub,
                       panel_tag = "") +
    labs(title = label) + label_theme
  k <- sum(!is.na(loo$estimate))
  panel(p, 7.4, 3.0 + 0.32 * (k + 1))
}

age_panels <- function(cmp, label_sub, label_reg) {
  sg  <- run_age_subgroup(clean_smd, "AB42", cmp)
  reg <- run_age_metareg(clean_smd, "AB42", cmp)
  p_sg <- subgroup_forest(sg) +
    plot_annotation(title = label_sub,
                    theme = annot_theme(cmp, sub_size = 10,
                                        small_caption = TRUE))
  p_rg <- bubble_plot(reg) + labs(title = label_reg) + label_theme
  list(sub = panel(p_sg, 6, max(5, 2 + 0.34 * nrow(sg$data))),
       reg = panel(p_rg, 7, 5.5))
}

# ---- Layout engines ---------------------------------------------------------
heights_of <- function(ps) vapply(ps, function(x) x$h, numeric(1))
widths_of  <- function(ps) vapply(ps, function(x) x$w, numeric(1))

# Pad a panel with empty space below it so it fills a cell of height h_cell.
pad_to <- function(pn, h_cell) {
  if (h_cell - pn$h < 1e-6) return(pn$p)
  wrap_plots(pn$p, plot_spacer(), ncol = 1, heights = c(pn$h, h_cell - pn$h))
}

# Stack panels vertically, padded at the bottom to total height H.
stack_col <- function(ps, H) {
  plots <- list(); hs <- numeric(0)
  for (i in seq_along(ps)) {
    if (i > 1) { plots <- c(plots, list(plot_spacer())); hs <- c(hs, ROW_GAP) }
    plots <- c(plots, list(ps[[i]]$p)); hs <- c(hs, ps[[i]]$h)
  }
  if (H - sum(hs) > 1e-6) {
    plots <- c(plots, list(plot_spacer())); hs <- c(hs, H - sum(hs))
  }
  wrap_plots(plots, ncol = 1, heights = hs)
}
col_height <- function(ps) sum(heights_of(ps)) + ROW_GAP * (length(ps) - 1)

# Two independent columns (Supplementary Figures S1, S2, S4, S5, S6).
compose_columns <- function(left, right) {
  H  <- max(col_height(left), col_height(right))
  wL <- max(widths_of(left)); wR <- max(widths_of(right))
  fig <- wrap_plots(stack_col(left, H), plot_spacer(), stack_col(right, H),
                    nrow = 1, widths = c(wL, GUTTER, wR))
  list(plot = fig, w = wL + GUTTER + wR, h = H)
}

# Rows of 1-2 panels side by side (Figure 3; Supplementary Figures S3,
# S7, S8, S9).
compose_rows <- function(rows) {
  nc   <- max(lengths(rows))
  colw <- vapply(seq_len(nc), function(j)
    max(vapply(rows, function(r) if (length(r) >= j) r[[j]]$w else 0,
               numeric(1))), numeric(1))
  plots <- list(); hs <- numeric(0)
  for (i in seq_along(rows)) {
    r  <- rows[[i]]
    rh <- max(heights_of(r))
    cells <- list(); ws <- numeric(0)
    for (j in seq_len(nc)) {
      if (j > 1) { cells <- c(cells, list(plot_spacer())); ws <- c(ws, GUTTER) }
      cells <- c(cells, list(if (j <= length(r)) pad_to(r[[j]], rh)
                             else plot_spacer()))
      ws <- c(ws, colw[j])
    }
    if (i > 1) { plots <- c(plots, list(plot_spacer())); hs <- c(hs, ROW_GAP) }
    plots <- c(plots, list(wrap_plots(cells, nrow = 1, widths = ws)))
    hs <- c(hs, rh)
  }
  list(plot = wrap_plots(plots, ncol = 1, heights = hs),
       w = sum(colw) + GUTTER * (nc - 1), h = sum(hs))
}

# ---- Save: TIFF at final size (uniform scaling) + vector PDF ----------------
save_composite <- function(comp, file) {
  s  <- FINAL_WIDTH / comp$w
  fh <- comp$h * s
  base <- file.path(OUT_DIR, file)
  ragg::agg_tiff(paste0(base, ".tiff"), width = FINAL_WIDTH, height = fh,
                 units = "in", res = RASTER_DPI, scaling = s,
                 compression = "lzw", background = "white")
  print(comp$plot)
  invisible(grDevices::dev.off())
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = comp$w, height = comp$h)
  print(comp$plot)
  invisible(grDevices::dev.off())
  message(sprintf("  saved: %s.tiff  (%.2f x %.2f in at %d dpi; scale %.2f)",
                  base, FINAL_WIDTH, fh, RASTER_DPI, s))
  if (fh > MAX_FINAL_H)
    message(sprintf("  ! %s is %.2f in tall (> %.1f in page limit)",
                    file, fh, MAX_FINAL_H))
}

safe_build <- function(name, expr) {
  tryCatch(expr, error = function(e)
    message("  ! ", name, " not built: ", conditionMessage(e)))
}

HC <- "iNPH_vs_HC"; AD <- "iNPH_vs_AD"
L  <- function(letter, bm) paste0(letter, ". ", biomarker_label(bm))

message("Building composite figures (final width ", FINAL_WIDTH, " in) ...")

# ---- Figure 3: diagnostic accuracy vs AD ------------------------------------
safe_build("Figure 3", {
  ab <- dr$reitsma[["ab42"]]
  p_sroc <- sroc_bist(ab$fit, dta = ab$data, title = "A. Aβ42",
                      auc = ab$summary$AUC, comp = AD, save = FALSE) +
    label_theme
  ps <- dplyr::filter(dr$per_study, biomarker == "p-tau181",
                      comparison == AD)
  p_dta <- paired_dta_forest(ps, dr$reitsma[["p-tau181"]]$summary,
                             comparison = AD) +
    plot_annotation(title = "B. p-tau181", theme = annot_theme(AD))
  save_composite(compose_rows(list(list(panel(p_sroc, 5.5, 5.5),
                                        panel(p_dta, 9, 5.5)))),
                 "Figure3_diagnostic_accuracy")
})

# ---- Supplementary Figures S1 / S4: principal biomarkers --------------------
safe_build("Supplementary Figure S1", save_composite(compose_columns(
  left  = list(smd_panel("AB42", HC, L("A", "AB42")),
               smd_panel("NFL", HC, L("C", "NFL"))),
  right = list(smd_panel("pTau181", HC, L("B", "pTau181")),
               smd_panel("sAPPalpha", HC, L("D", "sAPPalpha")))),
  "Supplementary_Figure_S1_HC_principal"))

safe_build("Supplementary Figure S4", save_composite(compose_columns(
  left  = list(smd_panel("pTau181", AD, L("A", "pTau181")),
               smd_panel("NFL", AD, L("C", "NFL"))),
  right = list(smd_panel("AB42", AD, L("B", "AB42")),
               smd_panel("sAPPalpha", AD, L("D", "sAPPalpha")))),
  "Supplementary_Figure_S4_AD_principal"))

# ---- Supplementary Figures S2 / S5: other pooled biomarkers -----------------
safe_build("Supplementary Figure S2", save_composite(compose_columns(
  left  = list(smd_panel("tTau", HC, L("A", "tTau"))),
  right = list(smd_panel("AB40", HC, L("B", "AB40")))),
  "Supplementary_Figure_S2_HC_other"))

safe_build("Supplementary Figure S5", save_composite(compose_columns(
  left  = list(smd_panel("tTau", AD, L("A", "tTau"))),
  right = list(smd_panel("AB40", AD, L("B", "AB40")))),
  "Supplementary_Figure_S5_AD_other"))

# ---- Supplementary Figure S3: ratios, HC (left) | AD (right) per row --------
safe_build("Supplementary Figure S3", save_composite(compose_rows(list(
  list(smd_panel("pTau/AB42", HC, L("A", "pTau/AB42")),
       smd_panel("pTau/AB42", AD, " ")),
  list(smd_panel("AB42/AB40", HC, L("B", "AB42/AB40")),
       smd_panel("AB42/AB40", AD, " ")),
  list(smd_panel("tTau/AB42", HC, L("C", "tTau/AB42")),
       smd_panel("tTau/AB42", AD, " ")))),
  "Supplementary_Figure_S3_ratios"))

# ---- Supplementary Figure S6: outcomes reported in fewer than MIN_K studies -
descr_keys <- function(cmp, exclude) {
  s <- do.call(rbind, lapply(smd_results, function(r) r$summary))
  s <- s[s$comparison == cmp & s$k < MIN_K & !(s$biomarker %in% exclude), ]
  s$biomarker[order(-s$k, s$biomarker)]
}
two_cols <- function(bms, cmp) {
  ps <- Map(function(bm, i) smd_panel(bm, cmp, L(LETTERS[i], bm)),
            bms, seq_along(bms))
  odd <- seq_along(ps) %% 2 == 1
  compose_columns(left = ps[odd], right = ps[!odd])
}
# HC ratios with k < 3 are already shown in Supplementary Figure S3, next to
# their AD panel.
safe_build("Supplementary Figure S6a", save_composite(
  two_cols(descr_keys(HC, exclude = c("pTau/AB42", "tTau/AB42")), HC),
  "Supplementary_Figure_S6a_descriptive_HC"))
safe_build("Supplementary Figure S6b", save_composite(
  two_cols(descr_keys(AD, exclude = character(0)), AD),
  "Supplementary_Figure_S6b_descriptive_AD"))

# ---- Supplementary Figures S7 / S8: leave-one-out ---------------------------
safe_build("Supplementary Figure S7", save_composite(compose_rows(list(list(
  loo_panel("AB42", "AB42", HC, "Aβ42"),
  loo_panel("AB42", "AB42", AD, " ")))),
  "Supplementary_Figure_S7_LOO_AB42"))

safe_build("Supplementary Figure S8", save_composite(compose_rows(list(list(
  loo_panel("PTAU181", "pTau181", HC, "p-tau181"),
  loo_panel("PTAU181", "pTau181", AD, " ")))),
  "Supplementary_Figure_S8_LOO_pTau181"))

# ---- Supplementary Figure S9: age (Abeta42) ---------------------------------
safe_build("Supplementary Figure S9", {
  a_hc <- age_panels(HC, "A. Subgroup forest plot", "B. Age meta-regression")
  a_ad <- age_panels(AD, " ", " ")
  save_composite(compose_rows(list(list(a_hc$sub, a_ad$sub),
                                   list(a_hc$reg, a_ad$reg))),
                 "Supplementary_Figure_S9_age")
})

message("\n11_composite_figures.R complete. See ", OUT_DIR, "/")
