# =============================================================================
# figure2_overview.R  —  Manuscript Figure 2 candidates
# -----------------------------------------------------------------------------
# The all-biomarker pooled-SMD overview, built TWO ways so they can be compared:
#   figure2_dotplot()  - point + 95% CI dot plot, two panels (vs HC | vs AD)
#   figure2_table()    - the same content as a formatted summary table
#
# Both reuse the project theme (theme_bist, comparison_style, biomarker_label,
# BIST_* colours) so whichever is chosen matches the supplementary forest plots.
#
# Data source: output/tables/smd_pooled_estimates.csv (written by scripts/02).
# Those estimates use whatever RE_METHOD was set in 00_setup.R; for the
# manuscript that must be REML. The DerSimonian-Laird reproduction is an
# internal audit and is deliberately NOT shown here.
# =============================================================================

# ---- Key biomarkers to highlight -------------------------------------------
# The four primary markers of the manuscript. In Figure 2 they get a bold
# label and a full-size point; all other biomarkers are de-emphasised in grey
# while keeping their true effect-size rank. These strings MUST match the
# output of biomarker_label() (the displayed labels, with Greek letters).
KEY_BIOMARKERS <- c("A\u03B242", "NfL", "p-tau181", "sAPP\u03B1")

# ---- Load + prepare the pooled-estimate table -------------------------------
load_overview <- function(path = file.path(PATHS$tables,
                                           "smd_pooled_estimates.csv")) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  # Figure 2 shows only biomarkers that met the pooling threshold (MIN_K from
  # 00_setup.R); biomarkers with fewer studies are not meta-analysed.
  n_before <- nrow(d)
  d <- d[d$k >= MIN_K, , drop = FALSE]
  message("  Figure 2: kept ", nrow(d), " of ", n_before,
          " biomarkers (k >= MIN_K = ", MIN_K, ").")
  d$label  <- biomarker_label(d$biomarker)
  d$sig    <- !(d$ci_lb < 0 & d$ci_ub > 0)   # TRUE = 95% CI excludes zero
  d$is_key <- d$label %in% KEY_BIOMARKERS    # one of the 4 highlighted markers
  d
}

# ---- Tidy data behind the table (also written to CSV) -----------------------
figure2_overview_data <- function(ov) {
  out <- data.frame(
    Comparison  = ifelse(ov$comparison == "iNPH_vs_HC",
                         "iNPH vs healthy controls",
                         "iNPH vs Alzheimer disease"),
    Biomarker   = ov$label,
    No_Studies  = ov$k,
    SMD         = round(ov$estimate, 2),
    CI_lower    = round(ov$ci_lb, 2),
    CI_upper    = round(ov$ci_ub, 2),
    SMD_95CI    = sprintf("%.2f (%.2f to %.2f)",
                          ov$estimate, ov$ci_lb, ov$ci_ub),
    I2_percent  = ifelse(is.na(ov$I2), NA_real_, round(ov$I2, 1)),
    Significant = ifelse(ov$sig, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  out[order(out$Comparison, out$SMD), ]
}

# =============================================================================
# CANDIDATE A — dot plot
# =============================================================================

# One comparison panel. xr is the shared x-range (so both panels align).
# Significance is shown by square fill (filled = CI excludes zero); the four
# KEY_BIOMARKERS are emphasised with a bold label and full-size point, all
# others de-emphasised in grey while keeping their true effect-size rank.
.overview_panel <- function(d, comp_key, xr) {
  cs <- comparison_style(comp_key)
  d  <- d[order(d$estimate), ]
  d$label   <- factor(d$label, levels = d$label)
  d$sig_chr <- as.character(d$sig)
  
  # Per-row emphasis for the key biomarkers.
  d$pt_size <- ifelse(d$is_key, 3.2, 2.2)
  d$edge    <- ifelse(d$is_key, cs$color, BIST_GREY)
  d$err_col <- ifelse(d$is_key, BIST_BLACK, BIST_GREY)
  ord       <- order(d$label)                       # bottom-to-top axis order
  axis_cols <- ifelse(d$is_key[ord], BIST_BLACK, BIST_GREY)
  axis_face <- ifelse(d$is_key[ord], "bold", "plain")
  
  ggplot(d, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 0, colour = BIST_GREY, linewidth = 0.4) +
    geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub, colour = err_col),
                   height = 0, linewidth = 0.5) +
    geom_point(aes(fill = sig_chr, size = pt_size, colour = edge),
               shape = 22, stroke = 0.6) +
    geom_text(aes(x = xr[2], label = paste0("k=", k)),
              hjust = 1, size = 2.5, colour = BIST_GREY) +
    scale_fill_manual(values = c("TRUE" = cs$color, "FALSE" = "white"),
                      guide = "none") +
    scale_colour_identity() +
    scale_size_identity() +
    scale_x_continuous(breaks = seq(-6, 6, 2)) +
    coord_cartesian(xlim = xr, clip = "off") +
    labs(x = "Standardized mean difference (95% CI)", title = cs$label,
         y = NULL) +                                # blank the "label" y-title
    theme_bist() +
    theme(plot.title         = element_text(colour = cs$color,
                                            face = "bold", size = 11),
          axis.title.y       = element_blank(),     # belt-and-braces
          axis.text.y        = element_text(colour = axis_cols,
                                            face = axis_face),
          panel.grid.major.x = element_line(colour = "grey90",
                                            linewidth = 0.3))
}

# Shared significance legend, drawn as its own tiny plot so it can sit once
# beneath both panels (filled square = significant, open square = not).
.significance_legend <- function() {
  ld <- data.frame(
    # Same x for both -> vertically stacked; y separates the two rows.
    x   = c(0.4, 0.4),
    y   = c(1.4, 0.5),
    lab = c("Significant (95% CI excludes 0)",
            "Not significant (95% CI includes 0)"),
    fil = c(BIST_BLACK, "white")
  )
  ggplot(ld, aes(x = x, y = y)) +
    geom_point(shape = 22, size = 3.0, stroke = 0.6,
               fill = ld$fil, colour = BIST_BLACK) +
    geom_text(aes(label = lab), hjust = 0, nudge_x = 0.018,
              size = 2.8, colour = BIST_BLACK) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(-0.6, 1.6)) +
    theme_void(base_family = BIST_FONT)
}

# Two-panel dot plot. ov = output of load_overview().
#
# The figure carries the SIGNIFICANCE KEY (filled vs open square) as a strip
# beneath the panels, but NOT a descriptive caption. By BIST convention the
# descriptive figure legend is submitted as separate manuscript text, not
# baked into the image - this also permanently removes the caption-clipping
# problem. The legend text to paste into the manuscript is in the constant
# FIGURE2_LEGEND below; figure2_dotplot() prints it when run.
figure2_dotplot <- function(ov) {
  # Shared x-range across both panels, with extra room on the right for the
  # "k=" labels. coord_cartesian (not scale limits) zooms the view WITHOUT
  # dropping any CI that runs past it.
  rng   <- range(c(ov$ci_lb, ov$ci_ub))
  pad_l <- diff(rng) * 0.04
  pad_r <- diff(rng) * 0.15            # room for k= labels
  xr    <- c(rng[1] - pad_l, rng[2] + pad_r)
  
  hc <- .overview_panel(ov[ov$comparison == "iNPH_vs_HC", ],
                        "iNPH_vs_HC", xr)
  ad <- .overview_panel(ov[ov$comparison == "iNPH_vs_AD", ],
                        "iNPH_vs_AD", xr)
  
  # Layout: two panels on top, the shared significance key as a thin strip
  # beneath. No caption - see FIGURE2_LEGEND for the manuscript legend text.
  panels <- patchwork::wrap_plots(hc, ad, nrow = 1)
  legend <- .significance_legend()
  
  message("\n  --- Figure 2 legend (paste into the manuscript) ---\n")
  message(FIGURE2_LEGEND, "\n")
  
  patchwork::wrap_plots(panels, legend, ncol = 1, heights = c(20, 2)) +
    patchwork::plot_annotation(
      theme = theme_bist() +
        theme(plot.margin = margin(12, 16, 10, 16))
    )
}

# ---- Manuscript figure legend (text, submitted separately from the image) ---
# Printed by figure2_dotplot() so it is easy to copy into the manuscript.
FIGURE2_LEGEND <- paste0(
  "Figure 2. Pooled standardized mean differences (SMD, Hedges' g) of CSF ",
  "biomarkers in idiopathic normal pressure hydrocephalus (iNPH) versus ",
  "healthy controls and versus Alzheimer disease. Each point is the ",
  "random-effects (REML) pooled SMD; error bars are 95% CIs. Filled squares ",
  "indicate a 95% CI that excludes zero (statistically significant); open ",
  "squares indicate a 95% CI that includes zero. k denotes the number of ",
  "contributing studies. Biomarkers in bold (A\u03B242, NfL, p-tau181, ",
  "sAPP\u03B1) are the primary markers of interest. A negative SMD indicates ",
  "a lower CSF concentration in iNPH than in the comparator. Only biomarkers ",
  "meta-analysed in at least ", if (exists("MIN_K")) MIN_K else 3,
  " studies are shown. HC indicates healthy controls; AD, Alzheimer disease."
)

# =============================================================================
# CANDIDATE B — summary table
# =============================================================================

# One comparison block of the table: Biomarker | No. of Studies | SMD (95% CI)
# | I-squared. Rendered as a text grid so it carries the project theme.
.table_block <- function(d, comp_key) {
  cs <- comparison_style(comp_key)
  d  <- d[order(d$estimate), ]
  n  <- nrow(d)
  d$row     <- n:1                                   # first row at the top
  d$smd_txt <- sprintf("%.2f (%.2f to %.2f)", d$estimate, d$ci_lb, d$ci_ub)
  d$i2_txt  <- ifelse(is.na(d$I2), "\u2014", sprintf("%.1f", d$I2))
  d$k_txt   <- as.character(d$k)
  
  x_bio <- 0.02; x_k <- 0.42; x_smd <- 0.50; x_i2 <- 0.985
  hdr_y <- n + 1.4
  
  ggplot() +
    # block title
    annotate("text", x = 0, y = n + 2.8, label = cs$label,
             hjust = 0, fontface = "bold", size = 3.4, colour = cs$color) +
    # column headers
    annotate("text", x = x_bio, y = hdr_y, label = "Biomarker",
             hjust = 0, fontface = "bold", size = 2.9, colour = BIST_BLACK) +
    annotate("text", x = x_k, y = hdr_y, label = "No. of\nStudies",
             hjust = 0.5, fontface = "bold", size = 2.9, colour = BIST_BLACK) +
    annotate("text", x = x_smd, y = hdr_y, label = "SMD (95% CI)",
             hjust = 0, fontface = "bold", size = 2.9, colour = BIST_BLACK) +
    annotate("text", x = x_i2, y = hdr_y, label = "I\u00b2, %",
             hjust = 1, fontface = "bold", size = 2.9, colour = BIST_BLACK) +
    # rules above and below the data rows
    annotate("segment", x = 0, xend = 1, y = n + 0.6, yend = n + 0.6,
             colour = BIST_BLACK, linewidth = 0.4) +
    annotate("segment", x = 0, xend = 1, y = 0.45, yend = 0.45,
             colour = BIST_BLACK, linewidth = 0.4) +
    # data rows
    geom_text(data = d, aes(y = row, label = label), x = x_bio,
              hjust = 0, size = 2.8, colour = BIST_BLACK) +
    geom_text(data = d, aes(y = row, label = k_txt), x = x_k,
              hjust = 0.5, size = 2.8, colour = BIST_BLACK) +
    geom_text(data = d, aes(y = row, label = smd_txt), x = x_smd,
              hjust = 0, size = 2.8, colour = BIST_BLACK) +
    geom_text(data = d, aes(y = row, label = i2_txt), x = x_i2,
              hjust = 1, size = 2.8, colour = BIST_BLACK) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.2, n + 3.4), expand = c(0, 0)) +
    theme_void(base_family = BIST_FONT) +
    theme(plot.margin = margin(6, 10, 6, 10))
}

# Full table: the two comparison blocks stacked. ov = output of load_overview().
# Like the dot plot, the table image carries no caption; the descriptive
# legend (TABLE2_LEGEND) is submitted as separate manuscript text.
figure2_table <- function(ov) {
  hc   <- .table_block(ov[ov$comparison == "iNPH_vs_HC", ], "iNPH_vs_HC")
  ad   <- .table_block(ov[ov$comparison == "iNPH_vs_AD", ], "iNPH_vs_AD")
  n_hc <- sum(ov$comparison == "iNPH_vs_HC")
  n_ad <- sum(ov$comparison == "iNPH_vs_AD")
  
  message("\n  --- Table legend (paste into the manuscript) ---\n")
  message(TABLE2_LEGEND, "\n")
  
  patchwork::wrap_plots(hc, ad, ncol = 1,
                        heights = c(n_hc + 3.4, n_ad + 3.4)) +
    patchwork::plot_annotation(
      theme = theme_bist() +
        theme(plot.margin = margin(10, 16, 10, 16))
    )
}

# ---- Manuscript table legend (text, submitted separately from the image) ----
TABLE2_LEGEND <- paste0(
  "Pooled standardized mean differences (SMD, Hedges' g) with 95% CIs of CSF ",
  "biomarkers in idiopathic normal pressure hydrocephalus (iNPH) versus ",
  "healthy controls and versus Alzheimer disease, estimated by ",
  "random-effects (REML) meta-analysis. I-squared is the percentage of total ",
  "variance attributable to between-study heterogeneity. A negative SMD ",
  "indicates a lower CSF concentration in iNPH than in the comparator. Only ",
  "biomarkers meta-analysed in at least ", if (exists("MIN_K")) MIN_K else 3,
  " studies are shown. HC indicates healthy controls; AD, Alzheimer disease."
)