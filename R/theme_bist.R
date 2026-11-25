# =============================================================================
# theme_bist.R  —  Shared ggplot2 theme for BIST-style figures
# -----------------------------------------------------------------------------
# BIST lab figure conventions: sans-serif type, no background grid,
# thin black axis rules, restrained colour, generous white space. Every figure
# in this project calls theme_bist() so the deck is typographically uniform.
#
# FONT: BIST style uses a sans-serif (Arial/Helvetica). On most systems "Arial" or
# "Helvetica" resolve; if you get a font warning, install the `extrafont`
# package or set BIST_FONT <- "" to fall back to the device default.
# =============================================================================

BIST_FONT  <- ""           # set to "Arial" or "Helvetica" if available
BIST_BLACK <- "#1A1A1A"
BIST_GREY  <- "#6E6E6E"
BIST_BLUE  <- "#00457C"    # accent for pooled diamonds / SROC curve
BIST_RED   <- "#A6093D"    # accent for "favours" annotations / warnings

# ---- Comparison styling -----------------------------------------------------
# Each comparison gets its OWN pretty label, accent colour, and short filename
# tag, so every forest plot is instantly identifiable as vs-HC or vs-AD - in
# the title, in the diamond colour, and in the saved file name.
COMPARISON_STYLE <- list(
  iNPH_vs_HC      = list(label = "iNPH vs Healthy Controls",
                         color = "#00457C", tag = "vs_HC"),   # blue
  iNPH_vs_AD      = list(label = "iNPH vs Alzheimer Disease",
                         color = "#9E4A06", tag = "vs_AD"),   # burnt amber
  `iNPH_vs_AD+HC` = list(label = "iNPH vs AD + Healthy Controls",
                         color = "#5B5B5B", tag = "vs_ADHC")  # grey
)
comparison_style <- function(comparison) {
  s <- COMPARISON_STYLE[[comparison]]
  if (is.null(s))
    s <- list(label = comparison, color = BIST_BLACK, tag = comparison)
  s
}

# ---- Biomarker display labels ----------------------------------------------
# Maps the cleaned tab name (left) to a publication-quality label (right).
# Tab names stay plain ASCII so they are valid Excel sheet names on every
# platform; the figures still show Greek letters. A biomarker not listed here
# simply uses its cleaned tab name verbatim.
BIOMARKER_LABELS <- c(
  "AB42"         = "A\u03B242",          # A-beta-42
  "AB40"         = "A\u03B240",
  "AB38"         = "A\u03B238",
  "tTau"         = "Total tau",
  "pTau181"      = "p-tau181",
  "sAPP"         = "sAPP",
  "sAPPalpha"    = "sAPP\u03B1",          # sAPP-alpha
  "sAPPbeta"     = "sAPP\u03B2",          # sAPP-beta
  "NFL"          = "NfL",
  "MCP-1"        = "MCP-1",
  "IL-8"         = "IL-8",
  "MBP"          = "MBP",
  "YKL-40"       = "YKL-40",
  "Sulfatide"    = "Sulfatide",
  "LRG-1"        = "LRG-1",
  "AB oligomers" = "A\u03B2 oligomers",
  "AB42/tTau"    = "A\u03B242/total tau",
  "AB42/AB40"    = "A\u03B242/A\u03B240",
  "pTau/tTau"    = "p-tau181/total tau",
  "pTau/AB42"    = "p-tau181/A\u03B242",
  "tTau/AB42"    = "Total tau/A\u03B242",      # revision: ratio reoriented
  "AB42/pTau"    = "A\u03B242/p-tau181",
  "ab42"         = "A\u03B242",              # diagnostic-workbook spelling
  # sheet-name variants where the tab spells the numerator "pTau181"
  "pTau181/tTau" = "p-tau181/total tau",
  "pTau181/AB42" = "p-tau181/A\u03B242"
)
biomarker_label <- function(x) {
  out <- BIOMARKER_LABELS[x]
  # unname(): ifelse() keeps the names of is.na(out), and an unmatched name is
  # NA - which breaks data.frame() row names downstream.
  unname(ifelse(is.na(out), x, unname(out)))
}

theme_bist <- function(base_size = 11, base_family = BIST_FONT) {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      axis.line.x      = element_line(colour = BIST_BLACK, linewidth = 0.4),
      axis.ticks.x     = element_line(colour = BIST_BLACK, linewidth = 0.4),
      axis.ticks.y     = element_blank(),
      axis.text        = element_text(colour = BIST_BLACK, size = base_size - 1),
      axis.title       = element_text(colour = BIST_BLACK, size = base_size),
      plot.title       = element_text(colour = BIST_BLACK, size = base_size + 1,
                                      face = "bold", hjust = 0,
                                      margin = margin(b = 6)),
      plot.subtitle    = element_text(colour = BIST_GREY, size = base_size - 1,
                                      hjust = 0, margin = margin(b = 8)),
      plot.caption     = element_text(colour = BIST_GREY, size = base_size - 2,
                                      hjust = 0, margin = margin(t = 8)),
      legend.position  = "none",
      plot.margin      = margin(12, 16, 12, 16)
    )
}

# Save a figure as both vector PDF (for submission) and PNG (for preview).
# ---- Raster resolution ------------------------------------------------------
# Submission rasters (line art / combination figures) are written at 600 dpi
# as LZW-compressed TIFF. THIS is the single value that sets the resolution of EVERY raster
# figure in the project (all 33 forest plots, the SROC plots, and Figure 2).
# Change it here and nowhere else.
RASTER_DPI <- 600

# Save a figure as vector PDF (proofing) and a 600-dpi LZW-compressed TIFF
# (submission).
save_bist <- function(plot, file, width = 9, height = 6, dir = PATHS$figures) {
  base <- file.path(dir, file)
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height,
         device = cairo_pdf)
  ggsave(paste0(base, ".tiff"), plot, width = width, height = height,
         dpi = RASTER_DPI, device = "tiff", compression = "lzw",
         bg = "white")
  message("  saved: ", base, ".pdf + .tiff (", RASTER_DPI, " dpi)")
  invisible(base)
}
