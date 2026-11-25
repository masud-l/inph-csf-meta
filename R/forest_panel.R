# =============================================================================
# forest_panel.R  —  Slim, grid-ready forest panel for multi-panel figures
# -----------------------------------------------------------------------------
# forest_bist() produces a STANDALONE figure (own title, caption, x-axis) - one
# biomarker per file. That is wrong for a supplementary figure that tiles
# several biomarkers in a 2x2 grid: it would stack redundant axes and captions.
#
# forest_panel() produces a PANEL instead - the same forest, but stripped for
# tiling: a compact panel header ("A  Aβ42" + comparison), no caption, and a
# heterogeneity line folded into the header. assemble_supp_figure() then composes
# four panels into one 2x2 supplementary figure with a single shared caption.
#
# Depends on helpers already defined in forest_bist.R (.forest_frame,
# .study_year) and on theme_bist.R - source those first.
# =============================================================================

# ---- One slim panel ---------------------------------------------------------
# result    : output of run_smd_one()
# panel_tag : the panel letter, e.g. "A" (NULL = none)
# order_by  : passed to .forest_frame()
forest_panel <- function(result, panel_tag = NULL, order_by = "effect") {

  biomarker  <- result$summary$biomarker
  comparison <- result$summary$comparison
  bm_label   <- biomarker_label(biomarker)
  cs         <- comparison_style(comparison)
  comp_col   <- cs$color

  df       <- .forest_frame(result, order_by)
  has_pool <- !is.null(result$fit)

  # Pooled estimate / diamond (row y = 1).
  if (has_pool) {
    pooled <- data.frame(est = as.numeric(result$fit$b),
                         lb  = result$fit$ci.lb, ub = result$fit$ci.ub)
  } else {
    pooled <- data.frame(est = df$estimate[1], lb = df$ci_lb[1],
                         ub = df$ci_ub[1])
  }
  diamond <- data.frame(
    x = c(pooled$lb, pooled$est, pooled$ub, pooled$est),
    y = c(1,         1 + 0.30,   1,         1 - 0.30)
  )

  # Header text: panel letter + biomarker, then comparison, then heterogeneity.
  head_main <- if (is.null(panel_tag)) bm_label
               else paste0(panel_tag, "  ", bm_label)
  head_sub  <- if (has_pool) {
    f <- result$fit
    sprintf("%s  \u00b7  k = %d  \u00b7  I\u00b2 = %.0f%%",
            cs$label, f$k, f$I2)
  } else {
    sprintf("%s  \u00b7  k = %d  \u00b7  descriptive, not pooled",
            cs$label, nrow(df))
  }

  # y scale: one row per study + the pooled row + headroom for the 2-line head.
  y_top <- max(df$y) + 2.4
  y_lim <- c(0.3, y_top)

  # x range padded around zero.
  x_rng <- range(c(df$ci_lb, df$ci_ub, diamond$x, 0))
  x_pad <- diff(x_rng) * 0.08
  # Reserve a label column to the LEFT of the data so study names never sit
  # on top of the confidence intervals; axis ticks cover the data range only.
  lab_w <- diff(x_rng) * 0.75
  x_lim <- c(x_rng[1] - x_pad - lab_w, x_rng[2] + x_pad)
  x_brk <- pretty(x_rng, 4)
  x_brk <- x_brk[x_brk >= x_rng[1] - x_pad & x_brk <= x_lim[2]]

  # Box size by weight (sqrt -> area tracks weight). Smaller than the
  # standalone figure because panels are physically smaller in a grid.
  df$box <- scales::rescale(sqrt(df$weight), to = c(1.6, 5))

  ggplot() +
    # header
    annotate("text", x = x_lim[1], y = y_top - 0.2, label = head_main,
             hjust = 0, fontface = "bold", size = 3.3, colour = BIST_BLACK) +
    annotate("text", x = x_lim[1], y = y_top - 1.1, label = head_sub,
             hjust = 0, size = 2.6, colour = comp_col) +
    # reference + pooled lines
    # Reference lines stop below the two-line header so they never cross it.
    annotate("segment", x = 0, xend = 0, y = y_lim[1], yend = max(df$y) + 0.6,
             colour = BIST_GREY, linewidth = 0.4) +
    { if (has_pool)
        annotate("segment", x = pooled$est, xend = pooled$est,
                 y = y_lim[1], yend = max(df$y) + 0.6, colour = comp_col,
                 linetype = "dotted", linewidth = 0.4) } +
    # study CIs + points
    geom_segment(data = df, aes(x = ci_lb, xend = ci_ub, y = y, yend = y),
                 colour = BIST_BLACK, linewidth = 0.45) +
    geom_point(data = df, aes(x = estimate, y = y, size = box),
               shape = 15, colour = BIST_BLACK) +
    scale_size_identity() +
    # study labels (left) and pooled label (row 1)
    geom_text(data = df, aes(x = x_lim[1], y = y, label = study),
              hjust = 0, size = 2.4, colour = BIST_BLACK) +
    { if (has_pool)
        annotate("text", x = x_lim[1], y = 1, label = "Pooled estimate",
                 hjust = 0, fontface = "bold", size = 2.4, colour = BIST_BLACK) } +
    # pooled diamond
    # (no diamond for outcomes that were not pooled, k < MIN_K)
    { if (has_pool)
        geom_polygon(data = diamond, aes(x = x, y = y),
                     fill = comp_col, colour = comp_col) } +
    scale_x_continuous(limits = x_lim, breaks = x_brk) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    labs(x = "Standardized mean difference (Hedges' g)") +
    theme_bist() +
    theme(axis.text.y  = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 8),
          plot.margin  = margin(6, 10, 6, 10))
}

# ---- Assemble one 2x2 supplementary figure ----------------------------------
# results   : list of run_smd_one() outputs (up to 4) for this figure
# tags      : panel letters, default A-D
# title     : figure title
# caption   : one shared caption for the whole figure
assemble_supp_figure <- function(results, title, caption,
                                 tags = c("A", "B", "C", "D")) {
  n <- length(results)
  if (n < 1 || n > 4)
    stop("assemble_supp_figure() takes 1-4 panels; got ", n)

  panels <- lapply(seq_len(n), function(i)
    forest_panel(results[[i]], panel_tag = tags[i]))

  # Always lay out on a 2-column grid; 1-3 panels simply leave a blank cell.
  patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title   = title,
      caption = stringr::str_wrap(caption, width = 190),
      theme   = theme_bist() +
        theme(plot.title   = element_text(face = "bold", size = 12),
              plot.caption = element_text(size = 8, colour = BIST_GREY,
                                          hjust = 0, lineheight = 1.4),
              plot.margin  = margin(12, 16, 16, 16))
    )
}
