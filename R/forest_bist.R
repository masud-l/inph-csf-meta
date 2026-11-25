# =============================================================================
# forest_bist.R  —  Shared, BIST-style SMD forest plot
# -----------------------------------------------------------------------------
# Strategy: metafor does the statistics; ggplot2 does the rendering. Every
# biomarker forest plot is produced by THIS ONE function, so all panels in the
# manuscript are typographically and structurally identical.
#
# Comparison encoding: the plot makes vs-HC vs. vs-AD unmistakable in THREE
# ways - the bold biomarker title, a colour-coded comparison subtitle, and the
# pooled diamond drawn in the comparison's accent colour.
#
# The plot is assembled with patchwork from three vertically-aligned panels:
#   [ study labels ] [ forest ] [ numeric SMD (95% CI) + weight ]
# Sharing one y-scale across panels guarantees the rows line up exactly.
# =============================================================================

# ---- Extract a 4-digit year from a study label, if present ------------------
.study_year <- function(s) {
  y <- stringr::str_extract(s, "(19|20)\\d{2}")
  suppressWarnings(as.numeric(y))
}

# ---- Build the per-row plotting frame from a run_smd_one() result -----------
.forest_frame <- function(result, order_by = "effect") {
  es <- result$es
  df <- data.frame(
    study    = es$study,
    estimate = as.numeric(es$yi),
    se       = sqrt(as.numeric(es$vi))
  )
  df$ci_lb <- df$estimate - 1.96 * df$se
  df$ci_ub <- df$estimate + 1.96 * df$se

  # Study weights from the fitted random-effects model (percentage).
  df$weight <- if (!is.null(result$fit)) as.numeric(weights(result$fit))
               else 100

  ord <- switch(order_by,
    effect = order(df$estimate),
    year   = order(.study_year(df$study), df$study),
    alpha  = order(df$study, decreasing = TRUE),
    none   = rev(seq_len(nrow(df)))
  )
  df <- df[ord, , drop = FALSE]
  df$y <- seq_len(nrow(df)) + 1L     # row 1 reserved for the pooled diamond
  df
}

# ---- Main entry point -------------------------------------------------------
# result     : output of run_smd_one()
# biomarker  : biomarker name; defaults to result$summary$biomarker
# comparison : "iNPH_vs_HC" / "iNPH_vs_AD"; defaults to result$summary$comparison
# xlab       : x-axis label
# order_by   : "effect" | "year" | "alpha" | "none"
forest_bist <- function(result,
                        biomarker  = NULL,
                        comparison = NULL,
                        xlab       = "Standardized mean difference (Hedges' g)",
                        order_by   = "effect") {

  # --- Resolve identity and styling -----------------------------------------
  if (is.null(biomarker))  biomarker  <- result$summary$biomarker
  if (is.null(comparison)) comparison <- result$summary$comparison
  bm_label <- biomarker_label(biomarker)
  cs       <- comparison_style(comparison)   # $label, $color, $tag
  comp_col <- cs$color

  df <- .forest_frame(result, order_by)
  has_pool <- !is.null(result$fit)

  # Pooled estimate (row y = 1).
  if (has_pool) {
    pooled <- data.frame(est = as.numeric(result$fit$b),
                         lb  = result$fit$ci.lb,
                         ub  = result$fit$ci.ub)
  } else {
    pooled <- data.frame(est = df$estimate[1], lb = df$ci_lb[1],
                         ub = df$ci_ub[1])
  }
  diamond <- data.frame(
    x = c(pooled$lb, pooled$est, pooled$ub, pooled$est),
    y = c(1,         1 + 0.32,   1,         1 - 0.32)
  )

  # Shared y scale across all three panels. The top limit leaves a full extra
  # row of headroom ABOVE the topmost study for the column-header row, plus a
  # margin so the header text is never clipped. Increase the "+ 1.8" if headers
  # still clip; decrease it to tighten the white space at the top.
  y_lim <- c(0.4, max(df$y) + 1.8)
  y_scale <- scale_y_continuous(limits = y_lim, expand = c(0, 0))

  # Box size scaled by weight (sqrt so area, not side, tracks weight).
  df$box <- scales::rescale(sqrt(df$weight), to = c(2, 6.5))

  # --- x-axis range: pad symmetrically around 0 ------------------------------
  x_rng <- range(c(df$ci_lb, df$ci_ub, diamond$x, 0))
  x_pad <- diff(x_rng) * 0.08
  x_lim <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)

  # Column-header row sits one row above the topmost study. Edit HEADER_Y to
  # move every header up/down together; it must stay below y_lim's top (which
  # is max(df$y) + 1.8) or the header text will clip again.
  header_y <- max(df$y) + 1.0

  # --- Panel 1: study labels -------------------------------------------------
  p_lab <- ggplot(df) +
    geom_text(aes(x = 0, y = y, label = study),
              hjust = 0, size = 3.1, colour = BIST_BLACK) +
    { if (has_pool)
        annotate("text", x = 0, y = 1, label = "Pooled estimate",
                 hjust = 0, fontface = "bold", size = 3.1, colour = BIST_BLACK) } +
    annotate("text", x = 0, y = header_y, label = "Study",
             hjust = 0, fontface = "bold", size = 3.0, colour = BIST_BLACK) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    y_scale +
    theme_void(base_family = BIST_FONT) +
    theme(plot.margin = margin(12, 2, 12, 8))

  # --- Panel 2: the forest (diamond + pooled line in comparison colour) ------
  p_for <- ggplot() +
    geom_vline(xintercept = 0, colour = BIST_GREY, linewidth = 0.4) +
    { if (has_pool)
        geom_vline(xintercept = pooled$est, colour = comp_col,
                   linetype = "dotted", linewidth = 0.4) } +
    geom_segment(data = df,
                 aes(x = ci_lb, xend = ci_ub, y = y, yend = y),
                 colour = BIST_BLACK, linewidth = 0.5) +
    geom_point(data = df, aes(x = estimate, y = y, size = box),
               shape = 15, colour = BIST_BLACK) +
    scale_size_identity() +
    { if (has_pool)
        geom_polygon(data = diamond, aes(x = x, y = y),
                     fill = comp_col, colour = comp_col) } +
    scale_x_continuous(limits = x_lim) +
    y_scale +
    labs(x = xlab) +
    theme_bist() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          plot.margin = margin(12, 6, 12, 6))

  # --- Panel 3: numeric column ----------------------------------------------
  df$lab_est <- sprintf("%.2f (%.2f to %.2f)",
                        df$estimate, df$ci_lb, df$ci_ub)
  df$lab_wt  <- sprintf("%.1f", df$weight)
  pooled_lab <- sprintf("%.2f (%.2f to %.2f)",
                        pooled$est, pooled$lb, pooled$ub)

  p_num <- ggplot(df) +
    geom_text(aes(x = 0, y = y, label = lab_est),
              hjust = 0, size = 3.0, colour = BIST_BLACK) +
    { if (has_pool)
        geom_text(aes(x = 1, y = y, label = lab_wt),
                  hjust = 1, size = 3.0, colour = BIST_GREY) } +
    { if (has_pool)
        annotate("text", x = 0, y = 1, label = pooled_lab,
                 hjust = 0, fontface = "bold", size = 3.0, colour = BIST_BLACK) } +
    annotate("text", x = 0, y = header_y,
             label = "SMD (95% CI)", hjust = 0, fontface = "bold",
             size = 3.0, colour = BIST_BLACK) +
    { if (has_pool)
        annotate("text", x = 1, y = header_y, label = "Wt, %",
                 hjust = 1, fontface = "bold", size = 3.0, colour = BIST_GREY) } +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    y_scale +
    theme_void(base_family = BIST_FONT) +
    theme(plot.margin = margin(12, 8, 12, 6))

  # --- Heterogeneity caption -------------------------------------------------
  cap <- if (has_pool) {
    f <- result$fit
    sprintf("Random-effects model (%s). k = %d studies. I-squared = %.0f%%, tau-squared = %.3f, Q-test P %s.",
            f$method, f$k, f$I2, f$tau2,
            ifelse(f$QEp < 0.001, "< .001", sprintf("= %.3f", f$QEp)))
  } else if (nrow(df) == 1) {
    "Single study - descriptive estimate only; not pooled."
  } else {
    sprintf("k = %d studies (fewer than %d) - individual estimates shown; not pooled.",
            nrow(df), MIN_K)
  }

  # --- Compose: title = biomarker, subtitle = comparison (colour-coded) ------
  combined <- patchwork::wrap_plots(p_lab, p_for, p_num,
                                    widths = c(2.6, 4.2, 2.4)) +
    patchwork::plot_annotation(
      title    = bm_label,
      subtitle = cs$label,
      caption  = cap,
      theme    = theme_bist() +
        theme(plot.subtitle = element_text(colour = comp_col,
                                           size = 11))
    )
  combined
}
