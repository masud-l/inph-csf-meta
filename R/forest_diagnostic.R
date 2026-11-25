# =============================================================================
# forest_diagnostic.R  —  Paired sensitivity/specificity forest + SROC
# -----------------------------------------------------------------------------
# Diagnostic accuracy is presented two ways, both standard in diagnostic
# test accuracy meta-analysis (Cochrane DTA reviews):
#   (1) paired_dta_forest()  - side-by-side sensitivity and specificity forest
#   (2) sroc_bist()          - summary ROC space with the bivariate operating
#                              point and confidence region
# =============================================================================

# ---- (1) Paired sensitivity / specificity forest ----------------------------
# per_study : output of diagnostic_per_study() for ONE biomarker
# pooled    : the $summary data frame from run_reitsma_one() (or NULL)
# comparison: comparison key; defaults to per_study$comparison
paired_dta_forest <- function(per_study, pooled = NULL, comparison = NULL) {

  if (is.null(comparison)) comparison <- per_study$comparison[1]
  cs       <- comparison_style(comparison)
  bm_label <- biomarker_label(per_study$biomarker[1])

  ps <- per_study[order(per_study$sens), ]
  ps$y <- seq_len(nrow(ps)) + 1L
  y_lim <- c(0.4, max(ps$y) + 0.8)

  has_pool <- !is.null(pooled)
  pooled_row <- if (has_pool) data.frame(
    sens = pooled$sensitivity, sens_lo = pooled$sens_lo, sens_hi = pooled$sens_hi,
    spec = pooled$specificity, spec_lo = pooled$spec_lo, spec_hi = pooled$spec_hi
  ) else NULL

  one_panel <- function(est, lo, hi, axis_title) {
    g <- ggplot(ps) +
      geom_segment(aes(x = .data[[lo]], xend = .data[[hi]], y = y, yend = y),
                   colour = BIST_BLACK, linewidth = 0.5) +
      geom_point(aes(x = .data[[est]], y = y), shape = 15, size = 2.6,
                 colour = BIST_BLACK)
    if (has_pool) {
      d <- data.frame(
        x = c(pooled_row[[lo]], pooled_row[[est]],
              pooled_row[[hi]], pooled_row[[est]]),
        y = c(1, 1.32, 1, 0.68))
      g <- g + geom_polygon(data = d, aes(x, y),
                            fill = cs$color, colour = cs$color)
    }
    g +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                         labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
      labs(x = axis_title) +
      theme_bist() +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
            plot.margin = margin(12, 8, 12, 8))
  }

  p_sens <- one_panel("sens", "sens_lo", "sens_hi", "Sensitivity")
  p_spec <- one_panel("spec", "spec_lo", "spec_hi", "Specificity")

  p_lab <- ggplot(ps) +
    geom_text(aes(x = 0, y = y, label = study), hjust = 0, size = 3.0) +
    { if (has_pool)
        annotate("text", x = 0, y = 1, label = "Pooled estimate",
                 hjust = 0, fontface = "bold", size = 3.0) } +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    theme_void(base_family = BIST_FONT) +
    theme(plot.margin = margin(12, 2, 12, 8))

  uni <- has_pool && "model" %in% names(pooled) &&
         grepl("^univariate", pooled$model[1])
  cap <- if (has_pool && uni)
    sprintf("Separate univariate random-effects logit models for sensitivity and specificity. k = %d studies.\nBivariate between-study correlation at its boundary: no summary ROC curve or AUC. 0.5 continuity correction applied to studies with a zero cell.",
            pooled$k)
  else if (has_pool)
    sprintf("Bivariate (Reitsma) random-effects model. k = %d studies. Pooled AUC = %.3f. 0.5 continuity correction applied to studies with a zero cell.",
            pooled$k, pooled$AUC)
  else
    "Single study - descriptive estimate only; no bivariate model fitted."

  patchwork::wrap_plots(p_lab, p_sens, p_spec, widths = c(2.6, 3.4, 3.4)) +
    patchwork::plot_annotation(
      title    = paste(bm_label, "- diagnostic accuracy"),
      subtitle = cs$label,
      caption  = cap,
      theme    = theme_bist() +
        theme(plot.subtitle = element_text(colour = cs$color,
                                           size = 11)))
}

# ---- (2) Summary ROC plot ---------------------------------------------------
# Built on mada's native plot method (correct ellipse / SROC math), restyled
# to project conventions. Saves directly to PDF + PNG.
# `dta` is the 2x2 data frame (TP/FN/FP/TN) returned by run_reitsma_one()$data;
# study points are computed from it directly so no reliance on fit internals.
sroc_bist <- function(reitsma_fit, dta = NULL, title = NULL,
                      auc = NULL, comp = "iNPH_vs_AD",
                      file = "sroc", dir = PATHS$figures, save = TRUE) {
  if (is.null(reitsma_fit)) {
    message("  sroc_bist: no fitted model (k too small) - skipped.")
    return(invisible(NULL))
  }

  cs      <- comparison_style(comp)
  acc_col <- cs$color

  # ---- Helper: build an ellipse in (FPR, Sens) space -----------------------
  # centre: c(mu_tfpr, mu_tsens) in logit space
  # Sigma:  2x2 covariance in logit space (vcov for CI, vcov+Psi for prediction)
  # level:  e.g. 0.95
  # Returns data frame with columns fpr, sens (back-transformed via plogis).
  make_ellipse <- function(centre, Sigma, level = 0.95, n = 200) {
    chi2 <- qchisq(level, df = 2)
    # Cholesky of Sigma
    L <- tryCatch(t(chol(Sigma)), error = function(e) NULL)
    if (is.null(L)) return(NULL)
    theta <- seq(0, 2 * pi, length.out = n + 1)
    pts   <- sqrt(chi2) * cbind(cos(theta), sin(theta))
    pts   <- sweep(pts %*% t(L), 2, centre, "+")   # rotate + translate
    data.frame(fpr  = plogis(pts[, 1]),
               sens = plogis(pts[, 2]))
  }

  # ---- Extract Reitsma model parameters ------------------------------------
  co  <- tryCatch(coef(reitsma_fit),  error = function(e) NULL)
  V   <- tryCatch(vcov(reitsma_fit),  error = function(e) NULL)
  Psi <- tryCatch(reitsma_fit$Psi,    error = function(e) NULL)
  if (is.null(co) || is.null(V)) {
    message("  sroc_bist: could not extract model parameters.")
    return(invisible(NULL))
  }

  # co[1] = logit-sensitivity, co[2] = logit-FPR (Reitsma parameterisation).
  # If names are present, use them; else fall back to position.
  nm <- names(co)
  i_s <- if (!is.null(nm)) grep("sens", nm, ignore.case = TRUE)[1] else 1L
  i_f <- if (!is.null(nm)) grep("fpr",  nm, ignore.case = TRUE)[1] else 2L
  if (is.na(i_s)) i_s <- 1L
  if (is.na(i_f)) i_f <- 2L

  mu_ts <- as.numeric(co[i_s])    # logit(mean sensitivity)
  mu_tf <- as.numeric(co[i_f])    # logit(mean FPR)

  # ---- 1. SROC curve --------------------------------------------------------
  # On the SROC, for a given logit-FPR value t:
  #   logit(sens) = mu_ts + (Psi_sf / Psi_ff) * (t - mu_tf)
  # If Psi is unavailable, use a horizontal line through the summary point.
  slope <- if (!is.null(Psi) && is.matrix(Psi) && all(dim(Psi) == 2)) {
    Psi[i_s, i_f] / Psi[i_f, i_f]
  } else 0

  tfpr_seq  <- seq(-5, 5, length.out = 300)
  tsens_seq <- mu_ts + slope * (tfpr_seq - mu_tf)
  sroc_pts  <- data.frame(fpr  = plogis(tfpr_seq),
                           sens = plogis(tsens_seq))

  # ---- 2. Confidence ellipse (uncertainty in the mean operating point) ------
  # Use the 2×2 submatrix of vcov corresponding to (tsens, tfpr).
  V_sub <- V[c(i_s, i_f), c(i_s, i_f)]
  # centre in (tfpr, tsens) order to match make_ellipse (x=FPR, y=sens)
  centre <- c(mu_tf, mu_ts)
  Sigma_ci <- V_sub[c(2, 1), c(2, 1)]    # reorder to (fpr, sens)
  conf_df  <- make_ellipse(centre, Sigma_ci, 0.95)

  # ---- 3. Prediction ellipse (where a new study might fall) ----------------
  pred_df <- if (!is.null(Psi) && is.matrix(Psi) && all(dim(Psi) == 2)) {
    Psi_sub     <- Psi[c(i_s, i_f), c(i_s, i_f)]
    Sigma_pred  <- V_sub[c(2, 1), c(2, 1)] + Psi_sub[c(2, 1), c(2, 1)]
    make_ellipse(centre, Sigma_pred, 0.95)
  } else NULL

  # ---- 4. Summary operating point ------------------------------------------
  sp <- data.frame(fpr  = plogis(mu_tf),
                   sens = plogis(mu_ts))

  # ---- 5. Individual study points ------------------------------------------
  study_df <- NULL
  if (!is.null(dta) && all(c("TP","FN","FP","TN") %in% names(dta))) {
    tp <- as.numeric(dta$TP); fn <- as.numeric(dta$FN)
    fp <- as.numeric(dta$FP); tn <- as.numeric(dta$TN)
    study_df <- data.frame(fpr  = fp / (fp + tn),
                           sens = tp / (tp + fn))
    study_df <- study_df[is.finite(study_df$fpr) & is.finite(study_df$sens), ]
  }

  # ---- 6. Labels inside the plot -------------------------------------------
  k_n_lab <- if (!is.null(study_df)) {
    sprintf("k\u00a0=\u00a0%d", nrow(study_df))
  } else ""
  auc_txt <- if (!is.null(auc)) sprintf("AUC\u00a0=\u00a0%.2f", auc) else ""

  # ---- 7. Build ggplot ------------------------------------------------------
  p <- ggplot() +
    # reference diagonal (AUC = 0.5)
    geom_abline(slope = 1, intercept = 0,
                colour = BIST_GREY, linewidth = 0.35, linetype = "dashed") +
    # prediction ellipse (dashed border, no fill)
    { if (!is.null(pred_df))
        geom_polygon(data = pred_df, aes(x = fpr, y = sens),
                     fill = NA, colour = acc_col,
                     linetype = "dashed", linewidth = 0.6) } +
    # confidence ellipse (shaded)
    { if (!is.null(conf_df))
        geom_polygon(data = conf_df, aes(x = fpr, y = sens),
                     fill = acc_col, alpha = 0.15,
                     colour = acc_col, linewidth = 0.8) } +
    # SROC curve
    geom_line(data = sroc_pts, aes(x = fpr, y = sens),
              colour = acc_col, linewidth = 1.3) +
    # individual study points
    { if (!is.null(study_df))
        geom_point(data = study_df, aes(x = fpr, y = sens),
                   shape = 1, size = 2.6, colour = BIST_BLACK, stroke = 0.7) } +
    # summary operating point
    geom_point(data = sp, aes(x = fpr, y = sens),
               shape = 21, size = 4.0,
               fill = acc_col, colour = BIST_BLACK, stroke = 0.9) +
    # k label top-left
    annotate("text", x = 0.02, y = 0.97, label = k_n_lab,
             hjust = 0, vjust = 1, size = 3.0, colour = BIST_GREY) +
    # AUC label bottom-right
    annotate("text", x = 0.98, y = 0.04, label = auc_txt,
             hjust = 1, vjust = 0, size = 3.2, fontface = "bold",
             colour = acc_col) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       expand = c(0.01, 0.01)) +
    labs(title    = if (!is.null(title)) title else "",
         subtitle = cs$label,
         x = "1-Specificity", y = "Sensitivity") +
    theme_bist() +
    theme(panel.grid.major = element_line(colour = "#EBEBEB", linewidth = 0.3),
          axis.title       = element_text(size = 9),
          plot.subtitle    = element_text(colour = acc_col, size = 10))

  # save = FALSE returns the plot without writing files (used by the
  # composite figures in scripts/11).
  if (save) save_bist(p, file = file, width = 5.5, height = 5.5, dir = dir)
  invisible(p)
}

# Null-coalescing helper.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
