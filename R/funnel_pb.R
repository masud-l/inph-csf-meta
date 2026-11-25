# =============================================================================
# funnel_pb.R  —  Funnel plots + Egger's test for small-study effects
# -----------------------------------------------------------------------------
# Publication-bias assessment for the biomarkers with enough studies to make it
# meaningful. Funnel-asymmetry tests are unreliable below ~10 studies, so this
# is applied ONLY to biomarker x comparison analyses with k >= PB_MIN_K.
#
# PRIMARY test (revision): the modified Egger-type test of Pustejovsky &
# Rodgers (Res Synth Methods 2019;10:57-71), with sqrt(1/n1 + 1/n2) as the
# predictor. For SMDs the conventional standard error depends partly on the
# observed effect size (vi contains a g^2 term), which inflates the false-
# positive rate of the conventional test when effects are large.
# SENSITIVITY: the conventional Egger test, metafor::regtest() with the
# standard error as predictor.
#
# IMPORTANT interpretation note (for the manuscript, not the code): with the
# high I-squared seen throughout this dataset, funnel asymmetry may reflect
# genuine between-study heterogeneity rather than publication bias. A
# significant Egger test should NOT be reported as "publication bias present"
# without that caveat.
#
# Depends on meta_smd.R (run_smd_one output) and theme_bist.R. metafor only.
# =============================================================================

PB_MIN_K <- 10        # minimum studies for a funnel-asymmetry test to be valid

# ---- Egger-type regression tests for one SMD meta-analysis -----------------
# result : a run_smd_one() output (needs $fit and $es with En / Cn columns).
# Returns a one-row data frame; all-NA if the tests cannot be run.
#   egger_z / egger_p   : standard Egger test (predictor = SE), regtest()
#   limit_est (+ CI)    : regtest() "limit estimate" - the predicted effect
#                         as SE -> 0 (NOT the regression intercept term)
#   mod_egger_z / _p    : modified test, predictor sqrt(1/n1 + 1/n2)
egger_test <- function(result, biomarker, comparison) {
  fit <- result$fit
  na_row <- data.frame(biomarker = biomarker, comparison = comparison,
                       k = if (is.null(fit)) NA_integer_ else fit$k,
                       egger_z = NA_real_, egger_p = NA_real_,
                       limit_est = NA_real_, limit_lb = NA_real_,
                       limit_ub = NA_real_,
                       mod_egger_z = NA_real_, mod_egger_p = NA_real_,
                       testable = FALSE)
  if (is.null(fit) || fit$k < PB_MIN_K) return(na_row)

  safe1 <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    as.numeric(x[1])
  }

  # Standard Egger test (regtest default: predictor = sei, same model/method).
  rt <- tryCatch(metafor::regtest(fit), error = function(e) NULL)
  if (is.null(rt)) return(na_row)

  # Modified test: same random-effects meta-regression, predictor free of g.
  es  <- result$es
  w   <- sqrt(1 / es$En + 1 / es$Cn)
  mod <- tryCatch(
    metafor::rma(yi = as.numeric(es$yi), vi = as.numeric(es$vi), mods = ~ w,
                 method = fit$method, test = fit$test),
    error = function(e) NULL)

  data.frame(
    biomarker   = biomarker,
    comparison  = comparison,
    k           = fit$k,
    egger_z     = safe1(rt$zval),
    egger_p     = safe1(rt$pval),
    limit_est   = safe1(rt$est),
    limit_lb    = safe1(rt$ci.lb),
    limit_ub    = safe1(rt$ci.ub),
    mod_egger_z = if (is.null(mod)) NA_real_ else as.numeric(mod$zval[2]),
    mod_egger_p = if (is.null(mod)) NA_real_ else as.numeric(mod$pval[2]),
    testable    = TRUE
  )
}

# ---- Funnel plot for one fitted model ---------------------------------------
# Saves a BIST-styled funnel as PDF + TIFF (RASTER_DPI). Built on metafor's
# funnel() method, which draws the correct pseudo-confidence contours.
funnel_bist <- function(fit, biomarker_label_txt, comparison_label,
                        file, dir = PATHS$figures, header_color = BIST_BLACK) {
  if (is.null(fit) || fit$k < PB_MIN_K) {
    message("  funnel: ", biomarker_label_txt, " (", comparison_label,
            ") k < ", PB_MIN_K, " - skipped.")
    return(invisible(NULL))
  }
  draw <- function() {
    # Extra top margin so the two-line header has room above the plot box.
    op <- graphics::par(mar = c(4.6, 4.4, 4.0, 1.6))
    on.exit(graphics::par(op), add = TRUE)

    metafor::funnel(fit,
                    level = c(90, 95, 99),          # pseudo-CI contours
                    shade = c("white", "gray90", "gray80"),
                    refline = 0,
                    back = "white",
                    pch = 19, cex = 0.8,
                    col = BIST_BLACK,
                    xlab = "Standardized mean difference (Hedges' g)",
                    main = "")
    # Two-line header matching the forest / LOO style: bold black biomarker
    # name on top, comparison label below in the comparison's accent colour.
    graphics::mtext(biomarker_label_txt, side = 3, line = 2.0, adj = 0,
                    font = 2, cex = 1.05, col = BIST_BLACK)
    graphics::mtext(comparison_label, side = 3, line = 0.9, adj = 0,
                    font = 1, cex = 0.9, col = header_color)
  }
  base <- file.path(dir, file)
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = 6, height = 5.2)
  draw(); grDevices::dev.off()
  grDevices::tiff(paste0(base, ".tiff"), width = 6, height = 5.2,
                  units = "in", res = RASTER_DPI, compression = "lzw",
                  bg = "white")
  draw(); grDevices::dev.off()
  message("  saved: ", base, ".pdf + .tiff (", RASTER_DPI, " dpi)")
  invisible(base)
}
