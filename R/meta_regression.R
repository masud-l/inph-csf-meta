# =============================================================================
# meta_regression.R  —  Age meta-regression (exploratory, study-level)
# -----------------------------------------------------------------------------
# Tests whether the per-study SMD varies with study-mean iNPH-group age, for a
# single biomarker x comparison (pre-specified: AB42, iNPH vs AD).
#
# This is a CONTINUOUS moderator analysis (meta-regression), NOT a subgroup /
# binning analysis - age enters as a continuous covariate, so there are no
# arbitrary low/moderate/high cutoffs.
#
# Hypothesis (pre-stated): older iNPH cohorts carry more incidental amyloid
# co-pathology, lowering their CSF AB42 toward AD levels, so the iNPH-vs-AD SMD
# should rise (become less negative / more positive) with increasing mean age.
#
# CAVEATS (state in the manuscript):
#   * study-level (ecological) association - no individual-patient inference;
#   * mean age confounded with publication era, assay generation, criteria;
#   * likely leaves most heterogeneity unexplained (a contributor, not a cause).
#
# Depends on meta_smd.R (escalc usage) and theme_bist.R. metafor + ggplot2.
# =============================================================================

# ---- Fit the meta-regression ------------------------------------------------
# dat        : clean per-study data (clean_data.rds$smd), must carry iNPH_age.
# biomarker  : raw biomarker name, e.g. "AB42".
# comparison : "iNPH_vs_AD" (the pre-specified comparison).
# knha       : Knapp-Hartung adjustment. Mirrors USE_KNHA for consistency with
#              the main analysis; note that test = "knha" is often recommended
#              specifically for meta-regression and can be switched on here.
# Returns a list: fit, data (with yi/vi/age), slope, ci, p, R2, k, n_dropped.
run_age_metareg <- function(dat, biomarker = "AB42",
                            comparison = "iNPH_vs_AD",
                            method = RE_METHOD, knha = USE_KNHA) {

  d <- dat[dat$biomarker == biomarker & dat$comparison == comparison, ,
           drop = FALSE]
  if (nrow(d) == 0) stop("no rows for ", biomarker, " / ", comparison)

  n_total <- nrow(d)
  # complete-case: drop studies missing age (NEVER impute a covariate)
  d <- d[!is.na(d$iNPH_age), , drop = FALSE]
  n_drop <- n_total - nrow(d)
  if (nrow(d) < 10)
    warning("only ", nrow(d), " studies with age - meta-regression unstable")

  # effect sizes (Hedges g); iNPH = group 1, comparator = group 2
  es <- metafor::escalc(measure = "SMD",
                        m1i = d$Em, sd1i = d$Esd, n1i = d$En,
                        m2i = d$Cm, sd2i = d$Csd, n2i = d$Cn,
                        data = d)

  fit <- metafor::rma(yi, vi, mods = ~ iNPH_age, data = es,
                      method = method,
                      test = if (knha) "knha" else "z")

  # slope is the 2nd coefficient (1st is intercept)
  b   <- as.numeric(fit$b)[2]
  ci  <- c(fit$ci.lb[2], fit$ci.ub[2])
  p   <- fit$pval[2]
  R2  <- if (!is.null(fit$R2)) fit$R2 else NA_real_

  list(fit = fit, data = es, slope = b, ci = ci, p = p, R2 = R2,
       k = nrow(es), n_dropped = n_drop,
       biomarker = biomarker, comparison = comparison)
}

# ---- Leave-one-out sensitivity on the slope ---------------------------------
# Re-fits the meta-regression dropping each study in turn; returns the range of
# slopes and flags whether the sign/significance is stable. This guards against
# a slope driven by one or two high-leverage (e.g. oldest-cohort) studies.
age_metareg_loo <- function(reg) {
  es <- reg$data
  rows <- lapply(seq_len(nrow(es)), function(i) {
    fit_i <- tryCatch(
      metafor::rma(yi, vi, mods = ~ iNPH_age, data = es[-i, ],
                   method = reg$fit$method,
                   test = reg$fit$test),
      error = function(e) NULL)
    if (is.null(fit_i)) return(NULL)
    data.frame(omitted   = es$study[i],
               age       = es$iNPH_age[i],
               slope     = as.numeric(fit_i$b)[2],
               p         = fit_i$pval[2])
  })
  out <- do.call(rbind, rows)
  out[order(-abs(out$slope - reg$slope)), ]   # most influential first
}

# ---- Subgroup analysis by age band ------------------------------------------
# Reproduces Dr. Nila's Stata decade-band subgroup analysis: studies are binned
# into age groups, each band is pooled separately (random-effects), and a test
# of group differences (Q_between) is reported. NOTE: included at the lead
# author's request; with very uneven band sizes (e.g. 2 / 24 / 3) the small
# bands are not reliably poolable and must be read with great caution -
# meta-regression (continuous age) is the statistically sounder tool.
#
# breaks/labels define the bands; default = decade bands 60-69, 70-79, 80-89.
run_age_subgroup <- function(dat, biomarker = "AB42",
                             comparison = "iNPH_vs_AD",
                             breaks = c(60, 70, 80, 90),
                             labels = c("60-69", "70-79", "80-89"),
                             method = RE_METHOD, knha = USE_KNHA) {

  d <- dat[dat$biomarker == biomarker & dat$comparison == comparison, ,
           drop = FALSE]
  d <- d[!is.na(d$iNPH_age), , drop = FALSE]
  if (nrow(d) == 0) stop("no rows with age for ", biomarker, "/", comparison)

  es <- metafor::escalc(measure = "SMD",
                        m1i = d$Em, sd1i = d$Esd, n1i = d$En,
                        m2i = d$Cm, sd2i = d$Csd, n2i = d$Cn,
                        data = d)
  # droplevels(): an empty band (e.g. no 80-89 study) would otherwise enter
  # the moderator model as a redundant predictor.
  es$age_band <- droplevels(cut(es$iNPH_age, breaks = breaks, labels = labels,
                                right = FALSE))

  test_arg <- if (knha) "knha" else "z"

  # per-band random-effects pool
  bands <- lapply(labels, function(L) {
    sub <- es[es$age_band == L & !is.na(es$age_band), ]
    if (nrow(sub) == 0) return(NULL)
    if (nrow(sub) == 1) {
      return(data.frame(band = L, k = 1,
                        estimate = sub$yi, ci_lb = NA, ci_ub = NA,
                        I2 = NA, tau2 = NA, note = "single study - not pooled"))
    }
    f <- metafor::rma(yi, vi, data = sub, method = method, test = test_arg)
    data.frame(band = L, k = nrow(sub),
               estimate = as.numeric(f$b), ci_lb = f$ci.lb, ci_ub = f$ci.ub,
               I2 = f$I2, tau2 = f$tau2, note = "")
  })
  band_tbl <- do.call(rbind, bands)

  # overall pool + test of group differences (Q_between)
  overall <- metafor::rma(yi, vi, data = es, method = method, test = test_arg)
  qbet <- tryCatch({
    fm <- metafor::rma(yi, vi, mods = ~ age_band, data = es,
                       method = method, test = test_arg)
    # Degrees of freedom of the omnibus moderator test = number of bands - 1.
    list(QM = fm$QM, QMp = fm$QMp,
         df = if (!is.null(fm$QMdf)) fm$QMdf[1] else fm$m)
  }, error = function(e) list(QM = NA, QMp = NA, df = NA))

  list(bands = band_tbl, overall = overall, qbet = qbet,
       data = es, biomarker = biomarker, comparison = comparison)
}

# ---- Subgroup forest plot ---------------------------------------------------
# Renders a BIST-style subgroup forest matching forest_bist(): studies grouped
# under age-band headers, a subtotal diamond per band, an overall diamond, and
# the comparison accent colour. Built from a run_age_subgroup() result.
subgroup_forest <- function(sg, title = NULL) {
  es  <- sg$data
  cs  <- comparison_style(sg$comparison)
  col <- cs$color
  labels <- levels(es$age_band)

  # Assemble plotting rows top-to-bottom: for each band -> header, studies,
  # subtotal diamond. Then overall diamond at the bottom. y decreases downward;
  # we build with increasing 'row' then flip.
  rows <- list(); r <- 0
  diamonds <- list()
  band_hdr <- list()
  for (L in labels) {
    sub <- es[!is.na(es$age_band) & es$age_band == L, ]
    if (nrow(sub) == 0) next
    r <- r + 1; band_hdr[[L]] <- r                      # band header row
    sub <- sub[order(sub$yi), ]
    for (i in seq_len(nrow(sub))) {
      r <- r + 1
      rows[[length(rows) + 1]] <- data.frame(
        row = r, study = sub$study[i], est = sub$yi[i],
        lb = sub$yi[i] - 1.96*sqrt(sub$vi[i]),
        ub = sub$yi[i] + 1.96*sqrt(sub$vi[i]),
        wt = 1/sub$vi[i], type = "study")
    }
    # subtotal
    bt <- sg$bands[sg$bands$band == L, ]
    if (nrow(bt) == 1 && !is.na(bt$ci_lb)) {
      r <- r + 1
      diamonds[[length(diamonds)+1]] <- data.frame(
        row = r, est = bt$estimate, lb = bt$ci_lb, ub = bt$ci_ub,
        label = sprintf("Subtotal (k=%d); I\u00b2=%s",
                        bt$k,
                        ifelse(is.na(bt$I2), "NA",
                               sprintf("%.0f%%", bt$I2))),
        type = "subtotal")
    }
  }
  # overall
  r <- r + 2
  ov <- sg$overall
  diamonds[[length(diamonds)+1]] <- data.frame(
    row = r, est = as.numeric(ov$b), lb = ov$ci.lb, ub = ov$ci.ub,
    label = sprintf("Overall; I\u00b2=%s",
                    ifelse(is.null(ov$I2) || is.na(ov$I2), "NA",
                           sprintf("%.0f%%", ov$I2))),
    type = "overall")

  study_df <- do.call(rbind, rows)
  diam_df  <- do.call(rbind, diamonds)
  rmax <- r
  # flip so row 1 is at top
  yof <- function(x) rmax - x + 1
  study_df$y <- yof(study_df$row)
  diam_df$y  <- yof(diam_df$row)
  band_y <- sapply(band_hdr, yof)

  study_df$box <- scales::rescale(sqrt(study_df$wt), to = c(2, 6))

  # diamond polygons
  poly <- do.call(rbind, lapply(seq_len(nrow(diam_df)), function(i) {
    d <- diam_df[i, ]; h <- if (d$type == "overall") 0.34 else 0.28
    data.frame(grp = i,
               x = c(d$lb, d$est, d$ub, d$est),
               y = c(d$y,  d$y+h, d$y,  d$y-h),
               type = d$type)
  }))

  x_rng <- range(c(study_df$lb, study_df$ub, diam_df$lb, diam_df$ub, 0))
  x_pad <- diff(x_rng) * 0.08
  x_lim <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)
  y_lim <- c(0.2, rmax + 1.6)

  # numeric labels
  study_df$lab <- sprintf("%.2f (%.2f to %.2f)",
                          study_df$est, study_df$lb, study_df$ub)
  diam_df$lab  <- ifelse(is.na(diam_df$lb), "",
                         sprintf("%.2f (%.2f to %.2f)",
                                 diam_df$est, diam_df$lb, diam_df$ub))
  header_y <- rmax + 1.0

  # Panel 1: labels (band headers bold, studies indented)
  p_lab <- ggplot() +
    annotate("text", x = 0, y = header_y, label = "Study / age band",
             hjust = 0, fontface = "bold", size = 3.0, colour = BIST_BLACK) +
    geom_text(data = study_df, aes(x = 0.04, y = y, label = study),
              hjust = 0, size = 3.0, colour = BIST_BLACK) +
    { if (length(band_y))
        annotate("text", x = 0, y = band_y, label = names(band_y),
                 hjust = 0, fontface = "bold", size = 3.0, colour = col) } +
    geom_text(data = diam_df, aes(x = 0.02, y = y, label = label),
              hjust = 0, fontface = "bold", size = 2.6, colour = BIST_BLACK) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    theme_void(base_family = BIST_FONT) +
    theme(plot.margin = margin(12, 2, 12, 8))

  # Panel 2: forest
  p_for <- ggplot() +
    geom_vline(xintercept = 0, colour = BIST_GREY, linewidth = 0.4) +
    geom_vline(xintercept = as.numeric(ov$b), colour = col,
               linetype = "dotted", linewidth = 0.4) +
    geom_segment(data = study_df, aes(x = lb, xend = ub, y = y, yend = y),
                 colour = BIST_BLACK, linewidth = 0.5) +
    geom_point(data = study_df, aes(x = est, y = y, size = box),
               shape = 15, colour = BIST_BLACK) +
    scale_size_identity() +
    geom_polygon(data = poly, aes(x = x, y = y, group = grp, fill = type,
                                  colour = type)) +
    scale_fill_manual(values = c(subtotal = col, overall = BIST_BLACK),
                      guide = "none") +
    scale_colour_manual(values = c(subtotal = col, overall = BIST_BLACK),
                        guide = "none") +
    annotate("text", x = 0, y = header_y, label = "", hjust = 0) +
    scale_x_continuous(limits = x_lim) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    labs(x = "Standardized mean difference (Hedges' g)") +
    theme_bist() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(12, 4, 12, 4))

  # Panel 3: numeric SMD (95% CI)
  p_num <- ggplot() +
    annotate("text", x = 0, y = header_y, label = "SMD (95% CI)",
             hjust = 0, fontface = "bold", size = 3.0, colour = BIST_BLACK) +
    geom_text(data = study_df, aes(x = 0, y = y, label = lab),
              hjust = 0, size = 3.0, colour = BIST_BLACK) +
    geom_text(data = diam_df, aes(x = 0, y = y, label = lab),
              hjust = 0, fontface = "bold", size = 2.9, colour = BIST_BLACK) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    theme_void(base_family = BIST_FONT) +
    theme(plot.margin = margin(12, 8, 12, 2))

  bm_lab <- biomarker_label(sg$biomarker)
  qtxt <- sprintf("Test for subgroup differences: Q=%.2f, df=%s, P=%s",
                  sg$qbet$QM, sg$qbet$df,
                  ifelse(is.na(sg$qbet$QMp), "NA",
                         ifelse(sg$qbet$QMp < .001, "<.001",
                                sub("^0","",sprintf("%.3f", sg$qbet$QMp)))))

  patchwork::wrap_plots(p_lab, p_for, p_num, widths = c(3.0, 4.0, 2.4)) +
    patchwork::plot_annotation(
      title = if (!is.null(title)) title else bm_lab,
      subtitle = cs$label,
      caption = qtxt,
      theme = theme_bist() +
        theme(plot.subtitle = element_text(colour = col, size = 10),
              plot.caption  = element_text(hjust = 0, size = 8,
                                           colour = BIST_GREY)))
}

# ---- Binned bubble plot (Dr. Nila's Stata-style) ----------------------------
# Companion to the SUBGROUP analysis: x-axis is the three discrete age bands
# (not continuous age), points jittered within each band column, bubble area
# ~ inverse-variance weight, with the subgroup linear prediction across bands.
# This visualises the binned data; the continuous meta-regression is the
# separate bubble_plot() below. Build from a run_age_subgroup() result.
bubble_plot_binned <- function(sg, title = NULL) {
  es  <- sg$data
  es  <- es[!is.na(es$age_band), ]
  cs  <- comparison_style(sg$comparison)
  col <- cs$color
  es$weight <- 1 / es$vi

  # numeric x position = band index (1,2,3); points jittered within the column
  es$xpos <- as.integer(es$age_band)
  bands   <- levels(es$age_band)

  # linear prediction across band index (mirrors a regression on band number,
  # matching the Stata "Linear prediction" line over the three categories)
  fit_lin <- lm(yi ~ xpos, data = es, weights = weight)
  xs <- seq(1, length(bands), length.out = 100)
  pr <- predict(fit_lin, newdata = data.frame(xpos = xs),
                interval = "confidence")
  line_df <- data.frame(x = xs, fit = pr[, "fit"],
                        lo = pr[, "lwr"], hi = pr[, "upr"])

  bm_lab <- biomarker_label(sg$biomarker)

  ggplot() +
    geom_hline(yintercept = 0, colour = BIST_GREY, linewidth = 0.3) +
    geom_ribbon(data = line_df, aes(x = x, ymin = lo, ymax = hi),
                fill = col, alpha = 0.12) +
    geom_line(data = line_df, aes(x = x, y = fit),
              colour = col, linewidth = 1.0) +
    geom_point(data = es,
               aes(x = xpos, y = yi, size = weight),
               position = position_jitter(width = 0.12, height = 0, seed = 1),
               shape = 21, fill = col, colour = BIST_BLACK,
               alpha = 0.6, stroke = 0.5) +
    scale_size_area(max_size = 9, guide = "none") +
    scale_x_continuous(breaks = seq_along(bands),
                       labels = paste0("Age ", bands),
                       limits = c(0.5, length(bands) + 0.5)) +
    labs(title = if (!is.null(title)) title else bm_lab,
         subtitle = cs$label,
         x = NULL,
         y = "Standardized mean difference (Hedges' g)",
         caption = "Weights: inverse-variance") +
    theme_bist() +
    theme(plot.subtitle = element_text(colour = col, size = 10),
          plot.caption  = element_text(hjust = 0, size = 8, colour = BIST_GREY),
          panel.grid.major.x = element_blank())
}

# ---- Bubble plot (continuous meta-regression) -------------------------------
# x = study-mean age, y = per-study SMD, point area ~ study weight (1/vi),
# with the fitted meta-regression line and 95% CI band.
bubble_plot <- function(reg, title = NULL) {
  es <- reg$data
  es$weight <- 1 / es$vi
  cs  <- comparison_style(reg$comparison)
  col <- cs$color

  # prediction line across the observed age range
  age_grid <- data.frame(
    iNPH_age = seq(min(es$iNPH_age), max(es$iNPH_age), length.out = 100))
  pr <- predict(reg$fit, newmods = age_grid$iNPH_age)
  line_df <- data.frame(age = age_grid$iNPH_age,
                        fit = pr$pred, lo = pr$ci.lb, hi = pr$ci.ub)

  bm_lab  <- biomarker_label(reg$biomarker)
  slope_txt <- sprintf("Slope = %.3f per year (95%% CI %.3f to %.3f); P = %s",
                       reg$slope, reg$ci[1], reg$ci[2],
                       ifelse(reg$p < 0.001, "<.001", sub("^0","",sprintf("%.3f", reg$p))))

  ggplot() +
    geom_hline(yintercept = 0, colour = BIST_GREY, linewidth = 0.35) +
    geom_ribbon(data = line_df, aes(x = age, ymin = lo, ymax = hi),
                fill = col, alpha = 0.12) +
    geom_line(data = line_df, aes(x = age, y = fit),
              colour = col, linewidth = 1.0) +
    geom_point(data = es, aes(x = iNPH_age, y = yi, size = weight),
               shape = 21, fill = col, colour = BIST_BLACK,
               alpha = 0.6, stroke = 0.5) +
    scale_size_area(max_size = 9, guide = "none") +
    labs(title = if (!is.null(title)) title else bm_lab,
         subtitle = cs$label,
         x = "Study mean iNPH age, years",
         y = "Standardized mean difference (Hedges' g)",
         caption = slope_txt) +
    theme_bist() +
    theme(plot.subtitle = element_text(colour = col, size = 10),
          plot.caption  = element_text(hjust = 0, size = 8, colour = BIST_GREY))
}
