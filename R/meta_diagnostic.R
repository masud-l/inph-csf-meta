# =============================================================================
# meta_diagnostic.R  —  Diagnostic test accuracy meta-analysis (mada)
# -----------------------------------------------------------------------------
# Bivariate (Reitsma) model for paired sensitivity / specificity, plus
# per-study exact (Clopper-Pearson) confidence intervals.
#
# SCOPE NOTE: only iNPH-vs-AD has enough studies to pool (ab42 k = 5,
# p-tau181 k = 5; revision counts, iNPH = positive condition). iNPH-vs-HC and
# sAPP-alpha each have a single study and are handled descriptively. The
# bivariate model needs >= 4 studies and is only stable at >= 10 - treat the
# pooled estimates as exploratory. When the between-study correlation reaches
# its boundary (|rho| > RHO_BOUNDARY), run_reitsma_one() reports separate
# univariate logit estimates and no AUC (p-tau181 in the revision).
#
# IMPLEMENTATION NOTE: the Reitsma model is parameterised on the logit scale
# in terms of sensitivity and the FALSE-POSITIVE RATE (FPR), not specificity.
# Specificity is therefore derived as 1 - FPR. Estimates are read from the
# fitted model's coefficients / vcov rather than from summary()'s printed
# coefficient matrix, whose row labels differ across mada versions.
# =============================================================================

# ---- Per-study sensitivity / specificity with exact CIs ---------------------
diagnostic_per_study <- function(d) {
  ci <- function(x, n) {
    bt <- binom.test(x, n)
    c(est = unname(bt$estimate), lo = bt$conf.int[1], hi = bt$conf.int[2])
  }
  sens <- t(mapply(ci, d$tp, d$tp + d$fn))
  spec <- t(mapply(ci, d$tn, d$tn + d$fp))
  data.frame(
    study      = d$study,
    biomarker  = d$biomarker,
    comparison = d$comparison,
    n_inph     = d$tp + d$fn,
    n_ctrl     = d$fp + d$tn,
    sens = sens[, "est"], sens_lo = sens[, "lo"], sens_hi = sens[, "hi"],
    spec = spec[, "est"], spec_lo = spec[, "lo"], spec_hi = spec[, "hi"],
    zero_cell = d$zero_cell
  )
}

# ---- AUC of the summary ROC curve, robust to mada version differences -------
.reitsma_auc <- function(fit) {
  tries <- list(
    function() { a <- summary(fit)$AUC; if (is.list(a)) a$AUC else a },
    function() { a <- summary(fit)$AUC; if (is.list(a)) a[["AUC"]] else a },
    function() { a <- mada::AUC(fit);   if (is.list(a)) a$AUC else a }
  )
  for (f in tries) {
    v <- tryCatch(as.numeric(f())[1], error = function(e) NA_real_)
    if (length(v) == 1 && !is.na(v)) return(v)
  }
  NA_real_
}

# ---- Back-transform pooled sensitivity / specificity from a reitsma fit -----
# Works from coef()/vcov(): the two fixed effects are logit-sensitivity and
# logit-FPR (default reitsma transformation, alphasens = alphafpr = 1). The
# intercept-only Reitsma model always has exactly these two coefficients, in
# the order (tsens, tfpr) - used as a positional fallback if names are absent.
.reitsma_estimates <- function(fit) {
  co <- tryCatch(coef(fit),  error = function(e) fit$coefficients)
  V  <- tryCatch(vcov(fit),  error = function(e) fit$vcov)
  if (is.null(co) || is.null(V))
    stop("could not obtain coefficients / vcov from the reitsma fit")
  nm  <- names(co)
  cov <- as.numeric(co)

  i_s <- if (!is.null(nm)) grep("sens", nm, ignore.case = TRUE)[1] else NA
  i_f <- if (!is.null(nm)) grep("fpr",  nm, ignore.case = TRUE)[1] else NA
  if (is.na(i_s)) i_s <- 1L          # positional fallback: 1st = tsens
  if (is.na(i_f)) i_f <- 2L          #                      2nd = tfpr

  z  <- qnorm(0.975)
  ts <- cov[i_s]; se_s <- sqrt(V[i_s, i_s])
  tf <- cov[i_f]; se_f <- sqrt(V[i_f, i_f])

  sens    <- plogis(ts)
  sens_ci <- plogis(c(ts - z * se_s, ts + z * se_s))
  fpr_ci  <- plogis(c(tf - z * se_f, tf + z * se_f))
  spec    <- 1 - plogis(tf)
  spec_ci <- sort(1 - fpr_ci)                 # bounds invert under 1 - x

  list(sensitivity = sens, sens_lo = sens_ci[1], sens_hi = sens_ci[2],
       specificity = spec, spec_lo = spec_ci[1], spec_hi = spec_ci[2])
}

# ---- Between-study correlation of a reitsma fit ------------------------------
# Correlation of the two random effects (logit sensitivity, logit FPR) from the
# fitted between-study covariance matrix Psi. NA if Psi is unavailable.
.reitsma_rho <- function(fit) {
  Psi <- tryCatch(fit$Psi, error = function(e) NULL)
  if (is.null(Psi) || !is.matrix(Psi) || any(dim(Psi) != 2)) return(NA_real_)
  den <- sqrt(Psi[1, 1] * Psi[2, 2])
  if (!is.finite(den) || den <= 0) return(NA_real_)
  as.numeric(Psi[1, 2] / den)
}

# |rho| above this value is treated as the parameter boundary (rho = +/-1):
# the bivariate summary curve and AUC are then not identifiable.
RHO_BOUNDARY <- 0.99

# ---- Separate univariate random-effects logit models (fallback) ------------
# Used when the bivariate correlation reaches its boundary. Same continuity
# rule as the bivariate fit: 0.5 added to every cell of a study table that
# contains a zero. Sensitivity = TP/(TP+FN), specificity = TN/(TN+FP);
# each pooled on the logit scale with REML and z-based CIs.
.univariate_logit <- function(d) {
  zero <- d$tp == 0 | d$fn == 0 | d$fp == 0 | d$tn == 0
  cc <- ifelse(zero, 0.5, 0)
  tp <- d$tp + cc; fn <- d$fn + cc; fp <- d$fp + cc; tn <- d$tn + cc
  pool <- function(x, y) {
    f <- metafor::rma(yi = log(x / y), vi = 1 / x + 1 / y,
                      method = "REML", test = "z")
    c(est = plogis(as.numeric(f$b)), lo = plogis(f$ci.lb),
      hi = plogis(f$ci.ub))
  }
  s <- pool(tp, fn)
  p <- pool(tn, fp)
  list(sensitivity = s[["est"]], sens_lo = s[["lo"]], sens_hi = s[["hi"]],
       specificity = p[["est"]], spec_lo = p[["lo"]], spec_hi = p[["hi"]])
}

# ---- Bivariate (Reitsma) pooled model for one biomarker ---------------------
# d : diagnostic rows for ONE biomarker x comparison.
# Returns the reitsma fit, the 2x2 data, and a tidy summary; fit NULL if k < 4.
run_reitsma_one <- function(d, min_k = 4) {
  k <- nrow(d)
  if (k < min_k) {
    message("  ", d$biomarker[1], " (", d$comparison[1], "): k = ", k,
            " < ", min_k, " - descriptive only, no bivariate model.")
    return(list(fit = NULL, k = k, summary = NULL, data = NULL))
  }
  # mada::reitsma expects columns TP, FN, FP, TN. A 0.5 continuity correction
  # is applied to any study containing a zero cell.
  dd  <- data.frame(TP = d$tp, FN = d$fn, FP = d$fp, TN = d$tn)
  fit <- mada::reitsma(dd, correction = 0.5, correction.control = "single")

  # Extract pooled sensitivity/specificity; on any unexpected object structure
  # dump a diagnostic and continue with NA rather than aborting the script.
  est <- tryCatch(.reitsma_estimates(fit), error = function(e) {
    message("  ! estimate extraction failed for ", d$biomarker[1],
            ": ", conditionMessage(e))
    message("    coef() names    : ",
            paste(names(stats::coef(fit)), collapse = ", "))
    message("    summary rownames: ",
            paste(tryCatch(rownames(summary(fit)$coefficients),
                           error = function(x) "<unavailable>"),
                  collapse = ", "))
    NULL
  })
  if (is.null(est))
    est <- list(sensitivity = NA, sens_lo = NA, sens_hi = NA,
                specificity = NA, spec_lo = NA, spec_hi = NA)

  # Boundary check. If the between-study correlation is at +/-1 the bivariate
  # model's summary ROC curve and AUC are not reliable: report separate
  # univariate logit estimates instead, no AUC, and return fit = NULL so that
  # no SROC curve is drawn (the bivariate fit is kept as fit_bivariate).
  rho      <- .reitsma_rho(fit)
  boundary <- !is.na(rho) && abs(rho) > RHO_BOUNDARY
  model    <- "bivariate (Reitsma)"
  auc      <- .reitsma_auc(fit)
  fit_out  <- fit
  if (is.na(rho))
    message("  ! ", d$biomarker[1], ": between-study correlation unavailable;",
            " boundary check not possible.")
  if (boundary) {
    message("  ", d$biomarker[1], ": bivariate correlation at boundary (rho = ",
            sprintf("%.3f", rho), ") - reporting separate univariate logit ",
            "models; no summary ROC curve or AUC.")
    est     <- .univariate_logit(d)
    model   <- "univariate logit (bivariate correlation at boundary)"
    auc     <- NA_real_
    fit_out <- NULL
  }

  summ <- data.frame(
    biomarker   = d$biomarker[1],
    comparison  = d$comparison[1],
    k           = k,
    sensitivity = est$sensitivity,
    sens_lo     = est$sens_lo,
    sens_hi     = est$sens_hi,
    specificity = est$specificity,
    spec_lo     = est$spec_lo,
    spec_hi     = est$spec_hi,
    AUC         = auc,
    model       = model,
    rho         = rho
  )
  list(fit = fit_out, fit_bivariate = fit, k = k, summary = summ, data = dd,
       boundary = boundary)
}

# ---- Heterogeneity of the univariate logit models ---------------------------
# tau^2 and I^2 for sensitivity and specificity pooled separately (same data
# handling as .univariate_logit). Used for the model-details table.
.univariate_heterogeneity <- function(d) {
  zero <- d$tp == 0 | d$fn == 0 | d$fp == 0 | d$tn == 0
  cc <- ifelse(zero, 0.5, 0)
  tp <- d$tp + cc; fn <- d$fn + cc; fp <- d$fp + cc; tn <- d$tn + cc
  het <- function(x, y) {
    f <- metafor::rma(yi = log(x / y), vi = 1 / x + 1 / y,
                      method = "REML", test = "z")
    c(tau2 = f$tau2, I2 = f$I2)
  }
  s <- het(tp, fn); p <- het(tn, fp)
  list(sens_tau2 = s[["tau2"]], sens_I2 = s[["I2"]],
       spec_tau2 = p[["tau2"]], spec_I2 = p[["I2"]])
}

# ---- Convergence flag of a reitsma fit ---------------------------------------
# mada fits the bivariate model with mvmeta-style code; the flag name differs
# between versions, so several are tried. NA if none is available.
.reitsma_converged <- function(fit) {
  tries <- list(function() fit$converged,
                function() fit$convergence == 0,
                function() fit$optim$convergence == 0)
  for (f in tries) {
    v <- tryCatch(as.logical(f())[1], error = function(e) NA)
    if (length(v) == 1 && !is.na(v)) return(v)
  }
  NA
}

# ---- Model estimation details for one pooled analysis (Reviewer 1.5) --------
# r : run_reitsma_one() output; d : the 2x2 rows it was fitted to.
reitsma_model_details <- function(r, d, analysis) {
  fb  <- r$fit_bivariate
  Psi <- tryCatch(fb$Psi, error = function(e) NULL)
  ok  <- is.matrix(Psi) && all(dim(Psi) == 2)
  het <- .univariate_heterogeneity(d)
  zero_studies <- d$study[d$tp == 0 | d$fn == 0 | d$fp == 0 | d$tn == 0]
  data.frame(
    analysis               = analysis,
    biomarker              = d$biomarker[1],
    k                      = nrow(d),
    estimator              = "REML",
    reported_model         = r$summary$model,
    bivariate_converged    = .reitsma_converged(fb),
    bivariate_rho          = r$summary$rho,
    rho_at_boundary        = isTRUE(r$boundary),
    bivariate_tau2_logit_sens = if (ok) Psi[1, 1] else NA_real_,
    bivariate_tau2_logit_fpr  = if (ok) Psi[2, 2] else NA_real_,
    univariate_tau2_sens   = het$sens_tau2,
    univariate_I2_sens     = het$sens_I2,
    univariate_tau2_spec   = het$spec_tau2,
    univariate_I2_spec     = het$spec_I2,
    continuity_correction  = if (length(zero_studies))
      paste0("0.5 added to every cell: ", paste(zero_studies, collapse = "; "))
      else "none needed",
    stringsAsFactors = FALSE)
}

# ---- Pool one sensitivity scenario --------------------------------------------
# k >= 4: run_reitsma_one() (bivariate, univariate fallback at the boundary);
# k == 3: separate univariate logit models (too few for the bivariate model);
# k < 3 : not pooled.
pool_dta_scenario <- function(d, scenario, biomarker,
                              comparison = "iNPH_vs_AD") {
  k <- nrow(d)
  out <- data.frame(scenario = scenario, biomarker = biomarker,
                    comparison = comparison, k = k,
                    studies = paste(d$study, collapse = "; "),
                    sensitivity = NA_real_, sens_lo = NA_real_,
                    sens_hi = NA_real_, specificity = NA_real_,
                    spec_lo = NA_real_, spec_hi = NA_real_, AUC = NA_real_,
                    model = NA_character_, rho = NA_real_,
                    stringsAsFactors = FALSE)
  if (k < 3) {
    out$model <- sprintf("not pooled (k < 3)")
    return(out)
  }
  if (k == 3) {
    est <- .univariate_logit(d)
    out[c("sensitivity", "sens_lo", "sens_hi",
          "specificity", "spec_lo", "spec_hi")] <- unlist(est)
    out$model <- "univariate logit (k = 3)"
    return(out)
  }
  r <- suppressMessages(run_reitsma_one(d))
  s <- r$summary
  out[c("sensitivity", "sens_lo", "sens_hi", "specificity", "spec_lo",
        "spec_hi", "AUC", "model", "rho")] <-
    s[c("sensitivity", "sens_lo", "sens_hi", "specificity", "spec_lo",
        "spec_hi", "AUC", "model", "rho")]
  out
}

# ---- Leave-one-out for an influential study (e.g. a perfect classifier) -----
reitsma_drop_study <- function(d, drop_pattern) {
  keep <- d[!grepl(drop_pattern, d$study, ignore.case = TRUE), , drop = FALSE]
  message("  leave-one-out: dropping ", nrow(d) - nrow(keep),
          " study/studies matching '", drop_pattern, "'")
  run_reitsma_one(keep)
}

# (Revision) smd_to_auc() removed - SMD-derived AUC is no longer reported.
