# =============================================================================
# meta_smd.R  —  Standardized mean difference meta-analysis (metafor)
# -----------------------------------------------------------------------------
# Wraps metafor::escalc + rma for one biomarker at a time, plus a batch runner
# that loops over every biomarker in a tidy data frame.
# =============================================================================

# ---- Single-biomarker meta-analysis -----------------------------------------
# df : tidy data frame for ONE biomarker x comparison (cols study,En,Em,Esd,
#      Cn,Cm,Csd). method / knha default to the globals set in 00_setup.R.
# Returns a list with the rma fit, the per-study escalc data, and a tidy summary
# row. Single-study biomarkers return est/CI from escalc with fit = NULL.
run_smd_one <- function(df, method = RE_METHOD, knha = USE_KNHA) {
  es <- metafor::escalc(
    measure = "SMD",
    m1i = df$Em, sd1i = df$Esd, n1i = df$En,   # iNPH  (group 1)
    m2i = df$Cm, sd2i = df$Csd, n2i = df$Cn,   # control (group 2)
    slab = df$study, data = df
  )
  k <- nrow(es)

  if (k < MIN_K) {
    # Not pooled. The summary row keeps the TRUE study count. With one study
    # the row carries that study's own estimate; with two or more there is no
    # single estimate, so estimate/CI are NA and the individual studies are
    # written to smd_descriptive_studies.csv by scripts/02 (and shown in the
    # forest plots).
    one <- k == 1
    se  <- if (one) sqrt(es$vi[1]) else NA_real_
    summ <- data.frame(
      biomarker = df$biomarker[1], comparison = df$comparison[1], k = k,
      estimate = if (one) es$yi[1] else NA_real_,
      ci_lb = if (one) es$yi[1] - 1.96 * se else NA_real_,
      ci_ub = if (one) es$yi[1] + 1.96 * se else NA_real_, se = se,
      I2 = NA_real_, tau2 = NA_real_, QEp = NA_real_, pval = NA_real_,
      note = sprintf("k = %d (< MIN_K = %d) - descriptive, not pooled%s",
                     k, MIN_K,
                     if (one) "" else "; see smd_descriptive_studies.csv"),
      pi_lb = NA_real_, pi_ub = NA_real_
    )
    return(list(fit = NULL, es = es, summary = summ))
  }

  fit <- metafor::rma(yi, vi, data = es, method = method,
                      test = if (knha) "knha" else "z")
  summ <- data.frame(
    biomarker = df$biomarker[1], comparison = df$comparison[1], k = k,
    estimate = as.numeric(fit$b), ci_lb = fit$ci.lb, ci_ub = fit$ci.ub,
    se = fit$se, I2 = fit$I2, tau2 = fit$tau2, QEp = fit$QEp,
    pval = fit$pval, note = ""
  )
  pi <- smd_prediction_interval(fit)
  summ$pi_lb <- pi[1]
  summ$pi_ub <- pi[2]
  list(fit = fit, es = es, summary = summ)
}

# ---- 95% prediction interval -------------------------------------------------
# Range expected to contain the effect of a new comparable study (Higgins et
# al. 2009; Cochrane Handbook): estimate +/- t(0.975, k - 2) * sqrt(tau^2 + SE^2).
# Added in the revision (Reviewer 1.7). Needs k >= 3.
smd_prediction_interval <- function(fit) {
  if (is.null(fit) || fit$k < 3) return(c(NA_real_, NA_real_))
  half <- qt(0.975, df = fit$k - 2) * sqrt(fit$tau2 + fit$se^2)
  as.numeric(fit$b) + c(-1, 1) * half
}

# ---- Batch runner over all biomarkers ---------------------------------------
# dat : tidy data frame possibly spanning many biomarkers (one comparison).
# Returns a list keyed by biomarker, each element the run_smd_one() output.
run_smd_all <- function(dat, method = RE_METHOD, knha = USE_KNHA) {
  split(dat, dat$biomarker) |>
    lapply(function(df) {
      out <- tryCatch(run_smd_one(df, method, knha),
                      error = function(e) {
                        message("  ! ", df$biomarker[1], ": ", conditionMessage(e))
                        NULL
                      })
      out
    }) |>
    purrr::compact()
}

# ---- Collect all summary rows into one table --------------------------------
smd_summary_table <- function(results) {
  purrr::map_dfr(results, "summary") |>
    arrange(comparison, desc(abs(estimate)))
}

# ---- Study-level estimates for outcomes that were not pooled ----------------
# One row per study for every biomarker x comparison with k < MIN_K.
smd_descriptive_table <- function(results) {
  rows <- lapply(results, function(r) {
    if (!is.null(r$fit)) return(NULL)
    es <- r$es
    se <- sqrt(as.numeric(es$vi))
    data.frame(biomarker = r$summary$biomarker,
               comparison = r$summary$comparison,
               k = r$summary$k,
               study = es$study,
               estimate = as.numeric(es$yi),
               ci_lb = as.numeric(es$yi) - 1.96 * se,
               ci_ub = as.numeric(es$yi) + 1.96 * se,
               se = se)
  })
  dplyr::bind_rows(rows)
}

# ---- Influence / leave-one-out diagnostics ----------------------------------
# Useful for the biomarkers driving the dual-purpose claim (Abeta42, NFL).
smd_leave1out <- function(result) {
  if (is.null(result$fit)) return(NULL)
  loo <- metafor::leave1out(result$fit)
  data.frame(omitted = result$es$study,
             estimate = loo$estimate, ci_lb = loo$ci.lb, ci_ub = loo$ci.ub,
             I2 = loo$I2)
}
