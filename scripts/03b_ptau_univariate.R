# =============================================================================
# Step 03b - Univariate pooling of sensitivity and specificity (revision)
# Script: 03b_ptau_univariate.R
# =============================================================================
#
# Why this script exists
# ----------------------
# For p-tau181 (iNPH vs AD) the bivariate Reitsma model reaches the boundary
# of the between-study correlation (rho = -1). Its summary ROC curve and AUC
# are therefore not reported (the earlier AUC of 0.96 was withdrawn). Instead,
# sensitivity and specificity are pooled in separate random-effects logit
# models (REML), as reported in Supplementary Table S5. Aβ42 is also pooled
# this way here as a cross-check; its primary estimate remains the bivariate
# model in 03.
#
# Conventions (identical to Supplementary Table S5)
# -------------------------------------------------
# - iNPH is the positive condition: tp + fn = iNPH, fp + tn = AD.
# - 0.5 is added to every cell of a study table that contains a zero.
# - Logit sensitivity = log(tp/fn), variance 1/tp + 1/fn.
#   Logit specificity = log(tn/fp), variance 1/tn + 1/fp.
#
# Inputs : data/diagnostic_performance.xlsx (Sheet1)
# Outputs: output/tables/diagnostic_univariate_pooled.csv
#          output/tables/diagnostic_univariate_loo.csv
# Run    : after 01_validate_data.R (needs R/00_setup.R for PATHS)
# =============================================================================

source("R/00_setup.R")
suppressPackageStartupMessages({ library(readxl); library(metafor) })

d <- readxl::read_excel(file.path(PATHS$data, "diagnostic_performance.xlsx"),
                        sheet = "Sheet1")
d <- as.data.frame(d)
d <- d[d$Comparison == "iNPH_vs_AD" & d$Biomarker %in% c("ab42", "p-tau181"), ]

# Full-table continuity correction for any study containing a zero cell
zero <- with(d, tp == 0 | fn == 0 | fp == 0 | tn == 0)
for (v in c("tp", "fn", "fp", "tn")) d[[v]] <- d[[v]] + ifelse(zero, 0.5, 0)

pool_logit <- function(x, y) {
  yi <- log(x / y)
  vi <- 1 / x + 1 / y
  fit <- metafor::rma(yi = yi, vi = vi, method = "REML", test = "z")
  c(k = fit$k,
    est = plogis(as.numeric(fit$b)),
    lo  = plogis(fit$ci.lb),
    hi  = plogis(fit$ci.ub),
    tau2 = fit$tau2,
    I2  = fit$I2)
}

run_one <- function(dd) {
  s <- pool_logit(dd$tp, dd$fn)   # sensitivity
  p <- pool_logit(dd$tn, dd$fp)   # specificity
  data.frame(k = s[["k"]],
             sensitivity = s[["est"]], sens_lo = s[["lo"]], sens_hi = s[["hi"]],
             sens_tau2 = s[["tau2"]], sens_I2 = s[["I2"]],
             specificity = p[["est"]], spec_lo = p[["lo"]], spec_hi = p[["hi"]],
             spec_tau2 = p[["tau2"]], spec_I2 = p[["I2"]])
}

res <- do.call(rbind, lapply(split(d, d$Biomarker), function(dd)
  cbind(biomarker = dd$Biomarker[1], comparison = "iNPH_vs_AD", run_one(dd))))

loo <- do.call(rbind, lapply(split(d, d$Biomarker), function(dd)
  do.call(rbind, lapply(seq_len(nrow(dd)), function(i)
    cbind(biomarker = dd$Biomarker[1], omitted = dd$Study[i],
          run_one(dd[-i, , drop = FALSE]))))))

pct <- function(df) {
  for (v in grep("^(sens|spec)", names(df), value = TRUE))
    if (!grepl("tau2|I2", v)) df[[v]] <- round(100 * df[[v]], 1)
  df
}
res <- pct(res); loo <- pct(loo)

write.csv(res, file.path(PATHS$tables, "diagnostic_univariate_pooled.csv"),
          row.names = FALSE)
write.csv(loo, file.path(PATHS$tables, "diagnostic_univariate_loo.csv"),
          row.names = FALSE)

print(res[, c("biomarker", "k", "sensitivity", "sens_lo", "sens_hi",
              "specificity", "spec_lo", "spec_hi")], row.names = FALSE)
message("\nExpected p-tau181 (Supplementary Table S5): sensitivity 90.8 (78.1-96.5), ",
        "specificity 87.6 (82.9-91.1).")
message("Expected Abeta42 univariate cross-check: sensitivity 52.9 (30.5-74.2), ",
        "specificity 83.3 (74.9-89.3). Primary Abeta42 = bivariate in script 03 ",
        "(52.4 / 82.1).")
