# =============================================================================
# Step 03 - Bivariate diagnostic test accuracy meta-analysis
# Script: 03_diagnostic_metaanalysis.R
# =============================================================================
#
# Project
# -------
# Cerebrospinal fluid biomarkers distinguishing idiopathic normal pressure
# hydrocephalus (iNPH) from healthy controls and Alzheimer disease:
# a systematic review and meta-analysis (41 studies).
#
# Purpose
# -------
# Pool empirical sensitivity, specificity, and AUC for biomarkers with
# sufficient studies reporting paired 2x2 data. In the included literature
# this was possible only for iNPH vs Alzheimer disease (AB42, k = 5;
# p-tau181, k = 5). iNPH is the positive condition (tp + fn = iNPH arm,
# fp + tn = AD arm; counts harmonised for the revision, Supplementary Table S5).
#
# Methods
# -------
# - Bivariate random-effects (Reitsma) model via mada::reitsma().
# - Continuity correction of 0.5 for studies with a zero cell.
# - Summary operating point, SROC curve, and pooled AUC.
# - If the between-study correlation reaches its boundary (|rho| > 0.99),
#   sensitivity and specificity are pooled in separate univariate random-
#   effects logit models and no SROC curve or AUC is reported (p-tau181).
# - Sensitivity analysis: p-tau181 with Kawamura et al. 2021 added back.
#   Kawamura is held out of the primary p-tau181 analysis because study
#   eligibility depended on p-tau181; its 2x2 row is on the "Not_primary"
#   sheet of diagnostic_performance.xlsx.
# - (Revision) The SMD-derived AUC, pnorm(|SMD|/sqrt(2)), was removed: a single
#   empirical comparison cannot validate it, and it is not reported.
# - Estimates from 5-6 studies are reported as preliminary.
#
# Inputs
# ------
# output/tables/clean_data.rds            (from Step 01; diagnostic data)
# output/tables/smd_pooled_estimates.csv  (from Step 02)
# R/meta_diagnostic.R
#
# Outputs
# -------
# output/tables/diagnostic_pooled.csv     (includes the add-Kawamura row)
# output/tables/diagnostic_per_study.csv
# output/tables/diagnostic_results.rds
# output/tables/diagnostic_sensitivity.csv      (revision: robustness)
# output/tables/diagnostic_model_details.csv    (revision: R1.5)
# output/tables/diagnostic_paired_comparison.csv (revision: R1.5)
#
# Run order
# ---------
# Requires Steps 01 and 02. Run before Steps 04 and 07.
# =============================================================================


source("R/00_setup.R")
source("R/load_data.R")
source("R/meta_diagnostic.R")

cd   <- readRDS(file.path(PATHS$tables, "clean_data.rds"))
diag <- cd$diagnostic

# ---- Per-study sensitivity / specificity (all rows) -------------------------
per_study <- diagnostic_per_study(diag)
write.csv(per_study, file.path(PATHS$tables, "diagnostic_per_study.csv"),
          row.names = FALSE)
message("Per-study sensitivity/specificity written to diagnostic_per_study.csv")

# ---- Bivariate models, iNPH vs AD only --------------------------------------
diag_ad <- dplyr::filter(diag, comparison == "iNPH_vs_AD")
reitsma_results <- list()
pooled_rows <- list()
for (bm in unique(diag_ad$biomarker)) {
  sub <- dplyr::filter(diag_ad, biomarker == bm)
  message("\nBiomarker: ", bm, " (iNPH vs AD), k = ", nrow(sub))
  r <- run_reitsma_one(sub)
  reitsma_results[[bm]] <- r
  if (!is.null(r$summary)) {
    pooled_rows[[bm]] <- r$summary
    print(r$summary, row.names = FALSE, digits = 3)
  }
}

# ---- Sensitivity analysis: add back Kawamura 2021 (p-tau181) ---------------
# Kawamura's p-tau181 result (sens = spec = 100%) is reserved for sensitivity
# analysis because eligibility depended on p-tau181. Read it from the
# "Not_primary" sheet and refit p-tau181 with it included.
cols_2x2 <- c("study", "biomarker", "comparison", "tp", "fn", "fp", "tn")
kaw_add  <- NULL               # p-tau181 rows with Kawamura added (reused below)
np <- tryCatch(
  suppressMessages(readxl::read_excel(
    file.path(PATHS$data, "diagnostic_performance.xlsx"),
    sheet = "Not_primary")),
  error = function(e) NULL)
if (!is.null(np) && "p-tau181" %in% names(reitsma_results)) {
  np  <- as.data.frame(np)
  kaw <- np[grepl("Kawamura", np[[1]]) & np[[2]] %in% "p-tau181" &
              grepl("^Reserved", np[[ncol(np)]]), 1:7, drop = FALSE]
  if (nrow(kaw) == 1) {
    names(kaw) <- cols_2x2
    kaw[c("tp", "fn", "fp", "tn")] <-
      lapply(kaw[c("tp", "fn", "fp", "tn")], as.numeric)
    ptau_ad <- dplyr::filter(diag_ad, biomarker == "p-tau181")
    add <- dplyr::bind_rows(as.data.frame(ptau_ad)[, cols_2x2], kaw)
    kaw_add <- add
    message("
Sensitivity: p-tau181 with Kawamura 2021 added back, k = ",
            nrow(add))
    sa <- run_reitsma_one(add)
    if (!is.null(sa$summary)) {
      sa$summary$biomarker <- "p-tau181 (add Kawamura)"
      print(sa$summary, row.names = FALSE, digits = 3)
      pooled_rows[["p-tau181_addKawamura"]] <- sa$summary
      reitsma_results[["p-tau181_addKawamura"]] <- sa
    }
  } else {
    message("\n(Kawamura p-tau181 sensitivity row not found on Not_primary)")
  }
}

# ---- Robustness of the pooled estimates (revision: R1.5, R2.1) -------------
# For each biomarker: leave-one-out; exclusion of count-source C tables
# (printed rates that cannot be reconciled with the stated group sizes);
# exclusion of studies in which CSF biomarkers contributed to the AD reference
# diagnosis (QUADAS-2 domain 3 high); the alternative denominators; and, for
# p-tau181, Kawamura 2021 added. Count source and QUADAS-2 domain 3 come from
# the "Notes" sheet, the alternative 2x2 counts from the "Alt_counts" sheet of
# diagnostic_performance.xlsx. Each scenario is pooled with the same rules as
# the primary analysis (pool_dta_scenario(): bivariate for k >= 4 with the
# univariate fallback at the boundary, univariate for k = 3, none for k < 3).
dx_path <- file.path(PATHS$data, "diagnostic_performance.xlsx")
read_dx <- function(sheet) tryCatch(
  as.data.frame(suppressMessages(readxl::read_excel(dx_path, sheet = sheet))),
  error = function(e) NULL)
dx_notes <- read_dx("Notes")
dx_alt   <- read_dx("Alt_counts")
note_of <- function(d, col) {
  if (is.null(dx_notes) || !col %in% names(dx_notes))
    return(rep(NA_character_, nrow(d)))
  dx_notes[[col]][match(paste(d$study, d$biomarker, d$comparison),
                        paste(dx_notes$Study, dx_notes$Biomarker,
                              dx_notes$Comparison))]
}

message("\nRobustness of the pooled diagnostic estimates ...")
sens_rows <- list(); detail_rows <- list()
for (bm in c("ab42", "p-tau181")) {
  d0 <- as.data.frame(dplyr::filter(diag_ad, biomarker == bm))
  if (nrow(d0) < 3 || is.null(reitsma_results[[bm]]$summary)) next
  add_s <- function(d, label) sens_rows[[length(sens_rows) + 1]] <<-
    pool_dta_scenario(d, label, bm)

  add_s(d0, "Primary analysis")
  for (i in seq_len(nrow(d0)))
    add_s(d0[-i, , drop = FALSE], paste("Omitting", d0$study[i]))

  src <- note_of(d0, "count_source")
  if (all(is.na(src))) {
    message("  (count_source not found on the Notes sheet - scenario skipped)")
  } else {
    add_s(d0[!src %in% "C", , drop = FALSE],
          "Excluding count source C tables")
  }
  d3 <- note_of(d0, "quadas2_D3")
  if (all(is.na(d3))) {
    message("  (quadas2_D3 not found on the Notes sheet - scenario skipped)")
  } else {
    add_s(d0[!d3 %in% "High", , drop = FALSE],
          "Excluding QUADAS-2 domain 3 high risk")
  }
  if (!is.null(dx_alt)) {
    a  <- dx_alt[dx_alt$Biomarker == bm & dx_alt$Comparison == "iNPH_vs_AD", ,
                 drop = FALSE]
    da <- d0
    for (j in seq_len(nrow(a))) {
      i <- match(a$Study[j], da$study)
      if (is.na(i)) {
        message("  ! Alt_counts row not matched: ", a$Study[j], " / ", bm)
        next
      }
      da[i, c("tp", "fn", "fp", "tn")] <-
        as.numeric(unlist(a[j, c("tp", "fn", "fp", "tn")]))
    }
    if (nrow(a)) {
      add_s(da, "Alternative denominators substituted")
      sens_rows[[length(sens_rows)]]$note <- paste(
        "alternative counts:", paste(intersect(a$Study, d0$study), collapse = "; "))
    }
  } else {
    message("  (Alt_counts sheet not found - scenario skipped)")
  }
  if (bm == "p-tau181" && !is.null(kaw_add))
    add_s(kaw_add, "Kawamura 2021 added")

  detail_rows[[length(detail_rows) + 1]] <-
    reitsma_model_details(reitsma_results[[bm]], d0, "Primary analysis")
}
if (!is.null(kaw_add) &&
    !is.null(reitsma_results[["p-tau181_addKawamura"]]))
  detail_rows[[length(detail_rows) + 1]] <-
    reitsma_model_details(reitsma_results[["p-tau181_addKawamura"]], kaw_add,
                          "Kawamura 2021 added")

dta_sens    <- dplyr::bind_rows(sens_rows)
dta_details <- dplyr::bind_rows(detail_rows)
write.csv(dta_sens, file.path(PATHS$tables, "diagnostic_sensitivity.csv"),
          row.names = FALSE)
write.csv(dta_details, file.path(PATHS$tables, "diagnostic_model_details.csv"),
          row.names = FALSE)
show <- dta_sens[, c("biomarker", "scenario", "k", "sensitivity",
                     "specificity", "model")]
show$sensitivity <- round(100 * show$sensitivity, 1)
show$specificity <- round(100 * show$specificity, 1)
print(show, row.names = FALSE)
message("Scenario results written to diagnostic_sensitivity.csv; model ",
        "details to diagnostic_model_details.csv")

# ---- Studies reporting both markers (Reviewer 1.5) -------------------------
# Descriptive only: too few studies to compare the two markers formally, and
# the evaluable participants are not always the same for both markers.
both_markers <- intersect(diag_ad$study[diag_ad$biomarker == "ab42"],
                          diag_ad$study[diag_ad$biomarker == "p-tau181"])
paired <- per_study[per_study$comparison == "iNPH_vs_AD" &
                      per_study$study %in% both_markers &
                      per_study$biomarker %in% c("ab42", "p-tau181"), ,
                    drop = FALSE]
paired <- paired[order(paired$study, paired$biomarker), , drop = FALSE]
same_n <- tapply(paste(paired$n_inph, paired$n_ctrl), paired$study,
                 function(x) length(unique(x)) == 1)
paired$same_group_sizes <- unname(same_n[paired$study])
write.csv(paired, file.path(PATHS$tables, "diagnostic_paired_comparison.csv"),
          row.names = FALSE)
message("Studies reporting both markers (", length(both_markers),
        "): diagnostic_paired_comparison.csv")

pooled_tbl <- dplyr::bind_rows(pooled_rows)
write.csv(pooled_tbl, file.path(PATHS$tables, "diagnostic_pooled.csv"),
          row.names = FALSE)
message("\nPooled diagnostic estimates written to diagnostic_pooled.csv")

saveRDS(list(per_study = per_study, reitsma = reitsma_results,
             pooled = pooled_tbl),
        file.path(PATHS$tables, "diagnostic_results.rds"))
message("\n03_diagnostic_metaanalysis.R complete.")
