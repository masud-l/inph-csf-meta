# =============================================================================
# Step 01 - Load, clean, and validate extracted study-level data
# Script: 01_validate_data.R
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
# Read the extracted study-level biomarker data from the source workbooks,
# standardise biomarker names, coerce numeric fields, and run validation
# checks before any pooling. Produces the single clean dataset used by all
# downstream scripts.
#
# Key behaviour
# -------------
# - Biomarker sheet names are normalised (e.g. "AB42 by AB40" -> "AB42/AB40").
# - Fixed column positions: En/Em/Esd (iNPH n/mean/SD), Cn/Cm/Csd
#   (comparator n/mean/SD), CSF_Type (column Q), iNPH_age (column R; AB42
#   sheets only). Missing optional columns are read as NA, never imputed.
# - Rows are checked for missing or non-positive n and SD.
#
# Inputs (not distributed; see data/README.md)
# ------
# data/CSF_iNPH_vs_HC.xlsx
# data/CSF_iNPH_vs_AD.xlsx
# R/00_setup.R, R/load_data.R
#
# Outputs
# -------
# output/tables/clean_data.rds
# output/tables/data_validation_report.csv
# output/tables/biomarker_study_counts.csv
#
# Run order
# ---------
# First script in the pipeline. Re-run whenever the source workbooks change.
# =============================================================================


source("R/00_setup.R")
source("R/load_data.R")

# ---- Load -------------------------------------------------------------------
message("Loading workbooks ...")
hc <- load_workbook(file.path(PATHS$data, "CSF_iNPH_vs_HC.xlsx"), "iNPH_vs_HC")
ad <- load_workbook(file.path(PATHS$data, "CSF_iNPH_vs_AD.xlsx"), "iNPH_vs_AD")
all_dat <- dplyr::bind_rows(hc, ad)

message(sprintf("  iNPH vs HC: %d rows across %d biomarkers",
                nrow(hc), dplyr::n_distinct(hc$biomarker)))
message(sprintf("  iNPH vs AD: %d rows across %d biomarkers",
                nrow(ad), dplyr::n_distinct(ad$biomarker)))

# ---- Validate ---------------------------------------------------------------
v <- validate_extraction(all_dat)
message(sprintf("\nValidation: %d rows in, %d clean, %d dropped.",
                v$n_in, v$n_clean, v$n_dropped))

if (nrow(v$report)) {
  message("Issues flagged (see report):")
  print(v$report, row.names = FALSE)
  write.csv(v$report,
            file.path(PATHS$tables, "data_validation_report.csv"),
            row.names = FALSE)
} else {
  message("No issues flagged.")
}

# ---- Per-biomarker study counts (what can / cannot be pooled) ---------------
counts <- v$clean |>
  dplyr::count(comparison, biomarker, name = "k") |>
  dplyr::mutate(poolable = k >= MIN_K) |>
  dplyr::arrange(comparison, dplyr::desc(k))
write.csv(counts, file.path(PATHS$tables, "biomarker_study_counts.csv"),
          row.names = FALSE)
message("\nBiomarker study counts written to ",
        file.path(PATHS$tables, "biomarker_study_counts.csv"))

# ---- Diagnostic data --------------------------------------------------------
diag <- load_diagnostic(file.path(PATHS$data, "diagnostic_performance.xlsx"))
diag_check <- diag |>
  dplyr::mutate(
    sens_gap = abs(sens_calc - sens_rep),
    spec_gap = abs(spec_calc - spec_rep)
  ) |>
  dplyr::select(study, biomarker, comparison, tp, fn, fp, tn,
                sens_calc, sens_rep, sens_gap,
                spec_calc, spec_rep, spec_gap, zero_cell)
write.csv(diag_check, file.path(PATHS$tables, "diagnostic_qc.csv"),
          row.names = FALSE)
message("Diagnostic QC (recomputed vs reported sens/spec, zero-cell flags) ",
        "written to ", file.path(PATHS$tables, "diagnostic_qc.csv"))

# ---- Persist the clean objects for downstream scripts -----------------------
saveRDS(list(smd = v$clean, diagnostic = diag),
        file.path(PATHS$tables, "clean_data.rds"))
message("\n01_validate_data.R complete.")
