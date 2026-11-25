# =============================================================================
# Step 07 - Build manuscript and supplementary tables
# Script: 07_tables.R
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
# Generate Word tables (via officer) and their underlying CSV data:
#   Table 2               - pooled SMD, principal biomarkers (AB42, NfL,
#                           p-tau181, sAPPa)
#   Supplementary Table S3 - pooled SMD, all biomarkers with k >= 3
#   Supplementary Table S4 - sample-size-weighted mean lumbar CSF AB42 and
#                           p-tau181 concentrations (pg/mL or ng/L rows only)
#   Supplementary Table S5 - empirical diagnostic accuracy vs Alzheimer disease
#                           (pooled estimates, Kawamura sensitivity row and
#                           study-level 2x2 counts)
#   Supplementary Table S7 - robustness of the diagnostic accuracy estimates
#                           (sensitivity analyses, model details, studies
#                           reporting both markers; from Step 03)
#   Numbering follows the manuscript. Supplementary Tables S1-S2 (JBI quality
#   appraisal) are prepared by hand.
#
# Key behaviour
# -------------
# - Biomarker names normalised (diagnostic vs SMD spelling).
# - (Revision) The former Supplementary Tables S6-S7 (SMD-derived AUC and its
#   comparison with the empirical AUC) were removed; the SMD-derived AUC is
#   not reported.
# - Supplementary Table S4 filters CSF_Type to lumbar (trimmed, case-insensitive) and reports
#   per-group concentration ranges, sample-size-weighted means, and the
#   between-group gap.
# - Final journal formatting was applied after export.
#
# Inputs
# ------
# output/tables/smd_pooled_estimates.csv     (from Step 02)
# output/tables/diagnostic_pooled.csv        (from Step 03)
# output/tables/clean_data.rds               (from Step 01)
# data/diagnostic_performance.xlsx           (2x2 counts and count source)
#
# Outputs
# -------
# output/tables/Table2.docx, Supplementary_Table_S3.docx,
#   Supplementary_Table_S5.docx, Supplementary_Table_S7_diagnostic_robustness.docx
#   (Supplementary_Table_S4.docx when unheld)
# output/tables/<table>_data.csv
#
# Run order
# ---------
# Requires Steps 01-03.
# =============================================================================

source("R/00_setup.R")
source("R/theme_bist.R")          # biomarker_label(), comparison_style()

if (!requireNamespace("officer", quietly = TRUE))
  install.packages("officer", repos = "https://cloud.r-project.org")
library(officer)

# ---- The four main biomarkers (raw CSV names) -------------------------------
MAIN_BIOMARKERS <- c("AB42", "NFL", "pTau181", "sAPPalpha")

# ---- Load pipeline tables ---------------------------------------------------
smd   <- read.csv(file.path(PATHS$tables, "smd_pooled_estimates.csv"),
                  stringsAsFactors = FALSE)
emp   <- tryCatch(read.csv(file.path(PATHS$tables, "diagnostic_pooled.csv"),
                           stringsAsFactors = FALSE),
                  error = function(e) data.frame())
# Per-study extraction data (group means/SDs/n + CSF_Type) - for Supplementary Table S4.
clean_smd <- readRDS(file.path(PATHS$tables, "clean_data.rds"))$smd

# ---- Name normalization -----------------------------------------------------
# diagnostic_pooled.csv labels biomarkers differently from the SMD pipeline
# ("ab42" vs "AB42", "p-tau181" vs "pTau181"). Map to the canonical SMD name so
# the empirical-AUC join works. The leave-one-out row "p-tau181 (excl.
# Kawamura)" / "p-tau181 (add Kawamura)" sensitivity rows are intentionally NOT
# mapped - they must not match a main biomarker.
canonical_name <- function(x) {
  m <- c("ab42" = "AB42", "p-tau181" = "pTau181", "ptau181" = "pTau181",
         "abeta42" = "AB42")
  out <- m[tolower(trimws(x))]
  ifelse(is.na(out), x, unname(out))
}
if (nrow(emp)) emp$biomarker_canon <- canonical_name(emp$biomarker)
# Prediction-interval columns come from Step 02 (revision); NA if 02 is older.
if (!"pi_lb" %in% names(smd)) { smd$pi_lb <- NA_real_; smd$pi_ub <- NA_real_ }

# ---- Formatting helpers -----------------------------------------------------
fmt_ci <- function(e, l, h) sprintf("%.2f (%.2f to %.2f)", e, l, h)
fmt_p  <- function(p) ifelse(is.na(p), "\u2014",
                      ifelse(p < 0.001, "<.001",
                             sub("^0", "", sprintf("%.3f", p))))
fmt_i2 <- function(i2) ifelse(is.na(i2), "\u2014", sprintf("%.1f", i2))
fmt_t2 <- function(t2) ifelse(is.na(t2), "\u2014", sprintf("%.3f", t2))
fmt_pc <- function(x)  sprintf("%.1f", 100 * x)          # proportion -> %, 1 dp
fmt_pi <- function(lb, ub, k)                            # 95% prediction interval
  ifelse(is.na(lb) | k < PI_MIN_K, "\u2014", sprintf("%.2f to %.2f", lb, ub))

# Empirical AUC lookup for one biomarker x comparison; NA if none.
empirical_auc <- function(bm, cmp) {
  if (!nrow(emp)) return(NA_real_)
  e <- emp[emp$biomarker_canon == bm & emp$comparison == cmp, , drop = FALSE]
  if (nrow(e) == 1) e$AUC[1] else NA_real_
}

# =============================================================================
# SMD TABLE BUILDER  (Table 2 and Supplementary Table S3)
# =============================================================================
# SMD only (the derived-AUC column was removed in the revision).
build_smd_table <- function(biomarkers = NULL, with_auc = FALSE) {
  d <- smd[smd$k >= MIN_K, , drop = FALSE]                 # only meta-analysed
  if (!is.null(biomarkers))
    d <- d[d$biomarker %in% biomarkers, , drop = FALSE]

  rows <- lapply(seq_len(nrow(d)), function(i) {
    r <- d[i, ]
    base <- data.frame(
      comparison = r$comparison,
      Biomarker  = biomarker_label(r$biomarker),
      Studies    = r$k,
      SMD        = fmt_ci(r$estimate, r$ci_lb, r$ci_ub),
      PI         = fmt_pi(r$pi_lb, r$pi_ub, r$k),
      I2         = fmt_i2(r$I2),
      tau2       = fmt_t2(r$tau2),
      P          = fmt_p(r$pval),
      est_sort   = abs(r$estimate),
      stringsAsFactors = FALSE
    )
    base
  })
  do.call(rbind, rows)
}

# Render an SMD table to .docx (blocks by comparison, largest effect first).
render_smd_docx <- function(tbl, title, file, with_auc) {
  cols <- c("Biomarker", "No. of Studies", "SMD (95% CI)", "95% PI",
            "I\u00B2, %", "\u03C4\u00B2", "P Value")

  doc <- read_docx()
  doc <- body_add_par(doc, title, style = "heading 2")

  for (cmp in c("iNPH_vs_HC", "iNPH_vs_AD")) {
    block <- tbl[tbl$comparison == cmp, , drop = FALSE]
    if (!nrow(block)) next
    block <- block[order(-block$est_sort), ]
    doc <- body_add_par(doc, comparison_style(cmp)$label, style = "Normal")

    mat <- data.frame(block$Biomarker, block$Studies, block$SMD, block$PI,
                       block$I2, block$tau2, block$P,
                       stringsAsFactors = FALSE)
    if (with_auc) mat <- cbind(mat, block$AUC)
    names(mat) <- cols
    doc <- body_add_table(doc, mat, style = "table_template")
    doc <- body_add_par(doc, "", style = "Normal")
  }

  fn <- c(
    paste0("Abbreviations: A\u03B2, amyloid-\u03B2; ",
           if (with_auc) "AUC, area under the receiver operating characteristic curve; " else "",
           "CI, confidence interval; CSF, cerebrospinal fluid; iNPH, idiopathic normal pressure hydrocephalus; NfL, neurofilament light chain; PI, prediction interval; sAPP\u03B1, soluble amyloid precursor protein \u03B1; SMD, standardized mean difference."),
    "SMDs are random-effects pooled estimates (Hedges g). A negative SMD indicates a lower CSF concentration in iNPH than in the comparator. I\u00B2 and \u03C4\u00B2 quantify between-study heterogeneity. Only biomarkers meta-analysed in at least 3 studies are shown.",
    paste0("The 95% prediction interval (PI) is the range expected to contain the effect in a new comparable study (t distribution with k \u2212 2 degrees of freedom). It is shown for analyses of at least ", PI_MIN_K, " studies; with fewer studies it is too wide to be informative.")
  )
  for (line in fn) doc <- body_add_par(doc, line, style = "Normal")

  print(doc, target = file)
  message("  wrote ", file)
}

# =============================================================================
# SUPPLEMENTARY TABLE S5  —  empirical diagnostic accuracy (iNPH vs AD)
# =============================================================================
# Layout and wording follow the lead author's Table S5 (Nov 2026):
#   Part 1 - pooled sensitivity / specificity / summary AUC per biomarker,
#            plus the labelled Kawamura 2021 p-tau181 sensitivity row.
#   Part 2 - study-level 2x2 counts used in the main analyses, with their
#            count source (B/C) from the "Notes" sheet.
# Pooled values come from diagnostic_pooled.csv (Step 03); counts are read
# from data/diagnostic_performance.xlsx, the same file Step 03 pools.
DIAG_XLSX     <- file.path(PATHS$data, "diagnostic_performance.xlsx")
S5_BIOMARKERS <- c("ab42", "p-tau181")          # display order: Abeta42 first
S5_KAWAMURA   <- "p-tau181 (add Kawamura)"      # row label written by Step 03
S5_MARK       <- "\u2020"                       # dagger for the sensitivity row

fmt_ci_pc   <- function(x, lo, hi)
  paste0(fmt_pc(x), " (", fmt_pc(lo), "\u2013", fmt_pc(hi), ")")
short_study <- function(x) gsub(" et al\\.,? ", " ", x)  # "Kapaki et al., 2007" -> "Kapaki 2007"
read_diag   <- function(sheet)
  as.data.frame(suppressMessages(readxl::read_excel(DIAG_XLSX, sheet = sheet)))

build_table_s5 <- function() {
  # --- Part 1: pooled estimates (main rows, then the Kawamura row) ----------
  e    <- emp[emp$comparison == "iNPH_vs_AD", , drop = FALSE]
  rows <- e[e$biomarker %in% c(S5_BIOMARKERS, S5_KAWAMURA), , drop = FALSE]
  rows <- rows[order(match(rows$biomarker, c(S5_BIOMARKERS, S5_KAWAMURA))), ,
               drop = FALSE]
  summary_tbl <- data.frame(
    Biomarker   = ifelse(rows$biomarker == S5_KAWAMURA,
                         paste0("p-tau181 with Kawamura 2021 (sensitivity analysis)",
                                S5_MARK),
                         biomarker_label(rows$biomarker_canon)),
    Studies     = rows$k,
    Sensitivity = fmt_ci_pc(rows$sensitivity, rows$sens_lo, rows$sens_hi),
    Specificity = fmt_ci_pc(rows$specificity, rows$spec_lo, rows$spec_hi),
    AUC         = ifelse(is.na(rows$AUC), "Not estimable",
                         sprintf("%.3f", rows$AUC)),
    stringsAsFactors = FALSE
  )

  # --- Part 2: study-level counts used in the main analyses -----------------
  s1 <- read_diag(1)
  names(s1)[1:7] <- c("study", "biomarker", "comparison", "tp", "fn", "fp", "tn")
  s1 <- s1[s1$comparison == "iNPH_vs_AD" & s1$biomarker %in% S5_BIOMARKERS, ,
           drop = FALSE]
  nt <- tryCatch(read_diag("Notes"), error = function(err) NULL)
  src <- if (is.null(nt)) rep(NA_character_, nrow(s1)) else
    nt$count_source[match(paste(s1$study, s1$biomarker, s1$comparison),
                          paste(nt$Study, nt$Biomarker, nt$Comparison))]
  ord <- order(match(s1$biomarker, S5_BIOMARKERS), s1$study)
  s1  <- s1[ord, , drop = FALSE]
  src <- src[ord]
  counts_tbl <- data.frame(
    Study     = short_study(s1$study),
    Biomarker = biomarker_label(canonical_name(s1$biomarker)),
    TP = as.integer(s1$tp), FP = as.integer(s1$fp),
    FN = as.integer(s1$fn), TN = as.integer(s1$tn),
    Source    = ifelse(is.na(src), "", src),
    stringsAsFactors = FALSE
  )

  # Consistency check: the counts must be the studies that were pooled.
  for (bm in S5_BIOMARKERS) {
    k_pool  <- e$k[e$biomarker == bm]
    k_count <- sum(s1$biomarker == bm)
    if (length(k_pool) == 1 && k_pool != k_count)
      warning("Supplementary Table S5: ", bm, " pooled k = ", k_pool,
              " but ", k_count, " count rows - re-run Step 03.", call. = FALSE)
  }

  # Kawamura p-tau181 counts for the sensitivity-row footnote (Not_primary).
  np  <- tryCatch(read_diag("Not_primary"), error = function(err) NULL)
  kaw <- if (is.null(np)) NULL else
    np[grepl("Kawamura", np[[1]]) & np[[2]] %in% "p-tau181" &
         grepl("^Reserved", np[[ncol(np)]]), , drop = FALSE]
  kaw_counts <- if (!is.null(kaw) && nrow(kaw) == 1)
    sprintf("TP %d, FP %d, FN %d, TN %d",
            as.integer(kaw$tp), as.integer(kaw$fp),
            as.integer(kaw$fn), as.integer(kaw$tn)) else NA_character_

  list(summary = summary_tbl, counts = counts_tbl,
       kawamura_counts = kaw_counts,
       has_kawamura = any(rows$biomarker == S5_KAWAMURA))
}

render_table_s5_docx <- function(parts, file) {
  doc <- read_docx()
  doc <- body_add_par(doc,
    "Supplementary Table S5. Diagnostic accuracy of CSF A\u03B242 and p-tau181 for differentiating iNPH from Alzheimer disease",
    style = "heading 2")
  doc <- body_add_par(doc,
    "iNPH was the positive diagnosis. Sensitivity is the proportion of iNPH cases classified as iNPH. Specificity is the proportion of AD cases classified as AD.",
    style = "Normal")

  m1 <- parts$summary
  names(m1) <- c("Biomarker", "Studies", "Sensitivity, % (95% CI)",
                 "Specificity, % (95% CI)", "Summary AUC")
  doc <- body_add_table(doc, m1, style = "table_template")
  doc <- body_add_par(doc, "", style = "Normal")

  doc <- body_add_par(doc, "Study-level counts used in the main analyses",
                      style = "Normal")
  doc <- body_add_table(doc, parts$counts, style = "table_template")
  doc <- body_add_par(doc, "", style = "Normal")

  fn <- c(
    "B denotes reconstructed counts consistent with the reported group sizes. C denotes a reconstruction for which the printed rates and stated group sizes could not be reconciled under conventional rounding. These counts can depend on the assumed ROC denominator or truncation of a printed percentage.",
    "The A\u03B242 summary uses a bivariate random-effects model. For p-tau181, the between-study correlation reached the model boundary. Sensitivity and specificity were pooled in separate random-effects logit models. For pooling, 0.5 was added to each cell of a study table containing a zero. The counts above precede this correction.",
    "Schirinzi\u2019s reconstructed denominators exceed enrollment. Kapaki\u2019s A\u03B242 threshold runs in the opposite direction to those of the other studies. Both were retained in the main analyses and examined by omission. Only the Kawamura derivation cohort contributes to the main A\u03B242 analysis. Its p-tau181 result was reserved for sensitivity analysis because eligibility depended on p-tau181."
  )
  if (isTRUE(parts$has_kawamura))
    fn <- c(fn, paste0(
      S5_MARK, " Sensitivity analysis: the main p-tau181 studies plus the Kawamura 2021 derivation cohort",
      if (!is.na(parts$kawamura_counts)) paste0(" (", parts$kawamura_counts, ")") else "",
      ", pooled in separate random-effects logit models. Not part of the main analysis because p-tau181 contributed to eligibility in that cohort."))
  fn <- c(fn,
    "AD, Alzheimer disease. AUC, area under the summary receiver operating characteristic curve. CI, confidence interval. CSF, cerebrospinal fluid. iNPH, idiopathic normal pressure hydrocephalus. TP, FP, FN and TN denote true positives, false positives, false negatives and true negatives.")
  for (line in fn) doc <- body_add_par(doc, line, style = "Normal")
  print(doc, target = file)
  message("  wrote ", file)
}

# =============================================================================
# SUPPLEMENTARY TABLE S7  —  robustness of the diagnostic accuracy estimates
# =============================================================================
# Revision (R1.5, R2.1). Built from the Step 03 outputs
# diagnostic_sensitivity.csv, diagnostic_model_details.csv and
# diagnostic_paired_comparison.csv. Proposed number: S7 became free when the
# empirical-vs-derived AUC table was dropped; change it here if needed.
DTA_ROBUST_LABEL <- "Supplementary Table S7"

short_model <- function(m) {
  ifelse(is.na(m), "",
  ifelse(grepl("^bivariate", m), "Bivariate",
  ifelse(grepl("boundary", m), "Separate logit (boundary)",
  ifelse(grepl("k = 3", m), "Separate logit (k = 3)",
  ifelse(grepl("^not pooled", m), "Not pooled", m)))))
}
fmt_num <- function(x, d = 2) ifelse(is.na(x), "\u2014", formatC(x, format = "f", digits = d))

render_dta_robustness_docx <- function(file) {
  rd <- function(f) {
    p <- file.path(PATHS$tables, f)
    if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE,
                                 fileEncoding = "UTF-8") else NULL
  }
  sens <- rd("diagnostic_sensitivity.csv")
  det  <- rd("diagnostic_model_details.csv")
  pair <- rd("diagnostic_paired_comparison.csv")
  if (is.null(sens)) {
    message("  ", DTA_ROBUST_LABEL, " skipped: diagnostic_sensitivity.csv ",
            "not found (re-run Step 03)")
    return(invisible(NULL))
  }
  lab_bm <- function(b) biomarker_label(canonical_name(b))
  pc_ci  <- function(x, lo, hi) ifelse(is.na(x), "\u2014", fmt_ci_pc(x, lo, hi))

  doc <- read_docx()
  doc <- body_add_par(doc, paste0(DTA_ROBUST_LABEL,
    ". Robustness of the pooled diagnostic accuracy of CSF A\u03B242 and p-tau181 (iNPH vs Alzheimer disease)"),
    style = "heading 2")

  # Part A - sensitivity analyses
  doc <- body_add_par(doc, "A. Sensitivity analyses", style = "Normal")
  sens <- sens[order(match(sens$biomarker, S5_BIOMARKERS)), , drop = FALSE]
  first <- !duplicated(sens$biomarker)
  ma <- data.frame(
    Biomarker = ifelse(first, lab_bm(sens$biomarker), ""),
    Analysis  = short_study(sens$scenario),
    Studies   = sens$k,
    Sens      = pc_ci(sens$sensitivity, sens$sens_lo, sens$sens_hi),
    Spec      = pc_ci(sens$specificity, sens$spec_lo, sens$spec_hi),
    Model     = short_model(sens$model),
    stringsAsFactors = FALSE)
  names(ma) <- c("Biomarker", "Analysis", "Studies", "Sensitivity, % (95% CI)",
                 "Specificity, % (95% CI)", "Model")
  doc <- body_add_table(doc, ma, style = "table_template")
  doc <- body_add_par(doc, "", style = "Normal")

  # Part B - model estimation details
  if (!is.null(det)) {
    doc <- body_add_par(doc, "B. Model estimation details", style = "Normal")
    mb <- data.frame(
      Biomarker = lab_bm(det$biomarker),
      Analysis  = det$analysis,
      Studies   = det$k,
      Model     = short_model(det$reported_model),
      Conv      = ifelse(is.na(det$bivariate_converged), "Not reported",
                         ifelse(det$bivariate_converged, "Yes", "No")),
      Rho       = fmt_num(det$bivariate_rho, 3),
      Tau2Biv   = paste0(fmt_num(det$bivariate_tau2_logit_sens), " / ",
                         fmt_num(det$bivariate_tau2_logit_fpr)),
      Tau2Uni   = paste0(fmt_num(det$univariate_tau2_sens), " (",
                         fmt_num(det$univariate_I2_sens, 1), ") / ",
                         fmt_num(det$univariate_tau2_spec), " (",
                         fmt_num(det$univariate_I2_spec, 1), ")"),
      CC        = short_study(det$continuity_correction),
      stringsAsFactors = FALSE)
    names(mb) <- c("Biomarker", "Analysis", "Studies", "Reported model",
                   "Bivariate fit converged", "\u03C1",
                   "Bivariate \u03C4\u00B2: logit sensitivity / logit FPR",
                   "Separate models \u03C4\u00B2 (I\u00B2, %): sensitivity / specificity",
                   "Continuity correction")
    doc <- body_add_table(doc, mb, style = "table_template")
    doc <- body_add_par(doc, "", style = "Normal")
  }

  # Part C - studies reporting both markers
  if (!is.null(pair) && nrow(pair)) {
    doc <- body_add_par(doc, "C. Studies reporting both markers", style = "Normal")
    first <- !duplicated(pair$study)
    mc <- data.frame(
      Study     = ifelse(first, short_study(pair$study), ""),
      Biomarker = lab_bm(pair$biomarker),
      nI        = pair$n_inph,
      nA        = pair$n_ctrl,
      Sens      = pc_ci(pair$sens, pair$sens_lo, pair$sens_hi),
      Spec      = pc_ci(pair$spec, pair$spec_lo, pair$spec_hi),
      stringsAsFactors = FALSE)
    names(mc) <- c("Study", "Biomarker", "iNPH, n", "AD, n",
                   "Sensitivity, % (95% CI)", "Specificity, % (95% CI)")
    doc <- body_add_table(doc, mc, style = "table_template")
    doc <- body_add_par(doc, "", style = "Normal")
  }

  alt_note <- sens$note[!is.na(sens$note) & nzchar(sens$note)]
  alt_studies <- unique(unlist(strsplit(sub("^alternative counts: ", "", alt_note), "; ")))
  diff_n <- if (!is.null(pair)) unique(pair$study[!pair$same_group_sizes]) else character(0)
  fn <- c(
    "iNPH was the positive diagnosis. Every analysis uses the rules of the main analysis: a bivariate random-effects (Reitsma) model when at least 4 studies contribute and the between-study correlation is not at its boundary; otherwise separate random-effects logit models for sensitivity and specificity; 0.5 added to each cell of a study table containing a zero. No pooled estimate is given for fewer than 3 studies.",
    "Count source C: the printed sensitivity or specificity cannot be reconciled with the stated group sizes (Supplementary Table S5).",
    "QUADAS-2 domain 3 high risk: the reference diagnosis was not independent of CSF biomarkers (Abu-Rumeileh 2019, AD classification used CSF profiles including A\u03B242; Kawamura 2021, p-tau181 determined iNPH eligibility).",
    if (length(alt_studies)) paste0("Alternative denominators: the count-source C tables were reconstructed with the alternative denominators documented in the extraction audit (", paste(short_study(alt_studies), collapse = ", "), ").") else NULL,
    "Kawamura 2021 added: the derivation-cohort p-tau181 result, excluded from the main analysis because p-tau181 contributed to eligibility.",
    "\u03C1, between-study correlation of logit sensitivity and logit false-positive rate (FPR); values at \u00B11 indicate the boundary at which the summary ROC curve is not identifiable. \u03C4\u00B2, between-study variance on the logit scale.",
    paste0("Studies reporting both markers are compared descriptively only; too few studies allowed a formal comparison",
           if (length(diff_n)) paste0(", and ", paste(short_study(diff_n), collapse = ", "), " reported different evaluable group sizes for the two markers") else "",
           ". Study confidence intervals are exact (Clopper-Pearson)."),
    "AD, Alzheimer disease; CI, confidence interval; CSF, cerebrospinal fluid; FPR, false-positive rate; iNPH, idiopathic normal pressure hydrocephalus; ROC, receiver operating characteristic."
  )
  for (line in fn) doc <- body_add_par(doc, line, style = "Normal")
  print(doc, target = file)
  message("  wrote ", file)
}

# =============================================================================
# SUPPLEMENTARY TABLE S4  —  reported CSF concentrations (lumbar only)
# =============================================================================
# Descriptive summary of the group-mean CSF concentrations reported by the
# studies of the primary meta-analysis - NOT a derived cutoff. For each
# biomarker x comparison and each group (iNPH, comparator): number of studies,
# total n, range of the study means and the sample-size-weighted mean
#     weighted mean = sum(n_i * mean_i) / sum(n_i)
# plus the between-group gap (iNPH minus comparator weighted mean).
# A row is combined only if its CSF is lumbar (CSF_Type, column Q) and its
# unit (column O) is pg/mL or ng/L (numerically identical). Rows failing either
# check are listed in the footnote and in Supplementary_Table_S4_excluded.csv,
# so the study counts always match what was combined (Reviewer 1 minor 4;
# lead author's check of the study counts, Nov 2026).
S4_UNITS <- c("pg/ml", "ng/l")

# Sample-size-weighted mean (weight = group n).
wmean <- function(x, w) sum(x * w) / sum(w)

build_table_s4 <- function(biomarkers = c("AB42", "pTau181")) {
  if (!"unit" %in% names(clean_smd))
    stop("clean_data.rds has no 'unit' column - re-run Step 01 first.")
  d <- clean_smd[clean_smd$biomarker %in% biomarkers, , drop = FALSE]

  site <- tolower(trimws(d$CSF_Type))
  unit <- tolower(gsub("\\s", "", d$unit))
  is_lumbar <- !is.na(site) & site == "lumbar"
  unit_ok   <- !is.na(unit) & unit %in% S4_UNITS
  keep      <- is_lumbar & unit_ok
  excluded <- data.frame(
    study = d$study, biomarker = d$biomarker, comparison = d$comparison,
    reason = ifelse(!is_lumbar,
      paste0("CSF site ", ifelse(is.na(d$CSF_Type), "not recorded", d$CSF_Type)),
      ifelse(is.na(d$unit) | trimws(d$unit) %in% c("_", "-", ""), "unit not recorded",
             paste0("unit ", d$unit, " (not pg/mL or ng/L)"))),
    stringsAsFactors = FALSE)[!keep, , drop = FALSE]
  d <- d[keep, , drop = FALSE]
  message("  Supplementary Table S4: kept ", nrow(d), " of ",
          nrow(d) + nrow(excluded), " rows (lumbar CSF, pg/mL or ng/L).")
  if (nrow(excluded)) print(excluded, row.names = FALSE)

  # How many of the combined study means were estimated from medians.
  conv_path <- file.path(PATHS$data, "conversion_status.csv")
  converted <- NA_integer_
  if (file.exists(conv_path)) {
    conv <- read.csv(conv_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    dm <- conv$direct_mean_sd[match(paste(d$comparison, d$study, d$sheet),
                                    paste(conv$comparison, conv$study, conv$sheet))]
    converted <- sum(dm %in% 0)
  }

  rows <- list()
  for (bm in biomarkers) {
    for (cmp in c("iNPH_vs_HC", "iNPH_vs_AD")) {
      s <- d[d$biomarker == bm & d$comparison == cmp, , drop = FALSE]
      if (!nrow(s)) next
      k_primary <- sum(clean_smd$biomarker == bm & clean_smd$comparison == cmp)
      if (nrow(s) != k_primary)
        message("  Table S4 ", bm, " ", cmp, ": ", nrow(s), " of ", k_primary,
                " primary studies combined (see footnote).")
      inph_w <- wmean(s$Em, s$En)            # iNPH group weighted mean
      comp_w <- wmean(s$Cm, s$Cn)            # comparator group weighted mean
      rows[[length(rows) + 1]] <- data.frame(
        biomarker  = biomarker_label(bm),
        comparison = comparison_style(cmp)$label,
        k          = nrow(s),
        k_primary  = k_primary,
        inph_n     = sum(s$En),
        inph_range = sprintf("%.0f\u2013%.0f", min(s$Em), max(s$Em)),
        inph_wmean = sprintf("%.0f", inph_w),
        comp_n     = sum(s$Cn),
        comp_range = sprintf("%.0f\u2013%.0f", min(s$Cm), max(s$Cm)),
        comp_wmean = sprintf("%.0f", comp_w),
        gap        = sprintf("%.0f", inph_w - comp_w),
        stringsAsFactors = FALSE
      )
    }
  }
  list(table = do.call(rbind, rows), excluded = excluded,
       converted = converted, n_rows = nrow(d))
}

render_table_s4_docx <- function(parts, file) {
  tbl <- parts$table
  doc <- read_docx()
  doc <- body_add_par(doc,
    "Supplementary Table S4. Sample-size-weighted mean lumbar CSF A\u03B242 and p-tau181 concentrations in iNPH compared with healthy controls and Alzheimer disease",
    style = "heading 2")

  # One block per biomarker.
  for (bm in unique(tbl$biomarker)) {
    block <- tbl[tbl$biomarker == bm, , drop = FALSE]
    doc <- body_add_par(doc, paste0("CSF ", bm), style = "Normal")
    m <- data.frame(
      block$comparison, block$k,
      block$inph_n, block$inph_range, block$inph_wmean,
      block$comp_n, block$comp_range, block$comp_wmean, block$gap,
      stringsAsFactors = FALSE, check.names = FALSE)
    names(m) <- c("Comparison", "No. of Studies",
                  "iNPH, n", "iNPH Range, pg/mL", "iNPH Weighted Mean, pg/mL",
                  "Comparator, n", "Comparator Range, pg/mL",
                  "Comparator Weighted Mean, pg/mL", "Between-Group Gap, pg/mL")
    doc <- body_add_table(doc, m, style = "table_template")
    doc <- body_add_par(doc, "", style = "Normal")
  }

  ex <- parts$excluded
  ex_txt <- if (nrow(ex)) {
    paste0("Concentrations reported in ng/L and pg/mL were combined (1 ng/L = 1 pg/mL). Not combined: ",
           paste0(short_study(ex$study), " (", biomarker_label(ex$biomarker),
                  ifelse(ex$comparison == "iNPH_vs_HC", " vs HC", " vs AD"),
                  "; ", ex$reason, ")", collapse = "; "),
           ". The number of studies is therefore lower than in Table 2 for ",
           if (length(unique(paste(ex$biomarker, ex$comparison))) == 1) "this comparison." else "these comparisons.")
  } else {
    "All records of the primary analysis used lumbar CSF and reported concentrations in pg/mL or ng/L, which were combined (1 ng/L = 1 pg/mL); the numbers of studies equal those in Table 2."
  }
  conv_txt <- if (!is.na(parts$converted) && parts$converted > 0)
    paste0(parts$converted, " of the ", parts$n_rows,
           " study means were estimated from reported medians and quartiles (Luo method); the others were reported as means.")
  else NULL

  fn <- c(
    "Abbreviations: A\u03B2, amyloid-\u03B2; AD, Alzheimer disease; CSF, cerebrospinal fluid; HC, healthy controls; iNPH, idiopathic normal pressure hydrocephalus.",
    "Values summarise the group-mean CSF concentrations reported by the studies of the primary meta-analysis; they are not derived diagnostic thresholds. Only records with lumbar CSF are included.",
    "For each group, the weighted mean is \u03A3(n \u00D7 m)/\u03A3n over the contributing studies, where m is the mean concentration reported for that group in a study and n the group's sample size; the n column gives \u03A3n. The range is the lowest to highest study mean; the between-group gap is the iNPH weighted mean minus the comparator weighted mean. No variance weighting, confidence intervals or statistical tests were applied.",
    ex_txt,
    conv_txt,
    "The studies used different assays and calibrations, so absolute concentrations are not directly comparable across studies; the values are descriptive and exploratory and should not be compared with assay-specific cut-offs."
  )
  for (line in fn) doc <- body_add_par(doc, line, style = "Normal")
  print(doc, target = file)
  message("  wrote ", file)
}


# =============================================================================
# BUILD ALL TABLES
# =============================================================================
message("Building manuscript tables ...")

# Table 2 - main biomarkers, SMD only
t2 <- build_smd_table(MAIN_BIOMARKERS, with_auc = FALSE)
render_smd_docx(t2,
  "Table 2. Pooled Standardized Mean Differences of Key CSF Biomarkers for iNPH Versus Healthy Controls and Alzheimer Disease",
  file.path(PATHS$tables, "Table2.docx"), with_auc = FALSE)

# Supplementary Table S3 - all meta-analysed biomarkers, SMD only
ts3 <- build_smd_table(NULL, with_auc = FALSE)
render_smd_docx(ts3,
  "Supplementary Table S3. Pooled standardized mean differences of all CSF biomarkers",
  file.path(PATHS$tables, "Supplementary_Table_S3.docx"), with_auc = FALSE)

# Supplementary Table S5 - empirical diagnostic accuracy
ts5 <- build_table_s5()
render_table_s5_docx(ts5, file.path(PATHS$tables, "Supplementary_Table_S5.docx"))

# Supplementary Table S7 - robustness of the diagnostic accuracy estimates
render_dta_robustness_docx(file.path(PATHS$tables,
                                     "Supplementary_Table_S7_diagnostic_robustness.docx"))

# Supplementary Table S4 - reported CSF concentrations, lumbar only.
# Reinstated in the revision at the lead author's request. Rows are combined
# only when the CSF is lumbar and the unit is pg/mL or ng/L; any other row is
# listed in the footnote. Set RENDER_TABLE_S4 <- FALSE to hold it again.
RENDER_TABLE_S4 <- TRUE
if (RENDER_TABLE_S4) {
  ts4 <- build_table_s4(c("AB42", "pTau181"))
  render_table_s4_docx(ts4, file.path(PATHS$tables, "Supplementary_Table_S4.docx"))
} else {
  message("  Supplementary Table S4 skipped (RENDER_TABLE_S4 = FALSE)")
}


# ---- Underlying data as CSV -------------------------------------------------
write.csv(t2[, setdiff(names(t2), "est_sort")],
          file.path(PATHS$tables, "Table2_data.csv"), row.names = FALSE)
write.csv(ts3[, setdiff(names(ts3), "est_sort")],
          file.path(PATHS$tables, "Supplementary_Table_S3_data.csv"), row.names = FALSE)
write.csv(ts5$summary, file.path(PATHS$tables, "Supplementary_Table_S5_data.csv"),
          row.names = FALSE)
write.csv(ts5$counts, file.path(PATHS$tables, "Supplementary_Table_S5_counts.csv"),
          row.names = FALSE)
if (RENDER_TABLE_S4) {
  write.csv(ts4$table, file.path(PATHS$tables, "Supplementary_Table_S4_data.csv"),
            row.names = FALSE)
  write.csv(ts4$excluded, file.path(PATHS$tables, "Supplementary_Table_S4_excluded.csv"),
            row.names = FALSE)
}

message("\n07_tables.R complete:")
message("  Table2.docx                  - principal biomarkers (SMD only)")
message("  Supplementary_Table_S3.docx  - all meta-analysed biomarkers (SMD)")
message("  Supplementary_Table_S4.docx  - weighted mean lumbar CSF concentrations (pg/mL)")
message("  Supplementary_Table_S5.docx  - empirical diagnostic accuracy (p-tau181, Abeta42 vs AD)")
message("  Supplementary_Table_S7_diagnostic_robustness.docx - sensitivity analyses, model details, paired markers")
