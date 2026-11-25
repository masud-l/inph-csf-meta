# =============================================================================
# Step 12 - Sensitivity analyses of the pooled SMDs (peer review)
# Script: 12_sensitivity_analyses.R
# =============================================================================
#
# Purpose
# -------
# Re-pool every biomarker x comparison that is meta-analysed in the primary
# analysis (k >= MIN_K) under the sensitivity rules agreed during the
# revision, so the reviewers' robustness questions are answered from one
# script:
#
#   Scenario                                   Reviewer comments
#   ---------------------------------------    -------------------------------
#   Directly reported means and SDs only       R1.7 converted data, R1.10, R3.3
#   All reports incl. possible cohort overlap  R1.1 participant overlap, R1.7
#   Yang L 2023 substituted for Yang F 2023    R1.1 overlap (same hospital)
#   Both Yang 2023 reports omitted             R1.1
#   de Geus 2023 substituted for 2025          R1.1 overlap (same biobank)
#   Said 2022 added (range-derived)            R1.7 converted data, R1.10
#   Torretta 2018 / 2021 added, one at a time  R1.7, R1.10 (range-derived)
#   Kamalian 2024 added (MCI-AD comparator)    R1.8 comparator definition
#   Clinical controls added (vs HC)            R1.7 control selection, R1.8,
#                                              R2.1
#   Studies with both an HC and an AD arm      R1.7 dual-comparator pattern
#
# It also writes the per-study HC and AD SMDs for studies with both arms, the
# age meta-regression with each Torretta report added (decision 4), and the
# R session information (package versions, R1.10). Prediction intervals
# (R1.7) are produced by Step 02 for the primary analysis and here for every
# scenario.
#
# Where the rows come from
# ------------------------
# - Primary rows: output/tables/clean_data.rds (Step 01).
# - Sensitivity rows: the "List nonprimary" sheet of each workbook. Rows are
#   selected by their Disposition text (e.g. "All-reports sensitivity",
#   "Substitute for Yang F") or by study name; rows on hold (lumbar n unknown,
#   overlapping arms, different index, non-AD comparator) are never used.
# - Direct vs converted means/SDs: data/conversion_status.csv, one row per
#   primary row, derived from Dr. Nila's conversion_audit.csv (24 Nov 2026).
#   direct_mean_sd = 0 marks means/SDs estimated from medians and quartiles
#   (Luo mean, Wan SD).
#
# Outputs
# -------
# output/tables/sensitivity_smd.csv
# output/tables/sensitivity_dual_comparator_studies.csv
# output/tables/sensitivity_age_addback.csv
# output/tables/Supplementary_Table_S6_sensitivity.docx
# output/tables/session_info.txt
#
# Run order
# ---------
# Requires Steps 01 and 02. Independent of Steps 03-11.
# =============================================================================

source("R/00_setup.R")
source("R/load_data.R")            # clean_biomarker()
source("R/meta_smd.R")             # run_smd_one(), smd_prediction_interval()
source("R/theme_bist.R")           # biomarker_label(), comparison_style()
source("R/meta_regression.R")      # run_age_metareg()

if (!requireNamespace("officer", quietly = TRUE))
  install.packages("officer", repos = "https://cloud.r-project.org")
library(officer)

# Proposed number: S6 became free when the SMD-derived AUC tables were
# dropped. Change here if the manuscript numbers it differently.
SENS_TABLE_LABEL <- "Supplementary Table S6"

# Mean iNPH age of the Torretta cohorts (as entered in column R before the
# rows moved to "List nonprimary"); used only for the age add-back.
TORRETTA_AGE <- c("Torretta et al., 2018" = 85, "Torretta et al., 2021" = 86)

COLS <- c("study", "biomarker", "comparison", "sheet",
          "En", "Em", "Esd", "Cn", "Cm", "Csd", "CSF_Type", "iNPH_age")
WB   <- c(iNPH_vs_HC = file.path(PATHS$data, "CSF_iNPH_vs_HC.xlsx"),
          iNPH_vs_AD = file.path(PATHS$data, "CSF_iNPH_vs_AD.xlsx"))

# ---- Primary rows -----------------------------------------------------------
prim <- readRDS(file.path(PATHS$tables, "clean_data.rds"))$smd[, COLS]

# ---- Non-primary rows (both workbooks) --------------------------------------
read_nonprimary <- function(path, comparison) {
  np <- as.data.frame(suppressMessages(readxl::read_excel(
    path, sheet = "List nonprimary", .name_repair = "minimal")))
  num <- function(v) suppressWarnings(as.numeric(trimws(as.character(v))))
  out <- data.frame(
    disposition = trimws(as.character(np[[1]])),
    study       = trimws(as.character(np[[2]])),
    sheet       = trimws(as.character(np[[3]])),
    En = num(np[[4]]), Em = num(np[[5]]), Esd = num(np[[6]]),
    Cn = num(np[[7]]), Cm = num(np[[8]]), Csd = num(np[[9]]),
    stringsAsFactors = FALSE)
  out <- out[!is.na(out$study) & out$study != "", , drop = FALSE]
  out$biomarker  <- clean_biomarker(out$sheet)
  out$comparison <- comparison
  out$CSF_Type   <- NA_character_
  out$iNPH_age   <- NA_real_
  usable <- with(out, !is.na(En) & !is.na(Em) & !is.na(Esd) & !is.na(Cn) &
                   !is.na(Cm) & !is.na(Csd) & Esd > 0 & Csd > 0)
  out[usable, , drop = FALSE]
}
np_all <- rbind(read_nonprimary(WB[["iNPH_vs_HC"]], "iNPH_vs_HC"),
                read_nonprimary(WB[["iNPH_vs_AD"]], "iNPH_vs_AD"))
message("Non-primary rows with complete data: ", nrow(np_all))

np_disp  <- function(prefix)                      # rows by Disposition text
  np_all[startsWith(np_all$disposition, prefix), COLS, drop = FALSE]
np_study <- function(pattern, comparison = NULL) { # rows by study name
  hit <- grepl(pattern, np_all$study)
  if (!is.null(comparison)) hit <- hit & np_all$comparison == comparison
  np_all[hit, COLS, drop = FALSE]
}
add_rows <- function(base, rows) rbind(base, rows[, COLS, drop = FALSE])
# Replace the rows of `drop_pattern` studies by `subs`, only for the
# biomarker x comparison combinations that `subs` covers.
swap_rows <- function(base, subs, drop_pattern) {
  covered <- paste(subs$comparison, subs$biomarker)
  drop <- grepl(drop_pattern, base$study) &
    paste(base$comparison, base$biomarker) %in% covered
  rbind(base[!drop, , drop = FALSE], subs[, COLS, drop = FALSE])
}
within_study <- function(d) {
  id   <- paste(d$study, d$biomarker)
  both <- intersect(id[d$comparison == "iNPH_vs_HC"],
                    id[d$comparison == "iNPH_vs_AD"])
  d[id %in% both, , drop = FALSE]
}

# ---- Direct vs converted means/SDs ------------------------------------------
conv <- read.csv(file.path(PATHS$data, "conversion_status.csv"),
                 stringsAsFactors = FALSE, fileEncoding = "UTF-8")
row_key <- function(d) enc2utf8(paste(d$comparison, d$study, d$sheet, sep = "|"))
direct  <- conv$direct_mean_sd[match(row_key(prim), row_key(conv))]
if (anyNA(direct)) {
  warning("conversion_status.csv has no entry for ", sum(is.na(direct)),
          " primary row(s); they are left out of the direct-only analysis:\n  ",
          paste(row_key(prim)[is.na(direct)], collapse = "\n  "), call. = FALSE)
}
message("Primary rows with directly reported means/SDs: ",
        sum(direct %in% 1), " of ", nrow(prim))

# ---- Study characteristics (coded by Dr. Nila, Nov 2026) --------------------
# Used for two scenarios added after the coding: AD comparisons in which CSF
# contributed to the AD diagnosis (incorporation; R1.8, R2.1), and Jurjevic
# 2017, whose AD comparator was possible iNPH with AD features identified
# clinically or by an AD-pattern CSF profile.
chars_path <- file.path(PATHS$data, "study_characteristics.xlsx")
csf_dx_studies <- character(0)
if (file.exists(chars_path)) {
  chars <- as.data.frame(suppressMessages(readxl::read_excel(chars_path,
                                                             sheet = "Coding")))
  csf_dx_studies <- chars$study[chars$AD_dx_used_CSF %in% "Yes"]
  message("AD comparisons with CSF-based AD diagnosis: ",
          paste(csf_dx_studies, collapse = "; "))
} else {
  message("(study_characteristics.xlsx not found - CSF-diagnosis scenario skipped)")
}
drop_ad <- function(d, studies)
  d[!(d$comparison == "iNPH_vs_AD" & d$study %in% studies), , drop = FALSE]

# ---- Scenarios --------------------------------------------------------------
clinical_ctrl <- rbind(np_disp("Clinical controls only"),
                       np_disp("Clinical CU only")[
                         grepl("^de Geus et al\\., 2025",
                               np_disp("Clinical CU only")$study), ])
SCENARIOS <- list(
  list(id = "primary", label = "Primary analysis", data = prim),
  list(id = "direct_only", label = "Directly reported means and SDs only",
       data = prim[direct %in% 1, , drop = FALSE]),
  list(id = "all_reports",
       label = "All reports, including possible cohort overlap",
       data = add_rows(prim, np_disp("All-reports"))),
  list(id = "yang_L", label = "Yang L 2023 substituted for Yang F 2023",
       data = swap_rows(prim, np_disp("Substitute for Yang F"),
                        "^Yang F et al")),
  list(id = "yang_none", label = "Both Yang 2023 reports omitted",
       data = prim[!grepl("^Yang [FL] et al", prim$study), , drop = FALSE]),
  list(id = "degeus_2023", label = "de Geus 2023 substituted for de Geus 2025",
       data = swap_rows(prim, np_disp("Substitute for de Geus 2025"),
                        "^de Geus et al\\., 2025")),
  list(id = "said", label = "Said 2022 added (range-derived)",
       data = add_rows(prim, np_study("^Said et al"))),
  list(id = "torretta_2018", label = "Torretta 2018 added (range-derived)",
       data = add_rows(prim, np_study("^Torretta et al\\., 2018"))),
  list(id = "torretta_2021", label = "Torretta 2021 added (range-derived)",
       data = add_rows(prim, np_study("^Torretta et al\\., 2021"))),
  list(id = "kamalian",
       label = "Kamalian 2024 added (biomarker-positive MCI as comparator)",
       data = add_rows(prim, np_study("^Kamalian", "iNPH_vs_AD"))),
  list(id = "clinical_controls",
       label = "Clinical controls added to healthy controls",
       data = add_rows(prim, clinical_ctrl)),
  list(id = "within_study",
       label = "Studies with both a healthy-control and an AD arm",
       data = within_study(prim)),
  list(id = "ad_csf_dx_omitted",
       label = "AD comparisons with CSF-based AD diagnosis omitted",
       data = drop_ad(prim, csf_dx_studies)),
  list(id = "jurjevic_ad_omitted",
       label = "Jurjevic 2017 AD comparison omitted",
       data = drop_ad(prim, prim$study[grepl("^Jurjevi", prim$study)]))
)

# Biomarker x comparison combinations pooled in the primary analysis.
k_prim  <- table(paste(prim$comparison, prim$biomarker, sep = "::"))
TARGETS <- names(k_prim)[k_prim >= MIN_K]

row_sig <- function(d) {
  d <- d[order(d$study), , drop = FALSE]
  paste(d$study, d$En, d$Em, d$Esd, d$Cn, d$Cm, d$Csd, sep = "|")
}

# ---- Run --------------------------------------------------------------------
message("Running ", length(SCENARIOS), " scenarios over ", length(TARGETS),
        " pooled biomarker x comparison combinations ...")
rows <- list()
for (s in SCENARIOS) {
  for (tg in TARGETS) {
    cmp <- strsplit(tg, "::")[[1]][1]; bm <- strsplit(tg, "::")[[1]][2]
    p0 <- prim[prim$comparison == cmp & prim$biomarker == bm, , drop = FALSE]
    d  <- s$data[s$data$comparison == cmp & s$data$biomarker == bm, ,
                 drop = FALSE]
    same <- identical(row_sig(d), row_sig(p0))
    if (s$id != "primary" && same) next          # scenario does not touch it
    if (anyDuplicated(d$study))
      stop("duplicate study rows in scenario '", s$id, "' for ", tg)
    added   <- setdiff(d$study, p0$study)
    removed <- setdiff(p0$study, d$study)
    base <- data.frame(
      scenario = s$id, analysis = s$label, comparison = cmp, biomarker = bm,
      k = nrow(d), estimate = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_,
      pval = NA_real_, I2 = NA_real_, tau2 = NA_real_,
      pi_lb = NA_real_, pi_ub = NA_real_,
      studies_added = paste(added, collapse = "; "),
      studies_removed = paste(removed, collapse = "; "),
      stringsAsFactors = FALSE)
    if (nrow(d) >= MIN_K) {
      r <- run_smd_one(d)$summary
      base[c("estimate", "ci_lb", "ci_ub", "pval", "I2", "tau2",
             "pi_lb", "pi_ub")] <-
        r[c("estimate", "ci_lb", "ci_ub", "pval", "I2", "tau2",
            "pi_lb", "pi_ub")]
    }
    rows[[length(rows) + 1]] <- base
  }
}
sens <- do.call(rbind, rows)

# Change from the primary estimate (same biomarker x comparison), and whether
# the 95% CI moves to the other side of zero (or onto it).
key     <- paste(sens$comparison, sens$biomarker)
is_prim <- sens$scenario == "primary"
is_sig  <- !is.na(sens$ci_lb) & (sens$ci_lb > 0 | sens$ci_ub < 0)
p_est   <- setNames(sens$estimate[is_prim], key[is_prim])
p_sig   <- setNames(is_sig[is_prim], key[is_prim])
sens$change_from_primary  <- sens$estimate - unname(p_est[key])
sens$significance_changed <- ifelse(is.na(sens$ci_lb), NA,
                                    is_sig != unname(p_sig[key]))

write.csv(sens, file.path(PATHS$tables, "sensitivity_smd.csv"),
          row.names = FALSE)
message("Sensitivity estimates written to sensitivity_smd.csv (",
        nrow(sens), " rows)")

flip <- sens[sens$significance_changed %in% TRUE, , drop = FALSE]
if (nrow(flip)) {
  message("\nScenarios in which the 95% CI changes side of zero:")
  print(flip[, c("analysis", "comparison", "biomarker", "k", "estimate",
                 "ci_lb", "ci_ub")], row.names = FALSE, digits = 3)
} else {
  message("\nNo scenario moves a 95% CI across zero.")
}

# ---- Dual-comparator pattern within studies ---------------------------------
ws <- within_study(prim)
es <- metafor::escalc(measure = "SMD", m1i = Em, sd1i = Esd, n1i = En,
                      m2i = Cm, sd2i = Csd, n2i = Cn, data = ws)
es$se <- sqrt(es$vi)
hc <- es[es$comparison == "iNPH_vs_HC", c("study", "biomarker", "yi", "se")]
ad <- es[es$comparison == "iNPH_vs_AD", c("study", "biomarker", "yi", "se")]
dual <- merge(hc, ad, by = c("study", "biomarker"), suffixes = c("_HC", "_AD"))
names(dual) <- c("study", "biomarker", "smd_vs_HC", "se_vs_HC",
                 "smd_vs_AD", "se_vs_AD")
dual$pattern <- ifelse(dual$smd_vs_HC < 0 & dual$smd_vs_AD > 0,
                       "lower than HC, higher than AD",
                ifelse(dual$smd_vs_HC < 0 & dual$smd_vs_AD < 0,
                       "lower than both",
                ifelse(dual$smd_vs_HC > 0 & dual$smd_vs_AD > 0,
                       "higher than both", "higher than HC, lower than AD")))
dual <- dual[order(dual$biomarker, dual$study), ]
write.csv(dual, file.path(PATHS$tables,
                          "sensitivity_dual_comparator_studies.csv"),
          row.names = FALSE)
message("\nStudies with both an HC and an AD arm (per-study SMDs in ",
        "sensitivity_dual_comparator_studies.csv):")
for (bm in intersect(c("AB42", "pTau181", "tTau", "NFL"), dual$biomarker)) {
  tab <- table(dual$pattern[dual$biomarker == bm])
  message(sprintf("  %-8s k = %2d: %s", bm, sum(tab),
                  paste(sprintf("%s %d", names(tab), tab), collapse = "; ")))
}

# ---- Age meta-regression with each Torretta report added (decision 4) -------
age_rows <- list()
for (cmp in c("iNPH_vs_HC", "iNPH_vs_AD")) {
  fit_row <- function(dat, label) {
    reg <- suppressWarnings(run_age_metareg(dat, "AB42", cmp))
    data.frame(comparison = cmp, analysis = label, k = reg$k,
               slope_per_year = reg$slope, ci_lb = reg$ci[1],
               ci_ub = reg$ci[2], p_value = reg$p, R2_percent = reg$R2)
  }
  age_rows[[length(age_rows) + 1]] <- fit_row(prim, "Primary analysis")
  for (st in names(TORRETTA_AGE)) {
    tr <- np_all[np_all$study == st & np_all$biomarker == "AB42" &
                   np_all$comparison == cmp, COLS, drop = FALSE]
    if (!nrow(tr)) next
    tr$iNPH_age <- TORRETTA_AGE[[st]]
    age_rows[[length(age_rows) + 1]] <-
      fit_row(rbind(prim, tr), sprintf("%s added (mean age %g)", st,
                                       TORRETTA_AGE[[st]]))
  }
}
age_tbl <- do.call(rbind, age_rows)
write.csv(age_tbl, file.path(PATHS$tables, "sensitivity_age_addback.csv"),
          row.names = FALSE)
message("\nAge meta-regression (Abeta42) with Torretta add-backs:")
print(transform(age_tbl, slope_per_year = round(slope_per_year, 4),
                ci_lb = round(ci_lb, 4), ci_ub = round(ci_ub, 4),
                p_value = signif(p_value, 3), R2_percent = round(R2_percent, 1)),
      row.names = FALSE)

# ---- Word table -------------------------------------------------------------
fmt_ci <- function(e, l, h) ifelse(is.na(e), "Not pooled",
                                   sprintf("%.2f (%.2f to %.2f)", e, l, h))
fmt_pi <- function(lb, ub, k) ifelse(is.na(lb) | k < PI_MIN_K, "\u2014",
                                     sprintf("%.2f to %.2f", lb, ub))
fmt_i2 <- function(x) ifelse(is.na(x), "\u2014", sprintf("%.1f", x))

doc <- read_docx()
doc <- body_add_par(doc, paste0(SENS_TABLE_LABEL,
  ". Sensitivity analyses of the pooled standardized mean differences"),
  style = "heading 2")
for (cmp in c("iNPH_vs_HC", "iNPH_vs_AD")) {
  blk <- sens[sens$comparison == cmp, , drop = FALSE]
  if (!nrow(blk)) next
  prim_blk <- blk[blk$scenario == "primary", , drop = FALSE]
  bm_order <- prim_blk$biomarker[order(-abs(prim_blk$estimate))]
  scen_order <- vapply(SCENARIOS, `[[`, "", "id")
  out <- list()
  for (bm in bm_order) {
    b <- blk[blk$biomarker == bm, , drop = FALSE]
    b <- b[order(match(b$scenario, scen_order)), , drop = FALSE]
    out[[bm]] <- data.frame(
      Biomarker = c(biomarker_label(bm), rep("", nrow(b) - 1)),
      Analysis  = b$analysis,
      Studies   = b$k,
      SMD       = fmt_ci(b$estimate, b$ci_lb, b$ci_ub),
      PI        = fmt_pi(b$pi_lb, b$pi_ub, b$k),
      I2        = fmt_i2(b$I2),
      stringsAsFactors = FALSE)
  }
  m <- do.call(rbind, out)
  names(m) <- c("Biomarker", "Analysis", "No. of Studies", "SMD (95% CI)",
                "95% PI", "I\u00B2, %")
  doc <- body_add_par(doc, comparison_style(cmp)$label, style = "Normal")
  doc <- body_add_table(doc, m, style = "table_template")
  doc <- body_add_par(doc, "", style = "Normal")
}
fn <- c(
  "Each analysis re-pools the biomarker with the same random-effects model as the primary analysis (REML, Hedges g, z-based 95% CI). Only analyses that change the set of contributing studies are listed under each biomarker.",
  "Directly reported means and SDs only: contrasts whose means and SDs were estimated from medians and interquartile ranges (Luo mean, Wan SD) are excluded.",
  "All reports, including possible cohort overlap: adds reports from the same centre or series that were excluded from the primary analysis because participants may overlap (Jeppsson 2013 and 2016; Santangelo 2019).",
  "Yang L 2023 and Yang F 2023 come from the same hospital, with possible overlap of iNPH cases; the primary analysis uses Yang F for shared outcomes. de Geus 2023 and 2025 come from the same biobank; the primary analysis uses the 2025 report.",
  "Said 2022 and Torretta 2018 and 2021 report medians and ranges, so their means were estimated; the Torretta reports are added one at a time because the independence of their cohorts is unresolved.",
  "Kamalian 2024: the comparator group consisted of biomarker-positive patients with mild cognitive impairment rather than AD dementia.",
  "Clinical controls added to healthy controls: control groups who underwent lumbar puncture for clinical reasons (Schirinzi 2015, suspected polyneuropathy; de Geus 2025, clinically cognitively unimpaired) are added to the healthy-control comparison.",
  "Studies with both a healthy-control and an AD arm: restricted to studies contributing both comparisons for the biomarker, so the two pooled estimates come from the same studies.",
  "AD comparisons with CSF-based AD diagnosis omitted: excludes the AD comparisons of studies in which CSF biomarkers explicitly contributed to the AD diagnosis (study-characteristics coding), to address incorporation bias.",
  "Jurjevic 2017 AD comparison omitted: its comparator group was patients with possible iNPH and AD features, identified clinically or by an AD-pattern CSF profile.",
  paste0("PI, 95% prediction interval (t distribution with k \u2212 2 degrees of freedom), shown for analyses of at least ", PI_MIN_K, " studies. Not pooled: fewer than ", MIN_K, " studies."),
  "AD, Alzheimer disease; CI, confidence interval; iNPH, idiopathic normal pressure hydrocephalus; MCI, mild cognitive impairment; SMD, standardized mean difference."
)
for (line in fn) doc <- body_add_par(doc, line, style = "Normal")
out_docx <- file.path(PATHS$tables, "Supplementary_Table_S6_sensitivity.docx")
print(doc, target = out_docx)
message("\n  wrote ", out_docx)

# ---- Session information (package versions, R1.10) --------------------------
writeLines(capture.output(sessionInfo()),
           file.path(PATHS$tables, "session_info.txt"))
message("  wrote ", file.path(PATHS$tables, "session_info.txt"))

message("\n12_sensitivity_analyses.R complete.")
