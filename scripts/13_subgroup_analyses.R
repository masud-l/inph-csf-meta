# =============================================================================
# Step 13 - Exploratory subgroup analyses by study characteristics
# Script: 13_subgroup_analyses.R
# =============================================================================
#
# Purpose
# -------
# Reviewer requests R1.7, R1.8, R2.1 and R3.3: do iNPH diagnostic criteria,
# control type, assay platform, the indication for CSF sampling, the AD
# criteria or restriction to shunt responders explain heterogeneity?
#
# Rules (fixed before looking at the results, agreed with the lead author)
# -----------------------------------------------------------------------
# - Characteristics come from data/study_characteristics.xlsx (sheet
#   "Coding"; coded by Dr. Nila, Nov 2026), joined to the primary rows by
#   study label.
# - "Not reported" is never merged with "No" or any other group; such studies
#   are counted but not tested.
# - A characteristic is tested for a biomarker x comparison only when at least
#   two documented groups each have >= SUB_MIN_K (3) studies. Smaller groups
#   are counted but not tested.
# - Two merges of closely related categories: all Japanese guideline editions
#   form one group, and assays are split into INNOTEST ELISA vs other
#   platforms.
# - Model: random-effects (REML, z-based) pooling within each group, and the
#   omnibus test of a mixed-effects meta-regression with the grouping as a
#   categorical moderator (as in the age subgroup analysis, Step 10).
# - Exploratory and not prespecified in the protocol; no multiplicity
#   adjustment. Report as hypothesis-generating.
# - Design, CSF use in the AD diagnosis, exclusion of AD-positive iNPH and iNPH
#   certainty cannot be tested (fewer than two documented groups of >= 3);
#   their study counts are reported descriptively.
#
# Outputs
# -------
# output/tables/subgroup_estimates.csv   (one row per group, tested or not)
# output/tables/subgroup_tests.csv       (one row per test)
# output/tables/subgroup_counts.csv      (study-level counts, all characteristics)
# output/tables/Supplementary_Table_S8_subgroups.docx
#
# Run order
# ---------
# Requires Steps 01 and 02.
# =============================================================================

source("R/00_setup.R")
source("R/load_data.R")
source("R/meta_smd.R")
source("R/theme_bist.R")           # biomarker_label(), comparison_style()

if (!requireNamespace("officer", quietly = TRUE))
  install.packages("officer", repos = "https://cloud.r-project.org")
library(officer)

# Proposed number (next free supplementary table); change here if needed.
SUBGROUP_TABLE_LABEL <- "Supplementary Table S8"
SUB_MIN_K <- 3
OUTCOMES  <- c("AB42", "pTau181", "tTau", "NFL")
ASSAY_COL <- c(AB42 = "assay_AB42", pTau181 = "assay_pTau181",
               tTau = "assay_tTau", NFL = "assay_NfL")
HC <- "iNPH_vs_HC"; AD <- "iNPH_vs_AD"

TESTED <- list(
  list(id = "iNPH_criteria", label = "iNPH diagnostic criteria", cmp = c(HC, AD)),
  list(id = "HC_control_type", label = "Healthy-control type", cmp = HC),
  list(id = "assay", label = "Assay platform", cmp = c(HC, AD)),
  list(id = "iNPH_CSF_indication", label = "Indication for CSF sampling in iNPH",
       cmp = c(HC, AD)),
  list(id = "AD_criteria", label = "AD diagnostic criteria", cmp = AD),
  list(id = "shunt_responders_only", label = "iNPH restricted to shunt responders",
       cmp = c(HC, AD))
)
DESCRIPTIVE <- c(design = "Recruitment design",
                 AD_dx_used_CSF = "CSF used in the AD diagnosis",
                 ADpos_iNPH_excluded = "AD-biomarker-positive iNPH excluded",
                 iNPH_certainty = "iNPH diagnostic certainty")

# ---- Data -------------------------------------------------------------------
chars <- as.data.frame(suppressMessages(readxl::read_excel(
  file.path(PATHS$data, "study_characteristics.xlsx"), sheet = "Coding")))
chars <- chars[!is.na(chars$study) & chars$study != "", , drop = FALSE]
prim  <- readRDS(file.path(PATHS$tables, "clean_data.rds"))$smd
missing_studies <- setdiff(unique(prim$study), chars$study)
if (length(missing_studies))
  stop("No coding for: ", paste(missing_studies, collapse = "; "))

recode_group <- function(var, x) {
  x <- ifelse(is.na(x) | trimws(x) == "", "Not reported", trimws(x))
  if (var == "iNPH_criteria")
    x <- ifelse(grepl("^Japanese", x), "Japanese guidelines (any edition)", x)
  if (var == "assay")
    x <- ifelse(grepl("^INNOTEST", x), "INNOTEST ELISA",
                ifelse(x == "Not reported", x, "Other platforms"))
  x
}
testable_group <- function(g) g != "Not reported" & !grepl("^Not applicable", g)

# ---- Run --------------------------------------------------------------------
est_rows <- list(); test_rows <- list()
for (v in TESTED) {
  for (bm in OUTCOMES) {
    for (cmp in v$cmp) {
      d <- prim[prim$biomarker == bm & prim$comparison == cmp, , drop = FALSE]
      if (!nrow(d)) next
      col <- if (v$id == "assay") ASSAY_COL[[bm]] else v$id
      d$group <- recode_group(v$id, chars[[col]][match(d$study, chars$study)])
      n_by <- table(d$group)
      ok <- names(n_by)[testable_group(names(n_by)) & n_by >= SUB_MIN_K]
      tested <- length(ok) >= 2

      es <- metafor::escalc(measure = "SMD", m1i = d$Em, sd1i = d$Esd,
                            n1i = d$En, m2i = d$Cm, sd2i = d$Csd, n2i = d$Cn,
                            data = d)
      for (g in names(n_by)) {
        r <- data.frame(variable = v$id, characteristic = v$label,
                        biomarker = bm, comparison = cmp, group = g,
                        k = as.integer(n_by[[g]]), in_test = tested && g %in% ok,
                        estimate = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_,
                        I2 = NA_real_, stringsAsFactors = FALSE)
        if (r$in_test) {
          f <- metafor::rma(yi, vi, data = es[es$group == g, ],
                            method = RE_METHOD, test = "z")
          r[c("estimate", "ci_lb", "ci_ub", "I2")] <-
            list(as.numeric(f$b), f$ci.lb, f$ci.ub, f$I2)
        }
        est_rows[[length(est_rows) + 1]] <- r
      }
      if (tested) {
        sub <- es[es$group %in% ok, ]
        sub$group <- factor(sub$group)
        fm <- metafor::rma(yi, vi, mods = ~ group, data = sub,
                           method = RE_METHOD, test = "z")
        test_rows[[length(test_rows) + 1]] <- data.frame(
          variable = v$id, characteristic = v$label, biomarker = bm,
          comparison = cmp, groups = paste(levels(sub$group), collapse = " | "),
          k_tested = nrow(sub),
          QM = fm$QM, df = if (!is.null(fm$QMdf)) fm$QMdf[1] else fm$m,
          p_between = fm$QMp, stringsAsFactors = FALSE)
      }
    }
  }
}
est   <- do.call(rbind, est_rows)
tests <- do.call(rbind, test_rows)
write.csv(est, file.path(PATHS$tables, "subgroup_estimates.csv"), row.names = FALSE)
write.csv(tests, file.path(PATHS$tables, "subgroup_tests.csv"), row.names = FALSE)

message("\nExploratory subgroup tests (", nrow(tests), "):")
print(transform(tests[, c("characteristic", "biomarker", "comparison",
                          "k_tested", "df", "QM", "p_between")],
                QM = round(QM, 2), p_between = signif(p_between, 3)),
      row.names = FALSE)

# ---- Study-level counts of every coded characteristic ----------------------
cnt_rows <- list()
all_vars <- c(setNames(vapply(TESTED, `[[`, "", "label"),
                       vapply(TESTED, `[[`, "", "id")), DESCRIPTIVE)
all_vars <- all_vars[names(all_vars) != "assay"]
has_ad <- unique(prim$study[prim$comparison == AD])
has_hc <- unique(prim$study[prim$comparison == HC])
for (id in names(all_vars)) {
  pool_set <- if (id %in% c("AD_criteria", "AD_dx_used_CSF")) has_ad
              else if (id == "HC_control_type") has_hc
              else unique(prim$study)
  x <- recode_group(id, chars[[id]][match(pool_set, chars$study)])
  tb <- table(x)
  cnt_rows[[length(cnt_rows) + 1]] <- data.frame(
    variable = id, characteristic = all_vars[[id]], group = names(tb),
    studies = as.integer(tb), of = length(pool_set), stringsAsFactors = FALSE)
}
counts <- do.call(rbind, cnt_rows)
write.csv(counts, file.path(PATHS$tables, "subgroup_counts.csv"), row.names = FALSE)

# ---- Word table -------------------------------------------------------------
fmt_est <- function(e, l, h) ifelse(is.na(e), "\u2014",
                                    sprintf("%.2f (%.2f to %.2f)", e, l, h))
fmt_p   <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<.001",
                              sub("^0", "", sprintf("%.3f", p))))
cmp_lab <- function(c) ifelse(c == HC, "vs HC", "vs AD")

doc <- read_docx()
doc <- body_add_par(doc, paste0(SUBGROUP_TABLE_LABEL,
  ". Exploratory subgroup analyses of the pooled standardized mean differences by study characteristics"),
  style = "heading 2")
doc <- body_add_par(doc, "A. Tested characteristics", style = "Normal")
tab <- list()
for (i in seq_len(nrow(tests))) {
  t <- tests[i, ]
  e <- est[est$variable == t$variable & est$biomarker == t$biomarker &
             est$comparison == t$comparison, , drop = FALSE]
  e <- e[order(!e$in_test, e$group), , drop = FALSE]
  first <- seq_len(nrow(e)) == 1
  tab[[i]] <- data.frame(
    Characteristic = ifelse(first, t$characteristic, ""),
    Outcome   = ifelse(first, paste(biomarker_label(t$biomarker),
                                    cmp_lab(t$comparison)), ""),
    Group     = paste0(e$group, ifelse(e$in_test, "", " (not tested)")),
    Studies   = e$k,
    SMD       = fmt_est(e$estimate, e$ci_lb, e$ci_ub),
    I2        = ifelse(is.na(e$I2), "\u2014", sprintf("%.1f", e$I2)),
    P         = ifelse(first, fmt_p(t$p_between), ""),
    stringsAsFactors = FALSE)
}
ma <- do.call(rbind, tab)
names(ma) <- c("Characteristic", "Outcome", "Group", "No. of Studies",
               "SMD (95% CI)", "I\u00B2, %", "P for subgroup difference")
doc <- body_add_table(doc, ma, style = "table_template")
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_par(doc, "B. Characteristics reported descriptively (number of studies)",
                    style = "Normal")
cd <- counts[counts$variable %in% names(DESCRIPTIVE), , drop = FALSE]
first <- !duplicated(cd$variable)
mb <- data.frame(Characteristic = ifelse(first, cd$characteristic, ""),
                 Group = cd$group,
                 Studies = paste0(cd$studies, " of ", cd$of),
                 stringsAsFactors = FALSE)
names(mb) <- c("Characteristic", "Category", "Studies")
doc <- body_add_table(doc, mb, style = "table_template")
doc <- body_add_par(doc, "", style = "Normal")

fn <- c(
  "Exploratory analyses, not prespecified in the protocol, without adjustment for multiple testing; they are hypothesis-generating.",
  paste0("Characteristics were coded from the primary publications. A characteristic was tested for an outcome only when at least two documented groups each contained at least ", SUB_MIN_K, " studies; groups marked \"not tested\" (fewer studies, or not reported) are shown for completeness. \"Not reported\" was kept separate from \"No\"."),
  "Japanese guideline editions were combined into one group, and assays were grouped as INNOTEST ELISA versus other platforms.",
  "Each group was pooled with a random-effects model (REML); the P value is the omnibus test of a mixed-effects meta-regression with the grouping as a categorical moderator.",
  "Characteristics in part B could not be tested because fewer than two documented categories contained at least 3 studies. Counts are of the studies in the primary analysis (AD-related characteristics: studies with an AD comparison).",
  "AD, Alzheimer disease; CI, confidence interval; HC, healthy controls; iNPH, idiopathic normal pressure hydrocephalus; SMD, standardized mean difference."
)
for (line in fn) doc <- body_add_par(doc, line, style = "Normal")
out_docx <- file.path(PATHS$tables, "Supplementary_Table_S8_subgroups.docx")
print(doc, target = out_docx)
message("\n  wrote ", out_docx)
message("\n13_subgroup_analyses.R complete.")
