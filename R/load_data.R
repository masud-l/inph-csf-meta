# =============================================================================
# load_data.R  —  Read and validate the CSF biomarker extraction workbooks
# -----------------------------------------------------------------------------
# Both workbooks share an identical column layout on every biomarker sheet:
#
#   col 1  : Study
#   col 9  : En   (iNPH group n)        col 12 : Cn   (control group n)
#   col 10 : Em   (iNPH group mean)     col 13 : Cm   (control group mean)
#   col 11 : Esd  (iNPH group SD)       col 14 : Csd  (control group SD)
#
# (1-indexed). One text sheet per workbook ("List of single biomarker(s)") is
# skipped. Every other sheet is one biomarker.
#
# IMPORTANT — sign convention: escalc(measure = "SMD", m1 = iNPH, m2 = control)
# yields iNPH minus control. A NEGATIVE SMD therefore means the biomarker is
# LOWER in iNPH than in the comparator. This matches the manuscript.
# =============================================================================

# Fixed 1-indexed column positions of the meta-analysis inputs. CSF_Type is
# column Q (17) on every biomarker sheet. iNPH_age is column R (18) and is
# present ONLY on the AB42 sheet of the vs-AD workbook (for the age
# meta-regression); every other sheet yields NA there, which is fine.
.COLPOS <- c(En = 9, Em = 10, Esd = 11, Cn = 12, Cm = 13, Csd = 14,
             unit = 15, CSF_Type = 17, iNPH_age = 18)

# ---- Clean a raw sheet name into a tidy biomarker label ---------------------
# Tab names are hand-cleaned to the biomarker only (see README). This function
# does the MINIMUM, predictable transformation:
#   * trims whitespace
#   * converts the ratio word "by" to "/" - handling BOTH the spaced form
#     ("AB42 by AB40") and the glued form ("AB42byAB40"), because "/" is an
#     illegal Excel sheet-name character so tabs must use "by".
# Whatever else you type in the tab name is preserved verbatim.
clean_biomarker <- function(sheet) {
  x <- str_trim(sheet)
  x <- str_replace_all(x, "(?i)\\s+by\\s+", "/")                       # spaced
  x <- str_replace_all(x, "(?i)([A-Za-z0-9])by([A-Za-z0-9])", "\\1/\\2") # glued
  str_squish(x)
}

# ---- Read one biomarker sheet into a tidy data frame ------------------------
read_biomarker_sheet <- function(path, sheet) {
  raw <- suppressMessages(
    readxl::read_excel(path, sheet = sheet, col_names = FALSE,
                       .name_repair = "minimal")
  )
  # Drop the header row, keep only rows with a study name in column 1.
  raw <- raw[-1, , drop = FALSE]
  raw <- raw[!is.na(raw[[1]]) & raw[[1]] != "", , drop = FALSE]
  if (nrow(raw) == 0) return(NULL)

  num <- function(v) suppressWarnings(as.numeric(trimws(as.character(v))))
  # CSF_Type is text, not numeric. A sheet narrower than column 17 (or a blank
  # cell) yields NA rather than an error.
  txt <- function(pos) {
    if (ncol(raw) < pos) return(rep(NA_character_, nrow(raw)))
    v <- trimws(as.character(raw[[pos]]))
    v[v == ""] <- NA_character_
    v
  }
  # Numeric positional reader: tolerates a missing column (returns NA), used
  # for iNPH_age which exists only on the AB42 vs-AD sheet.
  numpos <- function(pos) {
    if (ncol(raw) < pos) return(rep(NA_real_, nrow(raw)))
    suppressWarnings(as.numeric(trimws(as.character(raw[[pos]]))))
  }
  out <- data.frame(
    study = trimws(as.character(raw[[1]])),
    En  = num(raw[[.COLPOS["En"]]]),  Em  = num(raw[[.COLPOS["Em"]]]),
    Esd = num(raw[[.COLPOS["Esd"]]]), Cn  = num(raw[[.COLPOS["Cn"]]]),
    Cm  = num(raw[[.COLPOS["Cm"]]]),  Csd = num(raw[[.COLPOS["Csd"]]]),
    unit     = txt(.COLPOS["unit"]),       # column O, used by Table S4
    CSF_Type = txt(.COLPOS["CSF_Type"]),
    iNPH_age = numpos(.COLPOS["iNPH_age"]),
    stringsAsFactors = FALSE
  )
  out$biomarker  <- clean_biomarker(sheet)
  out$sheet      <- sheet
  out
}

# ---- Load an entire workbook ------------------------------------------------
# Returns a single tidy data frame for all biomarker sheets in one workbook.
# `comparison` is the label applied to every row ("iNPH_vs_HC" / "iNPH_vs_AD").
load_workbook <- function(path, comparison) {
  sheets <- readxl::excel_sheets(path)
  sheets <- sheets[!str_detect(tolower(sheets), "^\\s*list")]  # drop text sheet
  dat <- purrr::map_dfr(sheets, function(s) {
    df <- read_biomarker_sheet(path, s)
    if (is.null(df)) return(NULL)
    df$comparison <- comparison
    df
  })
  dat[, c("study", "biomarker", "comparison", "sheet",
          "En", "Em", "Esd", "Cn", "Cm", "Csd", "unit", "CSF_Type", "iNPH_age")]
}

# ---- Validation -------------------------------------------------------------
# Flags rows that cannot enter the meta-analysis and obviously implausible
# values. Returns a list: $clean (usable rows) and $report (a per-issue table).
validate_extraction <- function(dat) {
  issues <- list()

  bad_missing <- dat |>
    filter(if_any(c(En, Em, Esd, Cn, Cm, Csd), is.na))
  if (nrow(bad_missing))
    issues$missing <- bad_missing |>
      mutate(issue = "missing or non-numeric value in En/Em/Esd/Cn/Cm/Csd")

  bad_sd <- dat |> filter(!is.na(Esd) & Esd <= 0 | !is.na(Csd) & Csd <= 0)
  if (nrow(bad_sd))
    issues$sd <- bad_sd |> mutate(issue = "SD <= 0 (cannot compute SMD)")

  # An |SMD| this large almost always signals a unit / data-entry error.
  rough_smd <- with(dat, (Em - Cm) /
                      sqrt(((En - 1) * Esd^2 + (Cn - 1) * Csd^2) /
                             (En + Cn - 2)))
  bad_extreme <- dat[which(abs(rough_smd) > 5), , drop = FALSE]
  if (nrow(bad_extreme))
    issues$extreme <- bad_extreme |>
      mutate(issue = sprintf("implausible SMD ~ %.1f - check units/values",
                             rough_smd[which(abs(rough_smd) > 5)]))

  # Sheets misnamed relative to the workbook they live in.
  mismatch <- dat |>
    filter((comparison == "iNPH_vs_HC" & str_detect(sheet, "(?i)vsAD")) |
           (comparison == "iNPH_vs_AD" & str_detect(sheet, "(?i)vsHC"))) |>
    distinct(sheet, comparison)
  if (nrow(mismatch))
    issues$sheet_name <- mismatch |>
      mutate(study = NA, biomarker = NA,
             issue = "sheet name disagrees with workbook - verify comparison")

  report <- if (length(issues))
    dplyr::bind_rows(lapply(issues, function(x)
      x[, intersect(c("comparison", "sheet", "study", "biomarker", "issue"),
                    names(x))]))
  else
    data.frame()

  clean <- dat |> filter(if_all(c(En, Em, Esd, Cn, Cm, Csd), ~ !is.na(.))) |>
    filter(Esd > 0, Csd > 0)

  list(clean = clean, report = report,
       n_in = nrow(dat), n_clean = nrow(clean),
       n_dropped = nrow(dat) - nrow(clean))
}

# ---- Diagnostic 2x2 loader --------------------------------------------------
load_diagnostic <- function(path) {
  d <- suppressMessages(readxl::read_excel(path))
  names(d) <- c("study", "biomarker", "comparison",
                "tp", "fn", "fp", "tn",
                "sens_rep", "spec_rep", "cutoff", "auc", "En", "Cn")
  d |>
    mutate(across(c(tp, fn, fp, tn, En, Cn, auc), as.numeric),
           sens_calc = 100 * tp / (tp + fn),
           spec_calc = 100 * tn / (tn + fp),
           zero_cell = tp == 0 | fn == 0 | fp == 0 | tn == 0)
}
