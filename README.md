# inph-csf-meta

[![DOI](https://zenodo.org/badge/DOI/ZENODO_DOI_PLACEHOLDER.svg)](https://doi.org/ZENODO_DOI_PLACEHOLDER)

R pipeline and study-level data for the systematic review and meta-analysis

> **Cerebrospinal fluid biomarker profiles of idiopathic normal pressure hydrocephalus compared with healthy
> controls and Alzheimer's disease: a systematic review and meta-analysis.**
> Irin Sultana Nila, Md. Masudul Haque, Kyunghun Kang, Ki-Su Park. *Scientific Reports* (under review).

The pipeline pools standardized mean differences (SMDs) of CSF biomarkers between patients with idiopathic
normal pressure hydrocephalus (iNPH) and two comparator groups (healthy controls, Alzheimer's disease),
pools diagnostic sensitivity and specificity for iNPH versus Alzheimer's disease, and runs the leave-one-out,
sensitivity, exploratory subgroup, age meta-regression and small-study analyses. It renders Figure 3,
Table 2, Supplementary Figures S1–S9 and Supplementary Tables S3–S8 of the article.

## Repository layout

```
inph-csf-meta/
├── data/            extracted study-level data (see data/README.md)
├── R/               reusable modules: data loading, SMD and diagnostic models, meta-regression, plotting
├── scripts/         numbered run scripts 01–13 (run in order)
├── output/tables/   result tables (CSV and Word) from the published run
├── CITATION.cff     how to cite this archive
├── LICENSE          MIT (code)
└── LICENSE-DATA.md  CC BY 4.0 (data and result tables)
```

Figures are not stored; they are regenerated in `output/figures/` by the scripts.

## Requirements

R 4.6.1 with metafor 5.2-1, mada 0.5.12, mvmeta 1.0.3, readxl 1.6.0, dplyr 1.2.1, tidyr 1.3.2,
stringr 1.6.0, purrr 1.2.2, scales 1.4.0, ggplot2 4.0.4, patchwork 1.3.2 and officer 0.7.9 (the versions
used for the published results; the full session information is in `output/tables/session_info.txt`).
Script 01 installs missing packages on first run.

## How to run

Open `inph-csf-meta.Rproj` in RStudio and run the scripts from the project root, in order:

```r
source("scripts/01_validate_data.R")           # data checks; builds the analysis dataset
source("scripts/02_smd_metaanalysis.R")        # SMD meta-analysis, prediction intervals, leave-one-out
source("scripts/03_diagnostic_metaanalysis.R") # diagnostic accuracy vs AD + robustness analyses
source("scripts/03b_ptau_univariate.R")        # separate logit models for sensitivity and specificity
source("scripts/04_render_figures.R")          # single forest plots and summary ROC plots
source("scripts/05_figure2_overview.R")        # overview of pooled SMDs (Figure 2 data)
source("scripts/06_supplementary_preview.R")   # optional preview sheets
source("scripts/07_tables.R")                  # Table 2, Supplementary Tables S3, S4, S5 and S7 (Word)
source("scripts/08_publication_bias.R")        # funnel plots and Egger-type tests
source("scripts/09_sensitivity_loo.R")         # leave-one-out panels
source("scripts/10_age_metaregression.R")      # age subgroups and meta-regression
source("scripts/11_composite_figures.R")       # Figure 3 and Supplementary Figures S1–S9
source("scripts/12_sensitivity_analyses.R")    # SMD sensitivity analyses (Supplementary Table S6)
source("scripts/13_subgroup_analyses.R")       # exploratory subgroups (Supplementary Table S8)
```

Script 02 depends on 01; scripts 04–13 depend on 02, and 04, 07 and 11 also use the diagnostic outputs of 03.

## Methods implemented

| Setting (in `R/00_setup.R`) | Value | Meaning |
|---|---|---|
| `RE_METHOD` | `"REML"` | Random-effects estimator for SMDs (Hedges' g). |
| `USE_KNHA` | `FALSE` | Wald-type (z-based) confidence intervals. |
| `MIN_K` | 3 | Minimum number of independent studies for pooling; otherwise study-level estimates only. |
| `PI_MIN_K` | 5 | Minimum number of studies for reporting a 95% prediction interval (t distribution, k − 2 df). |

- **Independent comparisons.** Each participant group contributes once to each analysis. Two cohorts compared
  with the same comparator groups are combined; overlapping reports from the same centre are used only in
  sensitivity analyses (sheet `List nonprimary`).
- **Diagnostic accuracy (script 03).** iNPH is the positive condition. Bivariate (Reitsma) model when at least
  four studies contribute and the between-study correlation is not at its boundary; otherwise separate
  random-effects logit models without a summary AUC. A continuity correction of 0.5 is applied to tables with
  a zero cell; no estimate is pooled from fewer than three studies.
- **Sensitivity analyses (script 12).** Directly reported means and SDs only; all reports including possible
  cohort overlap; alternative reports from the same centre; median/range-derived data and a biomarker-positive
  MCI comparator added; clinical controls added; studies with both comparator arms; AD comparisons with
  CSF-based AD diagnosis omitted.
- **Subgroup analyses (script 13, exploratory).** iNPH criteria, control type, assay, CSF sampling setting,
  AD criteria and shunt-responder restriction; a characteristic is tested only when at least two documented
  categories each contain at least three studies, with the omnibus test of a mixed-effects meta-regression.
- **Small-study effects (script 08).** Funnel plots with the modified Egger test of Pustejovsky and Rodgers
  (primary) and the conventional Egger test (sensitivity) for biomarkers with at least 10 studies.

## Not produced by this pipeline

Figure 1 (PRISMA flow diagram), Figure 2 (drawn from `output/tables/Figure2_overview_data.csv`), Table 1,
Supplementary Tables S1–S2 (JBI appraisal), S9 (QUADAS-2) and S10 (study characteristics, from
`data/study_characteristics.xlsx`) were prepared outside the pipeline.

## Citation

Please cite the article and this archive (`CITATION.cff`; DOI: ZENODO_DOI_PLACEHOLDER). When reusing
the extracted values, please also cite the original studies.

## Licence

Code: MIT (`LICENSE`). Data in `data/` and tables in `output/tables/`: CC BY 4.0 (`LICENSE-DATA.md`).

## Contact

Ki-Su Park (corresponding author), Department of Neurosurgery, School of Medicine, Kyungpook National
University, Daegu, Republic of Korea — kisupark@knu.ac.kr
