# Input data

Study-level data extracted from the 41 studies included in the systematic review. Licence: CC BY 4.0
(see `../LICENSE-DATA.md`). The same values are provided as Supplementary Data 1 of the article.

| File | Content |
|------|---------|
| `CSF_iNPH_vs_HC.xlsx` | One sheet per biomarker, iNPH vs healthy controls: group sizes, means and SDs (`En, Em, Esd, Cn, Cm, Csd`), unit, assay, CSF sampling site (`CSF_Type`), mean age of the iNPH group. Sheet `List nonprimary`: rows used only in sensitivity analyses or not pooled, with the reason and source. Sheet `List decisions`: data decisions made during verification. |
| `CSF_iNPH_vs_AD.xlsx` | The same structure for iNPH vs Alzheimer's disease. |
| `diagnostic_performance.xlsx` | Two-by-two counts (TP, FP, FN, TN; iNPH = positive condition) for the diagnostic accuracy analysis, with the count source and QUADAS-2 reference-standard rating (`Notes`) and alternative denominators (`Alt_counts`). |
| `conversion_status.csv` | For every primary contrast: whether the mean and SD were reported directly or converted (Luo/Wan from median and IQR; two cohorts combined), and the source location in the publication. |
| `study_characteristics.xlsx` | Coding of iNPH criteria and certainty, shunt-responder restriction, CSF sampling setting, recruitment design, control type, AD criteria, use of CSF in the AD diagnosis, exclusion of AD biomarker-positive iNPH, and assays. |

## Sheet and column rules used by the loader (`R/load_data.R`)

1. Every sheet whose name does not start with `List` is read; the sheet name is the biomarker. Ratios are
   written with `by` (for example `AB42 by AB40`), because `/` is not allowed in Excel sheet names.
2. `Study` is column 1; `En, Em, Esd, Cn, Cm, Csd` are columns 9 to 14; the unit is column 15 and
   `CSF_Type` (`Lumbar`, `Ventricular` or `Brain biopsy`) column 17. Columns must not be moved.
