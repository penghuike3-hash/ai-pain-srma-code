# AI Pain Assessment SR-MA — Analysis Code

Reproducible R code for the systematic review and diagnostic test-accuracy meta-analysis of AI-based pain assessment in clinical patients (PROSPERO CRD420261354420).

## Quick start (one click)

The dataset is already included in `data/`. Run **`RUN_ME.R`** — any one of:
- **RStudio**: open `RUN_ME.R`, click **Source**.
- **Double-click** `RUN_ME.R` (Windows, if `.R` is associated with Rscript).
- **Console**: `source("RUN_ME.R")`

`RUN_ME.R` finds itself and the dataset automatically (no `setwd()` needed), installs any missing packages, runs every pool, prints each result next to the published value, and writes a full console log plus per-table CSVs to an `outputs/` folder beside this script.

It reproduces: pool sizes 11 / 62 / 31 / 22 (HSROC / AUROC / PCC / ICC; MAE is summarised descriptively, not pooled); HSROC summary sensitivity 0.882, specificity 0.867, AUC 0.929; pooled AUROC 0.887 [0.858, 0.911], I² 94.5%, 95% prediction interval 0.593–0.977; MAE median 0.50 (range 0.01–2.49), with 41 of 45 standardisable estimates below the MCID of 2; Pearson r 0.691; ICC 0.590; and all eight subgroup analyses (between-subgroup Q_M p ≥ 0.22).

## Repository contents

- `RUN_ME.R` — the reference implementation (self-contained; recommended). This is the script that produced every published value on the finalised dataset (N = 188).
- `R/01_data_prep.R … R/07_grade_publicationbias.R` — the same analysis split into modular, individually runnable steps (data prep, HSROC, AUROC, continuous pools, subgroups SG1–SG8, sensitivity/meta-regression, GRADE/publication bias). These mirror the methods in `RUN_ME.R`.
- `data/` — the extracted dataset (`Supplementary_Data_Extracted_Dataset_v6_final.xlsx`; `Studies` sheet, 188 rows × 99 columns). Also archived at the Zenodo data DOI below.
- `outputs/` — pre-computed reproduction summary (`reproduction_summary.txt`) and subgroup table (`SG1-SG8_results.csv`); the runner refreshes and extends this folder.
- `sessionInfo.txt` — the exact R session used for the published results.
- `LICENSE` — MIT.

## Method note (AUROC variance)

Per-study AUROC sampling variances are derived from each study's AUROC point estimate and sample size using the Hanley & McNeil (1982) formula (balanced-class N/2 split), mapped to the logit scale by the delta method. This is the standard ROC-meta-analysis approach and requires only AUROC and N — both present in the dataset — and reproduces the published pooled AUROC of 0.887 [0.858, 0.911] and I² of 94.5% exactly. Correlation and ICC pools use Fisher's z.

## Continuous outcomes (MAE)

The mean absolute error (MAE) is **not pooled**. Error metrics and pain scales are not comparable across studies, so MAE is summarised narratively (SWiM guidance): the median, range, and the proportion of standardisable estimates below the minimal clinically important difference (MCID = 2 NRS points).

## Environment

The published results were produced under **R 4.5.1** with **metafor 4.8-0**, **mada 0.5.12**, and **readxl 1.4.5** (full details in `sessionInfo.txt`). Packages auto-install on first run.

## Data and citation

- Dataset archive (Zenodo): DOI **10.5281/zenodo.20487779**
- Please cite the associated article: Ke P, Hu X, Wang L, Song Y, Xu G-H. *Diagnostic Performance of Artificial Intelligence Pain Assessment Systems in Clinical Patients: A Systematic Review and Meta-Analysis.* [Journal], [Year]. DOI **[article DOI]**.
