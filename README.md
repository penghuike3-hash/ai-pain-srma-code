# AI Pain Assessment SR-MA — Analysis Code

Reproducible R code for the systematic review and diagnostic test-accuracy meta-analysis of AI-based pain assessment in clinical patients (PROSPERO CRD420261354420).

## Quick start (one click)

1. Make sure the dataset is in the `data/` folder:
   ```
   data/Supplementary_Data_Extracted_Dataset.xlsx
   ```
   (Download it from the data DOI below and place it here.)
2. Run **`RUN_ME.R`** — any one of:
   - **RStudio**: open `RUN_ME.R`, click **Source**.
   - **Double-click** `RUN_ME.R` (Windows, if `.R` is associated with Rscript).
   - **Console**: `source("RUN_ME.R")`

`RUN_ME.R` finds itself and the dataset automatically (no `setwd()` needed), installs any missing packages, runs every pool, prints each result next to the published value, and writes `outputs/reproduction_summary.txt` and `outputs/SG1-SG8_results.csv`.

It reproduces: pool sizes 14/69/43/35/25; HSROC summary sensitivity 0.904, specificity 0.912, AUC 0.951; pooled AUROC 0.904 [0.879, 0.924], I² 94.5%, 95% prediction interval 0.618–0.982; MAE 0.79; Pearson r 0.688; ICC 0.638; and all eight subgroup analyses (Q_M p ≥ 0.25).

## Repository contents

- `RUN_ME.R` — one-click runner (self-contained; recommended).
- `R/01_data_prep.R … R/07_grade_publicationbias.R` — the same analysis split into modular, individually runnable steps (data prep, HSROC, AUROC, continuous pools, subgroups SG1–SG8, sensitivity/meta-regression, GRADE/publication bias).
- `data/` — place the extracted dataset here (download from the Zenodo data DOI; `Studies` sheet, 189 rows).
- `outputs/` — created at run time; holds the printed results.
- `sessionInfo.txt` — the exact R session used for the published results.
- `LICENSE` — MIT.

## Method note (AUROC variance)

Per-study AUROC sampling variances are derived from each study's AUROC point estimate and sample size using the Hanley & McNeil (1982) formula (balanced-class N/2 split), mapped to the logit scale by the delta method. This is the standard ROC-meta-analysis approach and requires only AUROC and N — both present in the dataset — and reproduces the published pooled AUROC and I² exactly. Correlation and ICC pools use Fisher's z; MAE uses SE = MAE/√N.

## Environment

The published results were produced under **R 4.5.1** with **metafor 4.8-0**, **mada 0.5.12**, and **readxl 1.4.5** (full details in `sessionInfo.txt`). Packages auto-install on first run.

## Data and citation

- Dataset archive (Zenodo): DOI **10.5281/zenodo.20487779**
- Please cite the associated article: Ke P, Hu X, Wang L, Song Y, Xu G-H. *Diagnostic Performance of Artificial Intelligence Pain Assessment Systems in Clinical Patients: A Systematic Review and Meta-Analysis.* [Journal], [Year]. DOI **[article DOI]**.
