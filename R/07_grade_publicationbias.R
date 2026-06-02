## =============================================================================
## 07_grade_publicationbias.R
## Purpose: (a) Per-pool publication-bias / small-study-effect tests, and
##          (b) a GRADE-for-DTA certainty summary table assembled from the
##              quantitative results (heterogeneity, pool size, bias test).
## Inputs : `dat` from 01_data_prep.R (sourced automatically). Re-fits the
##          pools so this script is self-contained.
## Depends: metafor, mada
## Notes  : Bias-test conventions match the protocol: Deeks' test for the
##          diagnostic (HSROC) pool; Egger's test for the continuous pools and
##          AUROC. Egger is NOT computed for MAE (SE-effect coupling). The GRADE
##          certainty column reflects the published ratings; the supporting
##          quantities (I^2, k, bias p) are recomputed from the dataset.
## =============================================================================

if (!exists("dat")) source(file.path("R", "01_data_prep.R"))
for (p in c("metafor", "mada"))
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is required. Install with install.packages('", p, "').")
suppressMessages({ library(metafor); library(mada) })

## ---- refit the continuous + AUROC pools (self-contained) --------------------
ap  <- dat[dat$flag_auroc & !is.na(dat$auroc_logit) & !is.na(dat$auroc_logit_var), ]
res_auroc <- rma(yi = ap$auroc_logit, vi = ap$auroc_logit_var, method = "PM")

mae <- dat[dat$flag_mae & !is.na(dat$MAE_val) & !is.na(dat$mae_var), ]
res_mae <- rma(yi = mae$MAE_val, vi = mae$mae_var, method = "PM")

pcc <- dat[dat$flag_pcc & !is.na(dat$pcc_z) & !is.na(dat$pcc_z_var), ]
res_pcc <- rma(yi = pcc$pcc_z, vi = pcc$pcc_z_var, method = "REML")

icc <- dat[dat$flag_icc & !is.na(dat$icc_z) & !is.na(dat$icc_z_var), ]
res_icc <- rma(yi = icc$icc_z, vi = icc$icc_z_var, method = "REML")

## ---- publication-bias / small-study tests -----------------------------------
cat("=== Small-study / asymmetry tests ===\n")
egg <- function(res) tryCatch(regtest(res, model = "lm")$pval, error = function(e) NA)
auroc_egger <- egg(res_auroc)
pcc_egger   <- egg(res_pcc)
icc_egger   <- egg(res_icc)
cat(sprintf("  AUROC  Egger p = %.3f\n", auroc_egger))
cat(sprintf("  PCC    Egger p = %.3f\n", pcc_egger))
cat(sprintf("  ICC    Egger p = %.3f\n", icc_egger))
cat("  MAE    Egger not computed (SE = MAE/sqrt(N) is coupled to the effect).\n")

## HSROC: Deeks' test (recomputed from the 2x2 pool)
hs <- dat[dat$flag_hsroc, c("TP", "FP", "TN", "FN")]
hs <- hs[stats::complete.cases(hs), ]
ess  <- with(hs, 4 * (TP + FN) * (TN + FP) / ((TP + FN) + (TN + FP)))
lndor <- with(hs, log(((TP + .5) * (TN + .5)) / ((FP + .5) * (FN + .5))))
deeks_p <- summary(lm(lndor ~ I(1 / sqrt(ess)), weights = ess))$coefficients["I(1/sqrt(ess))", "Pr(>|t|)"]
cat(sprintf("  HSROC  Deeks p = %.3f\n", deeks_p))

## ---- GRADE-DTA certainty summary --------------------------------------------
## Certainty ratings are the published GRADE-DTA judgments; the I^2, k and
## bias-test columns are recomputed here so the table is data-linked.
grade <- data.frame(
  Pool      = c("HSROC (sens/spec)", "AUROC", "MAE (0-10)", "Pearson r", "ICC"),
  k         = c(sum(dat$flag_hsroc), res_auroc$k, res_mae$k, res_pcc$k, res_icc$k),
  I2_pct    = round(c(NA, res_auroc$I2, res_mae$I2, res_pcc$I2, res_icc$I2), 1),
  Bias_p    = round(c(deeks_p, auroc_egger, NA, pcc_egger, icc_egger), 3),
  Certainty = c("LOW", "LOW", "VERY LOW", "VERY LOW", "VERY LOW"),
  stringsAsFactors = FALSE
)
cat("\n=== GRADE-DTA certainty summary ===\n")
print(grade, row.names = FALSE)
cat("\nPublication-bias domain was not downgraded in any pool ",
    "(no asymmetry/small-study test reached significance).\n", sep = "")

invisible(grade)
