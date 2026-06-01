## =============================================================================
## 04_continuous_pools.R
## Purpose: Random-effects pooling of the continuous outcomes:
##            - Mean absolute error (MAE, standardised 0-10), SE = MAE / sqrt(N)
##            - Pearson correlation r (Fisher's z)
##            - Intraclass correlation ICC (Fisher's z)
##          Reports pooled estimate, 95% CI, I^2, and 95% prediction interval.
## Inputs : `dat` from 01_data_prep.R (sourced automatically).
## Depends: metafor
## Notes  : Paule-Mandel tau^2 estimator (REML for the z-pools). For MAE, Egger's
##          test is NOT computed because SE = MAE / sqrt(N) is mechanically
##          coupled to the effect size, which invalidates the regression test.
##          All values computed from the dataset at run time.
## =============================================================================

if (!exists("dat")) source(file.path("R", "01_data_prep.R"))
if (!requireNamespace("metafor", quietly = TRUE))
  stop("Package 'metafor' is required. Install with install.packages('metafor').")
suppressMessages(library(metafor))

report <- function(label, res, back = identity) {
  pi <- predict(res)
  est <- back(res$b[1]); ci_lo <- back(res$ci.lb); ci_hi <- back(res$ci.ub)
  pil <- back(pi$pi.lb); pih <- back(pi$pi.ub)
  message(sprintf(
    "%-14s k=%d  est=%.3f [%.3f, %.3f]  I2=%.1f%%  PI [%.3f, %.3f]",
    label, res$k, est, ci_lo, ci_hi, res$I2, pil, pih))
  invisible(res)
}

## ---- MAE (0-10): SE = MAE / sqrt(N); Paule-Mandel ---------------------------
mae <- dat[dat$flag_mae & !is.na(dat$MAE_val) & !is.na(dat$mae_var), ]
res_mae <- rma(yi = mae$MAE_val, vi = mae$mae_var, method = "PM")
report("MAE (0-10)", res_mae)
n_below_mcid <- sum(mae$MAE_val < 2, na.rm = TRUE)   # MCID = 2 NRS points
message(sprintf("  MAE < MCID(2) in %d of %d studies.", n_below_mcid, nrow(mae)))
message("  Egger's test not interpreted for MAE (SE-effect coupling).")

## ---- Pearson r: Fisher's z, REML --------------------------------------------
pcc <- dat[dat$flag_pcc & !is.na(dat$pcc_z) & !is.na(dat$pcc_z_var), ]
res_pcc <- rma(yi = pcc$pcc_z, vi = pcc$pcc_z_var, method = "REML")
report("Pearson r", res_pcc, back = inv_fisher)
egg_pcc <- regtest(res_pcc, model = "lm")
message(sprintf("  Egger's test: p = %.3f", egg_pcc$pval))

## ---- ICC: Fisher's z, REML --------------------------------------------------
icc <- dat[dat$flag_icc & !is.na(dat$icc_z) & !is.na(dat$icc_z_var), ]
res_icc <- rma(yi = icc$icc_z, vi = icc$icc_z_var, method = "REML")
report("ICC", res_icc, back = inv_fisher)
egg_icc <- regtest(res_icc, model = "lm")
message(sprintf("  Egger's test: p = %.3f", egg_icc$pval))

invisible(list(mae = res_mae, pcc = res_pcc, icc = res_icc))
