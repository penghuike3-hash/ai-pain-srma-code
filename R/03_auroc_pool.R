## =============================================================================
## 03_auroc_pool.R
## Purpose: Random-effects meta-analysis of AUROC on the logit scale:
##          inverse-variance pooling (Paule-Mandel tau^2), I^2, 95% prediction
##          interval, and Egger's test for small-study effects.
## Inputs : `dat` from 01_data_prep.R (sourced automatically).
## Depends: metafor
## Method : Each study's AUROC sampling variance is computed from its AUROC point
##          estimate and sample size via the Hanley & McNeil (1982) formula
##          (set up in 01_data_prep.R) and mapped to the logit scale by the delta
##          method. This is the standard ROC-meta-analysis approach and requires
##          only AUROC and N -- both present in the dataset. It reproduces the
##          published pooled AUROC of 0.887 [0.858, 0.911] exactly.
## =============================================================================

if (!exists("dat")) source(file.path("R", "01_data_prep.R"))
if (!requireNamespace("metafor", quietly = TRUE))
  stop("Package 'metafor' is required. Install with install.packages('metafor').")
suppressMessages(library(metafor))

## ---- assemble the AUROC pool ------------------------------------------------
## Pool membership is defined by the dataset flag. One study reports no usable N
## and so cannot be weighted; it is in the pool but excluded from the model
## (reported as counted-not-weighted: 63 in pool, 62 weighted).
pool_all <- dat[dat$flag_auroc, ]
ap <- pool_all[!is.na(pool_all$auroc_logit) & !is.na(pool_all$auroc_logit_var), ]
message(sprintf("AUROC pool: %d studies in pool, %d weighted (with usable N).",
                nrow(pool_all), nrow(ap)))

## ---- inverse-variance random-effects model (Paule-Mandel) -------------------
res <- rma(yi = ap$auroc_logit, vi = ap$auroc_logit_var, method = "PM")

pi  <- predict(res)
est <- inv_logit(res$b[1]); ci_lo <- inv_logit(res$ci.lb); ci_hi <- inv_logit(res$ci.ub)
pi_lo <- inv_logit(pi$pi.lb); pi_hi <- inv_logit(pi$pi.ub)

message(sprintf("Pooled AUROC = %.3f [95%% CI %.3f, %.3f]", est, ci_lo, ci_hi))
message(sprintf("I^2 = %.1f%%   tau^2 = %.4f (logit scale)", res$I2, res$tau2))
message(sprintf("95%% prediction interval: [%.3f, %.3f]", pi_lo, pi_hi))

## ---- Egger's test for small-study effects -----------------------------------
egg <- regtest(res, model = "lm")
message(sprintf("Egger's test: z/t p = %.3f", egg$pval))

invisible(list(model = res, auroc = est, ci = c(ci_lo, ci_hi),
               pi = c(pi_lo, pi_hi), I2 = res$I2, egger_p = egg$pval))
