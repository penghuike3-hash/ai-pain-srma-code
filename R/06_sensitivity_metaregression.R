## =============================================================================
## 06_sensitivity_metaregression.R
## Purpose: The six pre-specified sensitivity analyses on the AUROC pool, plus a
##          dataset-silo moderator meta-regression.
##          Sensitivity analyses (re-pool AUROC after each restriction):
##            S1 exclude HIGH risk of bias (QUADAS-2)
##            S2 exclude N < 30
##            S3 external-tested studies only
##            S4 patient self-report reference standard only
##            S5 inter-rater reliability reported (strict) only
##            S6 exclude sole-public-benchmark (UNBC-only) studies
## Inputs : `dat` from 01_data_prep.R (sourced automatically).
## Depends: metafor
## Notes  : Pooling is on the AUROC pool (logit scale), consistent with
##          03_auroc_pool.R. All values computed from the dataset at run time.
## =============================================================================

if (!exists("dat")) source(file.path("R", "01_data_prep.R"))
if (!requireNamespace("metafor", quietly = TRUE))
  stop("Package 'metafor' is required. Install with install.packages('metafor').")
suppressMessages(library(metafor))

ap <- dat[dat$flag_auroc & !is.na(dat$auroc_logit) & !is.na(dat$auroc_logit_var), ]

pool_auroc <- function(df) {
  if (nrow(df) < 3) return(NULL)
  rma(yi = df$auroc_logit, vi = df$auroc_logit_var, method = "PM")
}
line <- function(label, df) {
  r <- pool_auroc(df)
  if (is.null(r)) { cat(sprintf("  %-34s k<3, not pooled\n", label)); return(invisible()) }
  cat(sprintf("  %-34s k=%2d  AUROC=%.3f [%.3f, %.3f]  I2=%.1f%%\n",
              label, r$k, inv_logit(r$b[1]), inv_logit(r$ci.lb),
              inv_logit(r$ci.ub), r$I2))
  invisible(r)
}

cat("=== AUROC pool: main and sensitivity analyses ===\n")
line("Main pool", ap)

## S1 exclude HIGH risk of bias
line("S1 exclude HIGH RoB", ap[ap$overall_RoB != "HIGH", ])

## S2 exclude N < 30 (uses the pre-computed SA2 flag where present, else N_val)
if ("in_AUROC_pool_SA2_N30" %in% names(dat)) {
  s2 <- dat[as_flag(dat$in_AUROC_pool_SA2_N30) &
              !is.na(dat$auroc_logit) & !is.na(dat$auroc_logit_var), ]
} else {
  s2 <- ap[!is.na(ap$N_val) & ap$N_val >= 30, ]
}
line("S2 exclude N<30", s2)

## S3 externally/independently tested (SG7 stratum independent OR external).
line("S3 externally/independently tested",
     ap[grepl("independent|external", ap$SG7_testing_strategy, ignore.case = TRUE), ])

## S4 patient self-report reference standard only
line("S4 self-report reference only",
     ap[grepl("self-report", tolower(ap$SG6_reference_standard)), ])

## S5 genuine reference inter-rater kappa >= 0.70. The kappa is parsed from the
## free-text reliability field; studies reporting model-vs-human agreement or a
## non-reference statistic are excluded (EXCL5), matching the manuscript SA5.
extract_kappa <- function(s) {
  s <- ifelse(is.na(s), "", as.character(s))
  pat <- "(?i).*?(?:cohen|fleiss|weighted)?\\s*(?:kappa|\u03ba)\\s*[()0-9,. ]*?=\\s*([01]?\\.[0-9]+).*"
  has <- grepl(pat, s, perl = TRUE)
  out <- rep(NA_real_, length(s)); out[has] <- as.numeric(sub(pat, "\\1", s[has], perl = TRUE)); out
}
.EXCL5 <- "Xu_2020|Aydin_2023|Wibowo_2023|Alkan_2025|Oznaneci_2025"
ap$kappa      <- extract_kappa(ap$inter_rater_reliability)
ap$is_excl5   <- grepl(.EXCL5, ap$study_id, ignore.case = TRUE)
ap$modelhuman <- grepl("model[ -]?(vs|versus|human|self)", ap$inter_rater_reliability, ignore.case = TRUE)
line("S5 IRR kappa>=0.70",
     ap[!is.na(ap$kappa) & ap$kappa >= 0.70 & !ap$is_excl5 & !ap$modelhuman, ])

## S6 exclude sole-public-benchmark (UNBC-only) studies (Cat6 flag).
line("S6 exclude sole-public-benchmark",
     ap[!as_flag(ap$Cat6_SoleRelianceOnPublicBenchmark), ])

## ---- silo-moderator meta-regression -----------------------------------------
cat("\n=== Dataset-silo moderator meta-regression (AUROC pool) ===\n")
ms <- ap[!is.na(ap$silo_4macro) & ap$silo_4macro != "", ]
ms$silo_4macro <- factor(ms$silo_4macro)
keep <- names(which(table(ms$silo_4macro) >= 3))
ms <- ms[ms$silo_4macro %in% keep, ]; ms$silo_4macro <- droplevels(ms$silo_4macro)
for (lev in levels(ms$silo_4macro)) {
  sub <- ms[ms$silo_4macro == lev, ]
  r <- rma(yi = sub$auroc_logit, vi = sub$auroc_logit_var, method = "REML")
  cat(sprintf("  %-16s k=%2d  AUROC=%.3f [%.3f, %.3f]\n",
              lev, r$k, inv_logit(r$b[1]), inv_logit(r$ci.lb), inv_logit(r$ci.ub)))
}
mr <- rma(yi = ms$auroc_logit, vi = ms$auroc_logit_var,
          mods = ~ ms$silo_4macro, method = "REML")
cat(sprintf("  Q_M = %.2f, df = %d, p = %.3f\n", mr$QM, mr$p - 1, mr$QMp))

invisible(TRUE)
