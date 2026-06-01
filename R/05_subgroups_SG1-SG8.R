## =============================================================================
## 05_subgroups_SG1-SG8.R
## Purpose: The eight pre-registered subgroup analyses (PROSPERO protocol).
##          For each axis: per-stratum pooled AUROC and the between-subgroup
##          moderator test (Q_M, df, p) from a mixed-effects meta-regression.
##            SG1 modality count        SG5 architecture
##            SG2 modality type         SG6 reference standard
##            SG3 clinical setting      SG7 testing strategy
##            SG4 task type             SG8 data-collection environment
## Inputs : `dat` from 01_data_prep.R (sourced automatically).
## Depends: metafor
## Notes  : Pooling is on the AUROC pool (logit scale), consistent with
##          03_auroc_pool.R. A stratum is pooled only if it has >= 3 studies
##          (protocol minimum). All values computed from the dataset at run time.
## =============================================================================

if (!exists("dat")) source(file.path("R", "01_data_prep.R"))
if (!requireNamespace("metafor", quietly = TRUE))
  stop("Package 'metafor' is required. Install with install.packages('metafor').")
suppressMessages(library(metafor))

## AUROC pool with usable logit value and variance.
ap <- dat[dat$flag_auroc & !is.na(dat$auroc_logit) & !is.na(dat$auroc_logit_var), ]
message(sprintf("AUROC pool for subgroup analyses: k = %d", nrow(ap)))

MIN_K <- 3   # protocol minimum studies per stratum

## Pool one subgroup axis: per-stratum estimates + moderator test.
run_subgroup <- function(df, factor_col, label) {
  df <- df[!is.na(df[[factor_col]]) & df[[factor_col]] != "", ]
  df[[factor_col]] <- factor(df[[factor_col]])
  ## keep only strata meeting the minimum-k rule
  keep <- names(which(table(df[[factor_col]]) >= MIN_K))
  df <- df[df[[factor_col]] %in% keep, ]
  df[[factor_col]] <- droplevels(df[[factor_col]])
  cat(sprintf("\n=== %s ===\n", label))
  if (nlevels(df[[factor_col]]) < 2) {
    cat("  < 2 strata meet the minimum-k rule; moderator test not estimable.\n")
    return(invisible(NULL))
  }
  ## per-stratum pooled AUROC (random-effects, REML), back-transformed
  for (lev in levels(df[[factor_col]])) {
    sub <- df[df[[factor_col]] == lev, ]
    r <- rma(yi = sub$auroc_logit, vi = sub$auroc_logit_var, method = "REML")
    cat(sprintf("  %-28s k=%2d  AUROC=%.3f [%.3f, %.3f]\n",
                lev, r$k, inv_logit(r$b[1]), inv_logit(r$ci.lb), inv_logit(r$ci.ub)))
  }
  ## moderator test: mixed-effects meta-regression with the factor as moderator
  mod <- rma(yi = df$auroc_logit, vi = df$auroc_logit_var,
             mods = ~ df[[factor_col]], method = "REML")
  cat(sprintf("  Q_M = %.2f, df = %d, p = %.3f\n",
              mod$QM, mod$p - 1, mod$QMp))
  invisible(mod)
}

run_subgroup(ap, "SG1_modality_count",      "SG1  Modality count")
run_subgroup(ap, "SG2_modality_type",       "SG2  Modality type")
run_subgroup(ap, "SG3_setting",             "SG3  Clinical setting")
run_subgroup(ap, "SG4_task_type",           "SG4  Task type")
run_subgroup(ap, "SG5_architecture",        "SG5  Architecture")
run_subgroup(ap, "SG6_reference_standard",  "SG6  Reference standard")
run_subgroup(ap, "SG7_testing_strategy",    "SG7  Testing strategy")
run_subgroup(ap, "SG8_environment",         "SG8  Data-collection environment")

invisible(TRUE)
