## =============================================================================
## 02_hsroc.R
## Purpose: Bivariate / HSROC diagnostic meta-analysis (Reitsma model) for the
##          studies with an extractable 2x2 table: summary sensitivity,
##          specificity, and AUC, plus Deeks' funnel-plot asymmetry test.
## Inputs : `dat` from 01_data_prep.R (sourced automatically).
## Outputs: prints the bivariate summary, the HSROC AUC, and Deeks' test.
## Depends: mada, metafor
## Notes  : Two zero-cell studies receive a 0.5 continuity correction (mada
##          default). All values are computed from the dataset at run time.
## =============================================================================

if (!exists("dat")) source(file.path("R", "01_data_prep.R"))
for (p in c("mada", "metafor"))
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is required. Install with install.packages('", p, "').")
suppressMessages({ library(mada); library(metafor) })

## ---- assemble the 2x2 data for the HSROC pool -------------------------------
hs <- dat[dat$flag_hsroc, c("study_id", "TP", "FP", "TN", "FN")]
hs <- hs[stats::complete.cases(hs[, c("TP", "FP", "TN", "FN")]), ]
message(sprintf("HSROC pool: %d studies with a complete 2x2 table.", nrow(hs)))

## mada applies a 0.5 continuity correction to studies containing a zero cell.
zero_cell <- with(hs, TP == 0 | FP == 0 | TN == 0 | FN == 0)
message(sprintf("Continuity correction (0.5) applies to %d study(ies): %s",
                sum(zero_cell), paste(hs$study_id[zero_cell], collapse = ", ")))

## ---- bivariate (Reitsma) model ----------------------------------------------
fit <- reitsma(hs, formula = cbind(tsens, tfpr) ~ 1,
               method = "ml", correction = 0.5)
sm  <- summary(fit)

## Summary sensitivity and specificity with 95% CI (back-transformed).
print(sm$coefficients[c("sensitivity", "specificity"), , drop = FALSE])

## ---- HSROC area under the curve ---------------------------------------------
auc_hsroc <- mada::AUC(fit)
message(sprintf("HSROC AUC = %.4f", auc_hsroc$AUC))

## ---- between-study heterogeneity (Zhou-Dendukuri style I^2 via mada) --------
## mada reports the correlation and variance components; print the full summary
## so the heterogeneity components are available alongside the point estimates.
print(sm)

## ---- Deeks' funnel-plot asymmetry test --------------------------------------
## Deeks' test regresses the log diagnostic odds ratio on 1/sqrt(ESS).
ess  <- with(hs, 4 * (TP + FN) * (TN + FP) / ((TP + FN) + (TN + FP)))   # effective sample size
lndor <- with(hs, log(((TP + 0.5) * (TN + 0.5)) / ((FP + 0.5) * (FN + 0.5))))
deeks <- summary(lm(lndor ~ I(1 / sqrt(ess))))
deeks_p <- coef(deeks)["I(1/sqrt(ess))", "Pr(>|t|)"]
message(sprintf("Deeks' test for funnel-plot asymmetry: p = %.3f", deeks_p))

invisible(list(fit = fit, auc = auc_hsroc$AUC, deeks_p = deeks_p))
