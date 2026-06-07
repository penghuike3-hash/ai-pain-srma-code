###############################################################################
# RUN_ME.R  —  AI Pain Assessment SR/MA : re-analysis on the finalised dataset
# -----------------------------------------------------------------------------
# Input  : Supplementary_Data_Extracted_Dataset_v6_final.xlsx  ("Studies" sheet)
# N      : 188 studies
# Pools  : AUROC (k=62) | HSROC (k=11) | PCC (k=31) | ICC (k=22)
#          MAE -> descriptive only (NOT pooled)
#
# HOW TO RUN
#   1. install.packages(c("metafor","mada","readxl"))     # one time
#   2. put this file next to the .xlsx (or edit INPUT_XLSX below)
#   3. source("RUN_ME.R")    # all results print to console + land in ./results/
#
# WHY R: the HSROC bivariate (Reitsma) model needs mada and has no Python
# equivalent; everything else is reproduced with the same methods used originally.
#
# SELF-CHECK -- AUTHORITATIVE values from this script's own R run (metafor PM /
# mada ML). Confirmed by independent Paule-Mandel root-finding. If a POINT
# ESTIMATE differs by > 0.01, a column was read wrong -- stop and check the
# flag/value columns; do NOT hand-edit numbers.
#     AUROC  k=62   0.887  [0.858, 0.911]   PI[0.593, 0.977]   I2 = 94.5%   (PM, tau^2=0.722)
#     PCC    k=31   0.691  [0.599, 0.765]                      I2 = 92.2%   (REML)
#     ICC    k=22   0.590  [0.492, 0.672]                      I2 = 86.4%   (REML)
#     HSROC  k=11   sens 0.882 [0.833,0.917]  spec 0.867 [0.790,0.919]  AUC 0.929   (Reitsma ML)
#     Subgroups: SG1-SG8, per-stratum REML, min k=3 (see section 7)
#
# METHODS match the deposited modular pipeline exactly: AUROC=PM (tau^2=0.722;
# an earlier hand-rolled PM root of 0.34 gave a wrong 0.878/89% -- ignore it),
# PCC/ICC=REML (04_continuous_pools.R), HSROC=reitsma ML + 0.5 correction
# (02_hsroc.R), subgroups=REML min-k 3 (05_subgroups_SG1-SG8.R).
###############################################################################

## ---- 0. CONFIG --------------------------------------------------------------
## --- self-locate the dataset so you can just Source / double-click (no setwd) -
.this_file <- function() {
  ca <- commandArgs(FALSE); m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]), mustWork = FALSE))
  for (i in rev(seq_len(sys.nframe()))) { of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, mustWork = FALSE)) }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(normalizePath(p, mustWork = FALSE)) }
  NA_character_
}
.find_input <- function(fname) {
  home <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else path.expand("~")
  sp <- .this_file(); sd <- if (!is.na(sp)) dirname(sp) else NULL
  cand <- unique(c(sd, file.path(sd, "data"), getwd(), file.path(getwd(), "data"),
            file.path(home, c("Downloads","Desktop","Documents",
                              "OneDrive/Downloads","OneDrive/Desktop","OneDrive/Documents"))))
  for (d in cand) { if (is.null(d) || is.na(d) || !nzchar(d) || !dir.exists(d)) next
    hit <- file.path(d, fname); if (file.exists(hit)) return(normalizePath(hit))
    g <- list.files(d, pattern = fname, full.names = TRUE)
    if (length(g)) return(normalizePath(g[1])) }
  if (interactive()) { message("Couldn't auto-find ", fname, " - please pick it in the dialog.")
                       return(file.choose()) }
  stop("Could not find ", fname, " in:\n  ", paste(cand[!is.na(cand)], collapse = "\n  "))
}
INPUT_XLSX <- .find_input("Supplementary_Data_Extracted_Dataset_v6_final.xlsx")
.SCRIPT_DIR <- tryCatch({ d <- dirname(.this_file()); if (is.na(d) || !nzchar(d)) getwd() else d },
                        error = function(e) getwd())
cat("Using dataset:\n  ", INPUT_XLSX, "\n", sep = "")
OUT_DIR        <- file.path(.SCRIPT_DIR, "outputs")  # outputs land next to THIS script (repo root), not inside data/
AUROC_METHOD   <- "PM"     # tau^2 estimator for AUROC + its SAs/subgroups (protocol = Paule-Mandel; "DL" reproduces the old draft's I2)
PCC_ICC_METHOD <- "REML"   # Fisher-z pools, matches original 04_continuous_pools.R (REML)
HSROC_METHOD   <- "ml"     # mada::reitsma estimation method (ML, per protocol)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
.need <- c("metafor", "mada", "readxl")
.miss <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss)) {
  message("Installing missing packages: ", paste(.miss, collapse = ", "), " ...")
  try(install.packages(.miss, repos = "https://cloud.r-project.org"))
  .miss <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
}
if (length(.miss)) stop("Still missing: ", paste(.miss, collapse = ", "),
                        "  -> install manually, then re-run.")
suppressMessages({ library(metafor); library(mada); library(readxl) })

## --- tee ALL console output into ONE file you send back to me ----------------
while (sink.number() > 0) sink()
.LOG <- file.path(OUT_DIR, "full_reproduction_log.txt")
sink(.LOG, split = TRUE)
cat("AI-Pain SR/MA -- full numbers run on v6_final  (", format(Sys.time()), ")\n", sep = "")

## ---- 1. LOAD + SHARED HELPERS ----------------------------------------------
S <- as.data.frame(read_excel(INPUT_XLSX, sheet = "Studies"))
cat(sprintf("Loaded 'Studies': N = %d rows\n", nrow(S)))     # expect 188

flag <- function(x) tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y")
n2   <- function(x) suppressWarnings(as.numeric(as.character(x)))

# robust prediction-interval extractor (metafor renamed cr.* -> pi.* across versions)
get_pi <- function(pp) {
  lb <- if (!is.null(pp$pi.lb)) pp$pi.lb else pp$cr.lb
  ub <- if (!is.null(pp$pi.ub)) pp$pi.ub else pp$cr.ub
  c(lb = lb, ub = ub)
}

cat("\nPool sizes (flagged):\n")
print(c(AUROC = sum(flag(S$in_AUROC_pool_main)),   # 63 flagged; 62 analysable (see below)
        HSROC = sum(flag(S$in_HSROC_pool)),         # 11
        PCC   = sum(flag(S$in_PCC_pool_main)),      # 31
        ICC   = sum(flag(S$in_ICC_pool_main)),      # 22
        MAE   = sum(flag(S$in_MAE_pool_main))))     # 0

## ============================================================================
## 2. AUROC POOL — logit-IV, Hanley & McNeil variance (prevalence 0.5),
##    Paule-Mandel tau^2.   k = 62
## ============================================================================
A <- S[flag(S$in_AUROC_pool_main), ]
A$auc <- n2(A$AUROC_numeric)
A$N   <- n2(A$N_numeric)

drop_auc <- A[is.na(A$auc) | is.na(A$N), c("study_id", "first_author", "AUROC_numeric")]
A <- A[!is.na(A$auc) & !is.na(A$N) & A$auc > 0 & A$auc < 1, ]
cat(sprintf("\nAUROC: flagged = %d | analysable (non-missing N) = %d | excluded = %d {%s}\n",
            sum(flag(S$in_AUROC_pool_main)), nrow(A), nrow(drop_auc),
            paste(drop_auc$study_id, collapse = ", ")))   # excluded: Cascella_2023 (no usable N)

hm_var <- function(auc, N) {                 # Hanley & McNeil (1982), n_pos = n_neg = N/2
  npos <- N / 2; nneg <- N / 2
  Q1 <- auc / (2 - auc); Q2 <- 2 * auc^2 / (1 + auc)
  (auc * (1 - auc) + (npos - 1) * (Q1 - auc^2) + (nneg - 1) * (Q2 - auc^2)) / (npos * nneg)
}
A$v_auc <- mapply(hm_var, A$auc, A$N)
A$yi <- qlogis(A$auc)                                   # logit scale
A$vi <- A$v_auc / (A$auc * (1 - A$auc))^2               # delta-method var on logit scale

m_auc  <- rma(yi = yi, vi = vi, data = A, method = AUROC_METHOD)
pp_auc <- predict(m_auc, transf = transf.ilogit)
pi_auc <- get_pi(pp_auc)
auroc_res <- data.frame(
  outcome = "AUROC", k = m_auc$k,
  est   = transf.ilogit(m_auc$b)[1],
  ci.lb = transf.ilogit(m_auc$ci.lb), ci.ub = transf.ilogit(m_auc$ci.ub),
  pi.lb = pi_auc["lb"], pi.ub = pi_auc["ub"],
  I2 = m_auc$I2, tau2 = m_auc$tau2)
cat("\n--- AUROC pooled ---\n"); print(auroc_res, row.names = FALSE, digits = 3)

## ============================================================================
## 3. PCC and ICC POOLS — Fisher-z random effects
## ============================================================================
fisher_pool <- function(df, valcol) {
  df$r <- n2(df[[valcol]]); df$N <- n2(df$N_numeric)
  df <- df[!is.na(df$r) & !is.na(df$N) & df$r > -1 & df$r < 1 & df$N > 3, ]
  es <- escalc(measure = "ZCOR", ri = r, ni = N, data = df)
  m  <- rma(yi, vi, data = es, method = PCC_ICC_METHOD)
  pp <- predict(m, transf = transf.ztor); pii <- get_pi(pp)
  data.frame(k = m$k,
             est   = transf.ztor(m$b)[1],
             ci.lb = transf.ztor(m$ci.lb), ci.ub = transf.ztor(m$ci.ub),
             pi.lb = pii["lb"], pi.ub = pii["ub"],
             I2 = m$I2, tau2 = m$tau2)
}
pcc_res <- fisher_pool(S[flag(S$in_PCC_pool_main), ], "pearson_r_numeric")
icc_res <- fisher_pool(S[flag(S$in_ICC_pool_main), ], "ICC_numeric")
cat("\n--- PCC pooled ---\n"); print(pcc_res, row.names = FALSE, digits = 3)
cat("\n--- ICC pooled ---\n"); print(icc_res, row.names = FALSE, digits = 3)

# Egger's test for PCC and ICC (matches original 04_continuous_pools.R)
egger_z <- function(df, valcol) {
  df$r <- n2(df[[valcol]]); df$N <- n2(df$N_numeric)
  df <- df[!is.na(df$r) & !is.na(df$N) & df$r > -1 & df$r < 1 & df$N > 3, ]
  es <- escalc(measure = "ZCOR", ri = r, ni = N, data = df)
  regtest(rma(yi, vi, data = es, method = PCC_ICC_METHOD), model = "lm")$pval
}
cat(sprintf("PCC Egger p = %.3f | ICC Egger p = %.3f\n",
            egger_z(S[flag(S$in_PCC_pool_main), ], "pearson_r_numeric"),
            egger_z(S[flag(S$in_ICC_pool_main), ], "ICC_numeric")))

## ============================================================================
## 4. HSROC POOL — bivariate (Reitsma) model, mada, method = ML.   k = 11
##    *** R-only: do not attempt to reproduce this in Python. ***
## ============================================================================
H <- S[flag(S$in_HSROC_pool), ]
H$TP <- n2(H$TP_num); H$FP <- n2(H$FP_num); H$TN <- n2(H$TN_num); H$FN <- n2(H$FN_num)
H <- H[complete.cases(H[, c("TP", "FP", "TN", "FN")]), ]
cat(sprintf("\nHSROC: k = %d\n", nrow(H)))

fit  <- reitsma(H, formula = cbind(tsens, tfpr) ~ 1, method = HSROC_METHOD, correction = 0.5)
sm   <- summary(fit)
cat("\n--- HSROC Reitsma summary (sensitivity / false-positive rate / AUC) ---\n")
print(sm)
hsroc_auc <- tryCatch(AUC(fit)$AUC, error = function(e) NA_real_)
cat(sprintf("HSROC summary AUC = %.3f\n", hsroc_auc))
writeLines(capture.output(print(sm), sprintf("AUC = %.4f", hsroc_auc)),
           file.path(OUT_DIR, "HSROC_reitsma_summary.txt"))
obs_sens <- H$TP / (H$TP + H$FN)         # observed points for the SROC plot
obs_fpr  <- H$FP / (H$FP + H$TN)

## ============================================================================
## 5. MAE — DESCRIPTIVE ONLY (Supplementary Table S4; no meta-analysis)
## ============================================================================
keepM <- !is.na(n2(S$MAE)) | !is.na(n2(S$MSE)) | !is.na(n2(S$MAE_standardized_0_10))
Mtab  <- S[keepM, c("study_id", "first_author", "year",
                    "MAE", "MSE", "MAE_scale", "MAE_standardized_0_10", "primary_task")]
write.csv(Mtab, file.path(OUT_DIR, "S4_MAE_descriptive.csv"), row.names = FALSE, na = "")
cat(sprintf("\nMAE descriptive: %d studies reporting MAE/MSE -> S4_MAE_descriptive.csv\n",
            nrow(Mtab)))

## ============================================================================
## 6. SENSITIVITY ANALYSES (AUROC primary outcome) -> Table S8
##    Masks taken verbatim from the deposited sensitivity_analyses_S8.R
##    (Main + SA1-SA6 + two Fig-cross-check rows). A already carries yi/vi, so
##    each row is just a re-pool of a subset of the analysable pool A.
##    NB: the deposited script's target comment (0.904, k=68) was the PRE-clean
##    pool; here the same masks run on the finalised pool (k=62) -> new numbers.
## ============================================================================
pool_yi <- function(df) {
  if (nrow(df) < 2)
    return(data.frame(k = nrow(df), est = NA, ci.lb = NA, ci.ub = NA, I2 = NA))
  m <- rma(yi, vi, data = df, method = AUROC_METHOD)
  data.frame(k = m$k, est = transf.ilogit(m$b)[1],
             ci.lb = transf.ilogit(m$ci.lb), ci.ub = transf.ilogit(m$ci.ub), I2 = m$I2)
}

# --- SA5: extract a GENUINE reference inter-rater kappa from free text -------
# matches Cohen/Fleiss/weighted/Greek kappa with a numeric value; ignores ICC,
# Pearson, percentages. Then drops (a) studies whose reported statistic is not a
# reference-standard kappa (EXCL5), and (b) model-vs-human agreement values.
extract_kappa <- function(s) {
  s   <- ifelse(is.na(s), "", as.character(s))
  pat <- "(?i).*?(?:cohen|fleiss|weighted)?\\s*(?:kappa|\u03ba)\\s*[()0-9,. ]*?=\\s*([01]?\\.[0-9]+).*"
  has <- grepl(pat, s, perl = TRUE)
  out <- rep(NA_real_, length(s)); out[has] <- as.numeric(sub(pat, "\\1", s[has], perl = TRUE)); out
}
EXCL5 <- "Xu_2020|Aydin_2023|Wibowo_2023|Alkan_2025|Oznaneci_2025"
A$kappa      <- extract_kappa(A$inter_rater_reliability)
A$is_excl5   <- grepl(EXCL5, A$study_id, ignore.case = TRUE)
A$modelhuman <- grepl("model[ -]?(vs|versus|human|self)", A$inter_rater_reliability, ignore.case = TRUE)

SA_masks <- list(
  "Main (all weighted)"          = rep(TRUE, nrow(A)),
  "SA1 exclude HIGH RoB"         = toupper(A$overall_RoB) != "HIGH",                                   # keeps LOW+UNCLEAR
  "SA2 N>=30"                    = A$N >= 30,
  "SA3 externally tested"        = grepl("independent|external", A$SG7_testing_stratum, ignore.case = TRUE),
  "SA4 self-report reference"    = grepl("self.?report", A$ref_std_macro, ignore.case = TRUE),
  "SA5 IRR kappa>=0.70"          = !is.na(A$kappa) & A$kappa >= 0.70 & !A$is_excl5 & !A$modelhuman,
  "SA6 exclude UNBC/public-only" = A$silo_4macro != "UNBC",
  "[chk] LOW RoB only (Fig6)"    = toupper(A$overall_RoB) == "LOW",
  "[chk] external only"          = grepl("external", A$SG7_testing_stratum, ignore.case = TRUE))

main_est <- pool_yi(A)$est
SA <- do.call(rbind, lapply(names(SA_masks), function(nm) {
  v <- pool_yi(A[SA_masks[[nm]], , drop = FALSE])
  data.frame(Analysis = nm, k = v$k, AUROC = round(v$est, 3),
             CI = ifelse(is.na(v$ci.lb), NA, sprintf("%.3f-%.3f", v$ci.lb, v$ci.ub)),
             I2 = round(v$I2, 1), Delta = round(v$est - main_est, 3), row.names = NULL)
}))
cat("\n--- Table S8: AUROC sensitivity analyses (Main + SA1-SA6 + 2 checks) ---\n")
print(SA, row.names = FALSE)
write.csv(SA, file.path(OUT_DIR, "S8_sensitivity_AUROC.csv"), row.names = FALSE)
cat("\nSA5 studies kept (genuine reference kappa >= 0.70):\n")
print(A[SA_masks[["SA5 IRR kappa>=0.70"]], c("study_id", "kappa", "inter_rater_reliability")], row.names = FALSE)

# Influence diagnostic (NOT one of the 6 registered SAs): leave-one-out on the pool
loo <- leave1out(m_auc)
loo_tab <- data.frame(omitted = A$study_id, est = transf.ilogit(loo$estimate), I2 = loo$I2)
write.csv(loo_tab, file.path(OUT_DIR, "influence_leave_one_out.csv"), row.names = FALSE)
cat(sprintf("\nLeave-one-out (influence, not a registered SA): AUROC %.3f - %.3f over %d omissions\n",
            min(loo_tab$est), max(loo_tab$est), nrow(loo_tab)))

## ============================================================================
## 7. SUBGROUPS SG1-SG8 (Table 4 / S7) -- derivations & method match the
##    deposited 01_data_prep.R + 05_subgroups_SG1-SG8.R EXACTLY:
##    per-stratum pooling = REML, minimum k = 3 per stratum, REML moderator test.
## ============================================================================
# --- SG factor derivations (verbatim logic from 01_data_prep.R) -------------
classify_modality <- function(ml) {
  s <- tolower(trimws(ml)); if (is.na(s) || s == "") return(NA_character_)
  has_face <- grepl("fac", s); n_types <- length(strsplit(s, ",")[[1]])
  has_physio <- grepl("physiolog|eeg|ecg|bioimped|skin", s)
  if (n_types >= 2 || (has_face && has_physio)) return("Multimodal")
  if (has_face) return("Facial"); if (has_physio) return("Physiological"); "Other"
}
classify_task <- function(t) {
  s <- tolower(trimws(t)); if (is.na(s) || s == "") return(NA_character_)
  if (grepl("binary", s)) return("Binary"); if (grepl("multi", s)) return("Multi-class")
  if (grepl("regress", s)) return("Regression"); "Other"
}
classify_arch <- function(a) {
  s <- tolower(trimws(a)); if (is.na(s) || s == "") return(NA_character_)
  dl <- grepl("cnn|convolutional|lstm|gru|rnn|transformer|vit|resnet|vgg|inception|densenet|deep|gan|autoencoder|u-net|efficientnet|alexnet|googlenet|mobilenet|bilstm|n-cnn|i3d|c3d|bert|gpt|gemini|belief network|neural network|mlp|perceptron", s)
  if (dl) "Deep learning" else "Traditional ML"
}
A$SG1_modality_count     <- ifelse(n2(A$modality_count) >= 2, "Multimodal", "Single-modality")
A$SG2_modality_type      <- vapply(A$modality_list, classify_modality, character(1))
A$SG3_setting            <- trimws(as.character(A$setting_macro))
A$SG4_task_type          <- vapply(A$task_type, classify_task, character(1))
A$SG5_architecture       <- vapply(A$architecture, classify_arch, character(1))
A$SG6_reference_standard <- trimws(as.character(A$ref_std_macro))
A$SG7_testing_strategy   <- trimws(as.character(A$SG7_testing_stratum))
A$SG8_environment        <- trimws(as.character(A$data_collection_env))

# --- PATCH (this run): collapse SG3->4 & SG6->3 levels, relabel SG8, to match
#     the manuscript's ORIGINAL S7 / Table-4 presentation. Pure relabel of the
#     macro columns; pooling method (per-stratum REML, min k=3) is unchanged. ---
.collapse <- function(x, map, default) {
  x <- trimws(as.character(x)); out <- unname(map[x]); out[is.na(out)] <- default
  out[is.na(x) | x == "" | x == "NA"] <- NA_character_; out
}
A$SG3_setting <- .collapse(A$SG3_setting,
  c("NICU" = "Neonatal", "ICU (adult)" = "ICU", "Postoperative" = "Postoperative"),
  default = "Chronic/other")        # Chronic pain / Other clinical / Oncology-SCD / Lab-benchmark -> Chronic/other
A$SG6_reference_standard <- .collapse(A$SG6_reference_standard,
  c("Self-report (NRS/VAS)" = "Self-report", "Other/mixed" = "Other"),
  default = "Expert observational") # PSPI-FACS + Neonatal-behavioral + Expert-observer + Adult-behavioral -> Expert observational
A$SG8_environment <- .collapse(A$SG8_environment,
  c("Clinical" = "Real bedside", "Lab" = "Lab", "Mixed" = "Mixed"), default = NA_character_)
cat("\n[collapsed for table] SG3 ->", paste(sort(unique(na.omit(A$SG3_setting))), collapse = " / "), "\n")
cat(  "[collapsed for table] SG6 ->", paste(sort(unique(na.omit(A$SG6_reference_standard))), collapse = " / "), "\n")

MIN_K <- 3
pool_reml <- function(df) {
  if (nrow(df) < 2) return(data.frame(k = nrow(df), est = NA, ci.lb = NA, ci.ub = NA))
  m <- rma(yi, vi, data = df, method = "REML")
  data.frame(k = m$k, est = transf.ilogit(m$b)[1],
             ci.lb = transf.ilogit(m$ci.lb), ci.ub = transf.ilogit(m$ci.ub))
}
run_sg <- function(col, label) {
  d <- A[!is.na(A[[col]]) & A[[col]] != "" & A[[col]] != "NA", ]
  keep <- names(which(table(d[[col]]) >= MIN_K)); d <- d[d[[col]] %in% keep, ]
  if (length(keep) < 2) { cat(sprintf("  [%s] <2 strata >= k%d\n", label, MIN_K)); return(NULL) }
  per <- do.call(rbind, lapply(sort(keep), function(L)
    cbind(SG = label, level = L, pool_reml(d[d[[col]] == L, ]))))
  qm <- tryCatch({ mm <- rma(yi, vi, mods = ~ factor(d[[col]]), data = d, method = "REML")
                   sprintf("Q_M=%.2f df=%d p=%.3f", mm$QM, mm$p - 1, mm$QMp) },
                 error = function(e) NA_character_)
  cat(sprintf("  [%s] %s\n", label, qm)); per
}
cat("\n--- AUROC subgroups SG1-SG8 (per-stratum REML, min k=3) ---\n")
SG <- rbind(
  run_sg("SG1_modality_count",     "SG1 Modality count"),
  run_sg("SG2_modality_type",      "SG2 Modality type"),
  run_sg("SG3_setting",            "SG3 Clinical setting"),
  run_sg("SG4_task_type",          "SG4 Task type"),
  run_sg("SG5_architecture",       "SG5 Architecture"),
  run_sg("SG6_reference_standard", "SG6 Reference standard"),
  run_sg("SG7_testing_strategy",   "SG7 Testing strategy"),
  run_sg("SG8_environment",        "SG8 Environment"))
print(SG, row.names = FALSE, digits = 3)
write.csv(SG, file.path(OUT_DIR, "S7_subgroups_SG1-SG8.csv"), row.names = FALSE)

# Silo meta-regression (SuppFig2 basis) -- 4-macro silo, REML
m_silo <- rma(yi, vi, mods = ~ factor(silo_4macro), data = A, method = "REML")
cat("\n--- Silo meta-regression (SuppFig2) ---\n")
cat(sprintf("Q_moderator = %.2f, df = %d, p = %.4f\n",
            m_silo$QM, m_silo$p - 1, m_silo$QMp))

# Publication bias: AUROC (Egger) and HSROC (Deeks, unweighted -- matches 02_hsroc.R)
cat("\n--- Egger's regression test (AUROC) ---\n")
print(regtest(m_auc, model = "lm"))

H$ess   <- 4 * (H$TP + H$FN) * (H$FP + H$TN) / ((H$TP + H$FN) + (H$FP + H$TN))
H$lnDOR <- log(((H$TP + 0.5) * (H$TN + 0.5)) / ((H$FP + 0.5) * (H$FN + 0.5)))
deeks <- summary(lm(lnDOR ~ I(1 / sqrt(ess)), data = H))   # Deeks' test (unweighted)
cat("\n--- Deeks' funnel asymmetry test (HSROC) ---\n")
print(deeks$coefficients)

## ============================================================================
## 8. SUMMARY TABLES regenerated for N = 188 (replace the stale 189 sheets)
## ============================================================================
quadas_cols <- c(`D1 Patient selection`  = "D1_PatientSelection_RoB",
                 `D2 Index test`         = "D2_IndexTest_RoB",
                 `D3 Reference standard`  = "D3_ReferenceStandard_RoB",
                 `D4 Flow & timing`      = "D4_FlowTiming_RoB",
                 `Overall`               = "overall_RoB")
qd <- sapply(quadas_cols, function(c)
  table(factor(toupper(trimws(S[[c]])), c("LOW", "UNCLEAR", "HIGH"))))
qd <- t(qd)   # rows = domains, cols = LOW / UNCLEAR / HIGH
write.csv(qd, file.path(OUT_DIR, "Summary_QUADAS2.csv"))

sx <- addmargins(table(Silo = S$silo_4macro, RefStd = S$ref_std_macro))
write.csv(sx, file.path(OUT_DIR, "Summary_Silo_x_RefStd.csv"))

cat_cols <- c(`Cat 1 Frame-level random sampling`            = "Cat1_FrameLevelRandom",
              `Cat 2 Sliding window w/o patient grouping`    = "Cat2_SlidingWindowNoGrouping",
              `Cat 3 Cross-validation w/o patient grouping`  = "Cat3_CrossValNoPatientGrouping",
              `Cat 5 IRR not reported or kappa < 0.70`       = "Cat5_IRR_NotReportedOrLow",
              `Cat 6 Sole reliance on public benchmark`      = "Cat6_SoleRelianceOnPublicBenchmark",
              `Cat 7 Class imbalance handling pre-split`     = "Cat7_ClassImbalanceHandlingPreSplit",
              `Cat 8 Hyperparameter tuning on test set`      = "Cat8_HyperparameterTuningOnTest")
lp <- data.frame(category = names(cat_cols),
                 n = sapply(cat_cols, function(c) sum(flag(S[[c]]))),
                 N = nrow(S), row.names = NULL)
lp$pct <- sprintf("%.1f%%", 100 * lp$n / lp$N)
write.csv(lp, file.path(OUT_DIR, "Summary_Leakage_Prevalence.csv"), row.names = FALSE)

lb <- as.data.frame(table(indicators = n2(S$n_leakage_indicators)))
lb$pct <- sprintf("%.1f%%", 100 * lb$Freq / nrow(S))
lb$cum <- sprintf("%.1f%%", 100 * cumsum(lb$Freq) / nrow(S))
names(lb)[2] <- "n"
write.csv(lb, file.path(OUT_DIR, "Summary_Leakage_Burden.csv"), row.names = FALSE)

cat("\nSummary tables written (QUADAS2 / Silo x RefStd / Leakage prevalence / burden).\n")
cat("NOTE: 'Cohort_Clusters' is a curated qualitative table; verify it no longer\n",
    "     references Cheng_2022 and keep its content otherwise unchanged.\n")

## ============================================================================
## 9. FIGURES — intentionally NOT produced here.
##    All submission figures are built separately in Python to match the
##    manuscript's existing house style (Fig3/Fig4/Fig5/Fig6/Fig7, SuppFig1-3/5).
##    This script is the NUMBERS engine only -> no plotting, nothing to crash.
## ============================================================================

## ============================================================================
## 8b. EXTRA DETERMINISTIC COUNTS  (Table 3 / Table 4 / abstract)
##     Pure counts & Cramer's V -- no estimator, so these match any tool exactly.
##     Each block is wrapped: a column mismatch warns but never stops the run.
## ============================================================================
cramer_v <- function(tab) {
  ch <- suppressWarnings(chisq.test(tab, correct = FALSE))
  n <- sum(tab); k <- min(nrow(tab), ncol(tab))
  sprintf("chi2=%.1f  df=%d  p=%.3g  V=%.3f",
          unname(ch$statistic), unname(ch$parameter), ch$p.value,
          sqrt(as.numeric(ch$statistic) / (n * (k - 1))))
}
try({
  cat("\n--- Studies contributing to >=1 of the 4 quantitative pools ---\n")
  inpool <- flag(S$in_AUROC_pool_main) | flag(S$in_HSROC_pool) |
            flag(S$in_PCC_pool_main)   | flag(S$in_ICC_pool_main)
  cat(sprintf("  %d of %d  (%.1f%%)\n", sum(inpool), nrow(S), 100 * sum(inpool) / nrow(S)))
})
try({
  cat("\n--- Per-pool RoB HIGH counts (GRADE / Table 4 / S10) ---\n")
  hi <- function(fl) { d <- S[flag(fl), ]
    sprintf("%d/%d", sum(toupper(trimws(d$overall_RoB)) == "HIGH"), nrow(d)) }
  cat(sprintf("  HSROC %s | AUROC %s | PCC %s | ICC %s\n",
      hi(S$in_HSROC_pool), hi(S$in_AUROC_pool_main), hi(S$in_PCC_pool_main), hi(S$in_ICC_pool_main)))
})
try({
  cat("\n--- Cramer's V structural couplings (Table 3, N=188) ---\n")
  cat("  S1 Silo x RefStd            ", cramer_v(table(S$silo_4macro, S$ref_std_macro)), "\n")
  cat("  S2 Silo x Setting           ", cramer_v(table(S$silo_4macro, S$setting_macro)), "\n")
  cat("  S3 Silo x Environment       ", cramer_v(table(S$silo_4macro, S$data_collection_env)), "\n")
  cat("  S4 Silo x (RefStd x Setting)", cramer_v(table(S$silo_4macro,
        interaction(S$ref_std_macro, S$setting_macro, drop = TRUE))), "\n")
})
try({
  cat("\n--- MAE descriptive stats (S4 / Fig4b) ---\n")
  mv <- n2(S$MAE_standardized_0_10); mv <- mv[!is.na(mv)]
  cat(sprintf("  k(standardized MAE)=%d | median=%.2f | range=%.2f-%.2f | below MCID(2)=%d/%d\n",
      length(mv), median(mv), min(mv), max(mv), sum(mv < 2), length(mv)))
})

## ---- 10. MASTER RESULTS TABLE ----------------------------------------------
table2 <- rbind(
  cbind(outcome = "AUROC", auroc_res[, c("k","est","ci.lb","ci.ub","pi.lb","pi.ub","I2")]),
  cbind(outcome = "PCC",   pcc_res [, c("k","est","ci.lb","ci.ub","pi.lb","pi.ub","I2")]),
  cbind(outcome = "ICC",   icc_res [, c("k","est","ci.lb","ci.ub","pi.lb","pi.ub","I2")])
)
write.csv(table2, file.path(OUT_DIR, "Table2_pooled_estimates.csv"), row.names = FALSE)

cat("\n=============================================================\n")
cat("DONE (numbers only). Outputs in ./", OUT_DIR, "/ :\n", sep = "")
cat("  Table2_pooled_estimates.csv  (AUROC/PCC/ICC; HSROC in the .txt)\n")
cat("  HSROC_reitsma_summary.txt    (sens/spec/AUC -- authoritative)\n")
cat("  S4_MAE_descriptive.csv\n")
cat("  S7_subgroups_SG1-SG8.csv\n")
cat("  S8_sensitivity_AUROC.csv  (Main + SA1-SA6 + 2 checks)\n")
cat("  influence_leave_one_out.csv  (diagnostic, not a registered SA)\n")
cat("  Summary_QUADAS2.csv / _Silo_x_RefStd.csv / _Leakage_Prevalence.csv / _Leakage_Burden.csv\n")
cat("  full_reproduction_log.txt  <-- SEND ME THIS ONE FILE\n")
cat("  (figures are generated in Python, not here)\n")
cat("=============================================================\n")
cat("\n>>> Full console log written to outputs/full_reproduction_log.txt\n")
while (sink.number() > 0) sink()
