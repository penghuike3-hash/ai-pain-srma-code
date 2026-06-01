###############################################################################
#  RUN_ME.R  —  AI Pain Assessment SR-MA: one-click reproduction
#
#  HOW TO RUN (any ONE of these — no setwd, no path editing needed):
#    • RStudio:  open this file, click "Source" (top-right of the editor)
#    • Double-click the file (Windows, if .R is associated with Rscript)
#    • Console:  source("RUN_ME.R")
#
#  It finds itself, finds the dataset automatically (looks next to this script,
#  in ./data/, on the Desktop, and in common folders), installs any missing
#  packages, runs every pool, prints the numbers next to the published values,
#  and writes results to an "outputs" folder beside this script.
#
#  Dataset needed (place it anywhere near this script, e.g. same folder or
#  a ./data subfolder):  Supplementary_Data_Extracted_Dataset.xlsx
###############################################################################

options(stringsAsFactors = FALSE, warn = 1)
cat("\n==================================================================\n")
cat(" AI pain assessment DTA meta-analysis — one-click reproduction\n")
cat("==================================================================\n\n")

## ---------------------------------------------------------------- 0. locate
# Find the folder this script lives in, robustly across run methods.
this_file <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f[1]) else NA_character_
}, error = function(e) NA_character_)
if (is.na(this_file) && requireNamespace("rstudioapi", quietly = TRUE)) {
  this_file <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                        error = function(e) NA_character_)
}
script_dir <- if (!is.na(this_file) && nzchar(this_file)) dirname(this_file) else getwd()
cat("Script folder:", script_dir, "\n")

## ---------------------------------------------------------------- 1. packages
need <- c("readxl", "metafor", "mada")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) {
  cat("Installing", p, "...\n")
  install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(readxl); library(metafor); library(mada) })

## ---------------------------------------------------------------- 2. find data
DATA_NAME <- "Supplementary_Data_Extracted_Dataset.xlsx"
search_dirs <- unique(c(
  script_dir,
  file.path(script_dir, "data"),
  getwd(), file.path(getwd(), "data"),
  file.path(path.expand("~"), "Desktop"),
  file.path(path.expand("~"), "Desktop"),
  file.path(path.expand("~"), "Desktop", "ai-pain-srma-code"),
  file.path(path.expand("~"), "Desktop", "ai-pain-srma-code", "data")
))
hits <- character(0)
# exact name first
for (d in search_dirs) { f <- file.path(d, DATA_NAME); if (file.exists(f)) hits <- c(hits, f) }
# then any .xlsx that looks like the dataset, searched recursively from script dir & desktop
if (!length(hits)) {
  roots <- unique(c(script_dir, getwd(), file.path(path.expand("~"), "Desktop")))
  for (r in roots) if (dir.exists(r)) {
    cand <- list.files(r, pattern = "(Extracted_Dataset|v17_final).*\\.xlsx$",
                       recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    hits <- c(hits, cand)
  }
}
hits <- hits[file.exists(hits)]
if (!length(hits)) {
  stop("\n>>> Could not find the dataset '", DATA_NAME, "'.\n",
       "    Put it in the SAME folder as this script (or a 'data' subfolder),\n",
       "    then run again. Looked in:\n      ",
       paste(search_dirs, collapse = "\n      "), "\n")
}
DATA_FILE <- hits[1]
cat("Using dataset:", DATA_FILE, "\n\n")

OUT_DIR <- file.path(script_dir, "outputs")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, showWarnings = FALSE)

## ---------------------------------------------------------------- 3. load + prep
d <- as.data.frame(read_excel(DATA_FILE, sheet = "Studies"))
flag <- function(x) toupper(trimws(as.character(x))) %in% c("1","TRUE","YES","Y","T")
num  <- function(x) suppressWarnings(as.numeric(as.character(x)))
logit <- function(p) log(p/(1-p)); inv_logit <- function(x) 1/(1+exp(-x))

d$f_hsroc <- flag(d$in_HSROC_pool); d$f_auroc <- flag(d$in_AUROC_pool_main)
d$f_mae <- flag(d$in_MAE_pool_main); d$f_pcc <- flag(d$in_PCC_pool_main); d$f_icc <- flag(d$in_ICC_pool_main)

cat(sprintf("Loaded %d studies, %d columns.\n", nrow(d), ncol(d)))
cat(sprintf("Pool sizes:  HSROC=%d  AUROC=%d  MAE=%d  PCC=%d  ICC=%d   (paper: 14/69/43/35/25)\n",
            sum(d$f_hsroc), sum(d$f_auroc), sum(d$f_mae), sum(d$f_pcc), sum(d$f_icc)))
cat(sprintf("In >=1 pool: %d of %d   (paper: 117 of 189)\n\n",
            sum(d$f_hsroc|d$f_auroc|d$f_mae|d$f_pcc|d$f_icc), nrow(d)))

results <- list()

## ---------------------------------------------------------------- 4. HSROC
cat("----- HSROC (2x2 pool) -----  paper: sens 0.904, spec 0.912, AUC 0.951\n")
tab <- data.frame(TP=num(d$TP_num), FP=num(d$FP_num), FN=num(d$FN_num), TN=num(d$TN_num))[d$f_hsroc, ]
tab <- tab[complete.cases(tab) & rowSums(tab) > 0, ]
fit <- reitsma(tab, method = "ml"); s <- summary(fit)
co <- s$coefficients
sens <- plogis(co["tsens.(Intercept)","Estimate"])
spec <- 1 - plogis(co["tfpr.(Intercept)","Estimate"])
sens_ci <- plogis(co["tsens.(Intercept)", c("CILB","CIUB")])
spec_ci <- rev(1 - plogis(co["tfpr.(Intercept)", c("CILB","CIUB")]))
auc  <- s$AUC$AUC
cat(sprintf("  k=%d   sens=%.3f [%.3f, %.3f]   spec=%.3f [%.3f, %.3f]   AUC=%.3f\n\n",
            nrow(tab), sens, sens_ci[1], sens_ci[2], spec, spec_ci[1], spec_ci[2], auc))
results$HSROC <- list(k=nrow(tab), sens=sens, spec=spec, auc=auc,
                      sens_ci=sens_ci, spec_ci=spec_ci)

## ---------------------------------------------------------------- 5. AUROC
cat("----- AUROC pool -----  paper: 0.904 [0.879, 0.924], I2 94.5%, PI 0.618-0.982\n")
ap <- d[d$f_auroc, ]; a <- num(ap$AUROC_numeric); n <- num(ap$N_numeric)
ok <- !is.na(a) & !is.na(n) & a>0 & a<1 & n>0
a <- a[ok]; n <- n[ok]
hanley <- function(au, np, nn){ q1<-au/(2-au); q2<-2*au^2/(1+au)
  (au*(1-au)+(np-1)*(q1-au^2)+(nn-1)*(q2-au^2))/(np*nn) }
vraw <- mapply(function(au,nn) hanley(au, nn/2, nn/2), a, n)
yi <- logit(a); vi <- vraw/((a*(1-a))^2)
m <- rma(yi=yi, vi=vi, method="PM"); pr <- predict(m)
cat(sprintf("  k=%d weighted   AUROC=%.3f [%.3f, %.3f]   I2=%.1f%%   PI [%.3f, %.3f]\n",
            length(a), inv_logit(m$b[1]), inv_logit(m$ci.lb), inv_logit(m$ci.ub),
            m$I2, inv_logit(pr$pi.lb), inv_logit(pr$pi.ub)))
egg <- regtest(m, model="lm")
cat(sprintf("  Egger p = %.3f   (paper 0.41)\n\n", egg$pval))
results$AUROC <- list(k=length(a), est=inv_logit(m$b[1]),
                      ci=c(inv_logit(m$ci.lb),inv_logit(m$ci.ub)), I2=m$I2,
                      pi=c(inv_logit(pr$pi.lb),inv_logit(pr$pi.ub)), egger=egg$pval)

## ---------------------------------------------------------------- 6. continuous
fisher_pool <- function(r, n, label, paper){
  ok <- !is.na(r) & !is.na(n) & abs(r)<1 & n>3
  z <- atanh(r[ok]); vi <- 1/(n[ok]-3)
  m <- rma(yi=z, vi=vi, method="PM"); pr <- predict(m)
  cat(sprintf("  %-4s k=%d   pooled=%.3f [%.3f, %.3f]   I2=%.1f%%   (paper %s)\n",
              label, sum(ok), tanh(m$b[1]), tanh(m$ci.lb), tanh(m$ci.ub), m$I2, paper))
  list(k=sum(ok), est=tanh(m$b[1]), ci=c(tanh(m$ci.lb),tanh(m$ci.ub)), I2=m$I2)
}
cat("----- Continuous pools -----\n")
# MAE: SE = MAE/sqrt(N)
mp <- d[d$f_mae, ]; mae <- num(mp$MAE_standardized_0_10); nm <- num(mp$N_numeric)
okm <- !is.na(mae)&!is.na(nm)&nm>0; vi <- (mae[okm]/sqrt(nm[okm]))^2
mm <- rma(yi=mae[okm], vi=vi, method="PM"); prm <- predict(mm)
cat(sprintf("  MAE  k=%d   pooled=%.2f [%.2f, %.2f]   I2=%.1f%%   (paper 0.79 [0.57,1.01], I2 100%%)\n",
            sum(okm), mm$b[1], mm$ci.lb, mm$ci.ub, mm$I2))
results$MAE <- list(k=sum(okm), est=as.numeric(mm$b[1]), I2=mm$I2)
results$PCC <- fisher_pool(num(d$pearson_r_numeric)[d$f_pcc], num(d$N_numeric)[d$f_pcc], "PCC", "0.688 [0.596,0.763], I2 93.2%")
results$ICC <- fisher_pool(num(d$ICC_numeric)[d$f_icc], num(d$N_numeric)[d$f_icc], "ICC", "0.638 [0.511,0.737], I2 94.3%")
cat("\n")

## ---------------------------------------------------------------- 7. subgroups SG1-SG8
cat("----- Subgroups SG1-SG8 (AUROC pool) -----  paper: all Q_M p >= 0.25 (no difference)\n")
ap2 <- d[d$f_auroc, ]; a2 <- num(ap2$AUROC_numeric); n2 <- num(ap2$N_numeric)
keep <- !is.na(a2)&!is.na(n2)&a2>0&a2<1&n2>0
ap2 <- ap2[keep,]; a2<-a2[keep]; n2<-n2[keep]
ap2$yi <- logit(a2); ap2$vi <- mapply(function(au,nn) hanley(au,nn/2,nn/2), a2, n2)/((a2*(1-a2))^2)

low <- function(x) tolower(ifelse(is.na(x),"",as.character(x)))
ap2$SG1 <- ifelse(num(ap2$modality_count)>=2,"Multimodal","Single")
ap2$SG2 <- sapply(ap2$modality_list, function(s){s<-low(s); if(s=="")NA else if(length(strsplit(s,",")[[1]])>=2)"Multimodal" else if(grepl("fac",s))"Facial" else if(grepl("physio|eeg|ecg|eda|emg",s))"Physiological" else if(grepl("voice|cry|audio",s))"Voice/cry" else "Other"})
ap2$SG3 <- sapply(ap2$clinical_setting, function(s){s<-low(s); if(grepl("nicu|neonat|infant",s))"Neonatal" else if(grepl("icu|critical",s))"ICU" else if(grepl("postop|surg|pacu",s))"Postoperative" else "Chronic/other"})
ap2$SG4 <- sapply(ap2$task_type, function(s){s<-low(s); if(startsWith(s,"binary"))"Binary" else if(grepl("multi",s))"Multi-class" else if(grepl("regress",s))"Regression" else "Other"})
ap2$SG5 <- sapply(ap2$architecture, function(s){s<-low(s); if(grepl("cnn|lstm|rnn|transformer|deep|resnet|vgg|gan|dnn|u-net|vit|bilstm",s))"Deep learning" else if(grepl("svm|forest|boost|knn|logistic|tree|bayes|lda",s))"Traditional ML" else "Other"})
ap2$SG6 <- sapply(ap2$ref_std_macro, function(s){if(grepl("Self-report",s))"Self-report" else if(grepl("Other|mixed",s))"Other" else "Expert observational"})
ap2$SG7 <- c("internal-only"="Internal CV","independent"="Independent","external"="External")[as.character(ap2$SG7_testing_stratum)]
ap2$SG8 <- sapply(ap2$data_collection_env, function(s){if(identical(s,"Lab"))"Lab" else if(identical(s,"Clinical"))"Bedside" else "Mixed"})

sg_names <- c(SG1="Modality count",SG2="Modality type",SG3="Setting",SG4="Task type",
              SG5="Architecture",SG6="Reference std",SG7="Testing strategy",SG8="Environment")
sg_rows <- list()
for (code in names(sg_names)) {
  g <- ap2[[code]]; tab <- table(g); keep <- names(tab)[tab>=3]
  sub <- ap2[g %in% keep, ]
  if (length(keep)>=2) {
    mm <- rma(yi=sub$yi, vi=sub$vi, mods=~factor(sub[[code]]), method="PM")
    cat(sprintf("  %-3s %-16s Q_M p=%.3f  (%d levels)\n", code, sg_names[code], mm$QMp, length(keep)))
    for (lv in names(tab)) if (tab[[lv]]>=3) {
      s1 <- ap2[g==lv,]; m1 <- rma(yi=s1$yi, vi=s1$vi, method="PM")
      sg_rows[[length(sg_rows)+1]] <- data.frame(SG=code, Axis=sg_names[code], Level=lv,
        k=as.integer(tab[[lv]]), AUROC=round(inv_logit(m1$b[1]),3), Q_M_p=round(mm$QMp,3))
    }
  }
}
cat("\n")

## ---------------------------------------------------------------- 8. save
sg_df <- do.call(rbind, sg_rows)
write.csv(sg_df, file.path(OUT_DIR, "SG1-SG8_results.csv"), row.names = FALSE)

summary_lines <- c(
  "AI pain assessment DTA meta-analysis — reproduction results",
  sprintf("Dataset: %s", basename(DATA_FILE)),
  sprintf("Pools: HSROC=%d AUROC=%d MAE=%d PCC=%d ICC=%d",
          results$HSROC$k, results$AUROC$k, results$MAE$k, results$PCC$k, results$ICC$k),
  sprintf("HSROC: sens=%.3f spec=%.3f AUC=%.3f", results$HSROC$sens, results$HSROC$spec, results$HSROC$auc),
  sprintf("AUROC: %.3f [%.3f, %.3f] I2=%.1f%% PI [%.3f, %.3f] Egger p=%.3f",
          results$AUROC$est, results$AUROC$ci[1], results$AUROC$ci[2], results$AUROC$I2,
          results$AUROC$pi[1], results$AUROC$pi[2], results$AUROC$egger),
  sprintf("MAE: %.2f I2=%.1f%%", results$MAE$est, results$MAE$I2),
  sprintf("PCC: %.3f [%.3f, %.3f] I2=%.1f%%", results$PCC$est, results$PCC$ci[1], results$PCC$ci[2], results$PCC$I2),
  sprintf("ICC: %.3f [%.3f, %.3f] I2=%.1f%%", results$ICC$est, results$ICC$ci[1], results$ICC$ci[2], results$ICC$I2),
  "",
  paste(capture.output(print(sg_df, row.names = FALSE)), collapse = "\n")
)
writeLines(summary_lines, file.path(OUT_DIR, "reproduction_summary.txt"))

cat("==================================================================\n")
cat(" DONE. Results saved to:\n")
cat("   ", file.path(OUT_DIR, "reproduction_summary.txt"), "\n")
cat("   ", file.path(OUT_DIR, "SG1-SG8_results.csv"), "\n")
cat("==================================================================\n")

if (!interactive()) { cat("\nPress [Enter] to close..."); tryCatch(readLines("stdin", 1), error=function(e) NULL) }
