## =============================================================================
## 01_data_prep.R
## Purpose: Load the extracted dataset and prepare all derived variables used by
##          the downstream analysis scripts: quantitative-pool membership flags,
##          the eight pre-registered subgroup factors, and the effect-size
##          transformations (logit, Fisher's z) with their sampling variances.
## Inputs : data/Supplementary_Data_Extracted_Dataset.xlsx  (sheet "Studies")
## Outputs: data frame `dat` (in the global environment) plus helper functions.
## Depends: readxl
## Notes  : Run from the repository root. All other scripts source this file.
##          No results are hard-coded; everything is derived from the dataset.
## =============================================================================

if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package 'readxl' is required. Install with install.packages('readxl').")

DATA_PATH <- file.path("data", "Supplementary_Data_Extracted_Dataset.xlsx")
if (!file.exists(DATA_PATH))
  stop("Dataset not found at ", DATA_PATH,
       ". Download it from the data archive (see README) and place it in ./data/.")

dat <- as.data.frame(readxl::read_excel(DATA_PATH, sheet = "Studies"))
message(sprintf("Loaded %d studies, %d columns.", nrow(dat), ncol(dat)))

## ---- small helpers ----------------------------------------------------------
## Coerce the project's TRUE/1/Yes-style flag columns to logical.
as_flag <- function(x) {
  toupper(trimws(as.character(x))) %in% c("1", "TRUE", "YES", "Y", "T")
}
## Safe numeric coercion (non-numeric -> NA).
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

## ---- quantitative-pool membership -------------------------------------------
## These flag columns are defined in the dataset (one column per pool).
dat$flag_hsroc <- as_flag(dat$in_HSROC_pool)
dat$flag_auroc <- as_flag(dat$in_AUROC_pool_main)
dat$flag_mae   <- as_flag(dat$in_MAE_pool_main)
dat$flag_pcc   <- as_flag(dat$in_PCC_pool_main)
dat$flag_icc   <- as_flag(dat$in_ICC_pool_main)

## ---- numeric effect sizes (use the pre-cleaned *_numeric columns) -----------
dat$AUROC_val <- num(dat$AUROC_numeric)
dat$sens_val  <- num(dat$sensitivity_numeric)
dat$spec_val  <- num(dat$specificity_numeric)
dat$MAE_val   <- num(dat$MAE_standardized_0_10)
dat$pcc_val   <- num(dat$pearson_r_numeric)
dat$icc_val   <- num(dat$ICC_numeric)
dat$N_val     <- num(dat$N_numeric)
dat$TP <- num(dat$TP_num); dat$FP <- num(dat$FP_num)
dat$TN <- num(dat$TN_num); dat$FN <- num(dat$FN_num)

## ---- effect-size transformations + sampling variances -----------------------
## Logit transform for proportions/AUROC.
logit     <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

## AUROC on the logit scale. Each study's AUROC sampling variance is derived
## from its AUROC point estimate and sample size using the Hanley & McNeil
## (1982) formula -- the standard approach in ROC meta-analysis, which does not
## require the source paper to report a confidence interval. The variance is
## then mapped to the logit scale by the delta method. This reproduces the
## published pooled AUROC of 0.904 [0.879, 0.924] exactly (see 03_auroc_pool.R).
##
## Hanley & McNeil variance of an AUROC `a` with `n_pos` positive and `n_neg`
## negative cases:
##   Q1 = a / (2 - a);  Q2 = 2 a^2 / (1 + a)
##   Var(a) = [ a(1-a) + (n_pos-1)(Q1 - a^2) + (n_neg-1)(Q2 - a^2) ] / (n_pos n_neg)
## When only the total N is available, positives and negatives are taken as N/2
## each (balanced-class assumption), as in the primary analysis.
hanley_var <- function(auc, n_pos, n_neg) {
  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)
  (auc * (1 - auc) +
     (n_pos - 1) * (q1 - auc^2) +
     (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg)
}

dat$auroc_logit <- ifelse(dat$AUROC_val > 0 & dat$AUROC_val < 1,
                          logit(dat$AUROC_val), NA_real_)
## Hanley-McNeil variance on the AUROC scale (balanced-class N/2 split),
## then delta-method transform to the logit scale: Var(logit a) = Var(a)/[a(1-a)]^2.
.auroc_var_raw <- ifelse(!is.na(dat$AUROC_val) & !is.na(dat$N_val) & dat$N_val > 0,
                         mapply(function(a, n) hanley_var(a, n / 2, n / 2),
                                dat$AUROC_val, dat$N_val),
                         NA_real_)
dat$auroc_logit_var <- ifelse(!is.na(.auroc_var_raw),
                              .auroc_var_raw / (dat$AUROC_val * (1 - dat$AUROC_val))^2,
                              NA_real_)

## Fisher's z for Pearson r and ICC: z = atanh(r), Var(z) = 1 / (N - 3).
fisher_z   <- function(r) atanh(r)
inv_fisher <- function(z) tanh(z)
dat$pcc_z      <- ifelse(abs(dat$pcc_val) < 1, fisher_z(dat$pcc_val), NA_real_)
dat$pcc_z_var  <- ifelse(!is.na(dat$pcc_z) & dat$N_val > 3, 1 / (dat$N_val - 3), NA_real_)
dat$icc_z      <- ifelse(abs(dat$icc_val) < 1, fisher_z(dat$icc_val), NA_real_)
dat$icc_z_var  <- ifelse(!is.na(dat$icc_z) & dat$N_val > 3, 1 / (dat$N_val - 3), NA_real_)

## MAE (standardised to a 0-10 scale) with SE = MAE / sqrt(N)
## (Cochrane Handbook approach for a mean with no reported dispersion).
dat$mae_se  <- ifelse(!is.na(dat$MAE_val) & !is.na(dat$N_val) & dat$N_val > 0,
                      dat$MAE_val / sqrt(dat$N_val), NA_real_)
dat$mae_var <- dat$mae_se^2

## ---- the eight pre-registered subgroup factors (SG1-SG8) --------------------
## Derived from the dataset's descriptive columns following the PROSPERO protocol.

## SG1 - modality count: single vs multimodal (>= 2 modalities).
dat$SG1_modality_count <- ifelse(num(dat$modality_count) >= 2, "Multimodal", "Single-modality")

## SG2 - modality type: facial / physiological / multimodal / other.
classify_modality <- function(ml) {
  s <- tolower(trimws(ml))
  if (is.na(s) || s == "") return(NA_character_)
  has_face   <- grepl("fac", s)
  n_types    <- length(strsplit(s, ",")[[1]])
  has_physio <- grepl("physiolog|eeg|ecg|bioimped|skin", s)
  if (n_types >= 2 || (has_face && has_physio)) return("Multimodal")
  if (has_face)   return("Facial")
  if (has_physio) return("Physiological")
  return("Other")
}
dat$SG2_modality_type <- vapply(dat$modality_list, classify_modality, character(1))

## SG3 - clinical setting (macro categories as coded in the dataset).
dat$SG3_setting <- trimws(as.character(dat$setting_macro))

## SG4 - task type: binary / multi-class / regression (collapse narrative-only).
classify_task <- function(t) {
  s <- tolower(trimws(t))
  if (is.na(s) || s == "") return(NA_character_)
  if (grepl("binary", s))      return("Binary")
  if (grepl("multi", s))       return("Multi-class")
  if (grepl("regress", s))     return("Regression")
  return("Other")
}
dat$SG4_task_type <- vapply(dat$task_type, classify_task, character(1))

## SG5 - architecture: traditional ML vs deep learning.
## Deep learning if any DL keyword appears; otherwise traditional ML.
classify_arch <- function(a) {
  s <- tolower(trimws(a))
  if (is.na(s) || s == "") return(NA_character_)
  dl <- grepl("cnn|convolutional|lstm|gru|rnn|transformer|vit|resnet|vgg|inception|densenet|deep|gan|autoencoder|u-net|efficientnet|alexnet|googlenet|mobilenet|bilstm|n-cnn|i3d|c3d|bert|gpt|gemini|belief network|neural network|mlp|perceptron", s)
  if (dl) "Deep learning" else "Traditional ML"
}
dat$SG5_architecture <- vapply(dat$architecture, classify_arch, character(1))

## SG6 - reference standard (macro categories as coded in the dataset).
dat$SG6_reference_standard <- trimws(as.character(dat$ref_std_macro))

## SG7 - testing strategy: pre-computed stratum column in the dataset.
dat$SG7_testing_strategy <- trimws(as.character(dat$SG7_testing_stratum))

## SG8 - data-collection environment (Lab / Clinical / Mixed).
dat$SG8_environment <- trimws(as.character(dat$data_collection_env))

## ---- overall risk of bias (QUADAS-2) ----------------------------------------
dat$overall_RoB <- toupper(trimws(as.character(dat$overall_RoB)))

## ---- console summary --------------------------------------------------------
message("Pool sizes: ",
        sprintf("HSROC=%d  AUROC=%d  MAE=%d  PCC=%d  ICC=%d",
                sum(dat$flag_hsroc), sum(dat$flag_auroc), sum(dat$flag_mae),
                sum(dat$flag_pcc), sum(dat$flag_icc)))
message("In >= 1 pool: ",
        sum(dat$flag_hsroc | dat$flag_auroc | dat$flag_mae |
              dat$flag_pcc | dat$flag_icc), " of ", nrow(dat))

invisible(TRUE)
