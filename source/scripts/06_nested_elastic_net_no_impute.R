#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop(paste(
    "Usage: 06_nested_elastic_net.R",
    "<glm_file> <geno_raw> <pheno> <pheno_name> <outdir>",
    "<alpha> <outer_folds> <inner_folds>",
    "[p_thresholds_csv] [label_mode] [seed]",
    "[missing_policy] [auto_impute_max_cell_rate] [auto_drop_max_sample_rate]"
  ))
}

glm_file   <- args[1]
geno_raw   <- args[2]
pheno_file <- args[3]
pheno_name <- args[4]
outdir     <- args[5]
alpha_val  <- as.numeric(args[6])
K_outer    <- as.integer(args[7])
K_inner    <- as.integer(args[8])

p_csv <- if (length(args) >= 9)  args[9]  else "0.2,0.3,0.5"
label_mode <- if (length(args) >= 10) args[10] else "percent"
seed <- if (length(args) >= 11) as.integer(args[11]) else 123L

missing_policy <- if (length(args) >= 12) args[12] else "error"
auto_cell_rate <- if (length(args) >= 13) as.numeric(args[13]) else 0.001
auto_sample_rate <- if (length(args) >= 14) as.numeric(args[14]) else 0.02

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

label_from_thresh <- function(t, mode=c("percent","legacy")) {
  mode <- match.arg(mode)
  if (mode == "legacy") {
    # legacy mapping: 0.2->p02, 0.5->p05 (ambiguous vs 0.05)
    s <- sub("^0\\.", "", format(t, scientific=FALSE))
    if (nchar(s) == 1) s <- paste0("0", s)
    return(paste0("p", s))
  } else {
    # percent mapping: 0.2->p20, 0.05->p05, 0.5->p50
    s <- sprintf("%02d", as.integer(round(t * 100)))
    return(paste0("p", s))
  }
}

# -------------------------
# Load association results
# -------------------------
glm_res <- fread(glm_file)
if ("TEST" %in% names(glm_res)) glm_res <- glm_res[TEST == "ADD"]
stopifnot(all(c("ID","P") %in% names(glm_res)))

# -------------------------
# Load genotype raw and allele map
# -------------------------
geno <- fread(geno_raw)
base_cols <- intersect(c("FID","IID","PAT","MAT","SEX","PHENOTYPE"), names(geno))
snp_cols <- setdiff(names(geno), base_cols)

allele_map <- data.table(
  raw_col = snp_cols,
  Variant = sub("_[A-Za-z]+$", "", snp_cols),
  Dosage_Allele = ifelse(grepl("_[A-Za-z]+$", snp_cols),
                         sub("^.*_([A-Za-z]+)$", "\\1", snp_cols),
                         NA_character_)
)
setnames(geno, snp_cols, allele_map$Variant)

fwrite(allele_map[, .(Variant, Dosage_Allele)],
       file.path(outdir, "allele_map_from_raw.tsv"), sep="\t")

# -------------------------
# Load phenotype and merge on IID
# -------------------------
ph <- fread(pheno_file)
if (!all(c("FID","IID", pheno_name) %in% names(ph))) {
  stop("Phenotype must contain FID IID and ", pheno_name)
}

setkey(geno, IID)
setkey(ph, IID)
dat <- ph[geno, on=.(IID), nomatch=0]

y <- dat[[pheno_name]]
keep_y <- !(is.na(y) | y %in% c(-9, "-9", "NA", "NaN", ""))
dat <- dat[keep_y]
y <- dat[[pheno_name]]

cat("Samples after phenotype filter:", nrow(dat), "\n")

# -------------------------
# Build variant sets
# -------------------------
thresholds <- as.numeric(strsplit(p_csv, ",")[[1]])
thresholds <- thresholds[!is.na(thresholds)]
thresholds <- unique(thresholds)

variant_sets <- list(all = unique(glm_res$ID))
for (t in thresholds) {
  nm <- label_from_thresh(t, label_mode)
  variant_sets[[nm]] <- unique(glm_res[P < t]$ID)
}

# -------------------------
# Missingness summary + policy selection
# -------------------------
all_vars_present <- intersect(variant_sets$all, names(dat))
if (length(all_vars_present) == 0) stop("No overlapping variants between GLM results and genotype matrix.")
X_all <- as.matrix(dat[, ..all_vars_present])

miss_cells <- sum(is.na(X_all))
total_cells <- length(X_all)
miss_cell_rate <- ifelse(total_cells > 0, miss_cells / total_cells, NA_real_)

miss_samples_n <- sum(!complete.cases(X_all))
miss_samples_rate <- miss_samples_n / nrow(X_all)

miss_vars_n <- sum(colSums(is.na(X_all)) > 0)
miss_report <- data.table(
  n_samples = nrow(X_all),
  n_variants = ncol(X_all),
  missing_cells = miss_cells,
  missing_cell_rate = miss_cell_rate,
  samples_with_any_missing = miss_samples_n,
  samples_with_any_missing_rate = miss_samples_rate,
  variants_with_any_missing = miss_vars_n
)
fwrite(miss_report, file.path(outdir, "missingness_report_all_variants.tsv"), sep="\t")

chosen_policy <- missing_policy
if (missing_policy == "auto") {
  if (miss_cells == 0) {
    chosen_policy <- "error"  # no missing => strict is fine
  } else if (!is.na(miss_cell_rate) && miss_cell_rate <= auto_cell_rate) {
    chosen_policy <- "mean_impute"
  } else if (miss_samples_rate <= auto_sample_rate) {
    chosen_policy <- "drop_samples"
  } else {
    stop(paste0(
      "AUTO policy could not choose a safe strategy. ",
      "missing_cell_rate=", signif(miss_cell_rate, 4),
      ", samples_with_any_missing_rate=", signif(miss_samples_rate, 4),
      ". Consider relaxing --geno/--mind or explicitly setting --missing-policy mean_impute."
    ))
  }
}

decision <- data.table(
  requested_policy = missing_policy,
  chosen_policy = chosen_policy,
  auto_impute_max_cell_rate = auto_cell_rate,
  auto_drop_max_sample_rate = auto_sample_rate
)
fwrite(decision, file.path(outdir, "missingness_policy_decision.tsv"), sep="\t")
cat("Missingness policy requested:", missing_policy, " | chosen:", chosen_policy, "\n")

# Apply drop_samples globally if chosen
if (chosen_policy == "drop_samples") {
  keep <- complete.cases(X_all)
  dat <- dat[keep]
  y <- y[keep]
  cat("After drop_samples on ALL variants, samples:", nrow(dat), "\n")
} else if (chosen_policy == "error") {
  if (anyNA(X_all)) stop("Found NA in genotype matrix; choose --missing-policy mean_impute or drop_samples, or enforce stricter --geno/--mind.")
} else if (chosen_policy == "mean_impute") {
  # do nothing here; imputation is performed fold-safely inside run_nested()
} else {
  stop("missing_policy must be one of: error, drop_samples, mean_impute, auto")
}

# Helper: fold-safe mean imputation (training-fold means)
impute_with_train_means <- function(X_train, X_test, fallback_means=NULL) {
  mu <- colMeans(X_train, na.rm = TRUE)

  # handle columns that are all-NA in training (rare)
  if (any(is.na(mu))) {
    if (!is.null(fallback_means)) {
      mu[is.na(mu)] <- fallback_means[is.na(mu)]
    }
    # still NA? set to 0 (very rare)
    mu[is.na(mu)] <- 0
  }

  for (j in seq_len(ncol(X_train))) {
    idx_tr <- is.na(X_train[, j])
    if (any(idx_tr)) X_train[idx_tr, j] <- mu[j]
    idx_te <- is.na(X_test[, j])
    if (any(idx_te)) X_test[idx_te, j] <- mu[j]
  }
  list(X_train = X_train, X_test = X_test, mu = mu)
}

# -------------------------
# Nested CV
# -------------------------
run_nested <- function(variant_ids, label) {
  variant_ids <- intersect(unique(variant_ids), names(dat))
  if (length(variant_ids) == 0) stop("No variants found for set: ", label)

  X <- as.matrix(dat[, ..variant_ids])

  # If chosen_policy is error or drop_samples, X should have no NA
  # If mean_impute, allow NAs and impute within folds
  if (chosen_policy %in% c("error", "drop_samples") && anyNA(X)) {
    stop("NA exists in X for set ", label, " but policy is ", chosen_policy)
  }

  # Outer folds
  set.seed(seed)
  fold_outer <- sample(rep(1:K_outer, length.out=nrow(X)))
  fwrite(data.table(IID=dat$IID, outer_fold=fold_outer),
         file.path(outdir, paste0("outer_folds_", label, ".tsv")), sep="\t")

  preds_all <- rep(NA_real_, length(y))
  r2_outer <- rep(NA_real_, K_outer)
  coef_list <- vector("list", K_outer)

  # fallback means (global) for any all-NA-in-training edge case
  fallback_means <- colMeans(X, na.rm = TRUE)
  fallback_means[is.na(fallback_means)] <- 0

  for (k in 1:K_outer) {
    test_idx <- which(fold_outer == k)
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)

    X_train <- X[train_idx, , drop=FALSE]
    y_train <- y[train_idx]
    X_test  <- X[test_idx, , drop=FALSE]
    y_test  <- y[test_idx]

    if (chosen_policy == "mean_impute") {
      imp <- impute_with_train_means(X_train, X_test, fallback_means)
      X_train <- imp$X_train
      X_test <- imp$X_test
    }

    if (anyNA(X_train) || anyNA(X_test)) stop("NA remains after preprocessing in fold ", k, " for set ", label)

    set.seed(seed + 1000 + k)
    cvfit <- cv.glmnet(X_train, y_train, alpha=alpha_val, nfolds=K_inner, standardize=TRUE)
    p <- as.numeric(predict(cvfit, X_test, s="lambda.min"))

    r2_outer[k] <- suppressWarnings(cor(y_test, p)^2)
    preds_all[test_idx] <- p

    b <- as.matrix(coef(cvfit, s="lambda.min"))
    b <- b[rownames(b)!="(Intercept)", , drop=FALSE]
    coef_vec <- as.numeric(b[,1]); names(coef_vec) <- rownames(b)
    coef_list[[k]] <- coef_vec

    cat(sprintf("[%s] outer fold %d/%d R2=%.4f\n", label, k, K_outer, r2_outer[k]))
  }

  overall_r2 <- suppressWarnings(cor(y, preds_all)^2)

   coef_mat <- do.call(cbind, coef_list)
  n_selected_per_fold <- colSums(coef_mat != 0)
  fwrite(data.table(
    Fold = seq_len(K_outer),
    n_selected = n_selected_per_fold
  ), file.path(outdir, paste0("n_selected_per_fold_", label, ".tsv")), sep="\t")
  colnames(coef_mat) <- paste0("Fold", seq_len(K_outer))
  sel_freq <- rowSums(coef_mat != 0)
  n_selected_at_least_5 <- sum(sel_freq >= 5)
  mean_coef <- rowMeans(coef_mat)

  stab <- data.table(
    Variant = rownames(coef_mat),
    Selection_Freq = sel_freq,
    Mean_Coefficient = mean_coef
  )[order(-Selection_Freq, -abs(Mean_Coefficient))]

  fwrite(data.table(Fold=seq_len(K_outer), R2=r2_outer),
         file.path(outdir, paste0("nested_outer_r2_", label, ".tsv")), sep="\t")
  fwrite(data.table(IID=dat$IID, Observed=y, Predicted=preds_all),
         file.path(outdir, paste0("nested_predictions_", label, ".tsv")), sep="\t")
  fwrite(stab, file.path(outdir, paste0("nested_stability_", label, ".tsv")), sep="\t")

  data.table(
    label=label,
    n_variants=length(variant_ids),
    n_selected_at_least_5=n_selected_at_least_5,
    mean_r2=mean(r2_outer, na.rm=TRUE),
    sd_r2=sd(r2_outer, na.rm=TRUE),
    overall_r2=overall_r2,
    missing_policy_used=chosen_policy
  )
}

summary_rows <- rbindlist(lapply(names(variant_sets), function(nm) run_nested(variant_sets[[nm]], nm)), fill=TRUE)
fwrite(summary_rows, file.path(outdir, "nested_EN_summary.tsv"), sep="\t")
print(summary_rows)
