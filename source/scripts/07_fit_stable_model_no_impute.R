#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop(paste(
    "Usage: 07_fit_stable_model.R",
    "<stability_tsv> <geno_raw> <pheno> <pheno_name> <outdir>",
    "<alpha> <min_folds> <label>",
    "[seed] [missing_policy] [auto_impute_max_cell_rate] [auto_drop_max_sample_rate]"
  ))
}

stab_file  <- args[1]
geno_raw   <- args[2]
pheno_file <- args[3]
pheno_name <- args[4]
outdir     <- args[5]
alpha_val  <- as.numeric(args[6])
min_folds  <- as.integer(args[7])
label      <- args[8]

seed <- if (length(args) >= 9) as.integer(args[9]) else 123L
missing_policy <- if (length(args) >= 10) args[10] else "error"
auto_cell_rate <- if (length(args) >= 11) as.numeric(args[11]) else 0.001
auto_sample_rate <- if (length(args) >= 12) as.numeric(args[12]) else 0.02

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

stab <- fread(stab_file)
stopifnot(all(c("Variant","Selection_Freq") %in% names(stab)))

stable_vars <- stab[Selection_Freq >= min_folds, Variant]
cat("Stable variants (>= ", min_folds, " folds): ", length(stable_vars), "\n", sep="")

# Load genotypes
geno <- fread(geno_raw)
base_cols <- intersect(c("FID","IID","PAT","MAT","SEX","PHENOTYPE"), names(geno))
snp_cols <- setdiff(names(geno), base_cols)
setnames(geno, snp_cols, sub("_[A-Za-z]+$", "", snp_cols))

# Load phenotype and merge
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

stable_vars <- intersect(stable_vars, names(dat))
if (length(stable_vars) == 0) stop("No stable variants present in genotype table.")

X <- as.matrix(dat[, ..stable_vars])

# Missingness summary for stable set
miss_cells <- sum(is.na(X))
total_cells <- length(X)
miss_cell_rate <- ifelse(total_cells > 0, miss_cells / total_cells, NA_real_)
miss_samples_n <- sum(!complete.cases(X))
miss_samples_rate <- miss_samples_n / nrow(X)

miss_report <- data.table(
  n_samples = nrow(X),
  n_variants = ncol(X),
  missing_cells = miss_cells,
  missing_cell_rate = miss_cell_rate,
  samples_with_any_missing = miss_samples_n,
  samples_with_any_missing_rate = miss_samples_rate
)
fwrite(miss_report, file.path(outdir, paste0("missingness_report_stable_set_", label, ".tsv")), sep="\t")

chosen_policy <- missing_policy
if (missing_policy == "auto") {
  if (miss_cells == 0) {
    chosen_policy <- "error"
  } else if (!is.na(miss_cell_rate) && miss_cell_rate <= auto_cell_rate) {
    chosen_policy <- "mean_impute"
  } else if (miss_samples_rate <= auto_sample_rate) {
    chosen_policy <- "drop_samples"
  } else {
    stop(paste0(
      "AUTO policy could not choose a safe strategy for stable-set. ",
      "missing_cell_rate=", signif(miss_cell_rate, 4),
      ", samples_with_any_missing_rate=", signif(miss_samples_rate, 4),
      ". Consider --missing-policy mean_impute."
    ))
  }
}

decision <- data.table(
  requested_policy = missing_policy,
  chosen_policy = chosen_policy,
  auto_impute_max_cell_rate = auto_cell_rate,
  auto_drop_max_sample_rate = auto_sample_rate
)
fwrite(decision, file.path(outdir, paste0("missingness_policy_decision_stable_set_", label, ".tsv")), sep="\t")
cat("Stable-set missingness policy requested:", missing_policy, " | chosen:", chosen_policy, "\n")

if (chosen_policy == "drop_samples") {
  keep <- complete.cases(X)
  X <- X[keep, , drop=FALSE]
  y <- y[keep]
  dat <- dat[keep]
  cat("After drop_samples on stable set, samples:", nrow(X), "\n")
} else if (chosen_policy == "mean_impute") {
  mu <- colMeans(X, na.rm=TRUE)
  mu[is.na(mu)] <- 0
  for (j in seq_len(ncol(X))) {
    idx <- is.na(X[, j])
    if (any(idx)) X[idx, j] <- mu[j]
  }
} else if (chosen_policy == "error") {
  if (anyNA(X)) stop("NA found in X; no imputation allowed under policy=error.")
} else {
  stop("missing_policy must be one of: error, drop_samples, mean_impute, auto")
}

if (anyNA(X)) stop("NA remains after preprocessing (stable model).")

# Fit final CV glmnet on all samples using stable variants
cvfit <- cv.glmnet(X, y, alpha=alpha_val, nfolds=10, standardize=TRUE)

b <- as.matrix(coef(cvfit, s="lambda.min"))
b <- b[rownames(b)!="(Intercept)", , drop=FALSE]
retained <- sum(b[,1] != 0)
cat("Final retained variants (non-zero):", retained, "\n")

# Save model + training-cohort scores
saveRDS(cvfit, file.path(outdir, paste0("final_EN_model_", label, "_stable.rds")))

pred <- as.numeric(predict(cvfit, X, s="lambda.min"))
scores <- data.table(IID=dat$IID, Observed=y, Genetic_Score=pred)
fwrite(scores, file.path(outdir, paste0("genetic_scores_", label, "_stable.tsv")), sep="\t")
