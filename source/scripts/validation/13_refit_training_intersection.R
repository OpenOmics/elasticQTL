#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop(paste(
    "Usage: 13_refit_training_intersection.R",
    "<train_raw> <train_pheno> <pheno_name> <variants_txt> <outdir>",
    "<alpha> <nfolds> <seed> <missing_policy> [auto_cell_rate] [auto_sample_rate]"
  ))
}

train_raw <- args[1]
train_pheno <- args[2]
pheno_name <- args[3]
variants_txt <- args[4]
outdir <- args[5]
alpha_val <- as.numeric(args[6])
nfolds <- as.integer(args[7])
seed <- as.integer(args[8])
missing_policy <- args[9]
auto_cell_rate <- if (length(args) >= 10) as.numeric(args[10]) else 0.001
auto_sample_rate <- if (length(args) >= 11) as.numeric(args[11]) else 0.02

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
set.seed(seed)

# Helpers
impute_mean <- function(X) {
  mu <- colMeans(X, na.rm=TRUE)
  mu[is.na(mu)] <- 0
  for (j in seq_len(ncol(X))) {
    idx <- is.na(X[, j])
    if (any(idx)) X[idx, j] <- mu[j]
  }
  X
}

# Load variant list
vars <- fread(variants_txt, header=FALSE)[[1]]
vars <- unique(as.character(vars))
if (length(vars) == 0) stop("No variants in variants_txt.")

# Load training raw
geno <- fread(train_raw)
std <- intersect(c("FID","IID","PAT","MAT","SEX","PHENOTYPE"), names(geno))
gt_cols <- setdiff(names(geno), std)

# Map raw columns -> Variant + counted allele (suffix after last underscore)
split_last <- function(x) {
  # returns list(id, allele)
  m <- regexpr("_[^_]+$", x)
  ifelse(m > 0, substr(x, 1, m-1), x)
}
get_allele <- function(x) {
  m <- regexpr("_[^_]+$", x)
  ifelse(m > 0, substr(x, m+1, nchar(x)), NA_character_)
}

allele_map <- data.table(
  raw_col = gt_cols,
  Variant = split_last(gt_cols),
  Train_Counted_Allele = toupper(get_allele(gt_cols))
)
# Rename genotype columns to Variant (strip suffix)
setnames(geno, gt_cols, allele_map$Variant)

# Keep only selected variants present
vars_use <- intersect(vars, names(geno))
if (length(vars_use) == 0) stop("None of the selected variants were present in train_raw after stripping suffix.")
fwrite(data.table(Variant=vars_use), file.path(outdir, "refit_variants_used.txt"), col.names=FALSE)

# Load phenotype and merge by IID
ph <- fread(train_pheno)
if (!all(c("FID","IID", pheno_name) %in% names(ph))) {
  stop("Training phenotype must include FID IID and ", pheno_name)
}

setkey(geno, IID)
setkey(ph, IID)
dat <- ph[geno, on=.(IID), nomatch=0]

y <- dat[[pheno_name]]
keep_y <- !(is.na(y) | y %in% c(-9, "-9", "NA", "NaN", ""))
dat <- dat[keep_y]
y <- dat[[pheno_name]]

X <- as.matrix(dat[, ..vars_use])

# Missingness metrics
miss_cells <- sum(is.na(X))
total_cells <- length(X)
miss_cell_rate <- miss_cells / total_cells
miss_samples_rate <- sum(!complete.cases(X)) / nrow(X)

fwrite(data.table(
  n_samples=nrow(X),
  n_variants=ncol(X),
  missing_cells=miss_cells,
  missing_cell_rate=miss_cell_rate,
  samples_with_any_missing_rate=miss_samples_rate
), file.path(outdir, "missingness_report_refit.tsv"), sep="\t")

chosen <- missing_policy
if (missing_policy == "auto") {
  if (miss_cells == 0) chosen <- "error"
  else if (miss_cell_rate <= auto_cell_rate) chosen <- "mean_impute"
  else if (miss_samples_rate <= auto_sample_rate) chosen <- "drop_samples"
  else stop("AUTO policy: missingness too high. Consider relaxing upstream QC or using mean_impute explicitly.")
}
fwrite(data.table(
  requested_policy=missing_policy,
  chosen_policy=chosen,
  auto_impute_max_cell_rate=auto_cell_rate,
  auto_drop_max_sample_rate=auto_sample_rate
), file.path(outdir, "missingness_policy_decision_refit.tsv"), sep="\t")

if (chosen == "drop_samples") {
  keep <- complete.cases(X)
  X <- X[keep, , drop=FALSE]
  y <- y[keep]
  dat <- dat[keep]
} else if (chosen == "mean_impute") {
  X <- impute_mean(X)
} else if (chosen == "error") {
  if (anyNA(X)) stop("NA found in X under policy=error.")
} else {
  stop("missing_policy must be error|drop_samples|mean_impute|auto")
}

set.seed(seed)
cvfit <- cv.glmnet(X, y, alpha=alpha_val, nfolds=nfolds, standardize=TRUE)
saveRDS(cvfit, file.path(outdir, "refit_model.rds"))

# Save weights
b <- as.matrix(coef(cvfit, s="lambda.min"))
wt <- data.table(Variant=rownames(b), Weight=as.numeric(b[,1]))
fwrite(wt, file.path(outdir, "refit_weights_all.tsv"), sep="\t")

wt_nz <- wt[Variant != "(Intercept)" & Weight != 0][order(-abs(Weight))]
fwrite(wt_nz, file.path(outdir, "refit_weights_nonzero.tsv"), sep="\t")

# Save training allele map for used variants
am_used <- unique(allele_map[Variant %in% vars_use, .(Variant, Train_Counted_Allele)])
fwrite(am_used, file.path(outdir, "train_counted_alleles_used.tsv"), sep="\t")

cat("Refit complete. Non-zero weights:", nrow(wt_nz), "\n")
cat("Saved:", file.path(outdir, "refit_model.rds"), "\n")
