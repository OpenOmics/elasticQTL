#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop(paste(
    "Usage: 14_score_validation_cohort.R",
    "<refit_model_rds> <train_raw> <cohort_raw> <cohort_mapping_tsv>",
    "<out_scores_csv> <out_qc_tsv> <ambiguous_policy> <missing_policy>"
  ))
}

refit_rds <- args[1]
train_raw <- args[2]
cohort_raw <- args[3]
map_tsv <- args[4]
out_scores <- args[5]
out_qc <- args[6]
ambig_policy <- args[7]   # drop|keep
missing_policy <- args[8] # error|drop_samples|mean_impute

# Helpers
is_snp_base <- function(x) toupper(x) %in% c("A","C","G","T")
comp_base <- function(x) {
  x <- toupper(x)
  out <- x
  out[x == "A"] <- "T"; out[x == "T"] <- "A"
  out[x == "C"] <- "G"; out[x == "G"] <- "C"
  out
}
impute_mean <- function(X) {
  mu <- colMeans(X, na.rm=TRUE)
  mu[is.na(mu)] <- 0
  for (j in seq_len(ncol(X))) {
    idx <- is.na(X[, j])
    if (any(idx)) X[idx, j] <- mu[j]
  }
  X
}
split_last <- function(x) {
  m <- regexpr("_[^_]+$", x)
  ifelse(m > 0, substr(x, 1, m-1), x)
}
get_suffix <- function(x) {
  m <- regexpr("_[^_]+$", x)
  ifelse(m > 0, substr(x, m+1, nchar(x)), NA_character_)
}

# Load refit model weights
cvfit <- readRDS(refit_rds)
b <- as.matrix(coef(cvfit, s="lambda.min"))
intercept <- as.numeric(b[1,1])
wt <- data.table(Variant=rownames(b), Weight_refit=as.numeric(b[,1]))
wt <- wt[Variant != "(Intercept)" & Weight_refit != 0]
if (nrow(wt) == 0) stop("No non-zero weights in refit model.")

# Training counted alleles map (from TRAIN raw header)
train_hdr <- names(fread(train_raw, nrows=0))
std <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE")
train_gt <- setdiff(train_hdr, std)

train_map <- data.table(
  raw_col = train_gt,
  Variant = split_last(train_gt),
  Train_Counted_Allele = toupper(get_suffix(train_gt))
)
train_map <- unique(train_map[, .(Variant, Train_Counted_Allele)], by="Variant")

# Cohort mapping (Variant -> cohort_id, strand_flip, ambiguous_pair)
mp <- fread(map_tsv)
need <- c("Variant","cohort_id","strand_flip","ambiguous_pair")
if (!all(need %in% names(mp))) stop("mapping_tsv missing required columns: ", paste(setdiff(need, names(mp)), collapse=", "))

# Restrict mapping to variants with non-zero weights
mp <- mp[Variant %in% wt$Variant]
if (nrow(mp) == 0) stop("No overlap between cohort mapping and refit weights.")

# Load cohort raw and parse header to map cohort_id -> raw_col + counted allele
geno <- fread(cohort_raw)
hdr <- names(geno)
gt_cols <- setdiff(hdr, std)
col_map <- data.table(
  raw_col = gt_cols,
  cohort_id = split_last(gt_cols),
  Val_Counted_Allele_Raw = toupper(get_suffix(gt_cols))
)
# Choose unique raw_col per cohort_id (if duplicates, keep first)
setorder(col_map, cohort_id)
col_map <- col_map[!duplicated(cohort_id)]

# Merge everything
dt <- merge(mp, wt, by="Variant", all.x=FALSE)
dt <- merge(dt, train_map, by="Variant", all.x=FALSE)
dt <- merge(dt, col_map, by="cohort_id", all.x=FALSE)

if (nrow(dt) == 0) stop("After merging mapping/weights/alleles, no variants remain.")

# Adjust validation counted allele for strand flips when allele is SNP base
dt[, Val_Counted_Allele_Physical := Val_Counted_Allele_Raw]
dt[strand_flip == TRUE & is_snp_base(Val_Counted_Allele_Raw),
   Val_Counted_Allele_Physical := comp_base(Val_Counted_Allele_Raw)]

# Identify ambiguous mismatches
dt[, ambiguous_mismatch := FALSE]
dt[ambiguous_pair == TRUE &
     is_snp_base(Train_Counted_Allele) & is_snp_base(Val_Counted_Allele_Physical) &
     Train_Counted_Allele != Val_Counted_Allele_Physical,
   ambiguous_mismatch := TRUE]

if (ambig_policy == "drop") {
  dt <- dt[ambiguous_mismatch == FALSE]
}
if (nrow(dt) == 0) stop("All variants dropped due to ambiguous mismatches.")

# Determine flips needed (after strand adjust)
dt[, flip := (Train_Counted_Allele != Val_Counted_Allele_Physical)]

# Build genotype matrix in the order of dt
X <- as.matrix(geno[, dt$raw_col, with=FALSE])
# Ensure numeric
mode(X) <- "numeric"

# Apply flips
flip_idx <- which(dt$flip)
if (length(flip_idx) > 0) {
  X[, flip_idx] <- 2 - X[, flip_idx]
}

# Missingness handling for scoring
if (missing_policy == "drop_samples") {
  keep <- complete.cases(X)
  X <- X[keep, , drop=FALSE]
  geno <- geno[keep]
} else if (missing_policy == "mean_impute") {
  X <- impute_mean(X)
} else if (missing_policy == "error") {
  if (anyNA(X)) stop("NA found in cohort genotype matrix under policy=error.")
} else {
  stop("missing_policy must be error|drop_samples|mean_impute")
}

# Score
score <- intercept + as.numeric(X %*% dt$Weight_refit)

scores_dt <- data.table(
  FID = geno$FID,
  IID = geno$IID,
  EN_score_refit = score,
  N_Variants_Used = ncol(X),
  N_Flipped = length(flip_idx),
  N_Dropped_Ambiguous = sum(mp$ambiguous_pair, na.rm=TRUE) - sum(dt$ambiguous_pair, na.rm=TRUE)
)
fwrite(scores_dt, out_scores)

# QC table
qc <- dt[, .(
  Variant, cohort_id, raw_col,
  Weight_refit,
  Train_Counted_Allele,
  Val_Counted_Allele_Raw,
  Val_Counted_Allele_Physical,
  strand_flip,
  ambiguous_pair,
  ambiguous_mismatch,
  flip
)]
fwrite(qc, out_qc, sep="\t")

cat("Scored cohort raw:", cohort_raw, "\n")
cat("Wrote scores:", out_scores, "\n")
cat("Wrote QC:", out_qc, "\n")
