#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: 09_add_maf_from_raw.R <final_table_csv> <geno_raw> <out_csv> <variant_colname>")

final_csv <- args[1]
geno_raw <- args[2]
out_csv <- args[3]
variant_col <- args[4]

tab <- fread(final_csv)
geno <- fread(geno_raw)

base_cols <- intersect(c("FID","IID","PAT","MAT","SEX","PHENOTYPE"), names(geno))
snp_cols <- setdiff(names(geno), base_cols)
setnames(geno, snp_cols, sub("_[A-Za-z]+$", "", snp_cols))

vars <- intersect(tab[[variant_col]], names(geno))

maf <- rbindlist(lapply(vars, function(v) {
  af <- mean(geno[[v]], na.rm=TRUE) / 2
  data.table(Variant=v, Allele_Freq=af, MAF=pmin(af, 1-af))
}))

setnames(maf, "Variant", variant_col)
out <- merge(tab, maf, by=variant_col, all.x=TRUE)
fwrite(out, out_csv)
cat("Wrote:", out_csv, "\n")
