#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: 00_make_keep_list.R <pheno> <pheno_name> <exclude_or_NONE> <out_keep> <out_report>")
}

pheno_file <- args[1]
pheno_name <- args[2]
exclude_file <- args[3]
out_keep <- args[4]
out_report <- args[5]

ph <- fread(pheno_file)
if (!all(c("FID","IID", pheno_name) %in% names(ph))) {
  stop("Phenotype file must contain columns: FID, IID, and ", pheno_name)
}

y <- ph[[pheno_name]]
is_missing <- is.na(y) | y %in% c(-9, "-9", "NA", "NaN", "")

ph_keep <- ph[!is_missing, .(FID, IID)]
setkey(ph_keep, FID, IID)

excluded_n <- 0L
if (!is.na(exclude_file) && exclude_file != "" && exclude_file != "NONE" && file.exists(exclude_file)) {
  ex <- fread(exclude_file, header = FALSE)
  if (ncol(ex) < 2) stop("Exclude file must have at least 2 columns: FID IID")
  setnames(ex, c("FID","IID"))
  setkey(ex, FID, IID)
  excluded_n <- nrow(ex)
  ph_keep <- ph_keep[!ex, on=.(FID, IID)]
}

fwrite(ph_keep, out_keep, sep="\t", col.names = FALSE)

report <- data.table(
  metric = c("pheno_total", "pheno_nonmissing", "excluded_ids", "kept_final"),
  value  = c(nrow(ph), nrow(ph[!is_missing]), excluded_n, nrow(ph_keep))
)
fwrite(report, out_report, sep="\t")

cat("Wrote keep list:", out_keep, "\n")
print(report)
