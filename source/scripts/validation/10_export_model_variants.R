#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: 10_export_model_variants.R <model_rds> <out_weights_tsv> <out_variants_txt>")
}

model_rds <- args[1]
out_weights <- args[2]
out_vars <- args[3]

cvfit <- readRDS(model_rds)
coefs <- as.matrix(coef(cvfit, s = "lambda.min"))

dt <- data.table(
  Variant = rownames(coefs),
  Weight = as.numeric(coefs[, 1])
)
dt <- dt[Variant != "(Intercept)" & Weight != 0]

# Parse chr:pos:ref:alt if present
dt[, c("chr_raw","pos","ref","alt") := tstrsplit(Variant, ":", fixed=TRUE)]
dt[, chr_clean := gsub("^chr", "", chr_raw)]
dt[, pos := suppressWarnings(as.integer(pos))]

# Keep only rows with parsed CHR+POS when possible; but still write all
fwrite(dt, out_weights, sep="\t")
fwrite(dt[, .(Variant)], out_vars, col.names = FALSE)

cat("Non-zero variants:", nrow(dt), "\n")
cat("Wrote weights:", out_weights, "\n")
cat("Wrote variant list:", out_vars, "\n")
