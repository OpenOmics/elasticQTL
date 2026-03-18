#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: 02_prepare_clump_inputs.R <glm_file> <meta_csv_or_NONE> <out_dir> <out_prefix>")
}

glm_file <- args[1]
meta_file <- args[2]
out_dir <- args[3]
out_prefix <- args[4]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gwas <- fread(glm_file)
if ("TEST" %in% names(gwas)) gwas <- gwas[TEST == "ADD"]

req <- c("ID","P","BETA")
if (!all(req %in% names(gwas))) stop("GLM missing required columns: ", paste(setdiff(req, names(gwas)), collapse=", "))

# De-duplicate: keep smallest P per ID
setorder(gwas, ID, P)
gwas <- gwas[!duplicated(ID)]

# CHR/BP from GLM if available, else parse ID
chr_col <- if ("#CHROM" %in% names(gwas)) "#CHROM" else if ("CHROM" %in% names(gwas)) "CHROM" else NA
pos_col <- if ("POS" %in% names(gwas)) "POS" else if ("BP" %in% names(gwas)) "BP" else NA

tmp <- copy(gwas)
if (!is.na(chr_col) && !is.na(pos_col)) {
  tmp[, CHR := as.character(get(chr_col))]
  tmp[, BP  := as.integer(get(pos_col))]
} else {
  # fallback parse for chr:pos:ref:alt
  tmp[, c("CHR","BP") := tstrsplit(ID, ":", fixed=TRUE)[1:2]]
  tmp[, CHR := gsub("^chr","", CHR)]
  tmp[, BP := as.integer(BP)]
}

tmp[, CHR := gsub("^chr","", CHR)]
tmp[, CHR := suppressWarnings(as.integer(CHR))]

clump_input <- tmp[, .(SNP = ID, CHR, BP, P, BETA)]
setorder(clump_input, CHR, BP)

snplist_file <- file.path(out_dir, paste0(out_prefix, ".consistent_variants.snplist"))
clump_file   <- file.path(out_dir, paste0(out_prefix, ".clump_input.tsv"))
annot_file   <- file.path(out_dir, paste0(out_prefix, ".qtl_annotated.tsv"))

fwrite(gwas[, .(ID)], snplist_file, col.names = FALSE)
fwrite(clump_input, clump_file, sep="\t")

# optional metadata merge
annotated <- copy(gwas)
if (!is.na(meta_file) && meta_file != "" && meta_file != "NONE" && file.exists(meta_file)) {
  meta <- fread(meta_file)
  if ("variant_id" %in% names(meta)) {
    meta[, variant_id := fifelse(grepl("^chr", variant_id), variant_id, paste0("chr", variant_id))]
    setkey(meta, variant_id)
    setkey(gwas, ID)
    annotated <- meta[gwas, on=.(variant_id = ID)]
  }
}
fwrite(annotated, annot_file, sep="\t")

cat("Wrote:\n", snplist_file, "\n", clump_file, "\n", annot_file, "\n")
cat("N variants:", nrow(gwas), "\n")
