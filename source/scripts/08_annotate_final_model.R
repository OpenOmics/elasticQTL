#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: 08_annotate_final_model.R <final_model_rds> <glm_assoc_file> <stability_tsv> <meta_csv_or_NONE> <allele_map_tsv> <out_csv>")
}

final_rds <- args[1]
assoc_file <- args[2]
stab_file <- args[3]
meta_file <- args[4]
allele_map_file <- args[5]
out_csv <- args[6]

cvfit <- readRDS(final_rds)

b <- as.matrix(coef(cvfit, s="lambda.min"))
b <- b[rownames(b)!="(Intercept)", , drop=FALSE]
final_weights <- data.table(Variant=rownames(b), EN_Final_Weight=as.numeric(b[,1]))[EN_Final_Weight != 0]

assoc <- fread(assoc_file)
if ("TEST" %in% names(assoc)) assoc <- assoc[TEST=="ADD"]
assoc <- assoc[, .(Variant=ID, QTL_Assoc_BETA=BETA, QTL_Assoc_P=P, Effect_Allele_A1=A1)]
assoc <- unique(assoc, by="Variant")

stab <- fread(stab_file)[, .(Variant, Selection_Freq)]
stab <- unique(stab, by="Variant")

amap <- fread(allele_map_file)
amap <- unique(amap, by="Variant")

final_tab <- merge(final_weights, assoc, by="Variant", all.x=TRUE)
final_tab <- merge(final_tab, stab, by="Variant", all.x=TRUE)
final_tab <- merge(final_tab, amap, by="Variant", all.x=TRUE)

# optional meta
if (!is.na(meta_file) && meta_file != "" && meta_file != "NONE" && file.exists(meta_file)) {
  meta <- fread(meta_file)
  if (all(c("variant_id","qtl_type","gene") %in% names(meta))) {
    meta[, Variant := fifelse(grepl("^chr", variant_id), variant_id, paste0("chr", variant_id))]
    meta <- unique(meta[, .(Variant, qtl_type, gene)], by="Variant")
    final_tab <- merge(final_tab, meta, by="Variant", all.x=TRUE)
  }
}

# Align sign to A1 when dosage allele available
final_tab[, EN_Weight_AlignedToA1 := EN_Final_Weight]
final_tab[!is.na(Dosage_Allele) & !is.na(Effect_Allele_A1) & Dosage_Allele != Effect_Allele_A1,
          EN_Weight_AlignedToA1 := -EN_Final_Weight]

setorder(final_tab, -abs(EN_Weight_AlignedToA1))
fwrite(final_tab, out_csv)
cat("Wrote:", out_csv, "\n")
