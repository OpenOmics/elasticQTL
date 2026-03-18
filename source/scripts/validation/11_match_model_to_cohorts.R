#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: 11_match_model_to_cohorts.R <model_weights_tsv> <cohorts_tsv> <outdir>")
}

model_weights <- args[1]
cohorts_tsv <- args[2]
outdir <- args[3]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Helpers
is_snp_base <- function(x) {
  x <- toupper(x)
  x %in% c("A","C","G","T")
}
comp_base <- function(x) {
  x <- toupper(x)
  out <- x
  out[x == "A"] <- "T"
  out[x == "T"] <- "A"
  out[x == "C"] <- "G"
  out[x == "G"] <- "C"
  out
}
is_ambiguous_pair <- function(a, b) {
  a <- toupper(a); b <- toupper(b)
  ( (a=="A" & b=="T") | (a=="T" & b=="A") | (a=="C" & b=="G") | (a=="G" & b=="C") )
}

read_bim <- function(bim_path) {
  bim <- fread(bim_path, header=FALSE)
  if (ncol(bim) < 6) stop("BIM file has <6 cols: ", bim_path)
  setnames(bim, c("chr","id","cm","pos","a1","a2"))
  bim[, chr_clean := gsub("^chr","", as.character(chr))]
  bim[, pos := as.integer(pos)]
  bim[, a1 := toupper(as.character(a1))]
  bim[, a2 := toupper(as.character(a2))]
  bim
}

# Load model weights
mw <- fread(model_weights)
if (!all(c("Variant","Weight") %in% names(mw))) stop("model_weights_tsv must contain Variant and Weight.")

# Parse CHR/POS/REF/ALT from Variant when possible
mw[, c("chr_raw","pos","ref","alt") := tstrsplit(Variant, ":", fixed=TRUE)]
mw[, chr_clean := gsub("^chr","", chr_raw)]
mw[, pos := suppressWarnings(as.integer(pos))]
mw[, ref := toupper(ref)]
mw[, alt := toupper(alt)]

# Keep only variants with chr+pos for matching
mw_pos <- mw[!is.na(chr_clean) & !is.na(pos)]
if (nrow(mw_pos) == 0) stop("No model variants had parseable chr:pos. Expected Variant IDs like chr#:pos:ref:alt.")

# Load cohort manifest
coh <- fread(cohorts_tsv)
if (!all(c("cohort","bfile") %in% names(coh))) stop("cohorts_tsv must have columns: cohort, bfile")

summary_list <- list()
available_sets <- list()

# Write per-cohort mapping + available sets
for (i in seq_len(nrow(coh))) {
  cohort <- as.character(coh$cohort[i])
  bfile <- as.character(coh$bfile[i])
  bim_path <- paste0(bfile, ".bim")
  if (!file.exists(bim_path)) stop("Missing BIM: ", bim_path)

  bim <- read_bim(bim_path)

  # Join by chr+pos (may create multiple matches per Variant)
  setkey(mw_pos, chr_clean, pos)
  setkey(bim, chr_clean, pos)
  j <- bim[mw_pos, on=.(chr_clean, pos), nomatch=0]

  # If no matches, write empty files
  if (nrow(j) == 0) {
    fwrite(data.table(), file.path(outdir, paste0(cohort, ".mapping.tsv")), sep="\t")
    fwrite(data.table(), file.path(outdir, paste0(cohort, ".available_model_variants.txt")), col.names=FALSE)
    summary_list[[cohort]] <- data.table(cohort=cohort, bfile=bfile, available_variants=0L)
    available_sets[[cohort]] <- character(0)
    next
  }

  # Determine allele-set match score where ref/alt exist and are SNP bases
  j[, model_has_alleles := is_snp_base(ref) & is_snp_base(alt) & !is.na(ref) & !is.na(alt)]
  j[, bim_has_alleles := is_snp_base(a1) & is_snp_base(a2) & !is.na(a1) & !is.na(a2)]

  j[, direct_match := FALSE]
  j[, complement_match := FALSE]
  j[, ambiguous_pair := FALSE]
  j[model_has_alleles & bim_has_alleles, `:=`(
    direct_match = ( (ref==a1 & alt==a2) | (ref==a2 & alt==a1) ),
    complement_match = ( (ref==comp_base(a1) & alt==comp_base(a2)) | (ref==comp_base(a2) & alt==comp_base(a1)) ),
    ambiguous_pair = is_ambiguous_pair(ref, alt)
  )]

  # Assign score: prefer direct, then complement; if alleles missing, allow (score 0)
  j[, match_score := fifelse(model_has_alleles & bim_has_alleles,
                             fifelse(direct_match, 2L, fifelse(complement_match, 1L, -9L)),
                             0L)]

  # Keep only viable matches
  j <- j[match_score >= 0]

  # Choose best match per model Variant (max score). If ties, keep first.
  j[, abs_weight := abs(Weight)]
  setorder(j, Variant, -match_score, -abs_weight)
  j <- j[!duplicated(Variant)]
  j[, abs_weight := NULL]

  # Strand flip is when complement_match used (and not ambiguous)
  j[, strand_flip := (match_score == 1L)]
  j[ambiguous_pair == TRUE, strand_flip := FALSE]  # ambiguous; cannot detect

  out_map <- j[, .(
    Variant,
    Model_Weight = Weight,
    chr_clean, pos,
    ref, alt,
    cohort_id = id,
    cohort_a1 = a1,
    cohort_a2 = a2,
    match_score,
    strand_flip,
    ambiguous_pair
  )]

  map_path <- file.path(outdir, paste0(cohort, ".mapping.tsv"))
  fwrite(out_map, map_path, sep="\t")

  # Available model variants in this cohort
  available_vars <- unique(out_map$Variant)
  available_sets[[cohort]] <- available_vars
  fwrite(data.table(Variant=available_vars),
         file.path(outdir, paste0(cohort, ".available_model_variants.txt")),
         col.names=FALSE)

  summary_list[[cohort]] <- data.table(
    cohort=cohort,
    bfile=bfile,
    available_variants=length(available_vars)
  )
}

# Availability summary
summary_dt <- rbindlist(summary_list, fill=TRUE)
fwrite(summary_dt, file.path(outdir, "availability_summary.tsv"), sep="\t")

# Intersection sets
cohort_names <- as.character(coh$cohort)
sets <- lapply(cohort_names, function(x) available_sets[[x]])
names(sets) <- cohort_names

intersection_all <- Reduce(intersect, sets)
intersection_any <- Reduce(union, sets)

fwrite(data.table(Variant=intersection_all),
       file.path(outdir, "intersection_all.model_variants.txt"),
       col.names=FALSE)
fwrite(data.table(Variant=intersection_any),
       file.path(outdir, "intersection_any.model_variants.txt"),
       col.names=FALSE)

# For each cohort, also write extract IDs corresponding to each intersection set
for (cohort in cohort_names) {
  map_path <- file.path(outdir, paste0(cohort, ".mapping.tsv"))
  mp <- fread(map_path)

  # Extract IDs for ALL intersection
  ids_all <- mp[Variant %in% intersection_all, unique(cohort_id)]
  fwrite(data.table(id=ids_all),
         file.path(outdir, paste0("extract_ids_all.", cohort, ".txt")),
         col.names=FALSE)

  # Extract IDs for ANY intersection (cohort-available subset)
  ids_any <- mp[Variant %in% intersection_any, unique(cohort_id)]
  fwrite(data.table(id=ids_any),
         file.path(outdir, paste0("extract_ids_any.", cohort, ".txt")),
         col.names=FALSE)

  # Extract IDs for cohort-only available set
  ids_av <- mp[, unique(cohort_id)]
  fwrite(data.table(id=ids_av),
         file.path(outdir, paste0("extract_ids_available.", cohort, ".txt")),
         col.names=FALSE)
}

cat("Wrote availability + intersections to:", outdir, "\n")
cat("Intersection ALL:", length(intersection_all), "\n")
cat("Intersection ANY:", length(intersection_any), "\n")
