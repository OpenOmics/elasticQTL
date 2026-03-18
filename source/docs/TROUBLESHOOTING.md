# Troubleshooting

## 1) Variant IDs differ across cohorts (rsIDs vs chr:pos:ref:alt)
Symptoms:
- Elastic net step says “0 variants found”
- Annotation/merges fail

Fix:
Standardize IDs in PLINK2 before running the pipeline. Example:

```bash
plink2 --bfile OLD_PREFIX   --set-all-var-ids @:#:$r:$a   --new-id-max-allele-len 50   --make-bed --out NEW_PREFIX
```

Then run the pipeline using `--bfile NEW_PREFIX`.

---

## 2) Allele sign mismatch between EN weights and PLINK2 A1
The pipeline writes `allele_map_from_raw.tsv` from the `.raw` export and aligns weights to A1 where possible.
If you still see unexpected sign flips, confirm:
- PLINK2 GLM `A1` allele is the allele counted in your genotype dosage coding
- `.raw` suffix allele indicates which allele the dosage is counting

---

## 3) Complete-case filters remove too many samples
Default is strict: `--geno 0 --mind 0`.

Options:
- relax to `--geno 0.01 --mind 0.01`
- set `--missing-policy drop_samples` to drop remaining NA rows in R (still no imputation)
- set `--missing-policy mean_impute` to impute missing genotypes with training-fold means (leakage-safe)
- set `--missing-policy auto` to choose a strategy based on missingness thresholds

---

## 4) PLINK2 GLM file has duplicates
PLINK2 output may include multiple TEST rows (ADD/DOM/etc).
The pipeline filters to `TEST == "ADD"` when present.

---

## 5) Running on HPC modules
Example:
```bash
module load plink/1.9
module load plink2/2.0
module load R/4.2
```

Then run:
```bash
bash scripts/run_pipeline.sh --config configs/my_study.env
```


## HPC modules / Biowulf

If you are on an HPC system using environment modules (e.g., Biowulf), and loading PLINK 2 unloads PLINK 1 (or vice versa),
enable module management in the pipeline by setting `USE_MODULES=1` and providing `MODULE_PLINK1` and `MODULE_PLINK2` in your config.
The driver script will run `module load` before each relevant step.
