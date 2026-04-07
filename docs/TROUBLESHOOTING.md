# Troubleshooting

This document covers common issues and their solutions for both the training pipeline (Step 1) and validation/scoring workflow (Step 2).

---

## General Issues

### Issue 1: Running on head node error

**Symptoms:**
```
############################################################
ERROR: Running this pipeline on a head node is not allowed.
Please submit this as a job to the cluster or run it from
an interactive node.
############################################################
```

**Cause:** You are trying to run the pipeline on the Skyline/OpenOmics head node, which is not allowed.

**Fix:** Grab an interactive node first:
```bash
srun -N 1 -n 1 --time=8:00:00 --mem=64gb -c 4 --pty bash
export PATH="/data/openomics/prod/elasticQTL/v0.1.0/bin:${PATH}"
pipeline --config config/study.env
```

---

### Issue 2: Pipeline not found in PATH

**Symptoms:**
```
bash: pipeline: command not found
```

**Cause:** The `bin/` directory is not in your PATH.

**Fix for Skyline/OpenOmics install:**
```bash
export PATH="/data/openomics/prod/elasticQTL/v0.1.0/bin:${PATH}"
```

**Fix for local clone:**
```bash
cd /path/to/elasticQTL
export PATH="$(pwd)/bin:${PATH}"
```

To make this permanent, add the export line to your `~/.bashrc` or `~/.bash_profile`.

---

### Issue 3: Permission denied when running pipeline

**Symptoms:**
```
bash: ./pipeline: Permission denied
```

**Cause:** The entrypoint scripts are not executable.

**Fix:**
```bash
chmod +x bin/pipeline bin/validate
```

---

### Issue 4: Config file not found

**Symptoms:**
```
ERROR: Config file not found: config/study.env
```

**Cause:** The config file path is incorrect or the file doesn't exist.

**Fix:**
- Check that the file exists: `ls -l config/study.env`
- Use an absolute path: `pipeline --config /full/path/to/config/study.env`
- Copy from template: `cp config/params.template.env config/study.env`

---

## Training Pipeline Issues (Step 1)

### Issue 5: Variant IDs differ across files (rsIDs vs chr:pos:ref:alt)

**Symptoms:**
- Elastic net step says "0 variants found"
- Annotation/merge steps fail
- Error: "No matching variants between files"

**Cause:** Inconsistent variant ID formats (e.g., `.bim` uses rsIDs but association results use chr:pos format).

**Fix:** Standardize variant IDs in PLINK2 before running the pipeline:
```bash
plink2 --bfile OLD_PREFIX \
  --set-all-var-ids @:#:\$r:\$a \
  --new-id-max-allele-len 50 \
  --make-bed \
  --out NEW_PREFIX
```

Then run the pipeline using `--bfile NEW_PREFIX`.

**Alternative format options:**
- `@:#` for chr:pos (no alleles)
- `@:#:\$r:\$a` for chr:pos:ref:alt (recommended)
- `@:#[hg38]\$r,\$a` for build-specific IDs

---

### Issue 6: Complete-case filters remove too many samples

**Symptoms:**
- Error: "Too few samples remaining after missingness filters"
- Warning: "XX% of samples dropped due to missing genotypes"
- Very low sample size in modeling step

**Cause:** Default strict filters (`--geno 0 --mind 0`) combined with `--missing-policy error` require perfect data.

**Fix (Option 1):** Relax genotype QC filters:
```bash
pipeline --config config/study.env \
  --geno 0.02 \
  --mind 0.05
```

**Fix (Option 2):** Allow sample dropping:
```bash
pipeline --config config/study.env \
  --missing-policy drop_samples
```

**Fix (Option 3):** Use mean imputation (if missingness is sparse):
```bash
pipeline --config config/study.env \
  --missing-policy mean_impute
```

**Fix (Option 4):** Let the pipeline decide automatically:
```bash
pipeline --config config/study.env \
  --missing-policy auto
```

See [`PIPELINE.md#missing-genotype-handling`](PIPELINE.md#missing-genotype-handling) for detailed policy explanations.

---

### Issue 7: PLINK2 GLM output has duplicate TEST entries

**Symptoms:**
- Multiple TEST rows per variant (ADD, DOMDEV, etc.)
- Unexpected number of variants in association results

**Cause:** PLINK2 outputs multiple test types by default.

**Fix:** The pipeline automatically filters to `TEST == "ADD"` when parsing GLM output. No action needed unless you want to change the test type.

**Manual verification:**
```bash
# Check what TEST types are present
cut -f12 01_glm_qtl/qtl_assoc.*.glm.linear | sort -u
```

---

### Issue 8: Allele sign mismatch between EN weights and PLINK2 A1

**Symptoms:**
- Coefficient signs seem reversed
- Positive weights for protective alleles (or vice versa)

**Cause:** The counted allele in the `.raw` export may differ from the A1 allele in PLINK2 association results.

**Fix:** The pipeline writes `allele_map_from_raw.tsv` in `06_en_nested/` to document which allele the dosage is counting. Use this to verify alignment.

**Check alignment:**
```bash
# View allele mapping
head -n 20 06_en_nested/allele_map_from_raw.tsv
```

The `.raw` suffix allele (e.g., `rs123_A`) indicates which allele the dosage counts.

---

### Issue 9: Out of memory errors during modeling

**Symptoms:**
```
Error: cannot allocate vector of size X Gb
```

**Cause:** Not enough memory allocated for the job, especially with large variant sets.

**Fix:** Request more memory when grabbing an interactive node:
```bash
srun -N 1 -n 1 --time=8:00:00 --mem=128gb -c 8 --pty bash
```

Or reduce the variant set size:
```bash
pipeline --config config/study.env \
  --clump-r2 0.1 \
  --p-thresholds 0.1,0.2
```

---

### Issue 10: Pipeline stops at specific step

**Symptoms:**
- Pipeline exits without error message
- Log shows step N completed but step N+1 never starts

**Cause:** Output file from previous step already exists and `--force` was not used.

**Fix:**
```bash
# Re-run with --force to overwrite existing outputs
pipeline --config config/study.env --force
```

Or manually remove the problematic output directory:
```bash
rm -rf /path/to/outdir/XX_stepname/
pipeline --config config/study.env
```

---

### Issue 11: No variants survive LD clumping

**Symptoms:**
- Step 4 log shows "0 variants after clumping"
- Downstream steps fail with empty variant lists

**Cause:** Clumping parameters are too stringent, or input variants are not in LD.

**Fix (Option 1):** Relax clumping parameters:
```bash
pipeline --config config/study.env \
  --clump-r2 0.5 \
  --clump-kb 1000
```

**Fix (Option 2):** Check that you have a reasonable number of input variants:
```bash
wc -l 01_glm_qtl/qtl_assoc.*.glm.linear
```

**Fix (Option 3):** Verify your candidate variant list is appropriate for QTL analysis.

---

### Issue 12: R package missing

**Symptoms:**
```
Error in library(glmnet) : there is no package called 'glmnet'
```

**Cause:** Required R package is not installed.

**Fix:** On Skyline/OpenOmics, packages should be pre-installed. If running elsewhere:
```r
# In R console
install.packages(c("data.table", "glmnet"))
```

---

## Validation/Scoring Issues (Step 2)

### Issue 13: Last cohort in manifest is skipped

**Symptoms:**
- Only N-1 cohorts are scored when manifest has N cohorts
- No error message, last cohort silently missing

**Cause:** Cohort manifest file is missing a final newline character.

**Fix:**
```bash
# Check if file ends with newline
tail -c 1 cohorts.tsv | od -An -tx1

# Add newline if needed (should show '0a')
echo "" >> cohorts.tsv
```

**Prevention:** Always ensure manifest files end with a newline.

---

### Issue 14: PLINK paths have `^M` characters (CRLF line endings)

**Symptoms:**
```
PLINK ERROR: --bfile prefix '/path/to/file^M' not found
```

**Cause:** Windows-style line endings (CRLF) in manifest or config files.

**Fix:**
```bash
# Convert to Unix line endings
dos2unix cohorts.tsv config/validation.env

# Or use sed
sed -i 's/\r$//' cohorts.tsv
sed -i 's/\r$//' config/validation.env
```

**Prevention:** Edit files on Linux/Mac or configure your editor to use LF line endings.

---

### Issue 15: Different variant counts across cohorts

**Symptoms:**
- `N_Variants_Used` differs across cohort score files
- Some cohorts have many fewer variants than others

**Cause:** Different cohorts have different genotyping platforms or coverage. This is expected behavior.

**Expected:** Each cohort score file reports its own `N_Variants_Used`. This can differ when:
- A variant is absent from a cohort (not in that cohort's `.bim`)
- A variant is dropped during harmonization (e.g., ambiguous SNP with allele disagreement)

**Fix (if you need identical variant sets):** Use strict intersection mode:
```bash
validate --config config/validation.env \
  --intersection-mode all
```

**Inspect per-cohort QC:**
```bash
# Review what happened to each variant
head -n 50 04_scores/wgs/wgs_qc_report.tsv
```

---

### Issue 16: Many variants dropped as ambiguous

**Symptoms:**
- High `N_Dropped_Ambiguous` in score files
- Warning: "XX ambiguous SNPs dropped"

**Cause:** Model includes many A/T or C/G SNPs, and validation cohorts code these alleles differently.

**Fix (Option 1 - Conservative):** Accept the reduction (recommended):
```bash
validate --config config/validation.env \
  --ambiguous-policy drop
```

**Fix (Option 2 - Permissive):** Keep ambiguous SNPs (use with caution):
```bash
validate --config config/validation.env \
  --ambiguous-policy keep
```

**Fix (Option 3 - Strict):** Error on ambiguous mismatches:
```bash
validate --config config/validation.env \
  --ambiguous-policy error
```

**Long-term fix:** Realign validation cohorts to a reference genome and re-export with consistent strand coding.

---

### Issue 17: Scores use Model_Weight instead of Weight_refit

**Symptoms:**
- Score files reference `Model_Weight` column
- No `refit_model.rds` in `03_refit/` directory

**Cause:** Refit step was skipped or failed.

**Fix (if intentional):** This is expected when using `--skip-refit`.

**Fix (if unintentional):** Check logs to see why refit failed:
```bash
cat logs/step_4_*.log
```

Then re-run without `--skip-refit`:
```bash
validate --config config/validation.env --force
```

---

### Issue 18: Training .raw file not found

**Symptoms:**
```
ERROR: Training .raw file not found: /path/to/file.raw
```

**Cause:** Path to training genotype matrix is incorrect, or training pipeline didn't complete successfully.

**Fix:** Verify the training pipeline completed Step 6:
```bash
ls -lh /path/to/training_outdir/05_genotypes/ld_variants_forEN.raw
```

Update `TRAIN_RAW` in your validation config to the correct path.

---

### Issue 19: Cohort bfile not found

**Symptoms:**
```
ERROR: Cohort bfile not found: /path/to/cohort_prefix
PLINK ERROR: Failed to open /path/to/cohort_prefix.bed
```

**Cause:** Bfile path in cohort manifest is incorrect or files don't exist.

**Fix:** Verify each bfile path in the manifest:
```bash
# Check each cohort
while read cohort bfile; do
  echo "Checking $cohort: $bfile"
  ls -lh ${bfile}.bed ${bfile}.bim ${bfile}.fam
done < cohorts.tsv
```

Update paths in the manifest to point to existing files.

---

### Issue 20: Allele flip counts seem wrong

**Symptoms:**
- `N_Flipped` is unexpectedly high or low
- Uncertainty about whether harmonization is correct

**Cause:** Validation cohort uses different allele coding or strand than training data.

**Fix:** Review the QC report for each cohort:
```bash
head -n 50 04_scores/wgs/wgs_qc_report.tsv
```

Look at these columns:
- `Strand_Flip_Detected` — Was a strand flip detected?
- `Dosage_Flipped` — Was the dosage flipped for this variant?
- `Training_Counted_Allele` vs `Cohort_Counted_Allele` — Which alleles are being counted?

See [`VALIDATION.md#understanding-allele-harmonization`](VALIDATION.md#understanding-allele-harmonization) for detailed explanation.

---

### Issue 21: No variants in intersection set

**Symptoms:**
```
ERROR: Intersection set is empty (0 variants)
```

**Cause:** No variants from the trained model are present in all cohorts (when using `--intersection-mode all`).

**Fix (Option 1):** Use a less strict intersection mode:
```bash
validate --config config/validation.env \
  --intersection-mode any
```

**Fix (Option 2):** Use a single cohort as reference:
```bash
validate --config config/validation.env \
  --intersection-mode cohort:wgs
```

**Fix (Option 3):** Check variant ID consistency:
```bash
# Compare variant IDs between training and validation
head -n 20 00_model/model_variants_chrpos.tsv
head -n 20 /path/to/validation_cohort.bim
```

May need to standardize IDs (see Issue 5).

---

## HPC and Environment Issues

### Issue 22: Module conflicts (PLINK 1.9 vs PLINK 2.0)

**Symptoms:**
- Loading PLINK2 module unloads PLINK 1.9 (or vice versa)
- "Command not found" for one PLINK version during pipeline run

**Cause:** HPC module system has conflicting PLINK modules.

**Fix:** Enable module management in the pipeline:

**In your config file:**
```bash
USE_MODULES=1
MODULE_PLINK1=plink/1.9
MODULE_PLINK2=plink/2.0
MODULE_R=R/4.2
MODULE_INIT=AUTO
```

The pipeline will automatically `module load` the correct version before each step.

---

### Issue 23: Module initialization script not found

**Symptoms:**
```
ERROR: Cannot find module initialization script
```

**Cause:** `MODULE_INIT=AUTO` detection failed.

**Fix:** Manually specify the module init script:
```bash
MODULE_INIT=/usr/share/modules/init/bash
```

Common paths:
- `/usr/share/modules/init/bash`
- `/etc/profile.d/modules.sh`
- `/usr/share/Modules/init/bash`

---

### Issue 24: Singularity/container issues

**Symptoms:**
- "Singularity not found"
- Container pull failures

**Cause:** Pipeline may be configured for containerized execution (not typical for Skyline/OpenOmics install).

**Fix:** Ensure PLINK and R are available in your PATH or via modules. The Skyline/OpenOmics install should have these pre-configured.

---

## Debugging Tips

### Tip 1: Use dry-run mode

Always test your command with `--dry-run` first:
```bash
pipeline --config config/study.env --dry-run
validate --config config/validation.env --dry-run
```

This shows what commands will be executed without actually running them.

---

### Tip 2: Check log files

Each step writes a log file to `logs/`:
```bash
# View most recent log
ls -lt logs/ | head

# Check for errors
grep -i error logs/*.log
grep -i warning logs/*.log
```

---

### Tip 3: Run partial pipelines for testing

Test data preparation without running expensive modeling:
```bash
# Training: run through genotype export only
pipeline --config config/study.env --to-step 6

# Validation: run matching only
validate --config config/validation.env --to-step 2
```

---

### Tip 4: Check manifest/params_used.txt

The training pipeline writes all parameters to `manifest/params_used.txt`:
```bash
cat manifest/params_used.txt
```

This is helpful for reproducing runs or debugging parameter issues.

---

### Tip 5: Verify input file formats

Common format issues:
```bash
# Check for CRLF line endings
file cohorts.tsv
# Should say "ASCII text", NOT "ASCII text, with CRLF line terminators"

# Check phenotype file has correct columns
head -n 2 phenotypes.tsv

# Check PLINK files are readable
plink2 --bfile PREFIX --freq
```

---

## Getting Help

If you encounter an issue not covered here:

1. **Check the logs** in `logs/` directory
2. **Review the parameters** in `manifest/params_used.txt`
3. **Try dry-run mode** to see what commands will execute
4. **Simplify** — test with a smaller dataset or fewer steps

When reporting a bug, please include:
- The config file used (redact sensitive paths if needed)
- Relevant log files from `logs/`
- `manifest/params_used.txt` (for training pipeline)
- Description of what you expected vs. what happened

---

## See Also

- [`PIPELINE.md`](PIPELINE.md) — Training workflow documentation
- [`VALIDATION.md`](VALIDATION.md) — Validation/scoring workflow documentation
- [`INPUT_FORMATS.md`](INPUT_FORMATS.md) — Input file format specifications
- [`METHODS.md`](METHODS.md) — Statistical methodology details
