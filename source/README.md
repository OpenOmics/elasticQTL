# QTL Elastic Net Analysis Pipeline

This repository runs a reproducible pipeline to:

1) run PLINK2 association on a candidate/QTL variant set,  
2) LD-clump variants using PLINK 1.9,  
3) export LD-independent genotypes,  
4) fit a nested cross-validated elastic net model (`glmnet`) **without genotype imputation**,  
5) fit a stability-filtered final model and export an annotated SNP table + MAF.

**No genotype data can be missing in the elastic net step** (glmnet does not accept NAs). This repo provides multiple defensible strategies:

- **Strict complete-case (default):** `--geno 0 --mind 0 --missing-policy error`
- **Drop samples with any remaining NA:** `--missing-policy drop_samples` (still no imputation)
- **Mean imputation:** `--missing-policy mean_impute` (leakage-safe: impute using *training-fold* means within the outer CV loop)
- **Auto:** `--missing-policy auto` chooses between `mean_impute`, `drop_samples`, or `error` based on missingness thresholds:
  - `--auto-impute-max-cell-rate` (default 0.001 = 0.1% missing genotype cells)
  - `--auto-drop-max-sample-rate` (default 0.02 = 2% of samples having ≥1 missing genotype)


---

## Repository layout

```text
configs/
  params.template.env
  example_study.env

scripts/
  run_pipeline.sh
  00_make_keep_list.R
  02_prepare_clump_inputs.R
  06_nested_elastic_net_no_impute.R
  07_fit_stable_model_no_impute.R
  08_annotate_final_model.R
  09_add_maf_from_raw.R
  99_manifest.sh

docs/
  INPUT_FORMATS.md
  METHODS.md
  TROUBLESHOOTING.md
```

---

## Requirements

### Software
- **PLINK 1.9** (for clumping)
- **PLINK 2.0** (for association testing and genotype export)
- **R ≥ 4.0**

### R packages
Required by the pipeline scripts:
- `data.table`
- `glmnet`

Optional (helpful for custom plotting / downstream wrangling, but not required for this pipeline):
- `dplyr`
- `tidyr`
- `ggplot2`
- `optparse`

---

## Input Files

### Required
1. **Genotype data**: PLINK binary format (`.bed/.bim/.fam`) via `--bfile`
2. **Phenotype file**: tab-delimited with columns `FID`, `IID`, and your phenotype column (provided as `--pheno-name`)

### Optional (but commonly used)
3. **Covariate file**: tab-delimited with `FID`, `IID`, and covariates
4. **QTL metadata**: CSV with variant annotations (used for final annotation table)
5. **Exclusion list**: individuals to exclude (two columns: `FID IID`)

See `docs/INPUT_FORMATS.md` for exact formatting expectations and examples.

---

## Biowulf / NIH HPC (modules)

On Biowulf, PLINK and R are typically provided via environment modules. Loading one version of a module can unload another version of the *same* module (e.g., PLINK 1.x vs PLINK 2.x), so this pipeline supports **module switching inside the driver script**.

### Recommended Biowulf settings
In your `configs/my_study.env`:

```bash
USE_MODULES="1"
MODULE_R="R"
MODULE_PLINK2="plink"        # provides plink2 on Biowulf
MODULE_PLINK1="plink/1.9"    # example; set to the PLINK 1.9 module available on your account
PLINK2="plink2"
PLINK1="plink"
RSCRIPT="Rscript"
```

If Biowulf does not provide a PLINK 1.9 module on your system, download a PLINK 1.9 binary to your project space and set:
```bash
MODULE_PLINK1="NONE"
PLINK1="/path/to/your/local/plink"
```

The pipeline will call `module load ...` before each step that needs PLINK1/PLINK2/R so the correct version is active.

---
## Quickstart

### 1) Install R packages
In R:

```r
install.packages(c("data.table", "glmnet"))
# optional:
install.packages(c("dplyr", "tidyr", "ggplot2", "optparse"))
```

### 2) Make scripts executable
From the repo root:

```bash
chmod +x scripts/run_pipeline.sh
chmod +x scripts/99_manifest.sh
chmod +x scripts/*.R
```

### 3) Create a config file
Copy the template:

```bash
cp configs/params.template.env configs/my_study.env
nano configs/my_study.env
```

### 4) Run the full pipeline

```bash
bash scripts/run_pipeline.sh --config configs/my_study.env
```

---

## Usage

### Quick Start (recommended)

```bash
bash scripts/run_pipeline.sh --config configs/my_study.env
```

### Run only some steps

```bash
# run through export only
bash scripts/run_pipeline.sh --config configs/my_study.env --to-step 6

# rerun modeling and downstream only
bash scripts/run_pipeline.sh --config configs/my_study.env --from-step 7 --to-step 9
```

### Dry run (prints commands, does not execute)

```bash
bash scripts/run_pipeline.sh --config configs/my_study.env --dry-run
```

### Overwrite outputs

```bash
bash scripts/run_pipeline.sh --config configs/my_study.env --force
```

---

## Detailed Workflow

The driver `scripts/run_pipeline.sh` runs the steps below.

### Step 1: Build analysis sample list
- removes missing phenotype values
- applies the exclusion list (if provided)
- writes `outdir/00_samples/analysis.keep`

### Step 2: Initial association testing (PLINK2)

```bash
plink2   --bfile <BFILE>   --keep outdir/00_samples/analysis.keep   --pheno <PHENO> --pheno-name <PHENO_NAME>   --covar <COVAR> --covar-variance-standardize   --glm hide-covar   --out outdir/01_glm_qtl/qtl_assoc
```

### Step 3: Prepare clumping input (R)
The script `scripts/02_prepare_clump_inputs.R`:
- filters to additive test if present (`TEST == "ADD"`)
- removes duplicate variant IDs (keeps smallest P)
- writes:
  - `outdir/02_clump_prep/qtl.clump_input.tsv`
  - `outdir/02_clump_prep/qtl.consistent_variants.snplist`

### Step 4: LD clumping (PLINK 1.9)

```bash
plink   --bfile <BFILE>   --keep outdir/00_samples/analysis.keep   --extract outdir/02_clump_prep/qtl.consistent_variants.snplist   --clump outdir/02_clump_prep/qtl.clump_input.tsv   --clump-p1 1.0 --clump-p2 1.0   --clump-r2 0.2 --clump-kb 500   --out outdir/03_clump/qtl_clumped
```

### Step 5: Association testing on LD-independent variants (PLINK2)
Runs PLINK2 `--glm` restricted to the clumped SNP list.

### Step 6: Export genotypes for elastic net (PLINK2) 
The pipeline creates a modeling dataset with missingness filters:

- `plink2 --geno <GENO> --mind <MIND> --make-bed ...`
- then exports dosages: `plink2 --export A`

### Step 7: Nested elastic net (R, `glmnet`)
The script `scripts/06_nested_elastic_net_no_impute.R`:
- runs nested CV (outer folds default 10, inner folds default 10)
- uses `alpha` (default 0.5)
- tests multiple p-value threshold sets (default thresholds: 0.2, 0.3, 0.5 + `all`)
- produces:
  - `nested_EN_summary.tsv`
  - `nested_predictions_<label>.tsv`
  - `nested_stability_<label>.tsv`

### Step 8: Stability-filtered final model (R)
The script `scripts/07_fit_stable_model_no_impute.R`:
- selects variants with selection frequency ≥ `--stability-min-folds` (default 5/10)
- fits a final CV elastic net model on all samples using only stable variants
- saves:
  - `final_EN_model_<label>_stable.rds`
  - `genetic_scores_<label>_stable.tsv`

### Step 9: Annotation + MAF (R)
The scripts `scripts/08_annotate_final_model.R` and `scripts/09_add_maf_from_raw.R`:
- extract non-zero final model weights
- merge association stats, stability frequency, optional QTL metadata
- align coefficient signs to PLINK2 `A1` when possible using `.raw` allele map
- compute MAF from the exported dosages
- write:
  - `final_model_annotated.csv`
  - `final_model_annotated_withMAF.csv`

---

## Output Files

### From PLINK
- `01_glm_qtl/qtl_assoc.*.glm.linear` — initial association results
- `03_clump/qtl_clumped.clumped` — LD-independent variants
- `04_glm_clumped/clumped_assoc.*.glm.linear` — association results for clumped variants
- `05_genotypes/ld_variants_forEN.raw` — genotype matrix for elastic net

### From the R pipeline
- `06_en_nested/nested_EN_summary.tsv` — summary of nested CV results across variant sets
- `06_en_nested/nested_predictions_<label>.tsv` — out-of-fold predictions vs observed
- `06_en_nested/nested_stability_<label>.tsv` — selection frequency across folds
- `07_en_stable/final_EN_model_<label>_stable.rds` — saved glmnet model
- `07_en_stable/genetic_scores_<label>_stable.tsv` — scores on the training cohort
- `08_annotation/final_model_annotated_withMAF.csv` — final annotated SNP table + MAF

### Manifest (reproducibility)
- `manifest/tool_versions.txt` — PLINK/R versions (best-effort)
- `manifest/params_used.txt` — exact parameters used

---

## Analysis options

### P-value thresholds
Controlled via:
- `--p-thresholds` (comma-separated), e.g. `0.2,0.3,0.5`
- `--label-mode percent` (recommended):
  - `0.20 → p20`, `0.50 → p50`
- `--label-mode legacy`:
  - `0.20 → p02`, `0.50 → p05` (**can be confused with 0.05**)

### LD clumping parameters
- `--clump-r2` (default 0.2)
- `--clump-kb` (default 500)

### Elastic net parameters
- `--alpha` (default 0.5)
- `--outer-folds` and `--inner-folds` (default 10/10)

### Missingness handling (no imputation)
- Strict complete-case (recommended):
  - `--geno 0 --mind 0 --missing-policy error`
- Allow some missingness upstream but still do not impute:
  - `--geno 0.01 --mind 0.01 --missing-policy drop_samples`
- Mean imputation (still reproducible; leakage-safe within outer CV):
  - `--geno 0.01 --mind 0.01 --missing-policy mean_impute`
- Auto mode (recommended if you want guardrails):
  - `--geno 0.01 --mind 0.01 --missing-policy auto --auto-impute-max-cell-rate 0.001 --auto-drop-max-sample-rate 0.02`

---

## Interpreting results

### Summary table (`nested_EN_summary.tsv`)
Columns:
- `label` — variant set (`all`, `p20`, `p30`, `p50`, etc.)
- `n_variants` — variants actually used
- `mean_r2` — mean outer-fold R²
- `sd_r2` — SD of outer-fold R²
- `overall_r2` — R² from concatenated out-of-fold predictions

### Stability analysis (`nested_stability_<label>.tsv`)
Columns:
- `Selection_Freq` — number of outer folds (0–10) with non-zero coefficient
- `Mean_Coefficient` — mean coefficient across outer folds

Variants with high selection frequency are more stable.

---

## Validation / refit workflow

This repo includes an optional validation pipeline for applying the trained model to external cohorts (WGS/WES/array)
and (optionally) refitting the model in the training cohort on a variant intersection set.

- Config template: `configs/validation.template.env`
- Cohort manifest example: `configs/validation_cohorts.example.tsv`
- Driver: `scripts/run_validation.sh`
- Documentation: `docs/VALIDATION.md`

Quick start:

```bash
cp configs/validation.template.env configs/my_validation.env
nano configs/my_validation.env
bash scripts/run_validation.sh --config configs/my_validation.env
```

The validation workflow:
1) exports non-zero variants from the trained model,
2) checks availability in each validation cohort (by BIM CHR/POS + allele set),
3) builds an intersection set (e.g., present in all cohorts),
4) refits EN in the training cohort using that set,
5) exports cohort genotype subsets and computes scores,
6) outputs QC tables (allele flips, strand checks, ambiguous SNP drops).

---


### Using the same variants across all validation cohorts

Each cohort score file reports `N_Variants_Used`. This number can differ across cohorts if a variant is:
- absent from a cohort (not present in that cohort’s `.bim`), or
- dropped during harmonization/QC (most commonly **strand-ambiguous A/T or C/G sites** when the counted allele cannot be safely aligned).

If your analysis requires **the exact same set of variants across WGS/WES/array** (recommended for strict comparability), define a **strict common set** as the intersection of variants actually used in scoring, then refit and rescore:

```bash
# Set this tag to match INTERSECTION_MODE in your run:
#   all -> all
#   any -> any
#   cohort:wgs -> cohort_wgs
TAG=all

# Collect variants actually used in each cohort (from QC tables)
awk 'NR>1{print $1}' 04_scores/wgs/wgs_qc_refit_${TAG}.tsv     | sort -u > /tmp/wgs.used
awk 'NR>1{print $1}' 04_scores/wes/wes_qc_refit_${TAG}.tsv     | sort -u > /tmp/wes.used
awk 'NR>1{print $1}' 04_scores/array/array_qc_refit_${TAG}.tsv | sort -u > /tmp/array.used

# Strict common set (variants used in ALL cohorts)
comm -12 /tmp/wgs.used /tmp/wes.used | comm -12 - /tmp/array.used   > 01_match/strict_common_used_variants.txt
wc -l 01_match/strict_common_used_variants.txt

# Refit EN in training cohort on strict common set
Rscript scripts/validation/13_refit_training_intersection.R   "$TRAIN_RAW" "$TRAIN_PHENO" "$PHENO_NAME"   01_match/strict_common_used_variants.txt   03_refit_strict_common   "$ALPHA" "$NFOLDS" "$SEED" "$MISSING_POLICY"   "$AUTO_IMPUTE_MAX_CELL_RATE" "$AUTO_DROP_MAX_SAMPLE_RATE"

# Rescore each cohort using the strict-common refit model
for cohort in wgs wes array; do
  Rscript scripts/validation/14_score_validation_cohort.R     03_refit_strict_common/refit_model.rds     "$TRAIN_RAW"     02_extract/${cohort}/${cohort}_extracted.raw     01_match/${cohort}.mapping.tsv     04_scores/${cohort}/${cohort}_scores_refit_strict_common.csv     04_scores/${cohort}/${cohort}_qc_refit_strict_common.tsv     "$AMBIGUOUS_POLICY" mean_impute
done
```

Alternative (less conservative): set `AMBIGUOUS_POLICY=keep` to avoid dropping strand-ambiguous SNPs, which can increase `N_Variants_Used` but may risk mis-orientation at A/T or C/G sites.

Tip: ensure your cohort manifest ends with a final newline; otherwise the last cohort line can be skipped in bash loops.

---
## Troubleshooting
See `docs/TROUBLESHOOTING.md`.

Common issues:
- variant IDs differ between cohorts (rsIDs vs chr:pos:ref:alt)
- allele counted in exported `.raw` differs from PLINK2 A1 allele
- complete-case filters remove too many samples (adjust `--geno/--mind`)

---

**Last updated**: March 3, 2026

