# Training Pipeline (Step 1)

The training workflow (`pipeline`) trains the elastic net model in your training cohort and estimates performance via nested cross-validation.

The pipeline executes the following stages in order: (1) PLINK2 association testing on candidate/QTL variant set, (2) LD clumping using PLINK 1.9 to identify independent variants, (3) genotype export of LD-independent variants, (4) nested cross-validated elastic net modeling (glmnet) with no genotype imputation, (5) stability filtering to select robust variants across CV folds, (6) final model training on the full dataset, and (7) variant annotation with MAF and metadata.

---

## Synopsis

```bash
$ pipeline [--help] [--config CONFIG] \
      --bfile BFILE \
      --pheno PHENO --pheno-name PHENO_NAME \
      --outdir OUTDIR \
      [--covar COVAR] [--exclude EXCLUDE] [--qtl-meta QTL_META] \
      [--plink1 PLINK1] [--plink2 PLINK2] [--rscript RSCRIPT] \
      [--geno GENO] [--mind MIND] \
      [--missing-policy {error,drop_samples,mean_impute,auto}] \
      [--auto-impute-max-cell-rate RATE] [--auto-drop-max-sample-rate RATE] \
      [--clump-r2 R2] [--clump-kb KB] [--clump-p1 P1] [--clump-p2 P2] \
      [--alpha ALPHA] [--outer-folds N] [--inner-folds N] [--seed SEED] \
      [--p-thresholds THRESHOLDS] [--label-mode {percent,legacy}] \
      [--stability-set LABEL] [--stability-min-folds N] \
      [--from-step STEP] [--to-step STEP] [--force] [--dry-run]
```

The synopsis for the command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user must provide PLINK genotype files via `--bfile`, phenotype data via `--pheno` and `--pheno-name`, and an output directory via `--outdir`.

You can always use the `-h` option for information on a specific command.

---

## Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

### `--bfile BFILE`

**PLINK genotype prefix.**  
**type:** path  

This should be the prefix for PLINK binary format files (.bed/.bim/.fam). The pipeline will look for `BFILE.bed`, `BFILE.bim`, and `BFILE.fam` files.

**Example:** `--bfile /data/$USER/genotypes/cohort1`

### `--pheno PHENO`

**Phenotype file.**  
**type:** file  

A phenotype table with at least three columns: FID (family ID), IID (individual ID), and the trait column specified by `--pheno-name`. The file should be tab-delimited or whitespace-delimited. Samples with missing phenotype values will be automatically excluded from analysis.

**Example:** `--pheno /data/$USER/phenotypes/traits.tsv`

### `--pheno-name PHENO_NAME`

**Phenotype column name.**  
**type:** string  

The column name in the phenotype file containing the trait to model. This should exactly match the column header in the file provided via `--pheno`.

**Example:** `--pheno-name pyroptosis_score`

### `--outdir OUTDIR`

**Output directory.**  
**type:** path  

This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically. All pipeline stages will write their outputs to subdirectories within this location.

**Example:** `--outdir /data/$USER/elasticQTL_out`

---

## Data Input Options

Each of the following arguments are optional, and do not need to be provided.

### `--config CONFIG`

**Configuration file.**  
**type:** file  
**default:** None

Path to a `.env` configuration file containing parameter assignments. This allows you to specify all parameters in a file rather than on the command line. Parameters specified on the command line will override values in the config file.

**Config file format:**
```bash
# Required
BFILE=/path/to/genotypes
PHENO=/path/to/phenotypes.txt
PHENO_NAME=trait_name
OUTDIR=/path/to/output

# Optional
COVAR=/path/to/covariates.txt
EXCLUDE=/path/to/exclude_samples.txt
```

**Example:** `--config config/study.env`

### `--covar COVAR`

**Covariate file.**  
**type:** file  
**default:** None (no covariates)

A covariate table with FID and IID columns plus one or more covariate columns. Covariates are highly recommended to control for population structure, batch effects, and other confounders. The file should be tab-delimited or whitespace-delimited and must include FID and IID columns matching those in the phenotype file.

**Example:** `--covar /data/$USER/covariates/pcs_age_sex.tsv`

### `--exclude EXCLUDE`

**Sample exclusion list.**  
**type:** file  
**default:** None (no exclusions)

A two-column file containing FID and IID pairs of samples to exclude from analysis. This is useful for removing outliers, duplicates, or samples failing QC. The file should be tab-delimited or whitespace-delimited with no header.

**Example:** `--exclude /data/$USER/qc/samples_to_remove.txt`

### `--qtl-meta QTL_META`

**QTL metadata file.**  
**type:** CSV file  
**default:** None (no metadata annotation)

A CSV file containing metadata for variant annotation in the final output. This file is used to annotate the final model output with additional information about each variant (e.g., gene names, functional annotations, tissue-specific expression).

**Example:** `--qtl-meta /data/$USER/annotations/qtl_metadata.csv`

---

## Software Path Options

Each of the following arguments are optional, and do not need to be provided.

### `--plink1 PLINK1`

**Path to PLINK 1.9 executable.**  
**type:** command or path  
**default:** `plink`

Specifies the command or full path to PLINK 1.9, which is used for LD clumping. If PLINK 1.9 is available in your PATH as `plink`, you do not need to specify this option.

**Example:** `--plink1 /usr/local/bin/plink1.9`

### `--plink2 PLINK2`

**Path to PLINK 2.0 executable.**  
**type:** command or path  
**default:** `plink2`

Specifies the command or full path to PLINK 2.0, which is used for association testing and genotype export. If PLINK 2.0 is available in your PATH as `plink2`, you do not need to specify this option.

**Example:** `--plink2 /usr/local/bin/plink2`

### `--rscript RSCRIPT`

**Path to Rscript executable.**  
**type:** command or path  
**default:** `Rscript`

Specifies the command or full path to the Rscript executable used for modeling and annotation steps. If Rscript is available in your PATH, you do not need to specify this option.

**Example:** `--rscript /usr/local/R/bin/Rscript`

---

## Missing Genotype Handling

> **IMPORTANT:** The elastic net step requires complete genotype data — glmnet does not accept NA values. The pipeline provides multiple strategies for ensuring no missing genotypes enter the modeling step.

### `--geno GENO`

**Variant-level missingness filter.**  
**type:** float (0.0–1.0)  
**default:** 0

Maximum missing rate per variant in the modeling set. Variants exceeding this threshold will be removed before elastic net modeling. A value of 0 means no missing genotypes are allowed per variant. A value of 0.02 means variants with >2% missing calls are removed.

**Example:** `--geno 0.02`

### `--mind MIND`

**Sample-level missingness filter.**  
**type:** float (0.0–1.0)  
**default:** 0

Maximum missing rate per sample across the modeling variant set. Samples exceeding this threshold will be removed before elastic net modeling. A value of 0 means no missing genotypes are allowed per sample. A value of 0.05 means samples with >5% missing calls are removed.

**Example:** `--mind 0.05`

### `--missing-policy {error,drop_samples,mean_impute,auto}`

**Strategy for handling any remaining missing genotypes.**  
**type:** enum  
**default:** `error`

After applying `--geno` and `--mind` filters, if any missing genotypes remain in the modeling matrix, this option determines how to proceed:

**`error`** (recommended, strict complete-case)  
The pipeline will error and exit if any NA values are detected. This is the strictest and most conservative approach. It is recommended to combine this with `--geno 0 --mind 0` to ensure complete data.

**`drop_samples`**  
Any samples with one or more missing genotypes across the variant set will be dropped from the analysis. No imputation is performed. This is useful when missingness is concentrated in a small number of samples.

**`mean_impute`**  
Missing genotypes are imputed using the mean dosage of the variant calculated from the training fold within each outer cross-validation loop. This is leakage-safe: held-out test samples do not contribute to the imputation mean. Use this when missingness is sparse and randomly distributed.

**`auto`**  
Automatically selects between `error`, `drop_samples`, and `mean_impute` based on the missingness pattern, using thresholds defined by `--auto-impute-max-cell-rate` and `--auto-drop-max-sample-rate` (see below).

**Example:** `--missing-policy drop_samples`

### `--auto-impute-max-cell-rate RATE`

**Cell-level missingness threshold for auto mode.**  
**type:** float  
**default:** 0.001 (0.1%)

When `--missing-policy auto` is used, this option sets the maximum proportion of missing genotype cells (missing calls / total cells) below which mean imputation is acceptable. If the overall missing rate is ≤0.1%, the auto mode will use mean imputation.

**Example:** `--auto-impute-max-cell-rate 0.002`

### `--auto-drop-max-sample-rate RATE`

**Sample-level missingness threshold for auto mode.**  
**type:** float  
**default:** 0.02 (2%)

When `--missing-policy auto` is used, this option sets the maximum proportion of samples with any missing genotypes below which dropping those samples is acceptable. If ≤2% of samples have at least one missing genotype, the auto mode will drop those samples.

**Example:** `--auto-drop-max-sample-rate 0.05`

---

## LD Clumping Parameters

Each of the following arguments are optional, and do not need to be provided.

### `--clump-r2 R2`

**LD r² threshold for clumping.**  
**type:** float  
**default:** 0.2

Variants in LD with an index variant at r² > this threshold will be clumped (removed). Lower values result in more stringent LD pruning and fewer variants retained.

**Example:** `--clump-r2 0.1`

### `--clump-kb KB`

**Clumping window size.**  
**type:** integer (kilobases)  
**default:** 500

The physical distance window (in kb) within which LD is calculated for clumping. Variants beyond this distance from the index variant are not considered for clumping.

**Example:** `--clump-kb 250`

### `--clump-p1 P1`

**Index variant p-value threshold.**  
**type:** float  
**default:** 1.0

Only variants with association p-value ≤ this threshold are considered as potential index (seed) variants for clumping. A value of 1.0 means all variants can serve as index variants.

**Example:** `--clump-p1 0.05`

### `--clump-p2 P2`

**Clumped variant p-value threshold.**  
**type:** float  
**default:** 1.0

Only variants with association p-value ≤ this threshold are considered for clumping around index variants. A value of 1.0 means all variants can be clumped.

**Example:** `--clump-p2 1.0`

---

## Elastic Net Parameters

Each of the following arguments are optional, and do not need to be provided.

### `--alpha ALPHA`

**Elastic net mixing parameter.**  
**type:** float (0.0–1.0)  
**default:** 0.5

The elastic net mixing parameter passed to glmnet. A value of 0 corresponds to ridge regression (L2 penalty only), 1 corresponds to lasso (L1 penalty only), and values in between use a mix of L1 and L2 penalties. The default of 0.5 balances feature selection (lasso) with handling correlated predictors (ridge).

**Example:** `--alpha 1.0`

### `--outer-folds N`

**Number of outer cross-validation folds.**  
**type:** integer  
**default:** 10

The number of outer folds for nested cross-validation. Performance estimates are computed by averaging across these held-out test sets. Typical values are 5 or 10.

**Example:** `--outer-folds 5`

### `--inner-folds N`

**Number of inner cross-validation folds.**  
**type:** integer  
**default:** 10

The number of inner folds used within each outer fold for hyperparameter tuning (lambda selection) via glmnet's built-in cross-validation. Typical values are 5 or 10.

**Example:** `--inner-folds 5`

### `--seed SEED`

**Random seed for reproducibility.**  
**type:** integer  
**default:** 123

Sets the random seed for cross-validation fold assignments and any stochastic steps. Using the same seed ensures reproducible results across runs.

**Example:** `--seed 42`

### `--p-thresholds THRESHOLDS`

**P-value thresholds for variant sets.**  
**type:** comma-separated floats  
**default:** `0.2,0.3,0.5`

Defines multiple p-value thresholds for creating variant sets to model. For each threshold P, a variant set is created containing all clumped variants with association p-value < P. The pipeline will fit models for each variant set and report performance metrics. An "all variants" set (no p-value filtering) is always included automatically.

**Example:** `--p-thresholds 0.1,0.2,0.3,0.5`

### `--label-mode {percent,legacy}`

**Labeling format for variant sets.**  
**type:** enum  
**default:** `percent`

Controls how variant sets are labeled in output files:

**`percent`**  
Uses percentile-style labels: `p20` for P<0.20, `p50` for P<0.50. This is clearer and less ambiguous.

**`legacy`**  
Uses decimal-style labels: `p02` for P<0.20, `p05` for P<0.50. This format can be ambiguous (p05 could mean 0.05 or 0.5).

**Example:** `--label-mode legacy`

---

## Stability Filter Parameters

Each of the following arguments are optional, and do not need to be provided.

### `--stability-set LABEL`

**Variant set for stability filtering.**  
**type:** string  
**default:** `p20` (when `--label-mode percent`) or `p02` (when `--label-mode legacy`)

Specifies which variant set to use for the stability filtering step. This set is used to identify variants that are consistently selected across outer cross-validation folds. The label should match one of the variant sets defined by `--p-thresholds`.

**Example:** `--stability-set p30`

### `--stability-min-folds N`

**Minimum folds for stability selection.**  
**type:** integer  
**default:** 5

A variant must be selected (have non-zero coefficient) in at least this many outer folds to be included in the final stability-filtered model. With 10 outer folds, a value of 5 means a variant must appear in ≥50% of folds.

**Example:** `--stability-min-folds 7`

---

## Run Control Options

Each of the following arguments are optional, and do not need to be provided.

### `--from-step STEP`

**Starting pipeline step.**  
**type:** integer (1–9)  
**default:** 1

The step number to start from. This is useful for resuming a partially completed run or re-running downstream steps after modifying parameters. Steps are:

1. Make keep list (non-missing phenotype, minus exclusions)
2. PLINK2 GLM on full (QTL candidate) dataset
3. Prepare clump inputs (from GLM)
4. PLINK1 clump
5. PLINK2 GLM on clumped variants
6. Make complete-case modeling bfile + export .raw
7. Nested CV elastic net (no imputation)
8. Stability-filtered final model fit (no imputation)
9. Annotation + MAF table

**Example:** `--from-step 7`

### `--to-step STEP`

**Ending pipeline step.**  
**type:** integer (1–9)  
**default:** 9

The step number to stop at. This is useful for running only the data preparation steps without modeling, or stopping before annotation.

**Example:** `--to-step 6`

### `--force`

**Overwrite existing outputs.**  
**type:** boolean flag  
**default:** off

When this option is provided, the pipeline will overwrite any existing output files from previous runs. Without this flag, the pipeline may skip steps if output files already exist.

**Example:** `--force`

### `--dry-run`

**Preview commands without executing.**  
**type:** boolean flag  
**default:** off

Displays what steps would be run and what commands would be executed. Does not actually execute anything. This is useful for validating your command and parameters before committing to a full run.

**Example:** `--dry-run`

---

## Miscellaneous Options

### `-h, --help`

**Display help message.**  
**type:** boolean flag

Shows the command's synopsis, help message, and an example command.

**Example:** `--help`

---

## Outputs

### Key Artifacts for Step 2 (Validation)

These files from Step 1 are needed for the external scoring/refit workflow:

| Path | Description |
|------|-------------|
| `07_en_stable/final_EN_model_<label>_stable.rds` | Final trained model |
| `05_genotypes/ld_variants_forEN.raw` | Genotype matrix used for elastic net |
| Your original phenotype file | Training phenotype (same file from `--pheno`) |
| Your `--pheno-name` value | Trait column name |

### Model Performance

| Path | Description |
|------|-------------|
| `06_en_nested/nested_EN_summary.tsv` | **Top-level nested CV summary** across all variant sets |
| `06_en_nested/nested_outer_r2_<label>.tsv` | Per-fold held-out R² for each variant set |
| `06_en_nested/nested_predictions_<label>.tsv` | Observed vs. predicted values across all outer folds |
| `06_en_nested/nested_stability_<label>.tsv` | Variant selection frequency and mean coefficient across folds |
| `06_en_nested/n_selected_per_fold_<label>.tsv` | Number of non-zero coefficients per outer fold |
| `06_en_nested/outer_folds_<label>.tsv` | Sample-to-outer-fold assignments |

### Final Model

| Path | Description |
|------|-------------|
| `08_annotation/final_model_annotated_withMAF.csv` | Annotated variant table with MAF |

### QC and Intermediate Files

| Path | Description |
|------|-------------|
| `01_glm_qtl/qtl_assoc.*.glm.linear` | Initial association results |
| `03_clump/qtl_clumped.clumped` | LD-independent variants |
| `06_en_nested/allele_map_from_raw.tsv` | Variant → dosage allele mapping from .raw header |
| `06_en_nested/missingness_report_all_variants.tsv` | Per-run missingness summary (cell rate, sample rate, variant count) |
| `06_en_nested/missingness_policy_decision.tsv` | Which missing-data policy was requested vs. actually applied |
| `manifest/params_used.txt` | Exact parameters used in this run |
| `logs/` | Step-level log files |

---

## Examples

### Example 1: Basic Run with Config File

```bash
# Step 1.) Create a config file
cat > config/study.env << EOF
BFILE=/data/$USER/genotypes/cohort1
PHENO=/data/$USER/phenotypes/traits.tsv
PHENO_NAME=pyroptosis_score
OUTDIR=/data/$USER/elasticQTL_out
COVAR=/data/$USER/covariates/pcs_age_sex.tsv
EOF

# Step 2.) Run the pipeline
pipeline --config config/study.env
```

### Example 2: Command-Line Parameters (No Config File)

```bash
pipeline \
  --bfile /data/$USER/genotypes/cohort1 \
  --pheno /data/$USER/phenotypes/traits.tsv \
  --pheno-name pyroptosis_score \
  --outdir /data/$USER/elasticQTL_out \
  --covar /data/$USER/covariates/pcs_age_sex.tsv \
  --exclude /data/$USER/qc/outliers.txt \
  --missing-policy error \
  --geno 0 --mind 0
```

### Example 3: Dry-Run to Preview Commands

```bash
# Preview what the pipeline will do without executing
pipeline --config config/study.env --dry-run
```

### Example 4: Relaxed Missingness Handling

```bash
# Allow some missingness, drop affected samples
pipeline --config config/study.env \
  --geno 0.02 \
  --mind 0.05 \
  --missing-policy drop_samples
```

### Example 5: Custom Elastic Net Parameters

```bash
# Use pure lasso and test multiple p-value thresholds
pipeline --config config/study.env \
  --alpha 1.0 \
  --p-thresholds 0.1,0.2,0.3,0.5 \
  --outer-folds 5 \
  --inner-folds 5
```

### Example 6: Run Only Data Preparation Steps

```bash
# Run through genotype export only (steps 1–6)
# Useful for validating inputs before running expensive modeling
pipeline --config config/study.env --to-step 6
```

### Example 7: Re-Run Modeling with Different Parameters

```bash
# Skip data prep, re-run modeling with different stability threshold
pipeline --config config/study.env \
  --from-step 7 \
  --stability-min-folds 7 \
  --force
```

### Example 8: Stringent LD Pruning and Stability Filtering

```bash
# More aggressive LD clumping and stability selection
pipeline --config config/study.env \
  --clump-r2 0.1 \
  --clump-kb 250 \
  --stability-min-folds 8
```

### Example 9: Complete Example on Skyline/OpenOmics

```bash
# Step 1.) Grab an interactive node
# Do not run on head node!
srun -N 1 -n 1 --time=8:00:00 --mem=64gb -c 4 --pty bash

# Step 2.) Add pipeline to PATH
export PATH="/data/openomics/prod/elasticQTL/v0.1.0/bin:${PATH}"

# Step 3.) Create config file
cat > config/my_study.env << EOF
BFILE=/data/$USER/genotypes/my_cohort
PHENO=/data/$USER/phenotypes/my_traits.tsv
PHENO_NAME=trait_score
OUTDIR=/data/$USER/results/elasticQTL
COVAR=/data/$USER/covariates/pcs_batch.tsv
MISSING_POLICY=error
GENO=0
MIND=0
EOF

# Step 4.) Dry-run to validate
pipeline --config config/my_study.env --dry-run

# Step 5.) Run the pipeline
pipeline --config config/my_study.env
```

---

## Next Steps

After training completes successfully:

1. **Review performance**: Check `06_en_nested/nested_EN_summary.tsv` for cross-validation R² across variant sets
2. **Inspect final model**: Review `08_annotation/final_model_annotated_withMAF.csv` for selected variants and their coefficients
3. **Proceed to validation**: Use the artifacts from Step 1 to score external cohorts — see [`VALIDATION.md`](VALIDATION.md)

---

## See Also

- [`VALIDATION.md`](VALIDATION.md) — external cohort scoring and refit workflow (Step 2)
- [`INPUT_FORMATS.md`](INPUT_FORMATS.md) — detailed input file format specifications
- [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) — common errors and fixes
