# Validation and External Scoring Workflow (Step 2)

The validation workflow (`validate`) applies a trained elastic net model to external cohorts (WGS/WES/array) and optionally refits the model on a shared marker set before scoring.

This workflow: (1) extracts non-zero variants from the trained model, (2) checks which variants are present in each validation cohort by CHR:POS and allele set, (3) builds an intersection variant set based on the specified strategy, (4) optionally refits the elastic net model in the training cohort using only the intersection set, (5) exports cohort-specific genotype subsets and computes scores per cohort, and (6) produces QC outputs describing allele flips, strand checks, and dropped ambiguous SNPs.

---

## Synopsis

```bash
$ validate [--help] [--config CONFIG] \
      --train-model-rds TRAIN_MODEL_RDS \
      --train-raw TRAIN_RAW \
      --train-pheno TRAIN_PHENO --pheno-name PHENO_NAME \
      --cohorts-tsv COHORTS_TSV \
      --outdir OUTDIR \
      [--intersection-mode {all,any,cohort:NAME}] \
      [--ambiguous-policy {drop,keep,error}] \
      [--skip-refit] \
      [--plink2 PLINK2] [--rscript RSCRIPT] \
      [--from-step STEP] [--to-step STEP] [--force] [--dry-run]
```

The synopsis for the command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user must provide the trained model artifacts from Step 1 (training workflow) via `--train-model-rds`, `--train-raw`, `--train-pheno`, and `--pheno-name`. External cohorts to score are specified via `--cohorts-tsv`, and results are written to `--outdir`.

You can always use the `-h` option for information on a specific command.

---

## Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

### `--train-model-rds TRAIN_MODEL_RDS`

**Trained elastic net model from Step 1.**  
**type:** file (.rds)

Path to the final trained model RDS file produced by the training workflow (Step 1). This file typically has the naming pattern `final_EN_model_<label>_stable.rds` and is found in the `07_en_stable/` directory of the training output.

The model object contains the fitted glmnet model with non-zero coefficients for selected variants. These variants will be extracted and used as the basis for scoring external cohorts.

**Example:** `--train-model-rds /data/$USER/training_out/07_en_stable/final_EN_model_p20_stable.rds`

### `--train-raw TRAIN_RAW`

**Training genotype matrix from Step 1.**  
**type:** file (.raw)

Path to the genotype matrix (`.raw` format) used for elastic net training in Step 1. This file is produced by PLINK2's `--export A` command and is typically found in the `05_genotypes/` directory of the training output with the name `ld_variants_forEN.raw`.

This file is required to:
- Determine the counted allele for each variant in the training data
- Match training variants to validation cohort variants
- Refit the model on the intersection set (if refitting is enabled)

**Example:** `--train-raw /data/$USER/training_out/05_genotypes/ld_variants_forEN.raw`

### `--train-pheno TRAIN_PHENO`

**Training phenotype file.**  
**type:** file

The same phenotype file used in Step 1 (training workflow). This file must contain at least FID, IID, and the trait column specified by `--pheno-name`. It is required for the optional refit step. If `--skip-refit` is provided, this file is not used but must still be specified.

**Example:** `--train-pheno /data/$USER/phenotypes/traits.tsv`

### `--pheno-name PHENO_NAME`

**Phenotype column name.**  
**type:** string

The same trait column name used in Step 1. This should match the column header in the file provided via `--train-pheno`. Required for the optional refit step.

**Example:** `--pheno-name pyroptosis_score`

### `--cohorts-tsv COHORTS_TSV`

**Cohort manifest file.**  
**type:** TSV or whitespace-delimited file

A plain-text file listing external cohorts to score. The file must contain at least two columns: `cohort` (cohort name) and `bfile` (PLINK prefix for that cohort's genotypes). The file should have a header row.

**File format:**
```
cohort    bfile
wgs       /data/WGS/genotypes_prefix
wes       /data/WES/genotypes_prefix
array     /data/ARRAY/genotypes_prefix
```

Each `bfile` should be a PLINK prefix pointing to `.bed`, `.bim`, and `.fam` files. The pipeline reads `<bfile>.bim` to assess variant availability and uses PLINK2 to extract and export genotypes.

> **IMPORTANT:** Ensure the manifest ends with a final newline character. Otherwise, the last cohort may be silently skipped in bash loops.

**Example:** `--cohorts-tsv /data/$USER/validation/cohorts.tsv`

### `--outdir OUTDIR`

**Output directory.**  
**type:** path

This location is where the validation workflow will create all of its output files. If the provided output directory does not exist, it will be created automatically. All workflow stages will write their outputs to subdirectories within this location.

**Example:** `--outdir /data/$USER/validation_out`

---

## Configuration Options

Each of the following arguments are optional, and do not need to be provided.

### `--config CONFIG`

**Configuration file.**  
**type:** file  
**default:** None

Path to a `.env` configuration file containing parameter assignments. This allows you to specify all parameters in a file rather than on the command line. Parameters specified on the command line will override values in the config file.

**Config file format:**
```bash
# Required
TRAIN_MODEL_RDS=/data/training_out/07_en_stable/final_EN_model_p20_stable.rds
TRAIN_RAW=/data/training_out/05_genotypes/ld_variants_forEN.raw
TRAIN_PHENO=/data/phenotypes/traits.tsv
PHENO_NAME=pyroptosis_score
COHORTS_TSV=/data/validation/cohorts.tsv
OUTDIR=/data/validation_out

# Optional
INTERSECTION_MODE=all
AMBIGUOUS_POLICY=drop
```

**Example:** `--config config/validation.env`

---

## Intersection Strategy Options

These options control how the variant set is defined across cohorts.

### `--intersection-mode {all,any,cohort:NAME}`

**Intersection strategy for multi-cohort variant sets.**  
**type:** enum  
**default:** `all`

Determines which variants are used for scoring across multiple cohorts:

**`all`** (strict common set)  
Uses only variants available in **all** listed cohorts. This is the most conservative approach and produces the most directly comparable scores across cohorts, as every cohort is scored using exactly the same variant set. However, this can substantially reduce the number of variants used if cohorts have different platforms (e.g., WGS vs. array).

**`any`** (union set)  
Uses the union of variants available across **any** cohort. Each cohort is scored using its available subset of model variants. This maximizes the number of variants used per cohort but means different cohorts may be scored on different variant sets, making cross-cohort comparisons less direct.

**`cohort:NAME`** (single cohort set)  
Uses only variants available in the specified cohort. Replace `NAME` with the cohort name from the `cohort` column in your manifest file (e.g., `cohort:wgs`). All cohorts are then scored using this single cohort's available variant set.

The chosen intersection set determines:
- Which variants are used for the optional refit step
- Which variants contribute to scores in each cohort (subject to per-cohort availability)

**Example:** `--intersection-mode all`  
**Example:** `--intersection-mode cohort:wgs`

---

## Allele Harmonization Options

These options control how allele mismatches are handled between training and validation cohorts.

### `--ambiguous-policy {drop,keep,error}`

**Strategy for handling strand-ambiguous SNPs.**  
**type:** enum  
**default:** `drop`

Controls how the pipeline handles strand-ambiguous SNP pairs (A/T or C/G) when the counted allele differs between training and validation cohorts.

**Background:** For ambiguous SNPs, strand flips cannot be reliably detected without a reference genome because the complement of A/T is T/A, and the complement of C/G is G/C. When the training counted allele is A and the validation counted allele is T, it's unclear whether this represents a strand flip or a genuine allele swap.

**`drop`** (conservative, recommended)  
Drops ambiguous variants only when the counted allele disagrees between training and validation after attempting strand flip detection. Ambiguous variants with matching counted alleles are retained. This prevents incorrect dosage flips from strand ambiguity.

**`keep`** (permissive, use with caution)  
Retains all ambiguous variants even when counted alleles disagree. The pipeline will attempt to harmonize by flipping dosages when needed, but this can introduce errors for strand-ambiguous sites.

**`error`** (strict)  
Errors and exits if any ambiguous variants are encountered with allele disagreements. Use this when you need absolute certainty about allele coding.

Per-cohort QC tables report which variants were dropped due to ambiguity.

**Example:** `--ambiguous-policy drop`

---

## Refit Options

These options control whether to refit the model on the intersection set.

### `--skip-refit`

**Skip the refit step.**  
**type:** boolean flag  
**default:** off (refit is performed by default)

When this option is provided, the pipeline will skip the model refit step and score cohorts using the original model weights from Step 1. By default, the pipeline refits the elastic net model in the training cohort using only the intersection variant set determined by `--intersection-mode`. 

If refit is skipped, scores are computed using the `Model_Weight` column from the original model. If refit is performed, scores use the `Weight_refit` column.

**Example:** `--skip-refit`

---

## Software Path Options

Each of the following arguments are optional, and do not need to be provided.

### `--plink2 PLINK2`

**Path to PLINK 2.0 executable.**  
**type:** command or path  
**default:** `plink2`

Specifies the command or full path to PLINK 2.0, which is used for variant extraction and genotype export. If PLINK 2.0 is available in your PATH as `plink2`, you do not need to specify this option.

**Example:** `--plink2 /usr/local/bin/plink2`

### `--rscript RSCRIPT`

**Path to Rscript executable.**  
**type:** command or path  
**default:** `Rscript`

Specifies the command or full path to the Rscript executable used for model parsing, matching, refitting, and scoring steps. If Rscript is available in your PATH, you do not need to specify this option.

**Example:** `--rscript /usr/local/R/bin/Rscript`

---

## Run Control Options

Each of the following arguments are optional, and do not need to be provided.

### `--from-step STEP`

**Starting workflow step.**  
**type:** integer (1–5)  
**default:** 1

The step number to start from. This is useful for resuming a partially completed run or re-running downstream steps after modifying parameters. Steps are:

1. Export model variants and parse training .raw alleles
2. Match model variants to cohorts and build intersection set
3. Extract cohort-specific genotype subsets
4. Refit model on intersection set (skipped if `--skip-refit`)
5. Score cohorts and produce QC tables

**Example:** `--from-step 3`

### `--to-step STEP`

**Ending workflow step.**  
**type:** integer (1–5)  
**default:** 5

The step number to stop at. This is useful for running only variant matching without scoring, or stopping before the refit step.

**Example:** `--to-step 2`

### `--force`

**Overwrite existing outputs.**  
**type:** boolean flag  
**default:** off

When this option is provided, the workflow will overwrite any existing output files from previous runs. Without this flag, the workflow may skip steps if output files already exist.

**Example:** `--force`

### `--dry-run`

**Preview commands without executing.**  
**type:** boolean flag  
**default:** off

Displays what steps would be run and what commands would be executed. Does not actually execute anything. This is useful for validating your command and parameters before committing to a full run.

**Example:** `--dry-run`

---

## Understanding Intersection Strategies

The intersection strategy (`--intersection-mode`) determines which variants from the trained model are used for scoring. This is critical when validation cohorts have different platforms or genotyping densities.

### Example Scenario

**Training model variants:** 100 SNPs  
**WGS cohort:** Has all 100 SNPs  
**WES cohort:** Has 80 SNPs (20 are outside exome targets)  
**Array cohort:** Has 60 SNPs (40 not on the array)

**Intersection mode outcomes:**

| Mode | Variants Used | WGS Scored On | WES Scored On | Array Scored On | Notes |
|------|---------------|---------------|---------------|-----------------|-------|
| `all` | 60 (common to all) | 60 SNPs | 60 SNPs | 60 SNPs | Most comparable, but smallest set |
| `any` | 100 (union) | 100 SNPs | 80 SNPs | 60 SNPs | Each cohort uses its available subset |
| `cohort:wgs` | 100 (WGS set) | 100 SNPs | 80 SNPs | 60 SNPs | Same as `any` if WGS has all variants |
| `cohort:array` | 60 (array set) | 60 SNPs | 60 SNPs | 60 SNPs | Same as `all` if array is the limiting platform |

### How refit interacts with intersection mode

If refit is enabled (default), the pipeline:
1. Identifies the intersection set based on `--intersection-mode`
2. Refits the elastic net model in the training cohort using only the intersection variants
3. Scores each cohort using the refit model weights

If refit is skipped (`--skip-refit`), the pipeline uses the original model weights but still respects the intersection strategy for determining which variants contribute to scores.

---

## Understanding Allele Harmonization

Scores are computed using the training counted allele (from the training `.raw` header) as the reference. For each cohort and variant, the pipeline:

### 1. Match by position
Variants are matched between training and validation cohorts using CHR:POS coordinates.

### 2. Confirm allele set agreement
The pipeline checks whether the two alleles in the validation `.bim` file match the two alleles from the training data (considering both possible orderings).

### 3. Detect strand flips
For non-ambiguous variants, the pipeline detects strand flips by checking if the validation allele set is the complement of the training allele set:
- A/G ↔ T/C (strand flip)
- A/C ↔ T/G (strand flip)
- But NOT: A/T ↔ A/T (ambiguous)
- And NOT: C/G ↔ C/G (ambiguous)

### 4. Flip genotype dosages when needed
When the validation counted allele (A1 in the `.bim`) differs from the training counted allele after strand adjustment, genotype dosages are flipped using: `dosage_flipped = 2 - dosage_original`

This ensures that the dosage always represents the same allele (the training counted allele) across all cohorts.

### Example

**Training data:**
- Variant: `rs123` at chr1:1000
- Alleles: A/G (counted allele = A)
- Dosages: 0, 1, 2 (counts of A allele)

**Validation cohort 1 (same strand, same coding):**
- `.bim`: chr1:1000 A G
- A1 (counted) = A → **No flip needed**
- Dosages used as-is: 0, 1, 2

**Validation cohort 2 (same strand, opposite coding):**
- `.bim`: chr1:1000 G A
- A1 (counted) = G → **Flip needed**
- Dosages flipped: 2, 1, 0 (now counts of A allele)

**Validation cohort 3 (strand flip, same coding after flip):**
- `.bim`: chr1:1000 T C (complement of A/G)
- Alleles T/C = complement of A/G → strand flip detected
- A1 (counted) = T, which is complement of A
- After accounting for strand flip, counted allele = A → **No flip needed**
- Dosages used as-is: 0, 1, 2

### Ambiguous SNPs

For A/T or C/G SNPs:
- Complement of A/T is T/A (same allele set)
- Complement of C/G is G/C (same allele set)
- Strand flips cannot be distinguished from allele coding differences

**Example of ambiguity:**
- Training: chr1:2000 A/T (counted = A)
- Validation: chr1:2000 T/A (counted = T)

Is this:
- (a) Same strand, opposite coding → flip dosages
- (b) Opposite strand, same coding → don't flip dosages

**Without a reference genome, we cannot tell.** This is why `--ambiguous-policy drop` is the default.

---

## Outputs

All outputs are written to subdirectories within `OUTDIR/`:

### Model Artifacts (`00_model/`)

| File | Description |
|------|-------------|
| `model_weights.tsv` | Parsed model weights with variant IDs and coefficients |
| `model_variants_chrpos.tsv` | Model variants with CHR:POS coordinates extracted |

### Matching and Intersection (`01_match/`)

| File | Description |
|------|-------------|
| `<cohort>_variant_mapping.tsv` | Per-cohort mapping of model variants to cohort variants |
| `<cohort>_availability_summary.tsv` | Summary of how many model variants are available in the cohort |
| `intersection_variants_<mode>.tsv` | Final intersection variant list based on `--intersection-mode` |

### Extracted Genotypes (`02_extract/`)

| File | Description |
|------|-------------|
| `<cohort>_extracted.raw` | Extracted genotype dosages for intersection variants |
| `<cohort>_extract.log` | PLINK2 log from extraction |

### Refit Model (`03_refit/`)

| File | Description |
|------|-------------|
| `refit_model.rds` | Elastic net model refit on intersection set |
| `refit_weights.tsv` | Refit model weights |

> **Note:** This directory is only created if refit is performed (default). If `--skip-refit` is used, this directory will not exist.

### Scores and QC (`04_scores/`)

| File | Description |
|------|-------------|
| `<cohort>/<cohort>_scores_refit.csv` | **Per-cohort prediction scores** (if refit performed) |
| `<cohort>/<cohort>_scores_original.csv` | Per-cohort scores using original model weights (if refit skipped) |
| `<cohort>/<cohort>_qc_report.tsv` | **Per-cohort QC table** listing variants used, flipped, and dropped |

**Score file columns:**
- `FID`, `IID` — Sample identifiers
- `EN_score_refit` — Elastic net score from refit model (or `EN_score_original` if refit skipped)
- `N_Variants_Used` — Number of variants contributing to this sample's score
- `N_Flipped` — Number of variants with flipped dosages
- `N_Dropped_Ambiguous` — Number of ambiguous variants dropped (if applicable)

**QC file columns:**
- `Variant_ID` — Variant identifier
- `CHR`, `POS` — Chromosome and position
- `Training_Counted_Allele` — Allele counted in training data
- `Cohort_Counted_Allele` — Allele counted in validation cohort
- `Strand_Flip_Detected` — Whether a strand flip was detected (TRUE/FALSE)
- `Dosage_Flipped` — Whether dosages were flipped for this variant (TRUE/FALSE)
- `Ambiguous` — Whether variant is strand-ambiguous (A/T or C/G)
- `Status` — USED, FLIPPED, DROPPED_AMBIGUOUS, or DROPPED_MISSING

### Logs (`logs/`)

| File | Description |
|------|-------------|
| `step_<N>_<name>.log` | Log file for each workflow step |

---

## Examples

### Example 1: Basic Run with Config File

```bash
# Step 1.) Create a validation config file
cat > config/validation.env << EOF
TRAIN_MODEL_RDS=/data/$USER/training_out/07_en_stable/final_EN_model_p20_stable.rds
TRAIN_RAW=/data/$USER/training_out/05_genotypes/ld_variants_forEN.raw
TRAIN_PHENO=/data/$USER/phenotypes/traits.tsv
PHENO_NAME=pyroptosis_score
COHORTS_TSV=/data/$USER/validation/cohorts.tsv
OUTDIR=/data/$USER/validation_out
EOF

# Step 2.) Run the validation workflow
validate --config config/validation.env
```

### Example 2: Command-Line Parameters (No Config File)

```bash
validate \
  --train-model-rds /data/$USER/training_out/07_en_stable/final_EN_model_p20_stable.rds \
  --train-raw /data/$USER/training_out/05_genotypes/ld_variants_forEN.raw \
  --train-pheno /data/$USER/phenotypes/traits.tsv \
  --pheno-name pyroptosis_score \
  --cohorts-tsv /data/$USER/validation/cohorts.tsv \
  --outdir /data/$USER/validation_out
```

### Example 3: Strict Common Variant Set (Comparable Scores)

```bash
# Use only variants present in ALL cohorts
# Best for cross-cohort comparisons
validate --config config/validation.env \
  --intersection-mode all \
  --ambiguous-policy drop
```

### Example 4: Maximize Variants Per Cohort

```bash
# Use all available variants in each cohort
# Best for within-cohort prediction accuracy
validate --config config/validation.env \
  --intersection-mode any
```

### Example 5: Use WGS Variant Set as Reference

```bash
# Score all cohorts using only WGS-available variants
validate --config config/validation.env \
  --intersection-mode cohort:wgs
```

### Example 6: Skip Refit, Use Original Model

```bash
# Score cohorts using original model weights
# No refit step performed
validate --config config/validation.env \
  --skip-refit
```

### Example 7: Permissive Ambiguous SNP Handling

```bash
# Keep ambiguous SNPs even when alleles disagree
# Use with caution - may introduce errors
validate --config config/validation.env \
  --ambiguous-policy keep
```

### Example 8: Dry-Run to Preview Commands

```bash
# Preview what the workflow will do
validate --config config/validation.env --dry-run
```

### Example 9: Run Matching Only (No Scoring)

```bash
# Run through intersection building, then stop
# Useful for inspecting variant overlap before scoring
validate --config config/validation.env --to-step 2
```

### Example 10: Re-Run Scoring with Different Parameters

```bash
# Skip variant matching, re-run scoring only
validate --config config/validation.env \
  --from-step 4 \
  --skip-refit \
  --force
```

### Example 11: Complete Example on Skyline/OpenOmics

```bash
# Step 1.) Grab an interactive node
# Do not run on head node!
srun -N 1 -n 1 --time=4:00:00 --mem=32gb -c 2 --pty bash

# Step 2.) Add to PATH
export PATH="/data/openomics/prod/elasticQTL/v0.1.0/bin:${PATH}"

# Step 3.) Create cohort manifest
cat > cohorts.tsv << EOF
cohort    bfile
wgs       /data/$USER/cohorts/wgs_genotypes
wes       /data/$USER/cohorts/wes_genotypes
array     /data/$USER/cohorts/array_genotypes
EOF

# Step 4.) Create validation config
cat > config/my_validation.env << EOF
TRAIN_MODEL_RDS=/data/$USER/training/07_en_stable/final_EN_model_p20_stable.rds
TRAIN_RAW=/data/$USER/training/05_genotypes/ld_variants_forEN.raw
TRAIN_PHENO=/data/$USER/phenotypes/traits.tsv
PHENO_NAME=trait_score
COHORTS_TSV=cohorts.tsv
OUTDIR=/data/$USER/validation_results
INTERSECTION_MODE=all
AMBIGUOUS_POLICY=drop
EOF

# Step 5.) Dry-run to validate
validate --config config/my_validation.env --dry-run

# Step 6.) Run the validation workflow
validate --config config/my_validation.env
```

---

## Common Issues and Troubleshooting

### Issue: Last cohort in manifest is skipped

**Cause:** The cohort manifest file is missing a final newline character.

**Solution:** Ensure your manifest file ends with a newline:
```bash
# Check if file ends with newline
tail -c 1 cohorts.tsv | od -An -tx1

# Add newline if needed
echo "" >> cohorts.tsv
```

### Issue: PLINK paths have `^M` characters

**Cause:** Windows-style line endings (CRLF) in the manifest or config file.

**Solution:** Convert to Unix line endings:
```bash
dos2unix cohorts.tsv
# or
sed -i 's/\r$//' cohorts.tsv
```

The validation driver automatically strips carriage returns, but it's best to fix the source files.

### Issue: Variant counts differ across cohorts

**Cause:** Different cohorts have different genotyping platforms or coverage.

**Expected behavior:** Each cohort score file reports `N_Variants_Used`. This value can differ across cohorts when:
- A variant is absent from a cohort (not in that cohort's `.bim`)
- A variant is dropped during harmonization (e.g., ambiguous SNP with allele disagreement)

**Solution:** Review per-cohort QC files in `04_scores/<cohort>/<cohort>_qc_report.tsv` to see which variants were used/dropped per cohort.

If you require exactly the same variant set across all cohorts, use `--intersection-mode all`.

### Issue: Scores use `Model_Weight` instead of `Weight_refit`

**Cause:** You skipped the refit step with `--skip-refit`, or the refit step failed.

**Solution:** 
- Check logs in `logs/` to confirm refit completed successfully
- If you want to use refit weights, re-run without `--skip-refit`
- If you intentionally skipped refit, this is expected behavior

### Issue: Many variants dropped as ambiguous

**Cause:** Your model includes many A/T or C/G SNPs, and validation cohorts code these alleles differently.

**Solutions:**
- Accept the reduction in variants (recommended with `--ambiguous-policy drop`)
- Use `--ambiguous-policy keep` if you're confident about strand consistency
- Realign validation cohorts to a reference genome and re-export with consistent strand coding

### Issue: Allele flip counts seem wrong

**Cause:** The validation cohort may use different allele coding or strand than the training data.

**Solution:** Review the QC report for each cohort (`04_scores/<cohort>/<cohort>_qc_report.tsv`). The `Strand_Flip_Detected` and `Dosage_Flipped` columns show exactly what harmonization was performed per variant.

---

## Requirements

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| PLINK | 2.0 | Genotype extraction and export |
| R | ≥ 4.0 | Model parsing, matching, refitting, scoring |

### R Packages

**Required:**
- data.table
- glmnet

**Optional:**
- dplyr
- tidyr

---

## Next Steps

After validation completes successfully:

1. **Review score distributions**: Check `04_scores/<cohort>/<cohort>_scores_*.csv` for prediction distributions
2. **Inspect QC reports**: Review `04_scores/<cohort>/<cohort>_qc_report.tsv` to understand variant harmonization
3. **Compare across cohorts**: If using `--intersection-mode all`, scores are directly comparable
4. **Downstream analysis**: Use scores for association testing, risk stratification, or other analyses

---

## See Also

- [`PIPELINE.md`](PIPELINE.md) — training workflow (Step 1) that produces the required model artifacts
- [`INPUT_FORMATS.md`](INPUT_FORMATS.md) — detailed cohort manifest format specifications
- [`METHODS.md`](METHODS.md) — statistical methodology for refitting and scoring
- [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) — additional common errors and fixes
