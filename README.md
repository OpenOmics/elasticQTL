# elasticQTL

> A reproducible QTL pipeline for PLINK-based association testing, LD clumping, nested elastic net modeling, stability filtering, and final variant annotation — plus an external scoring/refit workflow for WGS/WES/array cohorts.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Pipeline Components](#2-pipeline-components)
   - 2.1 [Training Workflow](#21-training-workflow-run_pipelinesh)
   - 2.2 [External Scoring/Refit Workflow](#22-external-scoringrefit-workflow-run_validationsh)
   - 2.3 [Missing Genotype Handling](#23-missing-genotype-handling)
3. [Repository Layout](#3-repository-layout)
4. [Requirements](#4-requirements)
5. [Inputs](#5-inputs)
   - 5.1 [Step 1 — Training Inputs](#51-step-1-training-inputs)
   - 5.2 [Step 2 — External Scoring/Refit Inputs](#52-step-2-external-scoringrefit-inputs)
6. [Getting Started](#6-getting-started)
   - 6.1 [Running on Skyline/OpenOmics](#61-running-on-skylineopenomics-install)
   - 6.2 [Clone the Repository](#62-clone-the-repository)
   - 6.3 [Quickstart](#63-quickstart)
7. [Usage](#7-usage)
   - 7.1 [Step 1 — Training](#71-step-1-training)
   - 7.2 [Step 2 — External Scoring/Refit](#72-step-2-external-scoringrefit)
8. [Outputs](#8-outputs)
9. [Validation Notes: Same Markers Across Cohorts](#9-validation-notes-same-markers-across-cohorts)
10. [Troubleshooting](#10-troubleshooting)
11. [Contributing](#11-contributing)

---

## 1. Overview

The pipeline executes the following stages in order:

1. Run **PLINK2 association** on a candidate/QTL variant set
2. **LD-clump** variants using PLINK 1.9
3. **Export** LD-independent genotypes
4. Fit a **nested cross-validated elastic net** model (`glmnet`) — no genotype imputation
5. Fit a **stability-filtered final model** and export an annotated SNP table with MAF

A second workflow is also included to **score external cohorts** (WGS/WES/array) and optionally **refit** the elastic net model on a marker intersection set before scoring.

---

## 2. Pipeline Components

### 2.1 Training Workflow (`run_pipeline.sh`)

**Purpose:** Train the model in the training cohort and estimate performance via nested cross-validation.

| Role | Parameter |
|------|-----------|
| Training genotypes | `BFILE` |
| Training phenotype | `PHENO` + `PHENO_NAME` |
| Optional covariates (highly recommended) | `COVAR` |
| Optional exclusions | `EXCLUDE` |

**Key outputs:**

- `07_en_stable/final_EN_model_<label>_stable.rds` — final trained model
- `05_genotypes/ld_variants_forEN.raw` — training genotype dosages used for EN

---

### 2.2 External Scoring/Refit Workflow (`run_validation.sh`)

**Purpose:** Apply the trained model to external cohorts (WGS/WES/array) to generate predictions/scores, and optionally refit on a shared marker set.

External cohorts to score are identified in `COHORTS_TSV`.

**Key outputs:**

- `04_scores/<cohort>/<cohort>_scores_refit_*.csv` — per-cohort predictions/scores
- `03_refit/refit_model.rds` — optional refit model
- QC tables per cohort describing allele flips and dropped ambiguous sites

---

### 2.3 Missing Genotype Handling

> **Important:** No genotype data can be missing in the elastic net step — `glmnet` does not accept `NA` values.

This repo provides multiple defensible strategies, configurable via `--missing-policy`:

| Policy | Flag | Description |
|--------|------|-------------|
| **Strict complete-case** *(default)* | `--missing-policy error` | Also requires `--geno 0 --mind 0`. Errors on any missingness. |
| **Drop samples** | `--missing-policy drop_samples` | Drops samples with any remaining NA. No imputation. |
| **Mean imputation** | `--missing-policy mean_impute` | Leakage-safe: imputes using *training-fold* means within the outer CV loop. |
| **Auto** | `--missing-policy auto` | Chooses between the above strategies based on missingness thresholds. |

**Auto-mode thresholds:**

- `--auto-impute-max-cell-rate` — default `0.001` (0.1% missing genotype cells → use mean imputation)
- `--auto-drop-max-sample-rate` — default `0.02` (2% of samples with ≥1 missing genotype → drop samples)

---

## 3. Repository Layout

```text
elasticQTL/
├── bin/
│   ├── pipeline          # Entrypoint for training workflow
│   └── validate          # Entrypoint for scoring/refit workflow
│
├── config/
│   └── *.env             # Example site-specific configs
│
├── source/
│   └── scripts/
│       ├── run_pipeline.sh
│       ├── run_validation.sh
│       ├── 00_make_keep_list.R
│       ├── 02_prepare_clump_inputs.R
│       ├── 06_nested_elastic_net_no_impute.R
│       ├── 07_fit_stable_model_no_impute.R
│       ├── 08_annotate_final_model.R
│       ├── 09_add_maf_from_raw.R
│       ├── 99_manifest.sh
│       └── validation/
│           ├── 10_export_model_variants.R
│           ├── 11_match_model_to_cohorts.R
│           ├── 13_refit_training_intersection.R
│           └── 14_score_validation_cohort.R
│
└── docs/
    ├── INPUT_FORMATS.md
    ├── METHODS.md
    ├── TROUBLESHOOTING.md
    └── VALIDATION.md
```

---

## 4. Requirements

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| PLINK | 1.9 | LD clumping |
| PLINK | 2.0 | Association testing and genotype export |
| R | ≥ 4.0 | Modeling and annotation |

### R Packages

**Required:**
- `data.table`
- `glmnet`

**Optional:**
- `dplyr`
- `tidyr`
- `ggplot2`
- `optparse`

> **Note:** On most HPC/container deployments, PLINK and R are provided by the container or are centrally installed. Manual package installation is typically not needed unless running outside a managed environment.

---

## 5. Inputs

### 5.1 Step 1 — Training Inputs

**Required:**

| Parameter | Description |
|-----------|-------------|
| `BFILE` | PLINK genotype prefix (`.bed`/`.bim`/`.fam`) |
| `PHENO` | Phenotype table with columns `FID`, `IID`, and the trait column |
| `PHENO_NAME` | Trait column name (e.g. `pyroptosis_score`) |
| `OUTDIR` | Output directory |

**Optional:**

| Parameter | Description |
|-----------|-------------|
| `COVAR` | Covariates table with `FID`, `IID` + covariate columns |
| `EXCLUDE` | `FID IID` pairs to remove |
| `QTL_META` | Metadata CSV used for annotation |

---

### 5.2 Step 2 — External Scoring/Refit Inputs

**Required:**

| Parameter | Description |
|-----------|-------------|
| `TRAIN_MODEL_RDS` | Model from Step 1 (e.g. `07_en_stable/final_EN_model_p20_stable.rds`) |
| `TRAIN_RAW` | Training `.raw` from Step 1 (e.g. `05_genotypes/ld_variants_forEN.raw`) |
| `TRAIN_PHENO` | Training phenotype file (needed for optional refit) |
| `PHENO_NAME` | Trait column name (same as Step 1) |
| `COHORTS_TSV` | Cohort manifest listing external cohorts to score |
| `OUTDIR` | Validation output directory |

**Cohort manifest format** (`COHORTS_TSV`):

A plain-text TSV (or whitespace-delimited) file with the following columns:

```tsv
cohort    bfile
wgs       /path/to/wgs_prefix
wes       /path/to/wes_prefix
array     /path/to/array_prefix
```

> **Tip:** Ensure the manifest ends with a final newline — otherwise the last cohort may be silently skipped in bash loops.

See [`docs/INPUT_FORMATS.md`](docs/INPUT_FORMATS.md) for exact formatting expectations and examples.

---

## 6. Getting Started


### 6.1 Running on Skyline/OpenOmics Install

The pipeline is installed at:

```
/data/openomics/prod/elasticQTL/v0.1.0/
```

> ⚠️ **Do not run on the head node.** The pipeline will error immediately if attempted.

```bash
# Grab an interactive node
srun -N 1 -n 1 --time=8:00:00 --mem=64gb -c 4 --pty bash

# Add executables to PATH
export PATH="/data/openomics/prod/elasticQTL/v0.1.0/bin:${PATH}"

# View help menus
pipeline -h
validate -h
```

If run on the head node, you will see:

```
############################################################
ERROR: Running this pipeline on a head node is not allowed.
Please submit this as a job to the cluster or run it from
an interactive node. You can grab an interactive node with:

  srun -N 1 -n 1 --time=8:00:00 --mem=64gb -c 4 --pty bash
############################################################
```

---
### 6.2 Clone the Repository

If you are setting up elasticQTL outside of a managed install, clone the repo directly:

```bash
# Clone
git clone https://github.com/your-org/elasticQTL.git
cd elasticQTL

# Make entrypoints executable
chmod +x bin/pipeline bin/validate

# Add to PATH (add this to your ~/.bashrc or ~/.bash_profile to persist)
export PATH="$(pwd)/bin:${PATH}"
```

Then verify the install:

```bash
pipeline -h
validate -h
```
---

### 6.3 Quickstart

```bash
# 1. Copy and edit a training config
cp config/study.template.env config/study.env
nano config/study.env

# 2. Run training
pipeline --config config/study.env

# 3. Copy and edit a validation config
cp config/validation.template.env config/validation.env
nano config/validation.env

# 4. Run external scoring/refit
validate --config config/validation.env
```

---

## 7. Usage

### 7.1 Step 1 — Training

```bash
pipeline --config config/study.env
```

**Common variations** (when calling the underlying script directly):

```bash
# Run through export only (steps 1–6)
bash source/scripts/run_pipeline.sh --config config/study.env --to-step 6

# Rerun modeling and downstream only (steps 7–9)
bash source/scripts/run_pipeline.sh --config config/study.env --from-step 7 --to-step 9

# Dry run (preview steps without executing)
bash source/scripts/run_pipeline.sh --config config/study.env --dry-run

# Overwrite existing outputs
bash source/scripts/run_pipeline.sh --config config/study.env --force
```

---

### 7.2 Step 2 — External Scoring/Refit

```bash
validate --config config/validation.env
```

**Where are predictions written?**

```
OUTDIR/04_scores/<cohort>/<cohort>_scores_refit_*.csv
```

Output columns include: `EN_score_refit`, `N_Variants_Used`, `N_Flipped`, `N_Dropped_Ambiguous`.

**Two-step artifact handoff** — what Step 2 needs from Step 1:

| Step 1 artifact | Used as in Step 2 |
|-----------------|-------------------|
| `07_en_stable/final_EN_model_<label>_stable.rds` | `TRAIN_MODEL_RDS` |
| `05_genotypes/ld_variants_forEN.raw` | `TRAIN_RAW` |
| `PHENO` + `PHENO_NAME` | `TRAIN_PHENO` + `PHENO_NAME` |

---

## 8. Outputs

### Step 1 — Training

| Path | Description |
|------|-------------|
| `01_glm_qtl/qtl_assoc.*.glm.linear` | Initial association results |
| `03_clump/qtl_clumped.clumped` | LD-independent variants |
| `05_genotypes/ld_variants_forEN.raw` | Genotype matrix for elastic net |
| `06_en_nested/allele_map_from_raw.tsv` | Variant → dosage allele mapping extracted from `.raw` header |
| `06_en_nested/missingness_report_all_variants.tsv` | Per-run missingness summary (cell rate, sample rate, variant count) |
| `06_en_nested/missingness_policy_decision.tsv` | Records which missing-data policy was requested vs. actually applied |
| `06_en_nested/outer_folds_<label>.tsv` | Sample-to-outer-fold assignments for each variant set |
| `06_en_nested/n_selected_per_fold_<label>.tsv` | Number of non-zero coefficients per outer fold |
| `06_en_nested/nested_outer_r2_<label>.tsv` | Per-fold held-out R² for each variant set |
| `06_en_nested/nested_predictions_<label>.tsv` | Observed vs. predicted values across all outer folds |
| `06_en_nested/nested_stability_<label>.tsv` | Variant selection frequency and mean coefficient across folds |
| `06_en_nested/nested_EN_summary.tsv` | **Top-level nested CV summary across all variant sets** |
| `07_en_stable/final_EN_model_<label>_stable.rds` | Final trained model |
| `08_annotation/final_model_annotated_withMAF.csv` | Annotated variant table with MAF |
| `manifest/params_used.txt` | Exact parameters used in this run |
| `logs/` | Step-level log files |

### Step 2 — External Scoring/Refit

| Path | Description |
|------|-------------|
| `03_refit/refit_model.rds` | Refit model (if enabled) |
| `04_scores/<cohort>/*scores*.csv` | Per-cohort predictions (`EN_score_refit`) |
| `04_scores/<cohort>/*qc*.tsv` | Per-cohort QC (flips, ambiguous drops) |
| `logs/` | Step-level log files |

---

## 9. Validation Notes: Same Markers Across Cohorts

Each cohort score file reports `N_Variants_Used`. This value can differ across cohorts when a variant is:

- **Absent** from a cohort (not present in that cohort's `.bim`), or
- **Dropped during harmonization/QC** — most commonly strand-ambiguous A/T or C/G sites where the counted allele cannot be safely aligned

**If you require the exact same variant set across WGS/WES/array**, use the *strict common set* procedure:

> Intersect variants actually used in scoring → refit → rescore

**Less conservative alternative:** Set `AMBIGUOUS_POLICY=keep` to skip dropping strand-ambiguous SNPs.

---

## 10. Troubleshooting

See [`docs/TROUBLESHOOTING.md`](docs/TROUBLESHOOTING.md) for full details.

**Common issues:**

- Variant IDs differ between cohorts (e.g. `rsIDs` vs `chr:pos:ref:alt`)
- Allele counted in exported `.raw` differs from PLINK2 A1 allele
- Complete-case filters remove too many samples (adjust `--geno` / `--mind`)
- Cohort manifest issues (Windows `CRLF` line endings, missing final newline)

---

## 11. Contributing

PRs and issues are welcome. When filing a bug report, please include:

- The config used *(redact paths if needed)*
- The relevant `logs/*.log` files
- `manifest/params_used.txt` (Step 1) or the validation logs (Step 2)
