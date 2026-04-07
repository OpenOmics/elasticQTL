# Input Formats

This document describes the expected format for all input files used by elasticQTL.

> **See working examples:** The repository includes example configuration files in `config/examples/` that you can copy and modify for your own analyses.

---

## Training Pipeline Inputs (Step 1)

### Phenotype File (`--pheno`)

**Format:** Tab-delimited or whitespace-delimited text file with header  
**Required columns:**
- `FID` — Family ID
- `IID` — Individual ID  
- `<trait>` — Your phenotype column (name specified by `--pheno-name`)

**Example:**
```
FID   IID   pyroptosis_score
1     1001  0.12
1     1002  -0.40
1     1003  0.85
1     1004  NA
2     2001  -0.15
```

**Missing values:**
- `NA`, `NaN`, blank cells, and `-9` are treated as missing
- Samples with missing phenotype values are automatically excluded from analysis
- No need to manually remove these samples

**Best practices:**
- Use continuous phenotypes (quantitative traits)
- Variance-standardize or rank-normalize phenotypes before running the pipeline
- Include FID and IID exactly as they appear in your PLINK `.fam` file

---

### Covariate File (`--covar`) [Optional]

**Format:** Tab-delimited or whitespace-delimited text file with header  
**Required columns:**
- `FID` — Family ID
- `IID` — Individual ID
- Additional columns are treated as covariates

**Example:**
```
FID   IID   age   sex   PC1      PC2      batch
1     1001  45    1     0.01     -0.02    A
1     1002  52    0     -0.03    0.05     A
1     1003  38    1     0.02     0.01     B
2     2001  61    0     0.00     -0.01    A
```

**Covariates used:**
- Typically: age, sex, principal components (PCs), batch indicators
- PLINK2 will variance-standardize covariates automatically when `--covar-variance-standardize` is used
- Both continuous and categorical (binary) covariates are supported

**Best practices:**
- **Always include covariates** to control for population structure and confounders
- Include at least 5-10 genetic PCs to account for ancestry
- Match FID/IID exactly with phenotype and genotype files

---

### Sample Exclusion List (`--exclude`) [Optional]

**Format:** Whitespace-delimited or tab-delimited text file **without header**  
**Required columns:**
- Column 1: `FID` — Family ID
- Column 2: `IID` — Individual ID

**Example:**
```
1  1009
1  1015
2  2003
2  2010
```

**Use cases:**
- Remove outliers identified in QC
- Exclude duplicates or related samples
- Remove samples failing quality control thresholds

**Best practices:**
- Generate this file from your QC pipeline
- One sample per line
- No header row

---

### Genotype Input (`--bfile`)

**Format:** PLINK binary format  
**Required files:**
- `<prefix>.bed` — Binary genotype data
- `<prefix>.bim` — Variant information (6 columns)
- `<prefix>.fam` — Sample information (6 columns)

**Variant ID format in `.bim` file:**

The 2nd column of the `.bim` file contains variant IDs. These can be:
- **rsIDs:** `rs123456`
- **chr:pos format:** `1:12345`
- **chr:pos:ref:alt format:** `1:12345:A:G` (recommended)

**Example `.bim` snippet:**
```
1  rs123456    0  12345   A  G
1  1:67890:C:T 0  67890   C  T
2  rs789012    0  23456   G  A
```

**Important:**
- Variant IDs must be **consistent across all pipeline steps**
- If validation cohorts use different ID formats, standardize them before running
- See [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) Issue 5 for ID standardization

**Best practices:**
- Use `chr:pos:ref:alt` format for maximum clarity and cross-cohort compatibility
- Ensure your `.bim` allele coding (A1/A2) is consistent
- Filter to QTL candidate variants before running the pipeline (if applicable)

---

### QTL Metadata File (`--qtl-meta`) [Optional]

**Format:** Comma-separated values (CSV) file with header  
**Purpose:** Annotate final model variants with additional metadata (e.g., gene names, QTL type, tissue)

**Example columns:**
```
variant_id,gene_name,qtl_type,tissue,chr,pos
rs123456,GENE1,eQTL,Whole_Blood,1,12345
1:67890:C:T,GENE2,pQTL,Plasma,1,67890
rs789012,GENE3,sQTL,Adipose,2,23456
```

**Required column:**
- A variant identifier column that matches your `.bim` variant IDs

**Optional columns:**
- `gene_name` — Associated gene
- `qtl_type` — Type of QTL (eQTL, pQTL, sQTL, mQTL, etc.)
- `tissue` — Tissue or cell type
- `chr`, `pos` — Genomic coordinates
- Any other metadata you want in the final annotated output

**Best practices:**
- Include all QTL variants you're testing, not just those in the final model
- Use consistent variant IDs with your `.bim` file
- The pipeline will join this metadata to final model output by variant ID

---

## Validation/Scoring Inputs (Step 2)

### Training Model Artifacts

These files are **outputs from Step 1** (training pipeline) that become inputs to Step 2 (validation):

#### `--train-model-rds`

**Format:** R data file (`.rds`)  
**Location:** `<training_outdir>/07_en_stable/final_EN_model_<label>_stable.rds`  
**Contents:** Fitted glmnet elastic net model object

**Example path:**
```
/data/user/training_out/07_en_stable/final_EN_model_p20_stable.rds
```

---

#### `--train-raw`

**Format:** PLINK `.raw` genotype matrix  
**Location:** `<training_outdir>/05_genotypes/ld_variants_forEN.raw`  
**Contents:** Additive dosage matrix (0/1/2) for LD-independent variants used in training

**Example path:**
```
/data/user/training_out/05_genotypes/ld_variants_forEN.raw
```

**Example snippet:**
```
FID IID PAT MAT SEX PHENOTYPE rs123456_A 1:67890:C:T_C rs789012_G
1   1001 0  0   1   1         2          1             0
1   1002 0  0   0   1         1          2             1
```

**Column naming:**
- Variant columns have suffix indicating counted allele (e.g., `rs123456_A` means dosage counts `A` allele)
- This allele information is used for harmonization with validation cohorts

---

#### `--train-pheno` and `--pheno-name`

**Format:** Same phenotype file used in Step 1  
**Purpose:** Required for the optional refit step

Use the exact same file and column name from your training run.

---

### Cohort Manifest (`--cohorts-tsv`)

**Format:** Tab-delimited or whitespace-delimited text file with header  
**Required columns:**
- `cohort` — Short cohort identifier (no spaces)
- `bfile` — PLINK prefix for that cohort's genotype files

> **See example:** `config/examples/validation_cohorts.example.txt` in the repository

**Example:**
```
cohort    bfile
wgs       /data/cohorts/wgs/genotypes_prefix
wes       /data/cohorts/wes/genotypes_prefix
array     /data/cohorts/array/genotypes_prefix
```

**Important:**
- Each `bfile` must point to valid `.bed`, `.bim`, `.fam` files
- Cohort names should be short, descriptive identifiers (used in output filenames)
- **File must end with a newline** (see [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) Issue 13)

**Common mistake:**
```
cohort    bfile
wgs       /data/wgs/genotypes_prefix
array     /data/array/genotypes_prefix[NO NEWLINE HERE - FILE ENDS]
```
This will cause the last cohort (`array`) to be silently skipped!

**Correct format:**
```
cohort    bfile
wgs       /data/wgs/genotypes_prefix
array     /data/array/genotypes_prefix
[NEWLINE HERE]
```

**Best practices:**
- Use absolute paths for `bfile` to avoid path resolution issues
- Test each bfile path before running: `ls ${bfile}.bed ${bfile}.bim ${bfile}.fam`
- Ensure all cohorts use the same genome build as training data
- Check for CRLF line endings (Windows) and convert to LF (Unix) if needed

---

## Common Format Issues

### Issue 1: Line Endings (CRLF vs LF)

**Problem:** Files created on Windows may have CRLF (`\r\n`) line endings, which can cause parsing errors.

**Check for CRLF:**
```bash
file myfile.tsv
# Bad: "ASCII text, with CRLF line terminators"
# Good: "ASCII text"
```

**Fix:**
```bash
dos2unix myfile.tsv
# or
sed -i 's/\r$//' myfile.tsv
```

---

### Issue 2: Missing Final Newline

**Problem:** Files that don't end with a newline can cause the last row to be skipped in bash loops.

**Check:**
```bash
tail -c 1 myfile.tsv | od -An -tx1
# Should show: 0a (hex for newline)
```

**Fix:**
```bash
echo "" >> myfile.tsv
```

---

### Issue 3: Inconsistent Delimiters

**Problem:** Mixing tabs and spaces as delimiters.

**Best practice:**
- Use **tabs** for all TSV files
- Be consistent within each file
- Don't mix tabs and multiple spaces

**Check tab consistency:**
```bash
cat -A myfile.tsv | head
# Tabs show as ^I
```

---

### Issue 4: Header Case Sensitivity

**Problem:** Column names are case-sensitive.

**Correct:**
```
FID   IID   pyroptosis_score
```

**Incorrect:**
```
fid   iid   pyroptosis_score  # lowercase FID/IID may not be recognized
```

**Best practice:** Always use uppercase `FID` and `IID` in PLINK-related files.

---

### Issue 5: Special Characters in Variant IDs

**Problem:** Variant IDs containing spaces, quotes, or special characters can cause parsing errors.

**Safe characters:**
- Alphanumeric: `A-Z`, `a-z`, `0-9`
- Punctuation: `_`, `-`, `:`, `.`

**Avoid:**
- Spaces, tabs, quotes (`"`, `'`), slashes (`/`, `\`)

**Example of safe IDs:**
```
rs123456
1:12345:A:G
chr1_12345_A_G
variant.123
```

---

## Quick Reference Table

| File Type | Format | Header | Delimiter | Required Columns | Optional |
|-----------|--------|--------|-----------|------------------|----------|
| Phenotype | Text | Yes | Tab/space | FID, IID, trait | - |
| Covariate | Text | Yes | Tab/space | FID, IID | + covariates |
| Exclude | Text | No | Tab/space | FID, IID | - |
| Genotype | PLINK binary | - | - | - | - |
| QTL metadata | CSV | Yes | Comma | variant_id | gene, qtl_type, etc. |
| Cohort manifest | Text | Yes | Tab/space | cohort, bfile | - |

---

## Example: Complete Input Setup

> **Tip:** See `config/examples/study.example.env` and `config/examples/validation.example.env` in the repository for complete working examples.

### Training Pipeline

```bash
# 1. Phenotype file: phenotypes.tsv
FID   IID   pyroptosis_score
1     1001  0.12
1     1002  -0.40
1     1003  0.85

# 2. Covariate file: covariates.tsv
FID   IID   age   sex   PC1      PC2
1     1001  45    1     0.01     -0.02
1     1002  52    0     -0.03    0.05
1     1003  38    1     0.02     0.01

# 3. Exclude file: exclude.txt (no header)
1  1009
1  1015

# 4. Genotype files (PLINK binary)
genotypes.bed
genotypes.bim
genotypes.fam

# 5. QTL metadata: qtl_meta.csv
variant_id,gene_name,qtl_type,tissue
rs123456,GENE1,eQTL,Blood
rs789012,GENE2,pQTL,Plasma

# 6. Config file: config/study.env
BFILE=genotypes
PHENO=phenotypes.tsv
PHENO_NAME=pyroptosis_score
COVAR=covariates.tsv
EXCLUDE=exclude.txt
QTL_META=qtl_meta.csv
OUTDIR=training_out
```

### Validation Pipeline

```bash
# 1. Cohort manifest: cohorts.tsv
cohort    bfile
wgs       /data/wgs/genotypes
array     /data/array/genotypes

# 2. Config file: config/validation.env
TRAIN_MODEL_RDS=training_out/07_en_stable/final_EN_model_p20_stable.rds
TRAIN_RAW=training_out/05_genotypes/ld_variants_forEN.raw
TRAIN_PHENO=phenotypes.tsv
PHENO_NAME=pyroptosis_score
COHORTS_TSV=cohorts.tsv
OUTDIR=validation_out
```

---

## Example Files in the Repository

The repository includes working example files that you can use as templates for your own analyses.

### Example Configuration Files

**Location:** `config/examples/`

| File | Purpose |
|------|---------|
| `study.example.env` | Complete training pipeline config with all parameters |
| `validation.example.env` | Complete validation pipeline config with all parameters |
| `validation_cohorts.example.txt` | Example cohort manifest for multi-cohort scoring |

### How to Use Example Files

**For training pipeline:**
```bash
# Copy the example config
cp config/examples/study.example.env config/my_study.env

# Edit with your paths
nano config/my_study.env

# Run the pipeline
pipeline --config config/my_study.env
```

**For validation pipeline:**
```bash
# Copy the example config
cp config/examples/validation.example.env config/my_validation.env

# Copy the cohort manifest example
cp config/examples/validation_cohorts.example.txt cohorts.txt

# Edit both files with your paths
nano config/my_validation.env
nano cohorts.txt

# Run validation
validate --config config/my_validation.env
```

### Template Files vs Example Files

**Template files** (`config/*.template.env`):
- Minimal examples showing required parameters only
- Good for quick starts

**Example files** (`config/examples/*.example.env`):
- Complete examples with optional parameters
- Include comments explaining each option
- Show realistic parameter values

**Recommendation:** Start with templates for simple use cases, use examples when you need to customize parameters.

---

## See Also

- [`PIPELINE.md`](PIPELINE.md) — Training workflow command-line options
- [`VALIDATION.md`](VALIDATION.md) — Validation/scoring workflow options
- [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) — Common format-related errors and fixes
