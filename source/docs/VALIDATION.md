# Validation and refit workflow

This repo includes an optional validation workflow that:

1) extracts the non-zero variants from a trained elastic net model (`.rds`),
2) checks which variants are present in each validation cohort (by CHR:POS and allele set),
3) builds an **intersection variant set** (e.g., present in all cohorts),
4) refits elastic net **in the training cohort** using only the intersection set,
5) exports cohort-specific genotype subsets and computes refit model scores per cohort,
6) produces QC outputs describing allele flips, strand checks, and dropped ambiguous SNPs.

## Inputs

### Training artifacts
- `TRAIN_MODEL_RDS`: final trained model (`final_EN_model_<label>_stable.rds`)
- `TRAIN_RAW`: training genotype `.raw` used for EN (from `--export A`)
- `TRAIN_PHENO` + `PHENO_NAME`: training phenotype

### Validation cohorts
- A manifest TSV (`COHORTS_TSV`) listing each cohort name and its PLINK prefix (`bfile`).

Example:

```text
cohort    bfile
wgs       /data/WGS/genotypes_prefix
wes       /data/WES/genotypes_prefix
array     /data/ARRAY/genotypes_prefix
```

The pipeline reads `<bfile>.bim` to assess availability and uses PLINK2 to extract/export `.raw`.

## Running on Biowulf (recommended)

1) Copy and edit config:

```bash
cp configs/validation.template.env configs/my_validation.env
nano configs/my_validation.env
```

2) Run:

```bash
bash scripts/run_validation.sh --config configs/my_validation.env
```

## Intersection strategies

- `INTERSECTION_MODE=all`: uses only variants available in *all* listed cohorts (most comparable scores)
- `INTERSECTION_MODE=any`: uses the union across cohorts; each cohort score uses its available subset
- `INTERSECTION_MODE=cohort:wgs`: uses variants available in a single cohort

## Allele harmonization notes

Scores are computed using the **training counted allele** (from the training `.raw` header) as the reference.
For each cohort and variant, the pipeline:

- matches variants by CHR+POS, then confirms allele set agreement using `.bim` alleles where possible
- detects non-ambiguous strand flips (via allele set complement)
- flips genotype dosage by `2 - dosage` when the validation counted allele differs from training (after strand adjustment)

Ambiguous SNP pairs (A/T or C/G) are problematic because strand flips cannot be reliably detected without a reference.
By default `AMBIGUOUS_POLICY=drop` will drop ambiguous variants *only when* the counted alleles disagree between
training and validation.

Outputs include per-cohort QC tables listing which variants were used/flipped/dropped.

## Outputs

Under `OUTDIR/`:

- `00_model/` model weights + parsed CHR/POS
- `01_match/` cohort mapping tables + availability summaries + intersection variant lists
- `02_extract/` cohort-specific extracted `.raw` exports
- `03_refit/` refit model `.rds` + refit weights
- `04_scores/` per-cohort score files + per-cohort QC tables


## Common issues

- **CRLF (Windows line endings)** in config/manifest files can add hidden `^M` characters to PLINK paths. The validation driver strips carriage returns and also supports plain text manifests.
- The cohort manifest can be `.tsv` or `.txt`; it is parsed as two whitespace-separated columns: `cohort` and `bfile`.
- Mapping outputs include `Model_Weight` (from the original model) and scoring uses `Weight_refit` (from the refit model).
