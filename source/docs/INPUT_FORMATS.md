# Input formats

## Phenotype file (--pheno)
Tab-delimited text with header.

Required columns:
- FID
- IID
- <trait> (your --pheno-name)

Example:
FID   IID   pyroptosis_score
1     1001  0.12
1     1002  -0.40

Missing values:
- NA, NaN, blank, -9 are treated as missing and excluded.

## Covariate file (--covar) [optional]
Tab-delimited with header.
Must include FID and IID. Remaining columns are covariates.

Example:
FID IID age sex PC1 PC2
1   1001 45  1   0.01 -0.02

## Exclude file (--exclude) [optional]
Whitespace or tab-delimited, no header.
Two columns: FID IID

Example:
1  1009
1  1015

## Genotype input (--bfile)
PLINK prefix: <prefix>.bed/.bim/.fam

Variant IDs:
- Can be rsIDs or chr:pos:ref:alt
- Pipeline will work with either *as long as IDs are consistent across steps*.

If cohorts differ in ID style, standardize IDs before running (see TROUBLESHOOTING).


## Missing genotypes

The pipeline can be run in strict complete-case mode (no missing genotypes allowed) or in modes that
drop samples with missing genotypes or perform mean imputation (see README). Missing genotypes in `.raw`
exports typically appear as `NA`.
