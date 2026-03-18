# Methods (manuscript-ready draft)

## Elastic net genetic score construction (nested cross-validation; no genotype imputation)

Genotype data were restricted to LD-independent QTL variants using LD clumping (PLINK v1.9 `--clump`)
with an LD threshold of r² = 0.2 within 500 kb (clump p1 = 1.0; clump p2 = 1.0).
Association statistics for the LD-independent set were obtained using PLINK2 linear regression (`--glm`)
under an additive genetic model with variance-standardized covariates (`--covar-variance-standardize`).
Individuals failing cohort-specific exclusions were removed prior to analysis.

For elastic net modeling, genotype dosages for the LD-independent variants were exported from PLINK2
(`--export A`) as additive allele counts (0/1/2). To avoid information leakage and to comply with a
no-imputation modeling protocol, we enforced a complete-case genotype matrix prior to model fitting by
excluding variants and individuals with any missing genotype values (PLINK2 `--geno 0 --mind 0` applied
to the modeling variant set). Individuals with missing phenotype values were also excluded. No missing
genotype imputation (e.g., mean imputation) was performed during training or cross-validation.

Predictive models were fit using elastic net regression as implemented in `glmnet`, with the mixing
parameter fixed at α = 0.5 and default feature standardization enabled (`standardize = TRUE`).
Model performance was assessed using nested cross-validation with 10 outer folds and 10 inner folds.
Within each outer training split, the elastic net penalty parameter λ was selected by 10-fold inner
cross-validation (`cv.glmnet`) using the λ_min criterion. Predictions were generated for held-out outer
test folds and aggregated across all outer folds to compute overall predictive performance. Predictive
accuracy was summarized using R² (squared correlation between observed and predicted phenotype) both per
outer fold and overall across out-of-fold predictions.

## Stability analysis and final model fitting

To evaluate feature stability, we recorded the non-zero coefficient set selected at λ_min in each outer
fold. Selection frequency for each variant was defined as the number of outer folds (out of 10) in which
the variant’s coefficient was non-zero. A stability-filtered variant panel was defined by retaining
variants selected in at least 5 of 10 outer folds. A final elastic net model was then fit on the full
(complete-case) cohort using only the stability-filtered variant panel, with λ again chosen by 10-fold
internal cross-validation. Final model coefficients were extracted at λ_min, and individual-level genetic
scores were generated as the model’s predicted values.

## Allele alignment, annotation, and MAF

Effect alleles for association statistics were taken from PLINK2 (`A1` column). For reporting, elastic
net coefficient signs were aligned to the association effect allele when possible by comparing the allele
counted in the exported dosage matrix (encoded in PLINK `.raw` column suffixes) to `A1`; coefficients
were sign-flipped when the counted allele differed from `A1`. Variants were optionally annotated with QTL
type and gene context using an external QTL metadata table. Minor allele frequencies (MAF) were computed
from the training cohort dosage matrix as MAF = min(p̂, 1 − p̂), where p̂ is the mean dosage divided by 2.
