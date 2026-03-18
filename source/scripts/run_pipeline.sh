\
#!/usr/bin/env bash
set -euo pipefail

# --------------------------
# Helper functions
# --------------------------
die() { echo "ERROR: $*" >&2; exit 1; }

run_cmd() {
  echo "+ $*"
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    eval "$@"
  fi
}

usage() {
  cat <<'EOF'
USAGE:
  bash scripts/run_pipeline.sh [--config configs/my.env] \
    --bfile <plink_prefix> \
    --pheno <phenotype.tsv> --pheno-name <trait_col> \
    --outdir <output_dir> \
    [options...]

REQUIRED:
  --bfile          PLINK prefix for bed/bim/fam
  --pheno          Phenotype file with at least FID IID <trait>
  --pheno-name     Column name of phenotype (e.g. pyroptosis_score)
  --outdir         Output directory for all results

COMMON OPTIONS:
  --covar          Covariate file (FID IID + covars). Optional.
  --exclude        FID IID list to remove. Optional; default NONE.
  --qtl-meta       QTL metadata CSV (optional). default NONE.
  --plink1         Path/command for plink1 (default: plink)
  --plink2         Path/command for plink2 (default: plink2)
  --rscript        Path/command for Rscript (default: Rscript)

MODULES (HPC; optional):
  --use-modules            If set, pipeline runs 'module load' before each step
  --module-init <path|AUTO> Modules init script path (or AUTO; default AUTO)
  --module-r <name>         e.g. R or R/4.3 (default NONE)
  --module-plink1 <name>    e.g. plink/1.9 for clumping (default NONE)
  --module-plink2 <name>    e.g. plink for plink2 steps (default NONE)

MISSINGNESS POLICY (NO IMPUTATION):
  --geno           Variant missingness filter for modeling set (default 0)
  --mind           Sample missingness filter for modeling set (default 0)
  --missing-policy error|drop_samples|mean_impute|auto  (default: error)
     error: stop if any NA in genotype matrix (recommended; strict complete-case)
     drop_samples: drop any sample with NA across the "all variants" set (still no imputation)
     mean_impute: leakage-safe mean imputation within outer folds (training-fold means)
     auto: choose mean_impute vs drop_samples vs error based on missingness thresholds

  --auto-impute-max-cell-rate   Used when --missing-policy auto (default 0.001)
  --auto-drop-max-sample-rate   Used when --missing-policy auto (default 0.02)

CLUMPING:
  --clump-r2       (default 0.2)
  --clump-kb       (default 500)
  --clump-p1       (default 1.0)
  --clump-p2       (default 1.0)

ELASTIC NET:
  --alpha          glmnet alpha (default 0.5)
  --outer-folds    default 10
  --inner-folds    default 10
  --seed           default 123
  --p-thresholds   comma-separated P thresholds for EN sets (default "0.2,0.3,0.5")
                  Always includes "all" set automatically.
  --label-mode     percent|legacy (default percent)
                  percent => p20 for P<0.20, p50 for P<0.50
                  legacy  => p02 for P<0.20, p05 for P<0.50 (can be ambiguous vs 0.05)

STABILITY FILTER:
  --stability-set        label of set to use for stability filtering (default p20 if label-mode percent; p02 if legacy)
  --stability-min-folds  selected in >= this many outer folds (default 5)

RUN CONTROL:
  --from-step      1..9  (default 1)
  --to-step        1..9  (default 9)
  --force          overwrite existing outputs (default off)
  --dry-run        print commands only

STEPS:
  1 = make keep list (non-missing phenotype, minus exclusions)
  2 = PLINK2 GLM on full (QTL candidate) dataset
  3 = prepare clump inputs (from GLM)
  4 = PLINK1 clump
  5 = PLINK2 GLM on clumped variants
  6 = make complete-case modeling bfile + export .raw
  7 = nested CV elastic net (NO IMPUTATION)
  8 = stability-filtered final model fit (NO IMPUTATION)
  9 = annotation + MAF table

EXAMPLE:
  bash scripts/run_pipeline.sh \
    --bfile /path/to/data_prefix \
    --pheno /path/to/pheno.tsv --pheno-name trait \
    --covar /path/to/covar.tsv \
    --exclude /path/to/exclude.tsv \
    --qtl-meta /path/to/qtl_meta.csv \
    --outdir /path/to/out \
    --p-thresholds 0.2,0.3,0.5 \
    --alpha 0.5 --outer-folds 10 --inner-folds 10 \
    --geno 0 --mind 0 \
    --missing-policy error
EOF
}

# --------------------------
# Defaults
# --------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONFIG=""
BFILE=""
PHENO=""
PHENO_NAME=""
COVAR=""
EXCLUDE="NONE"
QTL_META="NONE"
OUTDIR=""

PLINK1="plink"
PLINK2="plink2"
RSCRIPT="Rscript"

USE_MODULES="0"
MODULE_INIT="AUTO"
MODULE_R="NONE"
MODULE_PLINK1="NONE"
MODULE_PLINK2="NONE"

GENO="0"
MIND="0"
MISSING_POLICY="error"
AUTO_IMPUTE_MAX_CELL_RATE="0.001"
AUTO_DROP_MAX_SAMPLE_RATE="0.02"

CLUMP_R2="0.2"
CLUMP_KB="500"
CLUMP_P1="1.0"
CLUMP_P2="1.0"

ALPHA="0.5"
OUTER_FOLDS="10"
INNER_FOLDS="10"
SEED="123"
P_THRESHOLDS="0.2,0.3,0.5"
LABEL_MODE="percent"

STABILITY_SET=""
STABILITY_MIN_FOLDS="5"

FROM_STEP="1"
TO_STEP="9"
FORCE="0"
DRY_RUN="0"

# --------------------------
# Parse args
# --------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2; source "$CONFIG" ;;
    --bfile) BFILE="$2"; shift 2 ;;
    --pheno) PHENO="$2"; shift 2 ;;
    --pheno-name) PHENO_NAME="$2"; shift 2 ;;
    --covar) COVAR="$2"; shift 2 ;;
    --exclude) EXCLUDE="$2"; shift 2 ;;
    --qtl-meta) QTL_META="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --plink1) PLINK1="$2"; shift 2 ;;
    --plink2) PLINK2="$2"; shift 2 ;;
    --rscript) RSCRIPT="$2"; shift 2 ;;
    --geno) GENO="$2"; shift 2 ;;
    --mind) MIND="$2"; shift 2 ;;
    --missing-policy) MISSING_POLICY="$2"; shift 2 ;;
    --clump-r2) CLUMP_R2="$2"; shift 2 ;;
    --clump-kb) CLUMP_KB="$2"; shift 2 ;;
    --clump-p1) CLUMP_P1="$2"; shift 2 ;;
    --clump-p2) CLUMP_P2="$2"; shift 2 ;;
    --alpha) ALPHA="$2"; shift 2 ;;
    --outer-folds) OUTER_FOLDS="$2"; shift 2 ;;
    --inner-folds) INNER_FOLDS="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --p-thresholds) P_THRESHOLDS="$2"; shift 2 ;;
    --label-mode) LABEL_MODE="$2"; shift 2 ;;
    --stability-set) STABILITY_SET="$2"; shift 2 ;;
    --stability-min-folds) STABILITY_MIN_FOLDS="$2"; shift 2 ;;
    --from-step) FROM_STEP="$2"; shift 2 ;;
    --to-step) TO_STEP="$2"; shift 2 ;;
    --force) FORCE="1"; shift ;;
    --dry-run) DRY_RUN="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (use --help)" ;;
  esac
done

# --------------------------
# Validate required args
# --------------------------
[[ -n "$BFILE" ]] || die "--bfile is required"
[[ -n "$PHENO" ]] || die "--pheno is required"
[[ -n "$PHENO_NAME" ]] || die "--pheno-name is required"
[[ -n "$OUTDIR" ]] || die "--outdir is required"

# Stability set default depends on label mode
if [[ -z "${STABILITY_SET}" ]]; then
  if [[ "${LABEL_MODE}" == "legacy" ]]; then
    STABILITY_SET="p02"
  else
    STABILITY_SET="p20"
  fi
fi


# --------------------------
# Optional: module management (HPC)
# --------------------------
init_modules() {
  if [[ "${USE_MODULES}" -eq 0 ]]; then
    return 0
  fi

  if type module >/dev/null 2>&1; then
    return 0
  fi

  # Try to initialize environment modules if not present
  if [[ "${MODULE_INIT}" != "" && "${MODULE_INIT}" != "AUTO" && -f "${MODULE_INIT}" ]]; then
    # shellcheck source=/dev/null
    source "${MODULE_INIT}"
  elif [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck source=/dev/null
    source /etc/profile.d/modules.sh
  elif [[ -f /usr/share/Modules/init/bash ]]; then
    # shellcheck source=/dev/null
    source /usr/share/Modules/init/bash
  fi

  type module >/dev/null 2>&1 || die "USE_MODULES=1 but 'module' command not available. Run from a login shell (e.g. sbatch script with '#!/bin/bash -l'), or load modules outside the pipeline, or set PLINK1/PLINK2/RSCRIPT to absolute paths."
}

load_module() {
  local mod="$1"
  if [[ "${USE_MODULES}" -eq 0 ]]; then return 0; fi
  if [[ -z "${mod}" || "${mod}" == "NONE" ]]; then return 0; fi
  echo "+ module load ${mod}"
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    module load "${mod}"
  fi
}

ensure_R()      { load_module "${MODULE_R}"; }
ensure_plink1() { load_module "${MODULE_PLINK1}"; }
ensure_plink2() { load_module "${MODULE_PLINK2}"; }

# --------------------------
# Create directories
# --------------------------
init_modules
mkdir -p "${OUTDIR}"/{00_samples,01_glm_qtl,02_clump_prep,03_clump,04_glm_clumped,05_genotypes,06_en_nested,07_en_stable,08_annotation,logs,manifest}

# Write manifest header (tool versions + run time)
bash "${SCRIPT_DIR}/99_manifest.sh" "${OUTDIR}" "${PLINK1}" "${PLINK2}" "${RSCRIPT}" || true

# Record parameters for reproducibility
PARAMS_OUT="${OUTDIR}/manifest/params_used.txt"
{
  echo "DATE: $(date -Is)"
  echo "BFILE=${BFILE}"
  echo "PHENO=${PHENO}"
  echo "PHENO_NAME=${PHENO_NAME}"
  echo "COVAR=${COVAR}"
  echo "EXCLUDE=${EXCLUDE}"
  echo "QTL_META=${QTL_META}"
  echo "OUTDIR=${OUTDIR}"
  echo "PLINK1=${PLINK1}"
  echo "PLINK2=${PLINK2}"
  echo "RSCRIPT=${RSCRIPT}"
  echo "USE_MODULES=${USE_MODULES}"
  echo "MODULE_INIT=${MODULE_INIT}"
  echo "MODULE_R=${MODULE_R}"
  echo "MODULE_PLINK1=${MODULE_PLINK1}"
  echo "MODULE_PLINK2=${MODULE_PLINK2}"
  echo "GENO=${GENO}"
  echo "MIND=${MIND}"
  echo "MISSING_POLICY=${MISSING_POLICY}"
  echo "AUTO_IMPUTE_MAX_CELL_RATE=${AUTO_IMPUTE_MAX_CELL_RATE}"
  echo "AUTO_DROP_MAX_SAMPLE_RATE=${AUTO_DROP_MAX_SAMPLE_RATE}"
  echo "CLUMP_R2=${CLUMP_R2}"
  echo "CLUMP_KB=${CLUMP_KB}"
  echo "CLUMP_P1=${CLUMP_P1}"
  echo "CLUMP_P2=${CLUMP_P2}"
  echo "ALPHA=${ALPHA}"
  echo "OUTER_FOLDS=${OUTER_FOLDS}"
  echo "INNER_FOLDS=${INNER_FOLDS}"
  echo "SEED=${SEED}"
  echo "P_THRESHOLDS=${P_THRESHOLDS}"
  echo "LABEL_MODE=${LABEL_MODE}"
  echo "STABILITY_SET=${STABILITY_SET}"
  echo "STABILITY_MIN_FOLDS=${STABILITY_MIN_FOLDS}"
  echo "FROM_STEP=${FROM_STEP}"
  echo "TO_STEP=${TO_STEP}"
} > "${PARAMS_OUT}"

KEEP="${OUTDIR}/00_samples/analysis.keep"
KEEP_REPORT="${OUTDIR}/00_samples/analysis.keep.report.tsv"

GLM_QTL_PREFIX="${OUTDIR}/01_glm_qtl/qtl_assoc"
GLM_QTL_FILE="${GLM_QTL_PREFIX}.${PHENO_NAME}.glm.linear"

CLUMP_PREP_PREFIX="${OUTDIR}/02_clump_prep/qtl"
SNPLIST="${CLUMP_PREP_PREFIX}.consistent_variants.snplist"
CLUMPIN="${CLUMP_PREP_PREFIX}.clump_input.tsv"

CLUMP_PREFIX="${OUTDIR}/03_clump/qtl_clumped"
CLUMPED_SNPS="${OUTDIR}/03_clump/qtl_clumped.snplist"
CLUMPED_FILE="${CLUMP_PREFIX}.clumped"

GLM_CLUMP_PREFIX="${OUTDIR}/04_glm_clumped/clumped_assoc"
GLM_CLUMP_FILE="${GLM_CLUMP_PREFIX}.${PHENO_NAME}.glm.linear"

MODEL_BFILE="${OUTDIR}/05_genotypes/ld_clumped_completecase"
RAW_PREFIX="${OUTDIR}/05_genotypes/ld_variants_forEN"
RAW_FILE="${RAW_PREFIX}.raw"

NESTED_OUTDIR="${OUTDIR}/06_en_nested"
STABLE_OUTDIR="${OUTDIR}/07_en_stable"
ANNOT_OUTDIR="${OUTDIR}/08_annotation"

FINAL_RDS="${STABLE_OUTDIR}/final_EN_model_${STABILITY_SET}_stable.rds"
ALLELE_MAP="${NESTED_OUTDIR}/allele_map_from_raw.tsv"

FINAL_ANNOT="${ANNOT_OUTDIR}/final_model_annotated.csv"
FINAL_ANNOT_MAF="${ANNOT_OUTDIR}/final_model_annotated_withMAF.csv"

maybe_skip() {
  local outfile="$1"
  if [[ -f "$outfile" && "$FORCE" -eq 0 ]]; then
    echo "SKIP (exists): $outfile"
    return 0
  fi
  return 1
}

echo "Running steps ${FROM_STEP}..${TO_STEP}"
echo "Outdir: ${OUTDIR}"
echo "Label mode: ${LABEL_MODE}"
echo "Stability set: ${STABILITY_SET}, min folds: ${STABILITY_MIN_FOLDS}"
echo "Missing policy: ${MISSING_POLICY} (geno=${GENO}, mind=${MIND}, auto_cell=${AUTO_IMPUTE_MAX_CELL_RATE}, auto_sample=${AUTO_DROP_MAX_SAMPLE_RATE})"

# STEP 1: keep list
if [[ "${FROM_STEP}" -le 1 && "${TO_STEP}" -ge 1 ]]; then
  ensure_R
  if ! maybe_skip "${KEEP}"; then
    run_cmd "${RSCRIPT} '${SCRIPT_DIR}/00_make_keep_list.R' '${PHENO}' '${PHENO_NAME}' '${EXCLUDE}' '${KEEP}' '${KEEP_REPORT}' > '${OUTDIR}/logs/01_make_keep.log' 2>&1"
  fi
fi

# STEP 2: PLINK2 GLM QTL
if [[ "${FROM_STEP}" -le 2 && "${TO_STEP}" -ge 2 ]]; then
  ensure_plink2
  if ! maybe_skip "${GLM_QTL_FILE}"; then
    if [[ -n "${COVAR}" && "${COVAR}" != "NONE" ]]; then
      run_cmd "${PLINK2} --bfile '${BFILE}' --keep '${KEEP}' --pheno '${PHENO}' --pheno-name '${PHENO_NAME}' --covar '${COVAR}' --covar-variance-standardize --glm hide-covar --out '${GLM_QTL_PREFIX}' > '${OUTDIR}/logs/02_plink2_glm_qtl.log' 2>&1"
    else
      run_cmd "${PLINK2} --bfile '${BFILE}' --keep '${KEEP}' --pheno '${PHENO}' --pheno-name '${PHENO_NAME}' --glm hide-covar --out '${GLM_QTL_PREFIX}' > '${OUTDIR}/logs/02_plink2_glm_qtl.log' 2>&1"
    fi
  fi
fi

# STEP 3: prepare clump inputs
if [[ "${FROM_STEP}" -le 3 && "${TO_STEP}" -ge 3 ]]; then
  ensure_R
  if ! maybe_skip "${CLUMPIN}"; then
    run_cmd "${RSCRIPT} '${SCRIPT_DIR}/02_prepare_clump_inputs.R' '${GLM_QTL_FILE}' '${QTL_META}' '${OUTDIR}/02_clump_prep' 'qtl' > '${OUTDIR}/logs/03_prepare_clump_inputs.log' 2>&1"
  fi
fi

# STEP 4: PLINK1 clump
if [[ "${FROM_STEP}" -le 4 && "${TO_STEP}" -ge 4 ]]; then
  ensure_plink1
  if ! maybe_skip "${CLUMPED_SNPS}"; then
    run_cmd "${PLINK1} --bfile '${BFILE}' --keep '${KEEP}' --extract '${SNPLIST}' --clump '${CLUMPIN}' --clump-p1 '${CLUMP_P1}' --clump-p2 '${CLUMP_P2}' --clump-r2 '${CLUMP_R2}' --clump-kb '${CLUMP_KB}' --out '${CLUMP_PREFIX}' > '${OUTDIR}/logs/04_plink1_clump.log' 2>&1"
    run_cmd "awk 'NR==1{for(i=1;i<=NF;i++) if(\$i==\"SNP\") snp=i; next} NF>0 && snp>0 {print \$snp}' '${CLUMPED_FILE}' > '${CLUMPED_SNPS}'"
  fi
fi

# STEP 5: PLINK2 GLM on clumped set
if [[ "${FROM_STEP}" -le 5 && "${TO_STEP}" -ge 5 ]]; then
  ensure_plink2
  if ! maybe_skip "${GLM_CLUMP_FILE}"; then
    if [[ -n "${COVAR}" && "${COVAR}" != "NONE" ]]; then
      run_cmd "${PLINK2} --bfile '${BFILE}' --keep '${KEEP}' --extract '${CLUMPED_SNPS}' --pheno '${PHENO}' --pheno-name '${PHENO_NAME}' --covar '${COVAR}' --covar-variance-standardize --glm hide-covar --out '${GLM_CLUMP_PREFIX}' > '${OUTDIR}/logs/05_plink2_glm_clumped.log' 2>&1"
    else
      run_cmd "${PLINK2} --bfile '${BFILE}' --keep '${KEEP}' --extract '${CLUMPED_SNPS}' --pheno '${PHENO}' --pheno-name '${PHENO_NAME}' --glm hide-covar --out '${GLM_CLUMP_PREFIX}' > '${OUTDIR}/logs/05_plink2_glm_clumped.log' 2>&1"
    fi
  fi
fi

# STEP 6: complete-case modeling dataset + export raw
if [[ "${FROM_STEP}" -le 6 && "${TO_STEP}" -ge 6 ]]; then
  ensure_plink2
  if ! maybe_skip "${RAW_FILE}"; then
    run_cmd "${PLINK2} --bfile '${BFILE}' --keep '${KEEP}' --extract '${CLUMPED_SNPS}' --geno '${GENO}' --mind '${MIND}' --make-bed --out '${MODEL_BFILE}' > '${OUTDIR}/logs/06_make_completecase_bfile.log' 2>&1"
    run_cmd "${PLINK2} --bfile '${MODEL_BFILE}' --export A --out '${RAW_PREFIX}' > '${OUTDIR}/logs/06_export_raw.log' 2>&1"
  fi
fi

# STEP 7: nested EN
if [[ "${FROM_STEP}" -le 7 && "${TO_STEP}" -ge 7 ]]; then
  ensure_R
  if ! maybe_skip "${NESTED_OUTDIR}/nested_EN_summary.tsv"; then
    run_cmd "${RSCRIPT} '${SCRIPT_DIR}/06_nested_elastic_net_no_impute.R' '${GLM_CLUMP_FILE}' '${RAW_FILE}' '${PHENO}' '${PHENO_NAME}' '${NESTED_OUTDIR}' '${ALPHA}' '${OUTER_FOLDS}' '${INNER_FOLDS}' '${P_THRESHOLDS}' '${LABEL_MODE}' '${SEED}' '${MISSING_POLICY}' '${AUTO_IMPUTE_MAX_CELL_RATE}' '${AUTO_DROP_MAX_SAMPLE_RATE}' > '${OUTDIR}/logs/07_nested_en.log' 2>&1"
  fi
fi

# STEP 8: stability-filtered final model
if [[ "${FROM_STEP}" -le 8 && "${TO_STEP}" -ge 8 ]]; then
  ensure_R
  if ! maybe_skip "${FINAL_RDS}"; then
    run_cmd "${RSCRIPT} '${SCRIPT_DIR}/07_fit_stable_model_no_impute.R' '${NESTED_OUTDIR}/nested_stability_${STABILITY_SET}.tsv' '${RAW_FILE}' '${PHENO}' '${PHENO_NAME}' '${STABLE_OUTDIR}' '${ALPHA}' '${STABILITY_MIN_FOLDS}' '${STABILITY_SET}' '${SEED}' '${MISSING_POLICY}' '${AUTO_IMPUTE_MAX_CELL_RATE}' '${AUTO_DROP_MAX_SAMPLE_RATE}' > '${OUTDIR}/logs/08_fit_stable_model.log' 2>&1"
  fi
fi

# STEP 9: annotate + MAF
if [[ "${FROM_STEP}" -le 9 && "${TO_STEP}" -ge 9 ]]; then
  ensure_R
  if ! maybe_skip "${FINAL_ANNOT_MAF}"; then
    run_cmd "${RSCRIPT} '${SCRIPT_DIR}/08_annotate_final_model.R' '${FINAL_RDS}' '${GLM_CLUMP_FILE}' '${NESTED_OUTDIR}/nested_stability_${STABILITY_SET}.tsv' '${QTL_META}' '${ALLELE_MAP}' '${FINAL_ANNOT}' > '${OUTDIR}/logs/09_annotate.log' 2>&1"
    run_cmd "${RSCRIPT} '${SCRIPT_DIR}/09_add_maf_from_raw.R' '${FINAL_ANNOT}' '${RAW_FILE}' '${FINAL_ANNOT_MAF}' 'Variant' > '${OUTDIR}/logs/09_maf.log' 2>&1"
  fi
fi

echo "DONE."
echo "Final annotated + MAF table: ${FINAL_ANNOT_MAF}"
