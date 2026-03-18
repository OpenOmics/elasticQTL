\
#!/usr/bin/env bash
set -euo pipefail

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
  bash scripts/run_validation.sh --config configs/my_validation.env [--dry-run] [--force]

This validation workflow:
- exports non-zero variants from TRAIN_MODEL_RDS
- checks availability in each validation cohort (by BIM CHR/POS + allele set)
- chooses an intersection set (all/any/cohort:<name>)
- extracts those variants from each cohort bfile and exports .raw
- refits EN in the training cohort using the chosen intersection set
- scores each cohort using the refit model with allele harmonization

Required config variables (in --config):
  TRAIN_MODEL_RDS, TRAIN_RAW, TRAIN_PHENO, PHENO_NAME, COHORTS_TSV, OUTDIR

Optional / important:
  INTERSECTION_MODE (all|any|cohort:<name>)
  MISSING_POLICY (error|drop_samples|mean_impute|auto)
  AMBIGUOUS_POLICY (drop|keep)

HPC modules (optional):
  USE_MODULES=1, MODULE_R, MODULE_PLINK2, MODULE_INIT=AUTO
EOF
}

CONFIG=""
FORCE="0"
DRY_RUN="0"

# defaults (can be overridden by config)
PLINK2="plink2"
RSCRIPT="Rscript"
USE_MODULES="0"
MODULE_INIT="AUTO"
MODULE_R="NONE"
MODULE_PLINK2="NONE"

INTERSECTION_MODE="all"
ALPHA="0.5"
NFOLDS="10"
SEED="123"
MISSING_POLICY="mean_impute"
AUTO_IMPUTE_MAX_CELL_RATE="0.001"
AUTO_DROP_MAX_SAMPLE_RATE="0.02"
AMBIGUOUS_POLICY="drop"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2; source "$CONFIG" ;;
    --force) FORCE="1"; shift ;;
    --dry-run) DRY_RUN="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

[[ -n "${CONFIG}" ]] || die "--config is required"
[[ -n "${TRAIN_MODEL_RDS:-}" ]] || die "TRAIN_MODEL_RDS not set in config"
[[ -n "${TRAIN_RAW:-}" ]] || die "TRAIN_RAW not set in config"
[[ -n "${TRAIN_PHENO:-}" ]] || die "TRAIN_PHENO not set in config"
[[ -n "${PHENO_NAME:-}" ]] || die "PHENO_NAME not set in config"
[[ -n "${COHORTS_TSV:-}" ]] || die "COHORTS_TSV not set in config"
[[ -n "${OUTDIR:-}" ]] || die "OUTDIR not set in config"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

init_modules() {
  if [[ "${USE_MODULES}" -eq 0 ]]; then return 0; fi
  if type module >/dev/null 2>&1; then return 0; fi

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

  type module >/dev/null 2>&1 || die "USE_MODULES=1 but 'module' not available. Use a login shell '#!/bin/bash -l' or load modules outside."
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

ensure_R() { load_module "${MODULE_R}"; }
ensure_plink2() { load_module "${MODULE_PLINK2}"; }

init_modules

mkdir -p "${OUTDIR}"/{00_model,01_match,02_extract,03_refit,04_scores,logs}

MODEL_WEIGHTS="${OUTDIR}/00_model/model_weights.tsv"
MODEL_VARS="${OUTDIR}/00_model/model_variants.txt"

# 1) export model variants
ensure_R
run_cmd "${RSCRIPT} '${SCRIPT_DIR}/validation/10_export_model_variants.R' '${TRAIN_MODEL_RDS}' '${MODEL_WEIGHTS}' '${MODEL_VARS}' > '${OUTDIR}/logs/01_export_model_variants.log' 2>&1"

# 2) match to cohorts + create intersections and per-cohort extract lists
ensure_R
run_cmd "${RSCRIPT} '${SCRIPT_DIR}/validation/11_match_model_to_cohorts.R' '${MODEL_WEIGHTS}' '${COHORTS_TSV}' '${OUTDIR}/01_match' > '${OUTDIR}/logs/02_match_to_cohorts.log' 2>&1"

# Choose intersection variant list (model Variant IDs)
SELECTED_MODEL_SET=""
SELECTED_TAG=""
case "${INTERSECTION_MODE}" in
  all)
    SELECTED_MODEL_SET="${OUTDIR}/01_match/intersection_all.model_variants.txt"
    SELECTED_TAG="all"
    ;;
  any)
    SELECTED_MODEL_SET="${OUTDIR}/01_match/intersection_any.model_variants.txt"
    SELECTED_TAG="any"
    ;;
  cohort:*)
    cname="${INTERSECTION_MODE#cohort:}"
    SELECTED_MODEL_SET="${OUTDIR}/01_match/${cname}.available_model_variants.txt"
    SELECTED_TAG="cohort_${cname}"
    ;;
  *)
    die "Unknown INTERSECTION_MODE: ${INTERSECTION_MODE}"
    ;;
esac

[[ -f "${SELECTED_MODEL_SET}" ]] || die "Selected model set file missing: ${SELECTED_MODEL_SET}"

# 3) Refit EN in training cohort on selected set
ensure_R
run_cmd "${RSCRIPT} '${SCRIPT_DIR}/validation/13_refit_training_intersection.R' '${TRAIN_RAW}' '${TRAIN_PHENO}' '${PHENO_NAME}' '${SELECTED_MODEL_SET}' '${OUTDIR}/03_refit' '${ALPHA}' '${NFOLDS}' '${SEED}' '${MISSING_POLICY}' '${AUTO_IMPUTE_MAX_CELL_RATE}' '${AUTO_DROP_MAX_SAMPLE_RATE}' > '${OUTDIR}/logs/03_refit.log' 2>&1"

REFIT_RDS="${OUTDIR}/03_refit/refit_model.rds"

# 4) For each cohort: extract IDs and export raw, then score
ensure_plink2
while IFS=$'\t' read -r cohort bfile; do
  [[ "${cohort}" == "cohort" ]] && continue
  [[ -z "${cohort}" ]] && continue

  mkdir -p "${OUTDIR}/02_extract/${cohort}" "${OUTDIR}/04_scores/${cohort}"

  # Choose extract id list based on intersection selection
  EXTRACT_IDS=""
  if [[ "${INTERSECTION_MODE}" == "all" ]]; then
    EXTRACT_IDS="${OUTDIR}/01_match/extract_ids_all.${cohort}.txt"
  elif [[ "${INTERSECTION_MODE}" == "any" ]]; then
    EXTRACT_IDS="${OUTDIR}/01_match/extract_ids_any.${cohort}.txt"
  else
    EXTRACT_IDS="${OUTDIR}/01_match/extract_ids_available.${cohort}.txt"
  fi

  [[ -f "${EXTRACT_IDS}" ]] || die "Missing extract IDs for cohort ${cohort}: ${EXTRACT_IDS}"

  COHORT_OUT_PREFIX="${OUTDIR}/02_extract/${cohort}/${cohort}_extracted"
  COHORT_RAW="${COHORT_OUT_PREFIX}.raw"

  # Extract + export raw (skip if exists unless --force)
  if [[ -f "${COHORT_RAW}" && "${FORCE}" -eq 0 ]]; then
    echo "SKIP export (exists): ${COHORT_RAW}"
  else
    ensure_plink2
    run_cmd "${PLINK2} --bfile '${bfile}' --extract '${EXTRACT_IDS}' --export A --out '${COHORT_OUT_PREFIX}' > '${OUTDIR}/logs/04_export_${cohort}.log' 2>&1"
  fi

  MAP_TSV="${OUTDIR}/01_match/${cohort}.mapping.tsv"
  OUT_SCORES="${OUTDIR}/04_scores/${cohort}/${cohort}_scores_refit_${SELECTED_TAG}.csv"
  OUT_QC="${OUTDIR}/04_scores/${cohort}/${cohort}_qc_refit_${SELECTED_TAG}.tsv"

  ensure_R
  run_cmd "${RSCRIPT} '${SCRIPT_DIR}/validation/14_score_validation_cohort.R' '${REFIT_RDS}' '${TRAIN_RAW}' '${COHORT_RAW}' '${MAP_TSV}' '${OUT_SCORES}' '${OUT_QC}' '${AMBIGUOUS_POLICY}' 'mean_impute' > '${OUTDIR}/logs/05_score_${cohort}.log' 2>&1"

done < "${COHORTS_TSV}"

echo "DONE validation/refit workflow."
echo "Refit model: ${REFIT_RDS}"
echo "Scores under: ${OUTDIR}/04_scores/"
