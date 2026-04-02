#!/usr/bin/env bash
set -euo pipefail

SINGULARITY_SILENT_WARNINGS=1
err() { cat <<< "$@" 1>&2; }
die() {
	echo "ERROR: $*" >&2
	exit 1
}

get_current_server() {
    # Get the hostname
    local hn=$(hostname)

    # biowulf compute
    if [[ "$hn" =~ ^cn[0-9]{4}$ ]]; then
        echo "biowulf"
        return 0
    fi

    # skyline compute
    if [[ "$hn" =~ ^ai-hpcn[0-9]+ ]]; then
        echo "skyline"
        return 0
    fi

    # Check if user is attempting to run from
	# head node a biowulf or skyline head node
	if [[ "$hn" == "biowulf.nih.gov" || "$hn" =~ ^ai-hpcsubmit[0-9]+ ]]; then
	    err "############################################################"
	    err "ERROR: Running this pipeline on a head node is not allowed."
	    err "Please submit this as a job to the cluster or run it from"
	    err "an interactive node. You can grab an interactive node by"
		err "running the following command below:"
	    err "srun -N 1 -n 1 --time=8:00:00 --mem=64gb -c 4 --pty bash"
	    err "############################################################"
	    return 1
	fi

    # Unknown
    err "ERROR: Unknown host profile. Are you running this outside of biowulf or skyline?"
    return 1
}

server=$(get_current_server)

if [[ "$server" == "biowulf" ]]; then
	SINGULARITY_IMAGE="/data/OpenOmics/SIFs/elasticqtl_0.0.1.sif"
elif [[ "$server" == "skyline" ]]; then
	SINGULARITY_IMAGE="/data/openomics/SIFs/elasticqtl_0.0.1.sif"
else
	die "Unknown server, cannot determine singularity image path"
fi

run_cmd() {
	local cmd="$1"
	local extra_binds="${2:-}"

	# Collect all paths: PWD + any extras, deduplicated
	local -A seen
	local -a paths
	local deduped=""

	# Split extra_binds on comma into an array
	IFS=',' read -ra paths <<< "${PWD},${extra_binds}"

	for p in "${paths[@]}"; do
		[[ -z "$p" ]] && continue
		[[ -n "${seen[$p]+x}" ]] && continue
		seen["$p"]=1
		if [[ -z "$deduped" ]]; then
			deduped="$p"
		else
			deduped="${deduped},${p}"
		fi
	done

	echo "+ $cmd"
	if [[ "${DRY_RUN:-0}" -eq 0 ]]; then
		singularity exec \
			--pwd "${PWD}" \
			-c \
			-B "${deduped}" \
			"${SINGULARITY_IMAGE}" /bin/bash -c "$cmd"
	fi
}

parse_env() {
	local file="$1"
	local -n arr="$2" # nameref to the caller's associative array

	[[ -f "$file" ]] || {
		echo "ERROR: file not found: $file" >&2
		return 1
	}

	while IFS='=' read -r key value; do
		# skip blank lines and comments
		[[ -z "$key" || "$key" == \#* ]] && continue
		arr["$key"]="$value"
	done < "$file"
}

get_cohort_bfiles() {
    local file="$1"
    local col=-1

    [[ -f "$file" ]] || { echo "ERROR: file not found: $file" >&2; return 1; }

    # Read header to find the 'bfile' column index
    IFS=$'\t' read -ra header < "$file"
    for i in "${!header[@]}"; do
        if [[ "${header[$i]}" == "bfile" ]]; then
            col="$i"
            break
        fi
    done

    [[ "$col" -ge 0 ]] || { echo "ERROR: 'bfile' column not found in header" >&2; return 1; }

    # Print the bfile value from each data row
    local lineno=0
    while IFS=$'\t' read -ra fields; do
        lineno=$((lineno + 1))
        [[ "$lineno" -eq 1 ]] && continue       # skip header
        [[ -z "${fields[$col]:-}" ]] && continue # skip empty
        echo "${fields[$col]}"
    done < "$file"
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
PLINK2="plink2.3"
RSCRIPT="Rscript"

INTERSECTION_MODE="all"
ALPHA="0.5"
NFOLDS="10"
SEED="123"
MISSING_POLICY="mean_impute"
AUTO_IMPUTE_MAX_CELL_RATE="0.001"
AUTO_DROP_MAX_SAMPLE_RATE="0.02"
AMBIGUOUS_POLICY="drop"

HAS_CLI_FLAGS="0"
while [[ $# -gt 0 ]]; do
	case "$1" in
	--config)
		CONFIG="$2"
		shift 2
		source "$CONFIG"
		;;
	--force)
		FORCE="1"
		HAS_CLI_FLAGS="1"
		shift
		;;
	--dry-run)
		DRY_RUN="1"
		HAS_CLI_FLAGS="1"
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*) die "Unknown option: $1" ;;
	esac
done

# [[ -n "${CONFIG}" && "${HAS_CLI_FLAGS}" -eq 1 ]] && die "Either use --config or the other command line interface flags, do not use both."

[[ -n "${CONFIG}" ]] || die "--config is required"
[[ -n "${TRAIN_MODEL_RDS:-}" ]] || die "TRAIN_MODEL_RDS not set in config"
[[ -n "${TRAIN_RAW:-}" ]] || die "TRAIN_RAW not set in config"
[[ -n "${TRAIN_PHENO:-}" ]] || die "TRAIN_PHENO not set in config"
[[ -n "${PHENO_NAME:-}" ]] || die "PHENO_NAME not set in config"
[[ -n "${COHORTS_TSV:-}" ]] || die "COHORTS_TSV not set in config"
[[ -n "${OUTDIR:-}" ]] || die "OUTDIR not set in config"

declare -A config
parse_env "${CONFIG}" config

SCRIPT_DIR="$(dirname $(readlink -f "${BASH_SOURCE[0]}"))"
[[ "${SCRIPT_DIR}" == /vf/users/* ]] && SCRIPT_DIR="/data${SCRIPT_DIR#/vf/users}"

bfiles=( $(get_cohort_bfiles "$COHORTS_TSV") )

# path variables
# config["OUTDIR"]
# config["TRAIN_PHENO"]
# config["COHORTS_TSV"]
# config["TRAIN_RAW"]
# config["TRAIN_MODEL_RDS"]
declare -A bps
bps["TRAIN_MODEL_RDS"]="$(dirname "${TRAIN_MODEL_RDS}")"
bps["TRAIN_RAW"]="$(dirname "${TRAIN_RAW}")"
bps["TRAIN_PHENO"]="$(dirname "${TRAIN_PHENO}")"
bps["COHORTS_TSV"]="$(dirname "${COHORTS_TSV}")"
bps["OUTDIR"]="$(dirname "${OUTDIR}")"

echo ""
echo -e "\033[1;37m════════════════════════════════════════════\033[0m"
echo -e "\033[1;37m  elasticQTL Configuration\033[0m"
echo -e "\033[1;37m════════════════════════════════════════════\033[0m"
for key in "${!config[@]}"; do
	echo -e "  \033[1;36m${key}\033[0m = \033[0;33m${config[$key]}\033[0m"
done
echo -e "\033[1;37m════════════════════════════════════════════\033[0m"
echo ""

bindpaths="${SCRIPT_DIR},${bps["TRAIN_MODEL_RDS"]},${bps["TRAIN_RAW"]},${bps["TRAIN_PHENO"]},${bps["COHORTS_TSV"]},${bps["OUTDIR"]}"

for bf in "${bfiles[@]}"; do
	bindpaths="${bindpaths},$(dirname "$bf")"
done

echo ""
echo -e "\033[1;37m════════════════════════════════════════════\033[0m"
echo -e "\033[1;37m  Singularity Bind Paths\033[0m"
echo -e "\033[1;37m════════════════════════════════════════════\033[0m"
local_i=1
IFS=',' read -ra _bp_arr <<< "$bindpaths"
for p in "${_bp_arr[@]}"; do
	echo -e "  \033[1;36m[$local_i]\033[0m \033[0;33m${p}\033[0m"
	local_i=$((local_i + 1))
done
echo -e "\033[1;37m════════════════════════════════════════════\033[0m"
echo ""

mkdir -p "${OUTDIR}"/{00_model,01_match,02_extract,03_refit,04_scores,logs}

MODEL_WEIGHTS="${OUTDIR}/00_model/model_weights.tsv"
MODEL_VARS="${OUTDIR}/00_model/model_variants.txt"

# run_cmd "ls -al /data/OpenOmics/dev/katie_pipeline/from_katie/elasticQTL/" "$bindpaths"
# die ""

# 1) export model variants
run_cmd "${RSCRIPT} '${SCRIPT_DIR}/validation/10_export_model_variants.R' '${TRAIN_MODEL_RDS}' '${MODEL_WEIGHTS}' '${MODEL_VARS}' > '${OUTDIR}/logs/01_export_model_variants.log'" "$bindpaths"

# 2) match to cohorts + create intersections and per-cohort extract lists
run_cmd "${RSCRIPT} '${SCRIPT_DIR}/validation/11_match_model_to_cohorts.R' '${MODEL_WEIGHTS}' '${COHORTS_TSV}' '${OUTDIR}/01_match' > '${OUTDIR}/logs/02_match_to_cohorts.log'" "$bindpaths"

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
run_cmd "${RSCRIPT} '${SCRIPT_DIR}/validation/13_refit_training_intersection.R' '${TRAIN_RAW}' '${TRAIN_PHENO}' '${PHENO_NAME}' '${SELECTED_MODEL_SET}' '${OUTDIR}/03_refit' '${ALPHA}' '${NFOLDS}' '${SEED}' '${MISSING_POLICY}' '${AUTO_IMPUTE_MAX_CELL_RATE}' '${AUTO_DROP_MAX_SAMPLE_RATE}' > '${OUTDIR}/logs/03_refit.log'" "$bindpaths"
REFIT_RDS="${OUTDIR}/03_refit/refit_model.rds"

# # 4) For each cohort: extract IDs and export raw, then score
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
		run_cmd "${PLINK2} --bfile '${bfile}' --extract '${EXTRACT_IDS}' --export A --out '${COHORT_OUT_PREFIX}' > '${OUTDIR}/logs/04_export_${cohort}.log'" "$bindpaths"
	fi

	MAP_TSV="${OUTDIR}/01_match/${cohort}.mapping.tsv"
	OUT_SCORES="${OUTDIR}/04_scores/${cohort}/${cohort}_scores_refit_${SELECTED_TAG}.csv"
	OUT_QC="${OUTDIR}/04_scores/${cohort}/${cohort}_qc_refit_${SELECTED_TAG}.tsv"

	run_cmd "${RSCRIPT} '${SCRIPT_DIR}/validation/14_score_validation_cohort.R' '${REFIT_RDS}' '${TRAIN_RAW}' '${COHORT_RAW}' '${MAP_TSV}' '${OUT_SCORES}' '${OUT_QC}' '${AMBIGUOUS_POLICY}' 'mean_impute' > '${OUTDIR}/logs/05_score_${cohort}.log'" "$bindpaths"
done <"${COHORTS_TSV}"

echo ""
echo -e "\033[1;37m════════════════════════════════════════════\033[0m"
echo -e "\033[1;37m  DONE validation/refit workflow.\033[0m"
echo -e "\033[1;37m════════════════════════════════════════════\033[0m"
echo -e "  \033[1;36mRefit model\033[0m = \033[0;33m${REFIT_RDS}\033[0m"
echo -e "  \033[1;36mScores under\033[0m = \033[0;33m${OUTDIR}/04_scores/\033[0m"
echo -e "\033[1;37m════════════════════════════════════════════\033[0m"
echo ""