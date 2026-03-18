#!/usr/bin/env bash
set -euo pipefail

OUTDIR="$1"
PLINK1="${2:-plink}"
PLINK2="${3:-plink2}"
RSCRIPT="${4:-Rscript}"

mkdir -p "${OUTDIR}/manifest"

{
  echo "DATE: $(date -Is)"
  echo "HOST: $(hostname)"
  echo ""
  echo "PLINK1: $(${PLINK1} --version 2>/dev/null || echo 'unknown')"
  echo "PLINK2: $(${PLINK2} --version 2>/dev/null || echo 'unknown')"
  echo "R: $(${RSCRIPT} --version 2>/dev/null | head -n 1 || echo 'unknown')"
  echo ""
  if type module >/dev/null 2>&1; then
    echo "MODULE LIST (at manifest time):"
    module list 2>&1 || true
  fi
} > "${OUTDIR}/manifest/tool_versions.txt"
