#!/bin/bash
# ============================================================
# Sherlock-Lung Full Analysis Pipeline
# Environment: conda sherlock_lung
# Scripts: MAF → Expression → Integrated
# ============================================================

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${SHERLOCK_LUNG_BASE:-$(cd "$SCRIPT_DIR/.." && pwd)}"
LOGDIR="$BASE/logs"
FIGDIR="$BASE/figures"
STAMP="$(date +%Y%m%d_%H%M%S)"

mkdir -p "$LOGDIR" "$FIGDIR"
cd "$BASE" || exit 1

MAF_LOG="$LOGDIR/maf_run_${STAMP}.log"
EXPR_LOG="$LOGDIR/expression_run_${STAMP}.log"
INT_LOG="$LOGDIR/integrated_run_${STAMP}.log"
COMBINED_LOG="$LOGDIR/full_run_${STAMP}.log"

exec > >(tee -a "$COMBINED_LOG") 2>&1

run_step() {
  local label="$1"
  local script="$2"
  local log_file="$3"
  local required="$4"

  echo ""
  echo "============================================"
  echo "$label"
  echo "Start: $(date)"
  echo "Script: $script"
  echo "Log: $log_file"
  echo "============================================"

  Rscript "$script" 2>&1 | tee "$log_file"
  local status=${PIPESTATUS[0]}

  echo "$label exit status: $status — $(date)"
  if [ "$status" -ne 0 ]; then
    if [ "$required" = "required" ]; then
      echo "ERROR: $label failed. Check: $log_file"
      exit "$status"
    else
      echo "WARNING: $label failed. Check: $log_file"
    fi
  fi
}

# Activate conda only if conda is available in this shell.
if command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate sherlock_lung
fi

export SHERLOCK_LUNG_BASE="$BASE"

echo "============================================"
echo "Sherlock-Lung Full Analysis Pipeline"
echo "Start: $(date)"
echo "Base: $BASE"
echo "Rscript: $(which Rscript)"
echo "============================================"

# Quick input check.
for f in \
  "$BASE/data/TCGA_LUAD/TCGA_LUAD_clinical.tsv" \
  "$BASE/data/TCGA_LUAD/TCGA_LUAD_HiSeqV2.gz" \
  "$BASE/data/TCGA_LUAD/TCGA_LUAD_somatic.maf.gz"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: required input missing or empty: $f"
    exit 1
  fi
done

run_step "STEP 1: MAF Analysis" "scripts/02_MAF_analysis.R" "$MAF_LOG" "required"
run_step "STEP 2: Expression + Immune + GSEA + Survival" "scripts/01_TCGA_LUAD_analysis.R" "$EXPR_LOG" "required"
run_step "STEP 3: Integrated Driver-Expression + CNA" "scripts/03_integrated_analysis.R" "$INT_LOG" "optional"

echo ""
echo "============================================"
echo "Pipeline complete: $(date)"
echo "Figures generated ($(find "$FIGDIR" -maxdepth 1 -name '*.png' | wc -l) total):"
find "$FIGDIR" -maxdepth 1 -name '*.png' -printf '%s %p\n' 2>/dev/null | sort -n | awk '{printf "%8.1fK %s\n", $1/1024, $2}' | sed "s|$FIGDIR/||"
echo ""
echo "Logs:"
ls -lh "$LOGDIR"/*"${STAMP}"*.log
echo "============================================"
