#!/bin/bash
# ============================================================
# Run Sherlock-Lung Analysis
# Environment: conda sherlock_lung
# Logs every run separately and also creates a combined log
# ============================================================

set -o pipefail

BASE="/media/wrath/bioinfor_learning/sherlock_lung"
LOGDIR="$BASE/logs"
FIGDIR="$BASE/figures"
STAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p "$LOGDIR" "$FIGDIR"

MAF_LOG="$LOGDIR/maf_run_${STAMP}.log"
EXPR_LOG="$LOGDIR/expression_run_${STAMP}.log"
COMBINED_LOG="$LOGDIR/full_run_${STAMP}.log"

cd "$BASE" || exit 1

# Send everything printed by this shell script to combined log too
exec > >(tee -a "$COMBINED_LOG") 2>&1

echo "============================================"
echo "Sherlock-Lung Analysis"
echo "Start: $(date)"
echo "Base: $BASE"
echo "MAF log: $MAF_LOG"
echo "Expression log: $EXPR_LOG"
echo "Combined log: $COMBINED_LOG"
echo "============================================"

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sherlock_lung

echo ""
echo "============================================"
echo "Running MAF analysis..."
echo "Start MAF: $(date)"
echo "============================================"

Rscript scripts/02_MAF_analysis.R 2>&1 | tee "$MAF_LOG"
MAF_STATUS=${PIPESTATUS[0]}

echo "End MAF: $(date)"
echo "MAF exit status: $MAF_STATUS"

if [ "$MAF_STATUS" -ne 0 ]; then
  echo "ERROR: MAF analysis failed. Check: $MAF_LOG"
  exit "$MAF_STATUS"
fi

echo ""
echo "============================================"
echo "Running expression/ESTIMATE/GSEA/survival analysis..."
echo "Start expression: $(date)"
echo "============================================"

Rscript scripts/01_TCGA_LUAD_analysis.R 2>&1 | tee "$EXPR_LOG"
EXPR_STATUS=${PIPESTATUS[0]}

echo "End expression: $(date)"
echo "Expression exit status: $EXPR_STATUS"

if [ "$EXPR_STATUS" -ne 0 ]; then
  echo "ERROR: Expression analysis failed. Check: $EXPR_LOG"
  exit "$EXPR_STATUS"
fi

echo ""
echo "============================================"
echo "Complete: $(date)"
echo "Figures generated:"
ls -lh "$FIGDIR"/*.png 2>/dev/null || echo "No PNG figures found."
echo ""
echo "Logs generated:"
ls -lh "$LOGDIR"/*"${STAMP}"*.log
echo "============================================"
