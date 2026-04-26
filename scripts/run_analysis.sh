#!/bin/bash
# ============================================================
# Run Sherlock-Lung Analysis in Docker
# ============================================================

BASE="/media/wrath/bioinfor_learning/sherlock_lung"

echo "Starting Sherlock-Lung analysis container..."
echo "Data: $BASE/data"
echo "Figures will be saved to: $BASE/figures"

docker run --rm \
  --name sherlock_analysis \
  --memory="24g" \
  --cpus="14" \
  -v "$BASE/data:/analysis/data:ro" \
  -v "$BASE/scripts:/analysis/scripts:ro" \
  -v "$BASE/figures:/analysis/figures:rw" \
  sherlock_lung:v1.0 \
  Rscript /analysis/scripts/02_MAF_analysis.R

echo "Analysis complete. Check figures:"
ls -lh "$BASE/figures/"
