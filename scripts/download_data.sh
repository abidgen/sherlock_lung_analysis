#!/usr/bin/env bash
# ============================================================
# download_data.sh
# Download public TCGA LUAD data for Sherlock-Lung analysis
#
# Sources: GDC API + UCSC Xena TCGA Hub
# No authentication required for the required open-access files
#
# Usage:
#   bash scripts/download_data.sh
#   bash scripts/download_data.sh --skip-existing
#
# Notes:
#   - Required core files: MAF, expression, clinical
#   - Optional integrated-analysis file: PanCancer GISTIC2 CNA
#   - Telomere data are intentionally not downloaded because the
#     LUAD-specific Xena telomere URL is not available in this workflow.
# ============================================================

set -euo pipefail

# Resolve project root from this script location, unless explicitly supplied.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${SHERLOCK_LUNG_BASE:-$(cd "$SCRIPT_DIR/.." && pwd)}"
DATA="$BASE/data/TCGA_LUAD"
SKIP_EXISTING=false

for arg in "$@"; do
  case "$arg" in
    --skip-existing) SKIP_EXISTING=true ;;
    *) echo "Unknown argument: $arg"; exit 1 ;;
  esac
done

mkdir -p "$DATA"
cd "$DATA"

echo "============================================"
echo "Sherlock-Lung Data Download"
echo "Target directory: $DATA"
echo "Start: $(date)"
echo "============================================"

# ---- Helpers ----
file_ok() {
  local file="$1"
  [ -s "$file" ]
}

download_file() {
  local outfile="$1"
  local description="$2"
  local required="$3"
  shift 3
  local urls=("$@")

  if [ "$SKIP_EXISTING" = true ] && file_ok "$outfile"; then
    echo "[SKIP] $description — already exists ($(du -sh "$outfile" | cut -f1))"
    return 0
  fi

  echo ""
  echo "Downloading: $description"
  echo "  Output: $outfile"

  rm -f "${outfile}.tmp"
  local ok=false

  for url in "${urls[@]}"; do
    echo "  Trying: $url"
    if wget -c --show-progress --timeout=60 --tries=3 "$url" -O "${outfile}.tmp"; then
      if [ -s "${outfile}.tmp" ]; then
        mv "${outfile}.tmp" "$outfile"
        echo "  Done: $(du -sh "$outfile" | cut -f1)"
        ok=true
        break
      fi
    fi
    rm -f "${outfile}.tmp"
  done

  if [ "$ok" = false ]; then
    if [ "$required" = "required" ]; then
      echo "ERROR: failed to download required file: $description"
      return 1
    else
      echo "WARNING: failed to download optional file: $description"
      echo "         This is not fatal; integrated CNA figures will be skipped unless the file is added manually."
      return 0
    fi
  fi
}

verify_file() {
  local file="$1"
  local min_size_kb="$2"
  local description="$3"
  local required="$4"

  local actual_kb=0
  if [ -s "$file" ]; then
    actual_kb=$(du -k "$file" | cut -f1)
  fi

  if [ "$actual_kb" -lt "$min_size_kb" ]; then
    if [ "$required" = "required" ]; then
      echo "ERROR: $description missing/incomplete (${actual_kb}KB < ${min_size_kb}KB expected)"
      return 1
    else
      echo "OPTIONAL WARNING: $description missing/incomplete (${actual_kb}KB < ${min_size_kb}KB expected)"
      return 0
    fi
  fi

  echo "OK: $description ($(du -sh "$file" | cut -f1))"
  return 0
}

# ============================================================
# FILE 1: Somatic MAF — GDC API batch download
# Uses download_luad_maf.py to fetch all 618 per-sample MAFs
# Project: TCGA-LUAD | Data type: Masked Somatic Mutation
# Expected: 558 unique patients, GRCh38, ~60MB
# ============================================================
echo ""
echo "============================================"
echo "FILE 1: TCGA-LUAD Somatic MAF (GDC API)"
echo "============================================"

MAF_OUT="$DATA/TCGA_LUAD_somatic.maf.gz"

if [ "$SKIP_EXISTING" = true ] && [ -f "$MAF_OUT" ] && \
   [ "$(du -k "$MAF_OUT" | cut -f1)" -ge 40000 ]; then
  echo "[SKIP] Somatic MAF ($(du -sh "$MAF_OUT" | cut -f1))"
else
  if ! command -v python3 &>/dev/null; then
    echo "ERROR: python3 required. Install: conda install python"
    exit 1
  fi

  echo "Fetching 618 per-sample LUAD MAFs from GDC API..."
  echo "Expected time: 5-15 minutes"
  echo ""
  python3 "$BASE/scripts/download_luad_maf.py" --output "$MAF_OUT" --batch-size 100
  echo "MAF complete: $(du -sh "$MAF_OUT" | cut -f1)"
fi

# ============================================================
# Required file 2: Gene expression — UCSC Xena TCGA Hub
# Dataset: TCGA.LUAD.sampleMap/HiSeqV2.gz
# Values: normalized/log-transformed RNA-seq expression
# ============================================================
download_file \
  "TCGA_LUAD_HiSeqV2.gz" \
  "TCGA LUAD RNA-seq expression (UCSC Xena HiSeqV2)" \
  "required" \
  "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz" \
  "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/HiSeqV2.gz"

# ============================================================
# Required file 3: Clinical data — UCSC Xena TCGA Hub
# Dataset: TCGA.LUAD.sampleMap/LUAD_clinicalMatrix
# Key variables: tobacco_smoking_history, vital_status,
# days_to_death, days_to_last_followup
# ============================================================
download_file \
  "TCGA_LUAD_clinical.tsv" \
  "TCGA LUAD clinical matrix (UCSC Xena)" \
  "required" \
  "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FLUAD_clinicalMatrix" \
  "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix"

# ============================================================
# Optional file 4: Copy number alterations — UCSC Xena
# Dataset: PanCancer GISTIC2 thresholded by gene
# Values: -2 deep deletion, -1 deletion, 0 neutral, 1 gain, 2 amplification
# Local file expected by 03_integrated_analysis.R:
#   TCGA_PANCAN_CNA.tsv.gz
# ============================================================
download_file \
  "TCGA_PANCAN_CNA.tsv.gz" \
  "PanCancer TCGA GISTIC2 thresholded CNA by gene (optional)" \
  "optional" \
  "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz" \
  "https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"

# ============================================================
# Verification
# ============================================================
echo ""
echo "============================================"
echo "Verification"
echo "============================================"

ALL_OK=true
verify_file "TCGA_LUAD_somatic.maf.gz" 40000 "Somatic MAF" "required" || ALL_OK=false
verify_file "TCGA_LUAD_HiSeqV2.gz"     20000 "Expression HiSeqV2" "required" || ALL_OK=false
verify_file "TCGA_LUAD_clinical.tsv"   100   "Clinical matrix" "required" || ALL_OK=false
verify_file "TCGA_PANCAN_CNA.tsv.gz"   5000  "PanCancer GISTIC2 CNA" "optional" || true

echo ""
echo "Files present:"
find "$DATA" -maxdepth 1 -type f \( -name '*.gz' -o -name '*.tsv' -o -name '*.csv' -o -name '*.rds' \) -printf '%8s  %f\n' 2>/dev/null | sort || true

echo ""
echo "Total data directory size:"
du -sh "$DATA"

echo ""
echo "============================================"
echo "Complete: $(date)"
if [ "$ALL_OK" = true ]; then
  echo "Status: REQUIRED FILES OK"
  echo "Next step: bash scripts/run_analysis.sh"
else
  echo "Status: ERROR — one or more required files are missing or incomplete"
fi
echo "============================================"
