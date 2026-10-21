#!/usr/bin/env bash
# Stage 1: AF2-multimer_v3 screen of 267 GPCR-Galpha complexes (3 models, 3 recycles)
# Resumable: colabfold_batch skips any job whose .done.txt already exists.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${AF2_WORK:-$SCRIPT_DIR/af2_reproduction_run}"
IN="$SCRIPT_DIR/af2_original_queries/complexes_267.csv"
OUT="$WORK/stage1"
LOG="$WORK/logs"
mkdir -p "$OUT" "$LOG"

export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false

STAMP=$(date +%Y%m%d_%H%M%S)

colabfold_batch \
  --model-type alphafold2_multimer_v3 \
  --num-models 3 \
  --model-order 1,2,3 \
  --num-recycle 3 \
  --recycle-early-stop-tolerance 0.5 \
  --num-seeds 1 --random-seed 0 \
  --msa-mode mmseqs2_uniref_env \
  --pair-mode unpaired_paired --pair-strategy greedy \
  --rank multimer \
  --sort-queries-by length \
  "$IN" "$OUT" 2>&1 | tee -a "$LOG/stage1_${STAMP}.log"

echo "=== stage1: $(ls "$OUT"/*.done.txt 2>/dev/null | wc -l) / 267 complexes done ==="
