#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
cd "$ROOT/R2_computation"

"$PYTHON_BIN" analysis/reanalyse_computation.py
"$PYTHON_BIN" analysis/analyse_focal_seed.py

Rscript analysis/make_figure2_r2.R \
  inputs/atcc_gpcr_analysis/results/figure_source_data/Fig2a_WT_induction_atlas.tsv \
  inputs/atcc_gpcr_analysis/results/integration/candidate_evidence_by_gse.tsv \
  inputs/tree_annotation_metadata.tsv \
  inputs/selection_criterion.tsv \
  inputs/selection_regulator_effects.tsv \
  inputs/atcc_gpcr_analysis/results/design/series_design_manifest.tsv \
  analysis/figure2_9series/reproduced

echo "R2 tables and revised Figure 2 regenerated under R2_computation/analysis/."
