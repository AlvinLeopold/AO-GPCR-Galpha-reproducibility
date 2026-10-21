#!/usr/bin/env bash
set -euo pipefail

RUN_LOG="logs/run_all.log"
mkdir -p logs
{
  printf 'ATCC GPCR pipeline\n'
  printf 'started_utc: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  printf 'working_directory: %s\n' "$PWD"
  printf 'command: %s\n' "$0"
  Rscript scripts/01_prepare_inputs.R
  Rscript scripts/02_parse_expression.R
  Rscript scripts/03_analyze_candidates.R
  Rscript scripts/04_integrate_results.R
  Rscript scripts/05_write_report.R
  Rscript scripts/06_design_audit.R
  Rscript scripts/07_series_models.R
  Rscript scripts/08_time_response.R
  Rscript scripts/09_perturbation_consensus.R
  Rscript scripts/10_pathway_network.R
  Rscript scripts/11_condition_panels.R
  Rscript scripts/12_write_deep_report.R
  Rscript scripts/13_nature_figures.R
  printf 'finished_utc: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
} 2>&1 | tee "$RUN_LOG"
