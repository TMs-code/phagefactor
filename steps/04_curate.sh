#!/usr/bin/env bash
# =============================================================================
# 04_curate.sh — SLURM: curate annotations + build final output table
# Runs scripts/04_curate_annotations.py then scripts/05_build_output.py in one
# job. These always run back-to-back (curation feeds the output build
# directly, into the SAME 04_output/ folder) so there is no reason to split
# them into two SLURM submissions with a dependency between them.
#
# This is also the target of `SKIP_COMPARE=1 bash submit_all.sh` (see
# submit_all.sh usage note) — the one to resubmit when you've only edited a
# curation rule and want fresh output against an existing run's FoldSeek +
# comparison results, without re-running anything upstream.
# =============================================================================
#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -J pf_curate
#SBATCH -o logs/04_curate.out
#SBATCH -e logs/04_curate.err
#SBATCH --time=01:00:00

set -euo pipefail
# Resolve project root: SLURM_SUBMIT_DIR if sbatch-submitted (BASH_SOURCE points
# into the spool dir once queued), else this script's own parent dir.
PROJECT_DIR="${PHAGEFACTOR_HOME:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
source "${PROJECT_DIR}/config.sh"
activate_env

step_banner "Step 04a — Curate annotations"
python "${SCRIPTS_DIR}/04_curate_annotations.py"

step_banner "Step 04b — Build output (final table + GBK + review_suggested)"
python "${SCRIPTS_DIR}/05_build_output.py"

log "Done. Outputs in 04_output/. Optional next: sbatch steps/05_phynteny.sh"
