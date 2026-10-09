#!/usr/bin/env bash
# =============================================================================
# 03_compare.sh — SLURM: compare phold vs foldseek annotations
# Thin wrapper around scripts/03_compare_annotations.py (copied AS-IS)
# =============================================================================
#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -J pf_compare
#SBATCH -o logs/03_compare.out
#SBATCH -e logs/03_compare.err
#SBATCH --time=01:00:00

set -euo pipefail
# Resolve project root: SLURM_SUBMIT_DIR if sbatch-submitted (BASH_SOURCE points
# into the spool dir once queued), else this script's own parent dir.
PROJECT_DIR="${PHAGEFACTOR_HOME:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
source "${PROJECT_DIR}/config.sh"
activate_env
step_banner "Step 03 — Compare annotations"
python "${SCRIPTS_DIR}/03_compare_annotations.py"
log "Done. Next: sbatch steps/04_curate.sh"
