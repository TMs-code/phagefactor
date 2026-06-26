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
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
fi
source "${PROJECT_DIR}/config.sh"
activate_env
step_banner "Step 03 — Compare annotations"
python "${SCRIPTS_DIR}/03_compare_annotations.py"
log "Done. Next: sbatch steps/04_curate.sh"
