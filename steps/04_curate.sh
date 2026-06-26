#!/usr/bin/env bash
# =============================================================================
# 04_curate.sh — SLURM: automated curation of annotations
# Thin wrapper around scripts/04_curate_annotations.py (copied AS-IS, with DF fixes)
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
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
fi
source "${PROJECT_DIR}/config.sh"
activate_env
step_banner "Step 04 — Curate annotations"
python "${SCRIPTS_DIR}/04_curate_annotations.py"
log "Done. Next: sbatch steps/04_output.sh"
