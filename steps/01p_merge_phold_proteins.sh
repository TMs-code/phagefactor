#!/usr/bin/env bash
# =============================================================================
# 01p_merge_phold_proteins.sh — merge protein-mode phold batch outputs
#
# Thin sbatch wrapper around scripts/01p_merge_phold_proteins.py. Submitted
# with --dependency=afterok:<array_job_id> so it only runs once every batch
# in the 01p_phold_proteins_array.sh array has finished successfully.
#
# Lightweight (pandas I/O only — no model inference), so CPU + modest
# resources are sufficient regardless of PHOLD_PROTEINS_USE_GPU.
# =============================================================================

#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -J phold_proteins_merge
#SBATCH -o logs/01p_merge_phold_proteins.out
#SBATCH -e logs/01p_merge_phold_proteins.err
#SBATCH --time=00:30:00
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# Resolve project root: SLURM_SUBMIT_DIR if sbatch-submitted (BASH_SOURCE points
# into the spool dir once queued), else this script's own parent dir.
PROJECT_DIR="${PHAGEFACTOR_HOME:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
source "${PROJECT_DIR}/config.sh"

activate_env

step_banner "Merging protein-mode phold batches"
cd "${PROJECT_DIR}"
python3 scripts/01p_merge_phold_proteins.py

log "Done: $(date)"
log "Next: sbatch steps/02_foldseek_3di.sh  (auto-chained by submit_all.sh)"
