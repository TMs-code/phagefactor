#!/usr/bin/env bash
# =============================================================================
# 01c_merge.sh — SLURM: merge per-prophage phold outputs (TSV + 3Di/AA FASTA)
# Runs AFTER the phold array (01_phold_array.sh) completes, BEFORE FoldSeek (02).
# Wraps scripts/01c_merge_phold.py + steps/01d_merge_3di.sh so submit_all can
# chain 01 -> 01c_merge -> 02 automatically (no manual SKIP_PHOLD re-run).
# =============================================================================
#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -J phold_merge
#SBATCH -o logs/01c_merge.out
#SBATCH -e logs/01c_merge.err
#SBATCH --time=00:20:00

set -euo pipefail
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
fi
source "${PROJECT_DIR}/config.sh"
activate_env
step_banner "Step 01c — Merge phold outputs"

python "${SCRIPTS_DIR}/01c_merge_phold.py"     # per-prophage TSVs -> phold_all.tsv
bash   "${PROJECT_DIR}/steps/01d_merge_3di.sh"  # cat 3Di + AA FASTAs -> combined/

log "Done. FoldSeek (step 02) can now run."
