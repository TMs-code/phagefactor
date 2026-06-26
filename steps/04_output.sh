#!/usr/bin/env bash
# =============================================================================
# 04_output.sh — SLURM: build final output table + updated GenBank
# Thin wrapper around scripts/05_build_output.py.
# 2026-06: curation (04_curate.sh) + this build step now write to ONE
# general-output folder (04_output/). Phynteny + integration are the next step
# (05_phynteny.sh). Renamed from the old 05_output.sh.
# =============================================================================
#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -J pf_output
#SBATCH -o logs/04_output.out
#SBATCH -e logs/04_output.err
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
step_banner "Step 04 — Build output (final table + GBK + review_suggested)"
python "${SCRIPTS_DIR}/05_build_output.py"
log "Done. Outputs in 04_output/. Optional next: sbatch steps/05_phynteny.sh"
