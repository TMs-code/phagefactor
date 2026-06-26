#!/usr/bin/env bash
# =============================================================================
# 02_foldseek_3di.sh — SLURM: local FoldSeek 3di search (all 3 DBs)
#
# Wraps scripts/02d_foldseek_3di.py.
# Only config/config.yaml changes the DB paths, thread count, and output dirs.
#
# Submit:
#   sbatch steps/02_foldseek_3di.sh
#   (or as dependency: --dependency=afterok:<phold_merge_jobid>)
#
# Runtime estimate: ~4–8h for ~500k hypothetical sequences × 3 databases
# =============================================================================

#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G  # afdb50_ss.idx is 196 GB; loaded into RAM for the prefilter.
                     # (Do NOT lower + mmap: at <index-size RAM the index thrashes
                     #  from disk and the search goes from ~6min to >1h.)
                    # If common nodes don't have 256G, try bigmem partition or Option B (--prefilter-mode 1)
#SBATCH -J pf_foldseek
#SBATCH -o logs/02_foldseek.out
#SBATCH -e logs/02_foldseek.err
#SBATCH --time=08:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email@example.com

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
fi
source "${PROJECT_DIR}/config.sh"

# CRITICAL: cluster scratch disk for foldseek temp files
export TMPDIR="${SCRATCH_TMPDIR}"
mkdir -p "${TMPDIR}"

activate_env
require_tool foldseek

step_banner "Step 02 — FoldSeek 3di search"
log "SLURM job: ${SLURM_JOB_ID}  CPUs: ${SLURM_CPUS_PER_TASK}"

# Run the Python script (reads config/config.yaml for all paths / DB names).
#
# Optional FOLDSEEK_DBS env var lets you restrict to a subset of DBs at
# submit time without editing the Python.  Example:
#   FOLDSEEK_DBS="afdb-swissprot pdb100" sbatch steps/02_foldseek_3di.sh
# (Validates the pipeline on the two taxonomy-rich DBs before adding
#  baktfold-afdb.)
DBS_ARG=()
if [[ -n "${FOLDSEEK_DBS:-}" ]]; then
    log "Restricting to user-specified DBs: ${FOLDSEEK_DBS}"
    DBS_ARG=(--dbs ${FOLDSEEK_DBS})
fi

python "${SCRIPTS_DIR}/02d_foldseek_3di.py" \
    "${DBS_ARG[@]}" \
    2>&1

log "FoldSeek search complete: $(date)"
log "Next: sbatch steps/03_compare.sh"
