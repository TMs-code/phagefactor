#!/usr/bin/env bash
# =============================================================================
# 01_phold_array.sh — SLURM array job: phold on N prophages
#
# Rewritten 2026-05-27 for the new layout:
#   - reads pharokka GBK from 00c_pharokka/<NAME>/pharokka.gbk (not input/gbk/)
#   - uses GPU partition with --foldseek_gpu if PHOLD_USE_GPU=1 in config.sh
#   - passes --autotune so phold tunes the ProstT5 batch size to the hardware
#
# Submit:
#   cd phagefactor/
#   N=$(wc -l < input/prophage_list.txt)
#   sbatch --array=0-$((N-1))%500 steps/01_phold_array.sh
#
# After ALL array tasks complete:
#   python scripts/01c_merge_phold.py        # merge per-prophage TSVs → phold_all.tsv
#   bash steps/01d_merge_3di.sh              # cat 3di FASTA files → phold_3di.fasta
# =============================================================================

# --- The SBATCH header below uses CPU defaults; submit_all.sh overrides
#     --partition/--qos/--gres on the sbatch command line when PHOLD_USE_GPU=1.
#     Always-valid headers go here as a CPU fallback.
#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -J pf_phold_%a
#SBATCH -o logs/01_phold_%a.out
#SBATCH -e logs/01_phold_%a.err
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email@example.com

set -euo pipefail

# Resolve project root.  Under sbatch, SLURM copies this script to
# /var/spool/slurmd/jobXXXXX/ so BASH_SOURCE[0] no longer points to the repo.
# Prefer SLURM_SUBMIT_DIR (preserved by sbatch as the dir user submitted from);
# fall back to BASH_SOURCE when run interactively.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
fi
source "${PROJECT_DIR}/config.sh"

# CRITICAL: use cluster scratch disk, NOT /tmp (too small for phold temp files)
export TMPDIR="${SCRATCH_TMPDIR}"
mkdir -p "${TMPDIR}"

# Activate micromamba env inline (SLURM jobs don't source ~/.bashrc)
activate_env
require_tool phold

# ---------------------------------------------------------------------------
# Identify this task's prophage from the list file
# ---------------------------------------------------------------------------
PROPHAGE_NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${PROPHAGE_LIST}")
if [[ -z "${PROPHAGE_NAME}" ]]; then
    die "No prophage found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
fi

# Pharokka writes ${GBK_DIR}/<NAME>/<NAME>.gbk (the -p prefix is the basename).
GBK_FILE="${GBK_DIR}/${PROPHAGE_NAME}/${PROPHAGE_NAME}.gbk"
check_file "${GBK_FILE}"

OUT_DIR="${PHOLD_OUT_DIR}/${PROPHAGE_NAME}"
mkdir -p "${OUT_DIR}"

# ---------------------------------------------------------------------------
# Skip if already done (allows safe re-submission of failed tasks)
# ---------------------------------------------------------------------------
DONE_MARKER="${OUT_DIR}/.phold_done"
if [[ -f "${DONE_MARKER}" ]]; then
    log "Skipping ${PROPHAGE_NAME} — already done (${DONE_MARKER})"
    exit 0
fi

# ---------------------------------------------------------------------------
# Build phold args
# ---------------------------------------------------------------------------
PHOLD_ARGS=( -i  "${GBK_FILE}"
             -o  "${OUT_DIR}"
             -t  "${SLURM_CPUS_PER_TASK}"
             --prefix "${PROPHAGE_NAME}"
             -f )

# GPU vs CPU mode (controlled by config.sh PHOLD_USE_GPU)
if [[ "${PHOLD_USE_GPU:-0}" == "1" ]]; then
    PHOLD_ARGS+=( --foldseek_gpu )
    GPU_MODE_STR="GPU (--foldseek_gpu)"
    # --autotune requires a GPU (phold raises ERROR when run on CPU).  Only
    # add it when GPU mode is requested.  For CPU runs, phold uses batch_size=1
    # which is the right default.
    if [[ "${PHOLD_AUTOTUNE:-0}" == "1" ]]; then
        PHOLD_ARGS+=( --autotune )
    fi
else
    PHOLD_ARGS+=( --cpu )
    GPU_MODE_STR="CPU (--cpu)"
fi

# ---------------------------------------------------------------------------
# Run phold
# ---------------------------------------------------------------------------
step_banner "PHold: ${PROPHAGE_NAME}  [array task ${SLURM_ARRAY_TASK_ID}]"
log "Input:  ${GBK_FILE}"
log "Output: ${OUT_DIR}/"
log "CPUs:   ${SLURM_CPUS_PER_TASK}"
log "Mode:   ${GPU_MODE_STR}"

phold run "${PHOLD_ARGS[@]}"

# Touch done marker so re-submissions skip this prophage
touch "${DONE_MARKER}"

log "Done: ${PROPHAGE_NAME}  $(date)"

# One-line recap to a SINGLE shared file (instead of N per-task logs to scan).
N_PRED=$(ls "${OUT_DIR}"/*per_cds*.tsv 2>/dev/null | head -1 | xargs -r wc -l 2>/dev/null | awk '{print $1-1}')
printf "%s\t%s\tOK\t%s predictions\n" "$(date +%F_%T)" "${PROPHAGE_NAME}" "${N_PRED:-?}" \
    >> "${PROJECT_DIR}/logs/01_phold_recap.tsv"
