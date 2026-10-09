#!/usr/bin/env bash
# =============================================================================
# 01p_phold_proteins_array.sh — SLURM array job: phold proteins-{predict,compare}
#                               on fixed-size protein batches (PROTEIN MODE)
#
# Runs as a batch array rather than one single-job GPU submission. Input
# proteins are pre-split by scripts/00p_split_protein_batches.py into fixed
# PROTEIN_BATCH_SIZE-protein batches (default 50, reproducible across
# datasets — last batch may be partial). Each array task processes ONE batch
# FASTA independently — mirrors 01_phold_array.sh's per-prophage pattern.
#
# CPU-by-default (controlled by PHOLD_PROTEINS_USE_GPU in config.sh, NOT the
# genome-mode PHOLD_USE_GPU): decoupled from genome mode to avoid the same
# GPU-autotune failure mode described in 01_phold_array.sh.
#
# Submit (normally done by submit_all.sh, which computes N from the batch list):
#   N=$(wc -l < input/protein_batch_list.txt)
#   sbatch --array=0-$((N-1))%500 steps/01p_phold_proteins_array.sh
#
# After ALL array tasks complete:
#   python scripts/01p_merge_phold_proteins.py   # merge batches -> 01_phold/{combined,proteins}/...
# =============================================================================

# --- The SBATCH header below uses CPU defaults; submit_all.sh overrides
#     --partition/--qos/--gres on the sbatch command line only when
#     PHOLD_PROTEINS_USE_GPU=1. Always-valid headers go here as a CPU fallback
#     (mirrors 01_phold_array.sh).
#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -J phold_proteins_%a
#SBATCH -o logs/01p_phold_proteins_%a.out
#SBATCH -e logs/01p_phold_proteins_%a.err
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# Resolve project root: SLURM_SUBMIT_DIR if sbatch-submitted (BASH_SOURCE points
# into the spool dir once queued), else this script's own parent dir.
PROJECT_DIR="${PHAGEFACTOR_HOME:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
source "${PROJECT_DIR}/config.sh"

export TMPDIR="${SCRATCH_TMPDIR:-/local/scratch/tmp}/phold_proteins_${SLURM_ARRAY_TASK_ID:-0}_$$"
mkdir -p "${TMPDIR}"

activate_env
require_tool phold

# ---------------------------------------------------------------------------
# Identify this task's batch from the list file written by
# scripts/00p_split_protein_batches.py
# ---------------------------------------------------------------------------
BATCH_LIST="${PROTEIN_BATCH_LIST}"
check_file "${BATCH_LIST}"

BATCH_NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${BATCH_LIST}")
if [[ -z "${BATCH_NAME}" ]]; then
    die "No batch found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
fi

BATCH_FAA="${PROTEIN_BATCH_DIR}/${BATCH_NAME}.faa"
check_file "${BATCH_FAA}"

BATCH_OUT_DIR="${PHOLD_OUT_DIR}/proteins/batches/${BATCH_NAME}"
PREDICT_DIR="${BATCH_OUT_DIR}/predictions"
COMPARE_DIR="${BATCH_OUT_DIR}/compare"
mkdir -p "${PREDICT_DIR}" "${COMPARE_DIR}"

N_SEQS=$(grep -c "^>" "${BATCH_FAA}" || true)

# ---------------------------------------------------------------------------
# Skip if already done (allows safe re-submission of failed tasks)
# ---------------------------------------------------------------------------
DONE_MARKER="${BATCH_OUT_DIR}/.batch_done"
if [[ -f "${DONE_MARKER}" ]]; then
    log "Skipping ${BATCH_NAME} — already done (${DONE_MARKER})"
    exit 0
fi

step_banner "PHold proteins: ${BATCH_NAME}  [array task ${SLURM_ARRAY_TASK_ID}, ${N_SEQS} proteins]"
log "Input:  ${BATCH_FAA}"
log "Output: ${BATCH_OUT_DIR}/"
log "CPUs:   ${SLURM_CPUS_PER_TASK:-${THREADS}}"

# ---------------------------------------------------------------------------
# Step 1: phold proteins-predict (ProstT5 -> 3Di tokens)
#
# GPU/CPU selection is governed by PHOLD_PROTEINS_USE_GPU (config.sh), which
# is INTENTIONALLY decoupled from the genome-mode PHOLD_USE_GPU — see config.sh
# comment for the rationale (defaults to CPU=0).
# --autotune is GPU-only and fails immediately on CPU; we additionally probe
# PyTorch at runtime since SLURM may allocate a GPU node that PyTorch still
# can't see (the same fallback pattern as the old single-job 01p script).
# ---------------------------------------------------------------------------
step_banner "phold proteins-predict (${BATCH_NAME})"

PREDICT_DONE="${PREDICT_DIR}/.predict_done"
if [[ -f "${PREDICT_DONE}" ]]; then
    log "Already done — skipping (delete ${PREDICT_DONE} to re-run)"
else
    HAS_GPU=0
    if [[ "${PHOLD_PROTEINS_USE_GPU:-0}" == "1" ]]; then
        HAS_GPU=$(python3 -c "import torch; print(int(torch.cuda.is_available()))" 2>/dev/null || echo "0")
    fi

    if [[ "${HAS_GPU}" == "1" && "${PHOLD_AUTOTUNE:-0}" == "1" ]]; then
        GPU_ARGS="--autotune"
        log "GPU detected — using --autotune"
    elif [[ "${HAS_GPU}" == "1" ]]; then
        GPU_ARGS=""
        log "GPU detected — using default batch size"
    else
        GPU_ARGS="--batch_size 1 --cpu"
        log "Running on CPU (PHOLD_PROTEINS_USE_GPU=${PHOLD_PROTEINS_USE_GPU:-0})."
        log "  ~50 proteins on 8 CPUs takes a few minutes — fine for batch jobs."
    fi

    phold proteins-predict \
        -i  "${BATCH_FAA}" \
        -o  "${PREDICT_DIR}" \
        -t  "${SLURM_CPUS_PER_TASK:-${THREADS}}" \
        ${GPU_ARGS} \
        -f
    touch "${PREDICT_DONE}"
    log "proteins-predict done."
fi

# ---------------------------------------------------------------------------
# Step 2: phold proteins-compare (search against phrog/AFDB/sub-DBs)
# ---------------------------------------------------------------------------
step_banner "phold proteins-compare (${BATCH_NAME})"

COMPARE_DONE="${COMPARE_DIR}/.compare_done"
if [[ -f "${COMPARE_DONE}" ]]; then
    log "Already done — skipping (delete ${COMPARE_DONE} to re-run)"
else
    FS_GPU=""
    [[ "${PHOLD_PROTEINS_USE_GPU:-0}" == "1" ]] && FS_GPU="--foldseek_gpu"

    phold proteins-compare \
        -i  "${BATCH_FAA}" \
        --predictions_dir "${PREDICT_DIR}" \
        -o  "${COMPARE_DIR}" \
        -t  "${SLURM_CPUS_PER_TASK:-${THREADS}}" \
        ${FS_GPU} \
        -f
    touch "${COMPARE_DONE}"
    log "proteins-compare done."
fi

# Touch done marker so re-submissions skip this batch
touch "${DONE_MARKER}"

log "Done: ${BATCH_NAME}  (${N_SEQS} proteins)  $(date)"
log "Next (after ALL array tasks complete): python scripts/01p_merge_phold_proteins.py"
