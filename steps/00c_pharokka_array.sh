#!/usr/bin/env bash
# =============================================================================
# 00c_pharokka_array.sh — SLURM array job: pharokka on N prophages
#
# Runs pharokka.py with default Phanotate gene caller (catches more small phage
# ORFs than pyrodigal-gv).  Output GBKs feed into the next step (01_phold).
#
# Submit:
#   cd phagefactor/
#   N=$(wc -l < input/prophage_list.txt)
#   sbatch --array=0-$((N-1))%500 steps/00c_pharokka_array.sh
#
# Input  : input/fasta/<NAME>.fasta
# Output : 00c_pharokka/<NAME>/<NAME>.gbk
# =============================================================================

#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -J pf_pharokka_%a
#SBATCH -o logs/00c_pharokka_%a.out
#SBATCH -e logs/00c_pharokka_%a.err
#SBATCH --time=00:30:00
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

export TMPDIR="${SCRATCH_TMPDIR}"
mkdir -p "${TMPDIR}"

activate_env
require_tool pharokka.py

# ---------------------------------------------------------------------------
# Identify this task's prophage from the list file
# ---------------------------------------------------------------------------
# This is an ARRAY step: SLURM_ARRAY_TASK_ID must be set. The driver
# (submit_all.sh) submits it with --array=0-N; running it bare fails under set -u.
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    die "SLURM_ARRAY_TASK_ID not set — don't run this step directly. Use 'bash submit_all.sh', or for a manual run: sbatch --array=0-\$(( \$(wc -l < ${PROPHAGE_LIST}) - 1 )) steps/00c_pharokka_array.sh"
fi
PROPHAGE_NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${PROPHAGE_LIST}")
[[ -z "${PROPHAGE_NAME}" ]] && die "No prophage for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

FASTA_FILE="${FASTA_DIR}/${PROPHAGE_NAME}.fasta"
check_file "${FASTA_FILE}"
check_dir  "${PHAROKKA_DB}"

OUT_DIR="${PHAROKKA_OUT_DIR}/${PROPHAGE_NAME}"
mkdir -p "${OUT_DIR}"

# Skip if already done
DONE_MARKER="${OUT_DIR}/.pharokka_done"
if [[ -f "${DONE_MARKER}" ]]; then
    log "Skipping ${PROPHAGE_NAME} — already done"
    exit 0
fi

# ---------------------------------------------------------------------------
# Run pharokka  (Phanotate gene caller is the default — best for phages)
# ---------------------------------------------------------------------------
step_banner "Pharokka: ${PROPHAGE_NAME}  [array task ${SLURM_ARRAY_TASK_ID}]"
log "Input : ${FASTA_FILE}"
log "Output: ${OUT_DIR}/"
log "Threads: ${SLURM_CPUS_PER_TASK}"

# Run pharokka -- default Phanotate (NOT --fast which switches to pyrodigal-gv)
pharokka.py \
    -i "${FASTA_FILE}" \
    -o "${OUT_DIR}" \
    -d "${PHAROKKA_DB}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -p "${PROPHAGE_NAME}" \
    -f                              # force overwrite if dir exists

# Pharokka writes ${OUT_DIR}/${PROPHAGE_NAME}.gbk (the -p prefix becomes the
# basename).  01_phold_array.sh consumes it directly -- no symlink needed.
check_file "${OUT_DIR}/${PROPHAGE_NAME}.gbk"

touch "${DONE_MARKER}"
log "Done: ${PROPHAGE_NAME}  $(date)"

# One-line recap appended to a SINGLE shared file (so you don't have to open all
# N per-task logs). Per-task logs/01..00c_*.out remain only for debugging failures.
N_CDS=$(grep -c "^     CDS " "${OUT_DIR}/${PROPHAGE_NAME}.gbk" 2>/dev/null || echo "?")
printf "%s\t%s\tOK\t%s CDS\n" "$(date +%F_%T)" "${PROPHAGE_NAME}" "${N_CDS}" \
    >> "${PROJECT_DIR}/logs/00c_pharokka_recap.tsv"
