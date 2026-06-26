#!/usr/bin/env bash
# =============================================================================
# 01d_merge_3di.sh — Merge per-prophage phold 3di FASTA files
#
# Run AFTER the phold array (01_phold_array.sh) completes.
# Concatenates all per-prophage phold_3di.fasta and phold_aa.fasta files
# into single combined files for the FoldSeek 3di search step.
#
# The combined files are the input to:
#   scripts/02d_foldseek_3di.py  (PHOLD_3DI_FASTA, PHOLD_AA_FASTA in config/config.yaml)
#
# Run interactively (fast operation):
#   bash steps/01d_merge_3di.sh
# =============================================================================
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
fi
source "${PROJECT_DIR}/config.sh"

step_banner "Merging phold 3di and AA FASTA files"

COMBINED_DIR="${PHOLD_OUT_DIR}/combined"
mkdir -p "${COMBINED_DIR}"

# ---- 3di tokens FASTA -------------------------------------------------------
# Phold v1.2 with --prefix <NAME> writes to: 01_phold/<NAME>/<NAME>_3di.fasta
# (and <NAME>_aa.fasta).  We also check the older "phold_3di.fasta" name as a
# fallback in case phold was run without --prefix.
OUT_3DI="${COMBINED_DIR}/phold_3di.fasta"
log "Merging phold 3di FASTA files → ${OUT_3DI}"

> "${OUT_3DI}"   # empty/create
FOUND=0
while IFS= read -r PROPHAGE; do
    CANDIDATE_1="${PHOLD_OUT_DIR}/${PROPHAGE}/${PROPHAGE}_3di.fasta"
    CANDIDATE_2="${PHOLD_OUT_DIR}/${PROPHAGE}/${PROPHAGE}_phold_3di.fasta"
    CANDIDATE_3="${PHOLD_OUT_DIR}/${PROPHAGE}/phold_3di.fasta"
    if   [[ -f "${CANDIDATE_1}" ]]; then cat "${CANDIDATE_1}" >> "${OUT_3DI}"; FOUND=$((FOUND+1))
    elif [[ -f "${CANDIDATE_2}" ]]; then cat "${CANDIDATE_2}" >> "${OUT_3DI}"; FOUND=$((FOUND+1))
    elif [[ -f "${CANDIDATE_3}" ]]; then cat "${CANDIDATE_3}" >> "${OUT_3DI}"; FOUND=$((FOUND+1))
    else warn "3di FASTA not found for ${PROPHAGE} (checked ${CANDIDATE_1}, ${CANDIDATE_2}, ${CANDIDATE_3})"
    fi
done < "${PROPHAGE_LIST}"
log "3di FASTA: merged ${FOUND} prophages → $(grep -c '^>' "${OUT_3DI}") sequences"

# ---- AA FASTA ---------------------------------------------------------------
OUT_AA="${COMBINED_DIR}/phold_aa.fasta"
log "Merging phold AA FASTA files → ${OUT_AA}"

> "${OUT_AA}"
FOUND=0
while IFS= read -r PROPHAGE; do
    CANDIDATE_1="${PHOLD_OUT_DIR}/${PROPHAGE}/${PROPHAGE}_aa.fasta"
    CANDIDATE_2="${PHOLD_OUT_DIR}/${PROPHAGE}/${PROPHAGE}_phold_aa.fasta"
    CANDIDATE_3="${PHOLD_OUT_DIR}/${PROPHAGE}/phold_aa.fasta"
    if   [[ -f "${CANDIDATE_1}" ]]; then cat "${CANDIDATE_1}" >> "${OUT_AA}"; FOUND=$((FOUND+1))
    elif [[ -f "${CANDIDATE_2}" ]]; then cat "${CANDIDATE_2}" >> "${OUT_AA}"; FOUND=$((FOUND+1))
    elif [[ -f "${CANDIDATE_3}" ]]; then cat "${CANDIDATE_3}" >> "${OUT_AA}"; FOUND=$((FOUND+1))
    else warn "AA FASTA not found for ${PROPHAGE}"
    fi
done < "${PROPHAGE_LIST}"
log "AA FASTA: merged ${FOUND} prophages → $(grep -c '^>' "${OUT_AA}") sequences"

log "Done. Next: SKIP_PHOLD=1 bash submit_all.sh   (chains 02 FoldSeek → 03 compare → 04 curate)"
