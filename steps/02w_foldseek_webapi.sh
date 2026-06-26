#!/usr/bin/env bash
# =============================================================================
# 02w_foldseek_webapi.sh — FoldSeek search via the PUBLIC WEB SERVER (no local DB)
# Alternative to step 02 (local 02d) for users without the afdb50 DB + index.
#
# ⚠️ RUN ON THE LOGIN/SUBMIT NODE, NOT via sbatch: it makes HTTP calls to
#    search.foldseek.com, and cluster COMPUTE nodes usually have no internet.
#    Usage:  bash steps/02w_foldseek_webapi.sh        (after the phold merge)
#
# Submits the phold AA FASTA (server runs ProstT5 -> 3Di), then bridges the
# result to best_hit.csv / top3.csv using the SAME scoring path as local mode,
# so step 03 is identical for both search modes.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
source "${PROJECT_DIR}/config.sh"
activate_env

FS_DIR="${FOLDSEEK_DIR}/3di_tokens"
mkdir -p "${FS_DIR}"

if [[ ! -f "${PHOLD_AA_FASTA}" ]]; then
    echo "ERROR: ${PHOLD_AA_FASTA} not found — run phold + the merge (01c) first." >&2
    exit 1
fi

step_banner "Step 02w — FoldSeek WEB API search (no local DB)"
python "${SCRIPTS_DIR}/02w_foldseek_webapi.py" \
    --fasta    "${PHOLD_AA_FASTA}" \
    --out      "${FS_DIR}/webapi_hits.m8" \
    --best-hit "${FS_DIR}/best_hit.csv" \
    --top3     "${FS_DIR}/top3.csv" \
    --all-hits "${FS_DIR}/all_hits.csv"

log "Done. best_hit.csv / top3.csv written -> step 03 can run:"
log "  sbatch steps/03_compare.sh"
