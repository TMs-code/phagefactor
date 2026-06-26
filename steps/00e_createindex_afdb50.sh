#!/usr/bin/env bash
# =============================================================================
# 00e_createindex_afdb50.sh — Precompute FoldSeek k-mer index for afdb50
# =============================================================================
# Run ONCE via sbatch — requires a high-memory node (~200G+).
# Creates afdb50_ss.idx (3Di k-mer index) needed for low-RAM searches.
# Without this index, foldseek search loads ~180G into RAM at search time.
#
# IMPORTANT: Run AFTER deleting any incomplete index:
#   rm afdb50_ss.idx.0 afdb50_ss.idx.index.0  (if they exist)
#
# Usage:
#   sbatch steps/00e_createindex_afdb50.sh
#   # adjust --partition and --mem to match your cluster's fat-node queue
# =============================================================================

#SBATCH -N 1
#SBATCH --partition=common        # ← change to your fat/large-mem partition name
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G                # 3Di index fill needs ~180G; 200G gives headroom
#SBATCH -J afdb50_createindex
#SBATCH -o logs/00e_createindex_afdb50.out
#SBATCH -e logs/00e_createindex_afdb50.err
#SBATCH --time=04:00:00
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
activate_env

AFDB50_DB="${FOLDSEEK_DB_ROOT}/afdb50_db/afdb50"
TMPDIR="${SCRATCH_TMPDIR:-/local/scratch/tmp}"
mkdir -p "${TMPDIR}"

echo "[$(date)] Starting foldseek createindex for afdb50"
echo "  DB       : ${AFDB50_DB}"
echo "  Tmp      : ${TMPDIR}"
echo "  Threads  : 16"
echo "  RAM limit: 200G SLURM / ~180G estimated by foldseek"
echo ""

# Check for incomplete index from previous run and clean up
for f in "${AFDB50_DB}_ss.idx.0" "${AFDB50_DB}_ss.idx.index.0"; do
    if [[ -f "$f" ]]; then
        echo "  Removing incomplete index fragment: $f"
        rm -f "$f"
    fi
done

# --index-exclude 2 = skip Cα coordinates from the index
# Cα is only needed for TM-score and --sort-by-structure-bits scoring.
# Our ProstT5 query DB has no _ca, so excluding it saves ~100G of index.
foldseek createindex \
    "${AFDB50_DB}" \
    "${TMPDIR}/fsindex_${SLURM_JOB_ID}" \
    --index-exclude 2 \
    --threads "${SLURM_CPUS_PER_TASK:-16}"

echo ""
echo "[$(date)] createindex complete"
echo ""

# Verify expected output files
echo "Index files created:"
for sfx in ".idx" ".idx.dbtype" ".idx.index" "_ss.idx" "_ss.idx.dbtype" "_ss.idx.index"; do
    f="${AFDB50_DB}${sfx}"
    if [[ -f "$f" ]]; then
        printf "  OK   %-25s %s\n" "${f##*/}" "$(du -sh "$f" | cut -f1)"
    else
        printf "  MISS %s\n" "${f##*/}"
    fi
done

# Clean tmp
rm -rf "${TMPDIR}/fsindex_${SLURM_JOB_ID}" 2>/dev/null || true

echo ""
echo "Next: sbatch steps/02_foldseek_3di.sh"
