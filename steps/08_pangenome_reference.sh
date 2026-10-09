#!/usr/bin/env bash
# =============================================================================
# 08_pangenome_reference.sh — PAN-GENOME sub-mode post-processing (after curate).
#
# Rejoins final_annotations_table.csv (keyed on the clean pang_NNNNN ids) back to
# the biological group names + protein sequences via pangenome_id_map.csv, and
# emits the homology-reference set (CSV + XLSX + a DIAMOND/MMseqs-ready .faa)
# under 04_output/. See scripts/08_pangenome_reference.py.
#
# Runs ONLY when an id map is present (i.e. the input came through
# 00n_translate_pangenome.py). submit_all.sh submits it automatically after
# 04_curate for pan-genome runs; it can also be run by hand:
#   sbatch --chdir=<run_dir> --dependency=afterok:<curate_jobid> steps/08_pangenome_reference.sh
# =============================================================================

#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -J pangenome_reference
#SBATCH -o logs/08_pangenome_reference.out
#SBATCH -e logs/08_pangenome_reference.err
#SBATCH --time=00:30:00
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# Resolve project root: PHAGEFACTOR_HOME (exported by submit_all) if present,
# else SLURM_SUBMIT_DIR, else this script's own parent dir.
PROJECT_DIR="${PHAGEFACTOR_HOME:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
source "${PROJECT_DIR}/config.sh"

activate_env

FINAL_TABLE="${OUTPUT_DIR}/final_annotations_table.csv"
ID_MAP="${PANGENOME_ID_MAP:-${INPUT_DIR}/pangenome_id_map.csv}"
# Clean-id protein FASTA. The run-local combined faa (headers 'PROPH:pang_NNNNN')
# is always present after the split; 08.py matches the id map on the ':' suffix.
REF_FAA="${PROTEIN_COMBINED_FAA}"
OUT_PREFIX="${OUTPUT_DIR}/pangenome_reference"

if [[ ! -f "${ID_MAP}" ]]; then
    log "[SKIP] No pangenome id map at ${ID_MAP} — not a pan-genome run; nothing to do."
    exit 0
fi
if [[ ! -f "${FINAL_TABLE}" ]]; then
    log "[ERROR] ${FINAL_TABLE} not found — did 04_curate finish?"
    exit 1
fi

step_banner "Step 08 — pan-genome homology-reference (group + sequence + annotation)"
log "final=${FINAL_TABLE}"
log "id_map=${ID_MAP}"
log "faa=${REF_FAA}"

python3 "${SCRIPTS_DIR}/08_pangenome_reference.py" \
    --final    "${FINAL_TABLE}" \
    --id-map   "${ID_MAP}" \
    --faa      "${REF_FAA}" \
    --out-prefix "${OUT_PREFIX}"

log "Done -> ${OUT_PREFIX}.{csv,xlsx,faa}"
log "  Build a homology DB:   diamond makedb --in ${OUT_PREFIX}.faa -d pangenome_ref"
log "  Transfer annotations:  diamond blastp -q <your_proteins.faa> -d pangenome_ref ... "
log "    (each hit's sseqid == group_id; join to ${OUT_PREFIX}.csv for the annotation)"
