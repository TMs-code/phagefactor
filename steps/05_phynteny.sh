#!/usr/bin/env bash
# =============================================================================
# 05_phynteny.sh — SLURM: Phynteny + synteny integration (merged step)
#
# Phynteny and integration are ONE step writing to ONE folder
# (05_phynteny/). It (1) runs Phynteny on the curated GenBank from the general
# output folder (04_output/updated_prophages.gb), then (2) integrates the
# phynteny categories + C1/Cro synteny hints into the final table.
#
#   Part 1 (phynteny env, python<3.11):  04_output/updated_prophages.gb
#       -> 05_phynteny/phynteny.tsv (+ GBK with phynteny qualifier)
#   Part 2 (pharokka env):  merge into 05_phynteny/final_annotations_integrated.csv
#       (+ 05_phynteny/final_with_synteny.gb)
#
# Optional-but-recommended: the final_annotations_table.csv from steps/04_curate.sh
# (05_build_output.py) is already complete; this step adds phynteny categories +
# synteny naming.
#
# Submit: sbatch steps/05_phynteny.sh
# =============================================================================

#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -J pf_phynteny
#SBATCH -o logs/05_phynteny.out
#SBATCH -e logs/05_phynteny.err
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
# --mail-user omitted on purpose: submit_all.sh passes it (SLURM_EMAIL). For a
# manual run, add --mail-user=you@your.institute (some clusters reject others).

set -euo pipefail

# Resolve project root: SLURM_SUBMIT_DIR if sbatch-submitted (BASH_SOURCE points
# into the spool dir once queued), else this script's own parent dir.
PROJECT_DIR="${PHAGEFACTOR_HOME:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
source "${PROJECT_DIR}/config.sh"

export TMPDIR="${SCRATCH_TMPDIR}"
mkdir -p "${TMPDIR}" "${PHYNTENY_DIR}"

INPUT_GBK="${OUTPUT_DIR}/updated_prophages.gb"

# Graceful skip if the curated GenBank is missing (e.g. 05_build_output.py could not build it).
if [[ ! -f "${INPUT_GBK}" ]]; then
    log "[SKIP] Required file not found: ${INPUT_GBK}"
    log "  Expected if 05_build_output.py (steps/04_curate.sh) could not locate/build the GenBank."
    log "  final_annotations_table.csv from steps/04_curate.sh is still valid and complete."
    exit 0
fi

# ---- Part 1: Phynteny (own env, python=3.10) --------------------------------
export MAMBA_ROOT_PREFIX
eval "$("${MAMBA_ROOT_PREFIX}/bin/micromamba" shell hook --shell bash)"
micromamba activate phynteny

step_banner "Step 05 (a) — Phynteny"
log "SLURM job: ${SLURM_JOB_ID:-none (local run)}  CPUs: ${SLURM_CPUS_PER_TASK}"

python "${SCRIPTS_DIR}/06_phynteny.py" \
    --input  "${INPUT_GBK}" \
    --output "${PHYNTENY_DIR}" \
    --threshold 0.8 \
    2>&1

micromamba deactivate || true
log "Phynteny complete: outputs in ${PHYNTENY_DIR}/"

# ---- Part 2: Integrate phynteny + synteny into the final table --------------
activate_env
step_banner "Step 05 (b) — Integrate phynteny + synteny"

PHY_TSV="${PHYNTENY_DIR}/phynteny.tsv"
if [[ ! -f "${PHY_TSV}" ]]; then
    log "[SKIP] ${PHY_TSV} not found — phynteny produced no table."
    log "  final_annotations_table.csv from steps/04_curate.sh is still valid."
    exit 0
fi

python "${SCRIPTS_DIR}/07_integrate.py" \
    --final    "${OUTPUT_DIR}/final_annotations_table.csv" \
    --gbk      "${INPUT_GBK}" \
    --phynteny "${PHY_TSV}" \
    --out      "${PHYNTENY_DIR}/final_annotations_integrated.csv" \
    --threshold 0.8

log "Done."
log "  Integrated table -> ${PHYNTENY_DIR}/final_annotations_integrated.csv"
log "  GBK with notes   -> ${PHYNTENY_DIR}/final_with_synteny.gb"
