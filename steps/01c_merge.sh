#!/usr/bin/env bash
# =============================================================================
# 01c_merge.sh — SLURM: merge per-prophage phold outputs (TSV + 3Di/AA FASTA)
# Runs AFTER the phold array (01_phold_array.sh) completes, BEFORE FoldSeek (02).
# Wraps scripts/01c_merge_phold.py, then concatenates the per-prophage 3Di/AA
# FASTAs itself, so submit_all can chain 01 -> 01c_merge -> 02 automatically
# (no manual SKIP_PHOLD re-run).
# =============================================================================
#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -J phold_merge
#SBATCH -o logs/01c_merge.out
#SBATCH -e logs/01c_merge.err
#SBATCH --time=00:20:00

set -euo pipefail
# Resolve project root: SLURM_SUBMIT_DIR if sbatch-submitted (BASH_SOURCE points
# into the spool dir once queued), else this script's own parent dir.
PROJECT_DIR="${PHAGEFACTOR_HOME:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
source "${PROJECT_DIR}/config.sh"
activate_env
step_banner "Step 01c — Merge phold outputs"

python "${SCRIPTS_DIR}/01c_merge_phold.py"     # per-prophage TSVs -> phold_all.tsv

# ---- Merge per-prophage phold 3Di + AA FASTAs into combined/ ---------------
# Phold v1.2 with --prefix <NAME> writes to: 01_phold/<NAME>/<NAME>_3di.fasta
# (and <NAME>_aa.fasta). Also checks older naming as a fallback in case phold
# was run without --prefix. Combined files feed scripts/02d_foldseek_3di.py
# (PHOLD_3DI_FASTA, PHOLD_AA_FASTA).
COMBINED_DIR="${PHOLD_OUT_DIR}/combined"
mkdir -p "${COMBINED_DIR}"

_merge_phold_fasta() {
    local suffix="$1" out="$2" found=0 prophage c1 c2 c3
    log "Merging phold ${suffix} FASTA files → ${out}"
    > "${out}"
    while IFS= read -r prophage; do
        c1="${PHOLD_OUT_DIR}/${prophage}/${prophage}_${suffix}.fasta"
        c2="${PHOLD_OUT_DIR}/${prophage}/${prophage}_phold_${suffix}.fasta"
        c3="${PHOLD_OUT_DIR}/${prophage}/phold_${suffix}.fasta"
        if   [[ -f "${c1}" ]]; then cat "${c1}" >> "${out}"; found=$((found + 1))
        elif [[ -f "${c2}" ]]; then cat "${c2}" >> "${out}"; found=$((found + 1))
        elif [[ -f "${c3}" ]]; then cat "${c3}" >> "${out}"; found=$((found + 1))
        else warn "${suffix} FASTA not found for ${prophage} (checked ${c1}, ${c2}, ${c3})"
        fi
    done < "${PROPHAGE_LIST}"
    log "${suffix} FASTA: merged ${found} prophages → $(grep -c '^>' "${out}") sequences"
}

_merge_phold_fasta "3di" "${COMBINED_DIR}/phold_3di.fasta"
_merge_phold_fasta "aa"  "${COMBINED_DIR}/phold_aa.fasta"

log "Done. FoldSeek (step 02) can now run."
