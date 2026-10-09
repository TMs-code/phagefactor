#!/usr/bin/env bash
# =============================================================================
# 02w_then_curate.sh — LOGIN-NODE one-shot for WEBAPI mode.
#
# In webAPI mode, submit_all.sh only submits pharokka -> phold -> merge, then
# stops, because the FoldSeek web search (02w) needs INTERNET and must run on the
# login node (compute nodes have none). This script does the rest in one go:
#   1. run the web search (02w)              -> best_hit.csv / top3.csv
#   2. verify the search actually produced results
#   3. submit compare -> curate -> (optional) phynteny as SLURM jobs
#
# Run it on the LOGIN node AFTER the phold merge job has finished:
#     bash steps/02w_then_curate.sh
# It follows the most recent run automatically (RUN_BASE/.current_run); to target
# a specific run instead:  PHAGEFACTOR_RUN_DIR=/path/to/run bash steps/02w_then_curate.sh
# Skip the optional phynteny step with:  RUN_PHYNTENY=0 bash steps/02w_then_curate.sh
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
source "${PROJECT_DIR}/config.sh"    # resolves RUN_DIR / FOLDSEEK_DIR / OUTPUT_DIR

echo "Run dir : ${RUN_DIR}"
echo "Search  : webAPI (search.foldseek.com, login node)"

# 1) the web search (writes best_hit.csv / top3.csv under 02_foldseek/3di_tokens/)
bash "${SCRIPT_DIR}/02w_foldseek_webapi.sh"

# 2) gate: do NOT chain 03/04 if the search produced nothing (e.g. server outage)
BEST="${FOLDSEEK_DIR}/3di_tokens/best_hit.csv"
if [[ ! -s "${BEST}" ]]; then
    echo "ERROR: ${BEST} is missing/empty — the web search did not complete" >&2
    echo "       (server outage, or all proteins failed). Not submitting 03/04." >&2
    echo "       Try SEARCH_MODE=local with on-disk databases instead." >&2
    exit 1
fi
echo "OK: ${BEST} present -> chaining compare/curate."

# 3) compare -> curate -> optional phynteny (SLURM; --chdir keeps logs in the run)
J03=$(sbatch --parsable --chdir="${RUN_DIR}" "${SCRIPT_DIR}/03_compare.sh")
echo "[SUBMITTED] 03_compare   ${J03}"
J04=$(sbatch --parsable --chdir="${RUN_DIR}" --dependency=afterok:${J03} "${SCRIPT_DIR}/04_curate.sh")
echo "[SUBMITTED] 04_curate    ${J04}  (after ${J03})"
if [[ "${RUN_PHYNTENY:-1}" == "1" ]]; then
    J05=$(sbatch --parsable --chdir="${RUN_DIR}" --dependency=afterok:${J04} "${SCRIPT_DIR}/05_phynteny.sh")
    echo "[SUBMITTED] 05_phynteny  ${J05}  (after ${J04}, optional)"
fi
echo "Done. Final outputs will land in: ${OUTPUT_DIR}/"
