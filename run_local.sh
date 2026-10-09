#!/usr/bin/env bash
# run_local.sh — run the whole pipeline on ONE machine, no SLURM.
#
# submit_all.sh is a SLURM driver: every step is an sbatch call, and the two
# heavy steps are job ARRAYS. That makes the pipeline unrunnable on a laptop —
# and it is also why WebAPI search mode was unusable, since HPC compute nodes
# have no outbound internet while login nodes have no queue.
#
# This runner calls the SAME steps/*.sh scripts with plain bash, in dependency
# order, replacing each job array with a sequential loop over SLURM_ARRAY_TASK_ID.
# Nothing in steps/ needed changing: the #SBATCH lines are comments to bash.
#
# USAGE
#   bash run_local.sh                          # auto-detect mode, local FoldSeek
#   SEARCH_MODE=webapi bash run_local.sh       # no local DB (needs internet)
#   SEARCH_MODE=none   bash run_local.sh       # no FoldSeek at all (pharokka+phold)
#   INPUT_MODE=protein bash run_local.sh       # force protein mode
#   SKIP_PHAROKKA=1 bash run_local.sh          # resume from phold
#   PHAGEFACTOR_RUN_DIR=/path bash run_local.sh   # reuse an exact run
#   RESUME=1 bash run_local.sh                 # reuse the most recent run
#
# Sequential by design: on one machine the arrays are the bottleneck, and
# running them in parallel would oversubscribe the CPU that phold already wants.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PHAGEFACTOR_HOME="${SCRIPT_DIR}"

[[ -z "${SLURM_JOB_ID:-}" ]] || {
    echo "ERROR: do not sbatch run_local.sh — it IS the runner." >&2
    echo "       On a cluster use submit_all.sh instead." >&2; exit 1; }

SEARCH_MODE="${SEARCH_MODE:-local}"
SKIP_PHAROKKA="${SKIP_PHAROKKA:-0}"
SKIP_PHOLD="${SKIP_PHOLD:-0}"
SKIP_FOLDSEEK="${SKIP_FOLDSEEK:-0}"


# A run dir must exist before config.sh is sourced (it reads PHAGEFACTOR_RUN_DIR).
_RUN_BASE="${SCRIPT_DIR}/runs"
if [[ -z "${PHAGEFACTOR_RUN_DIR:-}" ]]; then
    if [[ "${RESUME:-0}" == "1" && -f "${_RUN_BASE}/.current_run" ]]; then
        export PHAGEFACTOR_RUN_DIR="$(cat "${_RUN_BASE}/.current_run")"
    else
        export PHAGEFACTOR_RUN_DIR="${_RUN_BASE}/local_$(date +%Y%m%d_%H%M%S)"
    fi
fi
mkdir -p "${PHAGEFACTOR_RUN_DIR}" "${_RUN_BASE}"
# Same pointer submit_all.sh writes: steps/05_phynteny.sh and the smoke test
# both follow it to find the run that just finished.
printf '%s\n' "${PHAGEFACTOR_RUN_DIR}" > "${_RUN_BASE}/.current_run"
source "${SCRIPT_DIR}/config.sh"

# CLR_ prefix on purpose: a bare `N` for the reset code collided with the `N`
# holding the task count, so every banner printed "[local]2" instead of resetting.
CLR_C='\033[0;36m'; CLR_G='\033[0;32m'; CLR_Y='\033[1;33m'
CLR_R='\033[0;31m'; CLR_0='\033[0m'
say()  { echo -e "${CLR_C}[local]${CLR_0} $*"; }
ok()   { echo -e "${CLR_G}[ ok ]${CLR_0}  $*"; }
warn() { echo -e "${CLR_Y}[warn]${CLR_0}  $*"; }
die()  { echo -e "${CLR_R}[FAIL]${CLR_0}  $*" >&2; exit 1; }

run_step() {   # run_step <label> <script> [extra env assignments...]
    local label="$1" script="$2"; shift 2
    local t0=$SECONDS
    say "${label}"
    env "$@" bash "${SCRIPT_DIR}/steps/${script}" \
        || die "${label} failed (steps/${script}). Log above."
    ok "${label}  [$((SECONDS - t0))s]"
}

run_array() {  # run_array <label> <script> <n_tasks>
    local label="$1" script="$2" n="$3" i t0=$SECONDS
    say "${label} — ${n} task(s), sequential"
    for (( i=0; i<n; i++ )); do
        printf '      task %d/%d\n' "$((i+1))" "${n}"
        SLURM_ARRAY_TASK_ID="${i}" bash "${SCRIPT_DIR}/steps/${script}" \
            || die "${label}: task ${i} failed (steps/${script})"
    done
    ok "${label}  [$((SECONDS - t0))s]"
}

echo
say "run dir     : ${PHAGEFACTOR_RUN_DIR}"
say "search mode : ${SEARCH_MODE}"

# --- input mode -------------------------------------------------------------
# Same call submit_all.sh makes: the flag is --fasta-dir (NOT --input-dir), and
# stderr is deliberately NOT swallowed -- the detector explains WHY it failed
# (no FASTAs, mixed genome/protein content), which a generic message hides.
if [[ -z "${INPUT_MODE:-}" ]]; then
    DETECT_PY="${SCRIPTS_DIR}/lib/detect_input_type.py"
    DETECT_ARGS=()
    [[ "${REDETECT:-0}" == "1" ]] && DETECT_ARGS+=(--redetect)
    INPUT_MODE=$(python3 "${DETECT_PY}" \
        --fasta-dir "${FASTA_DIR}" "${DETECT_ARGS[@]}") \
        || die "auto-detection failed (see error above). Debug with:
  python3 ${DETECT_PY} --fasta-dir ${FASTA_DIR} --verbose
or override: INPUT_MODE=genome bash run_local.sh"
fi
say "input mode  : ${INPUT_MODE}"

if [[ "${SEARCH_MODE}" == "webapi" ]]; then
    # Fail early rather than after Pharokka + phold have burned an hour.
    curl -fsS --max-time 15 -o /dev/null https://search.foldseek.com 2>/dev/null \
        || warn "search.foldseek.com is not reachable. WebAPI mode will fail at
       step 02; the server also rate-limits and has had outages. Use
       SEARCH_MODE=local if you have the DBs."
fi

# =============================================================================
# PROTEIN MODE
# =============================================================================
if [[ "${INPUT_MODE}" == "pangenome" ]]; then
    die "pangenome mode needs the translate pre-step; run it via submit_all.sh,
       or translate first and re-run with INPUT_MODE=protein."
fi

if [[ "${INPUT_MODE}" == "protein" ]]; then
    if [[ "${SKIP_PHOLD}" != "1" ]]; then
        say "splitting protein input into batches"
        ( cd "${SCRIPT_DIR}"
          export PROTEIN_BATCH_SIZE="${PROTEIN_BATCH_SIZE}"
          [[ -n "${PROTEIN_FASTA_DIR:-}" ]] && export PROTEIN_FASTA_DIR
          python3 scripts/00p_split_protein_batches.py ) \
            || die "protein split failed"
        [[ -f "${PROTEIN_BATCH_LIST}" ]] || die "no ${PROTEIN_BATCH_LIST}"
        N_TASKS=$(wc -l < "${PROTEIN_BATCH_LIST}")
        (( N_TASKS > 0 )) || die "0 batches in ${PROTEIN_BATCH_LIST}"
        run_array "phold (proteins)" "01p_phold_proteins_array.sh" "${N_TASKS}"
        run_step  "merge phold proteins" "01p_merge_phold_proteins.sh"
    fi
else
    # =========================================================================
    # GENOME MODE
    # =========================================================================
    [[ -f "${PROPHAGE_LIST}" ]] || die "missing ${PROPHAGE_LIST}"
    N_TASKS=$(wc -l < "${PROPHAGE_LIST}")
    (( N_TASKS > 0 )) || die "prophage_list.txt is empty: ${PROPHAGE_LIST}"

    if [[ "${SKIP_PHAROKKA}" != "1" && "${SKIP_PHOLD}" != "1" ]]; then
        run_array "pharokka (FASTA -> GBK)" "00c_pharokka_array.sh" "${N_TASKS}"
    fi
    if [[ "${SKIP_PHOLD}" != "1" ]]; then
        run_array "phold (genome)" "01_phold_array.sh" "${N_TASKS}"
        run_step  "merge phold"    "01c_merge.sh"
    fi
fi

# =============================================================================
# FoldSeek -> compare -> curate  (shared by both modes)
# =============================================================================
if [[ "${SKIP_FOLDSEEK}" != "1" ]]; then
    case "${SEARCH_MODE}" in
      webapi) run_step "FoldSeek (web API)" "02w_foldseek_webapi.sh" ;;
      none)
        # No custom-FoldSeek evidence at all: pharokka + phold only. 03 REQUIRES
        # best_hit.csv and exits 1 if it is absent, so write a well-formed EMPTY
        # one rather than skipping the step. Useful when you have neither the
        # local DBs nor a working web API -- you still get the pharokka trusted
        # layer and phold's structural calls through the curation rules.
        say "SEARCH_MODE=none — no FoldSeek evidence (pharokka + phold only)"
        mkdir -p "${FOLDSEEK_DIR}/3di_tokens"
        python3 - "${FOLDSEEK_DIR}/3di_tokens" <<'PYEOF'
import sys, pandas as pd
from pathlib import Path
d = Path(sys.argv[1])
cols = ["gene", "accession", "description", "taxname", "pident", "evalue",
        "score", "composite_score", "qcov_aa", "same_host", "defense_flag",
        "informative_hit_found", "aa_length", "foldseek_subdb"]
pd.DataFrame(columns=cols).to_csv(d / "best_hit.csv", index=False)
pd.DataFrame(columns=["gene", "rank", "description", "accession",
                      "evalue", "score"]).to_csv(d / "top3.csv", index=False)
print(f"  wrote empty best_hit.csv / top3.csv -> {d}")
PYEOF
        ;;
      *)      run_step "FoldSeek (local DBs)" "02_foldseek_3di.sh" ;;
    esac
    run_step "compare (03)" "03_compare.sh"
fi
# 04 runs curation AND output-building in one script
run_step "curate + build output (04)" "04_curate.sh"

echo
ok "pipeline finished"
echo "    results : ${RUN_DIR}/04_output/final_annotations_table.csv"
echo "    review  : ${RUN_DIR}/04_output/curation/review_suggested.csv"
echo
echo "    optional synteny/phynteny step (separate env):"
echo "      bash steps/05_phynteny.sh"
