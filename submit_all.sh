#!/usr/bin/env bash
# =============================================================================
# submit_all.sh — Orchestrated SLURM submission for phageFACTor
#
# AUTO-DETECTS input type from input/fasta/ and selects the correct pipeline:
#
#   GENOME MODE  (sequences > 5kb, nucleotide):
#     00c (pharokka, FASTA → GBK) → 01  (phold genome)
#     → [manual merge] → 02–06
#
#   PROTEIN MODE (sequences are AA, < 5kb):
#     00p split (fixed PROTEIN_BATCH_SIZE-protein batches, inline)
#     → 01p phold-proteins ARRAY (one batch per task, CPU-default)
#     → 01p merge (afterok) → [automatic] → 02–06
#     All input proteins treated as hypothetical.
#
# Usage:
#   cd phagefactor/
#   bash submit_all.sh                          # auto-detect input type
#   INPUT_MODE=genome bash submit_all.sh        # force genome mode
#   INPUT_MODE=protein bash submit_all.sh       # force protein mode
#   SKIP_PHAROKKA=1 bash submit_all.sh          # genome mode, skip pharokka
#   SKIP_PHOLD=1   bash submit_all.sh           # continue from step 02
#   REDETECT=1 bash submit_all.sh               # ignore cached .input_type
# =============================================================================

set -euo pipefail

# submit_all.sh is a DRIVER: run it with `bash submit_all.sh` from the pipeline
# dir (it issues the sbatch calls itself). If sbatch'd, SLURM copies it to
# /var/spool/slurmd/<job>/ and `source config.sh` below cannot find config.sh.
if [[ "${BASH_SOURCE[0]}" == /var/spool/* || -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: do not 'sbatch submit_all.sh' — it is a driver." >&2
    echo "       Run:  micromamba activate pharokka && bash submit_all.sh" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

SKIP_PHAROKKA="${SKIP_PHAROKKA:-0}"
SKIP_PHOLD="${SKIP_PHOLD:-0}"
REDETECT="${REDETECT:-0}"
# SEARCH_MODE=local (default, step 02 local FoldSeek/afdb50) | webapi (no local DB,
# search via search.foldseek.com on the LOGIN node — see step-02 section below).
SEARCH_MODE="${SEARCH_MODE:-local}"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[SUBMITTED]${NC}  $*"; }
info() { echo -e "              $*"; }
mode() { echo -e "${CYAN}[MODE]${NC}       $*"; }

mkdir -p "${LOGS_DIR}"

step_banner "phageFACTor — submitting SLURM jobs"

# ---------------------------------------------------------------------------
# Auto-detect input type (genome vs protein)
# Can be overridden with INPUT_MODE=genome or INPUT_MODE=protein env var.
# ---------------------------------------------------------------------------
if [[ -n "${INPUT_MODE:-}" ]]; then
    mode "INPUT_MODE forced to '${INPUT_MODE}' by environment variable"
elif [[ "${SKIP_PHOLD}" == "1" ]]; then
    # When continuing from step 02, check cache to know which mode was used
    CACHE="${INPUT_DIR}/.input_type"
    if [[ -f "${CACHE}" ]]; then
        INPUT_MODE=$(cat "${CACHE}")
        mode "Resuming from step 02 — cached input type: ${INPUT_MODE}"
    else
        INPUT_MODE="genome"
        warn "No .input_type cache found; assuming genome mode for step 02+ resume"
    fi
else
    DETECT_ARGS=()
    [[ "${REDETECT}" == "1" ]] && DETECT_ARGS+=(--redetect)

    # NOTE 2026-06: was calling bare `python` with stderr silenced — on
    # some submit nodes only `python3` exists (no python->python3 symlink),
    # so this died with exit 127 *silently* (2>/dev/null hid "command not
    # found", and `set -e` then killed the script before any [MODE] line).
    # Fixed: use python3 explicitly, and surface failures with a clear message
    # instead of swallowing stderr.
    INPUT_MODE=$(python3 "${SCRIPTS_DIR}/detect_input_type.py" \
        --fasta-dir "${FASTA_DIR}" \
        "${DETECT_ARGS[@]}") \
        || die "Auto-detection failed (see error above). Debug with:
  python3 ${SCRIPTS_DIR}/detect_input_type.py --fasta-dir ${FASTA_DIR} --verbose
or override explicitly: INPUT_MODE=protein bash submit_all.sh"

    mode "Auto-detected input type: ${INPUT_MODE}"
fi

echo ""
if [[ "${INPUT_MODE}" == "protein" ]]; then
    echo "  ┌─────────────────────────────────────────────────────────┐"
    echo "  │  PROTEIN MODE: skipping Pharokka; running phold proteins│"
    echo "  │  All input sequences treated as hypothetical candidates  │"
    echo "  └─────────────────────────────────────────────────────────┘"
else
    echo "  ┌─────────────────────────────────────────────────────────┐"
    echo "  │  GENOME MODE: Pharokka → phold genome → FoldSeek        │"
    echo "  └─────────────────────────────────────────────────────────┘"
fi
echo ""

# ---------------------------------------------------------------------------
# Common sbatch overrides
# ---------------------------------------------------------------------------
MAIL_OVERRIDES=( --mail-user="${EMAIL}" --mail-type="${MAIL_TYPE:-BEGIN,END,FAIL}" )

# ═══════════════════════════════════════════════════════════════════════════════
# PROTEIN MODE: split → phold-proteins array (50-protein batches) → merge →
#               then steps 02–06 chained
#
# Rewritten 2026-06: replaces the old single-job GPU submission
# (01p_phold_proteins.sh) with a fixed-size batch array, mirroring the
# genome-mode split → 01_phold_array → 01c_merge_phold chain. User decisions:
#   - fixed PROTEIN_BATCH_SIZE-protein batches (default 50, config.sh —
#     "dont divide by 4 ... reproducible for next datasets ... ok if last
#     batch is not full")
#   - CPU by default for the initial submission (PHOLD_PROTEINS_USE_GPU=0,
#     decoupled from genome-mode PHOLD_USE_GPU — "if genome mode used CPU,
#     then perhaps we should have the initial submission as CPU too ...
#     would avoid GPU bug")
# ═══════════════════════════════════════════════════════════════════════════════
if [[ "${INPUT_MODE}" == "protein" ]]; then

    DEP_FS=""
    if [[ "${SKIP_PHOLD}" != "1" ]]; then

        # --- Step 0: split input proteins into fixed-size batches (inline,
        #     not sbatch — pure I/O, needed up front to compute array size) ---
        step_banner "Splitting protein-mode input into ${PROTEIN_BATCH_SIZE}-protein batches"
        # PROTEIN_BATCH_SIZE / PROTEIN_FASTA_DIR are bash-only config.sh values
        # (no python-config equivalent) — export so the python split step
        # sees them; PROTEIN_FASTA_DIR may be unset (script then defaults to
        # ${INPUT_DIR}/fasta, same as the old 01p_phold_proteins.sh default).
        (
            cd "${SCRIPT_DIR}"
            export PROTEIN_BATCH_SIZE="${PROTEIN_BATCH_SIZE}"
            if [[ -n "${PROTEIN_FASTA_DIR:-}" ]]; then
                export PROTEIN_FASTA_DIR
            fi
            python3 scripts/00p_split_protein_batches.py
        )

        BATCH_LIST="${INPUT_DIR}/protein_batch_list.txt"
        [[ -f "${BATCH_LIST}" ]] || { echo "ERROR: ${BATCH_LIST} not written by split step" >&2; exit 1; }
        N_BATCHES=$(wc -l < "${BATCH_LIST}")
        [[ "${N_BATCHES}" -gt 0 ]] || { echo "ERROR: 0 batches in ${BATCH_LIST}" >&2; exit 1; }
        info "Protein input split into ${N_BATCHES} batch(es) of <= ${PROTEIN_BATCH_SIZE} proteins"

        # --- Step 1: SLURM array job, one task per batch ---
        PHOLD_P_ARGS=( --parsable --array=0-$((N_BATCHES - 1))%500 "${MAIL_OVERRIDES[@]}" )

        if [[ "${PHOLD_PROTEINS_USE_GPU:-0}" == "1" ]]; then
            PHOLD_P_ARGS+=(
                --partition="${PHOLD_GPU_PARTITION}"
                --qos="${PHOLD_GPU_QOS}"
                --gres="${PHOLD_GPU_GRES}"
                --time="${TIME_PHOLD}"
            )
            info "phold proteins array will run on GPU (${PHOLD_GPU_PARTITION})"
        else
            info "phold proteins array will run on CPU (PHOLD_PROTEINS_USE_GPU=0 — default;"
            info "  decoupled from genome-mode PHOLD_USE_GPU, see config.sh comment)"
        fi

        JOB_PHOLD_P=$(sbatch "${PHOLD_P_ARGS[@]}" "${SCRIPT_DIR}/steps/01p_phold_proteins_array.sh")
        ok "01p_phold_proteins_array.sh  JOBID=${JOB_PHOLD_P}  (array 0-$((N_BATCHES - 1)), ${N_BATCHES} batches)"

        # --- Step 2: merge job, runs after ALL array tasks succeed ---
        JOB_MERGE=$(sbatch --parsable \
            --dependency="afterok:${JOB_PHOLD_P}" \
            "${MAIL_OVERRIDES[@]}" \
            "${SCRIPT_DIR}/steps/01p_merge_phold_proteins.sh")
        ok "01p_merge_phold_proteins.sh  JOBID=${JOB_MERGE}  (after array ${JOB_PHOLD_P})"

        DEP_FS="afterok:${JOB_MERGE}"
    else
        info "SKIP_PHOLD=1 — skipping phold proteins, going straight to step 02"
        DEP_FS=""
    fi

    # Steps 02–06 (same as genome mode, chained from the merge job)
    FS_ARGS=( --parsable "${MAIL_OVERRIDES[@]}" )
    [[ -n "${DEP_FS}" ]] && FS_ARGS+=( --dependency="${DEP_FS}" )

    JOB_FS=$(sbatch "${FS_ARGS[@]}" "${SCRIPT_DIR}/steps/02_foldseek_3di.sh")
    ok "02_foldseek_3di.sh    JOBID=${JOB_FS}"

    JOB_COMPARE=$(sbatch --parsable \
        --dependency=afterok:${JOB_FS} \
        "${MAIL_OVERRIDES[@]}" \
        "${SCRIPT_DIR}/steps/03_compare.sh")
    ok "03_compare.sh         JOBID=${JOB_COMPARE}  (after ${JOB_FS})"

    JOB_CURATE=$(sbatch --parsable \
        --dependency=afterok:${JOB_COMPARE} \
        "${MAIL_OVERRIDES[@]}" \
        "${SCRIPT_DIR}/steps/04_curate.sh")
    ok "04_curate.sh          JOBID=${JOB_CURATE}  (after ${JOB_COMPARE})"

    # 2026-06: curation is now near-fully automatic (few review rows), so chain
    # the build-output step automatically. review_suggested.csv is emitted in
    # 04_output/curation/ for any genes that still want a human look -- no blocking gate.
    JOB_OUTPUT=$(sbatch --parsable \
        --dependency=afterok:${JOB_CURATE} \
        "${MAIL_OVERRIDES[@]}" \
        "${SCRIPT_DIR}/steps/04_output.sh")
    ok "04_output.sh          JOBID=${JOB_OUTPUT}  (after ${JOB_CURATE})"

    echo ""
    echo -e "${YELLOW}[OPTIONAL, recommended]${NC} phynteny + synteny integration (one merged step):"
    echo "    sbatch --dependency=afterok:${JOB_OUTPUT} steps/05_phynteny.sh"
    echo "  (04_output/final_annotations_table.csv is already complete after step 04;"
    echo "   05_phynteny adds phynteny categories + C1/Cro synteny naming.)"
    echo ""
    step_banner "Protein mode submission complete — squeue -u \$USER"
    exit 0
fi

# ═══════════════════════════════════════════════════════════════════════════════
# GENOME MODE (original pipeline)
# ═══════════════════════════════════════════════════════════════════════════════

# ---------------------------------------------------------------------------
# Determine array size (used for both pharokka and phold arrays)
# ---------------------------------------------------------------------------
N_PROPHAGES=$(wc -l < "${PROPHAGE_LIST}")
[[ "${N_PROPHAGES}" -gt 0 ]] || die "prophage_list.txt is empty: ${PROPHAGE_LIST}"
ARRAY_SPEC="0-$((N_PROPHAGES - 1))%500"

# ---------------------------------------------------------------------------
# Step 00c — pharokka array (FASTA -> GBK)
# ---------------------------------------------------------------------------
JOB_PHAROKKA=""
if [[ "${SKIP_PHAROKKA}" != "1" && "${SKIP_PHOLD}" != "1" ]]; then
    JOB_PHAROKKA=$(sbatch \
        --array="${ARRAY_SPEC}" \
        --parsable \
        "${MAIL_OVERRIDES[@]}" \
        "${SCRIPT_DIR}/steps/00c_pharokka_array.sh")
    ok "00c_pharokka_array.sh  JOBID=${JOB_PHAROKKA}  (${N_PROPHAGES} tasks)"
fi

# ---------------------------------------------------------------------------
# Step 01 — phold array (GBK -> 3Di + annotations)
# ---------------------------------------------------------------------------
if [[ "${SKIP_PHOLD}" != "1" ]]; then
    PHOLD_SBATCH_ARGS=( --array="${ARRAY_SPEC}" --parsable "${MAIL_OVERRIDES[@]}" )

    if [[ "${PHOLD_USE_GPU:-0}" == "1" ]]; then
        PHOLD_SBATCH_ARGS+=(
            --partition="${PHOLD_GPU_PARTITION}"
            --qos="${PHOLD_GPU_QOS}"
            --gres="${PHOLD_GPU_GRES}"
            --time="${TIME_PHOLD}"
        )
        info "Phold will run on GPU partition (${PHOLD_GPU_PARTITION}, gres=${PHOLD_GPU_GRES})"
    fi

    if [[ -n "${JOB_PHAROKKA}" ]]; then
        PHOLD_SBATCH_ARGS+=( --dependency="aftercorr:${JOB_PHAROKKA}" )
        info "Phold depends on pharokka job ${JOB_PHAROKKA}"
    fi

    JOB_PHOLD=$(sbatch "${PHOLD_SBATCH_ARGS[@]}" "${SCRIPT_DIR}/steps/01_phold_array.sh")
    ok "01_phold_array.sh      JOBID=${JOB_PHOLD}  (${N_PROPHAGES} tasks)"

    # Auto-chain the phold merge (01c+01d) after ALL phold array tasks finish, so
    # the user no longer needs a manual merge + 'SKIP_PHOLD=1' re-run.
    # aftercorr would be per-task; we need the WHOLE array done -> afterok:<arrayjob>.
    JOB_MERGE=$(sbatch --parsable \
        --dependency="afterok:${JOB_PHOLD}" \
        "${MAIL_OVERRIDES[@]}" \
        "${SCRIPT_DIR}/steps/01c_merge.sh")
    ok "01c_merge.sh           JOBID=${JOB_MERGE}  (after phold array ${JOB_PHOLD})"
    DEP_FS="afterok:${JOB_MERGE}"   # step 02 below waits on the merge
fi

# ---------------------------------------------------------------------------
# Steps 02–04 — FoldSeek, compare, curate (auto-chained)
# Reached either by falling through from the phold block (DEP_FS = the merge job)
# or directly via SKIP_PHOLD=1 (DEP_FS empty -> step 02 starts immediately).
# ---------------------------------------------------------------------------
echo ""

# --- WEBAPI search mode: login-node step, not part of the sbatch chain --------
if [[ "${SEARCH_MODE}" == "webapi" ]]; then
    echo -e "${YELLOW}[WEBAPI MODE]${NC} Step 02 runs against search.foldseek.com (no local DB)."
    echo "  It must run on the LOGIN node (compute nodes have no internet) AFTER the"
    echo "  phold merge${JOB_MERGE:+ (job ${JOB_MERGE})} completes:"
    echo "    bash steps/02w_foldseek_webapi.sh        # writes best_hit.csv/top3.csv"
    echo "  then chain compare -> curate:"
    echo "    J03=\$(sbatch --parsable steps/03_compare.sh)"
    echo "    J04=\$(sbatch --parsable --dependency=afterok:\$J03 steps/04_curate.sh)"
    echo "    J04b=\$(sbatch --parsable --dependency=afterok:\$J04 steps/04_output.sh)"
    echo "    sbatch --dependency=afterok:\$J04b steps/05_phynteny.sh   # optional"
    echo ""
    step_banner "Submitted through phold merge. WebAPI search 02w is a manual login-node step."
    exit 0
fi

info "Submitting steps 02–04 with dependency chain..."
echo ""

FS_DEP_ARGS=()
[[ -n "${DEP_FS:-}" ]] && FS_DEP_ARGS+=( --dependency="${DEP_FS}" )
JOB_FS=$(sbatch --parsable "${FS_DEP_ARGS[@]}" "${MAIL_OVERRIDES[@]}" "${SCRIPT_DIR}/steps/02_foldseek_3di.sh")
ok "02_foldseek_3di.sh  JOBID=${JOB_FS}${DEP_FS:+  (after ${DEP_FS#afterok:})}"

JOB_COMPARE=$(sbatch --parsable \
    --dependency=afterok:${JOB_FS} \
    "${MAIL_OVERRIDES[@]}" \
    "${SCRIPT_DIR}/steps/03_compare.sh")
ok "03_compare.sh       JOBID=${JOB_COMPARE}  (after ${JOB_FS})"

JOB_CURATE=$(sbatch --parsable \
    --dependency=afterok:${JOB_COMPARE} \
    "${MAIL_OVERRIDES[@]}" \
    "${SCRIPT_DIR}/steps/04_curate.sh")
ok "04_curate.sh        JOBID=${JOB_CURATE}  (after ${JOB_COMPARE})"

# 2026-06: the build-output step is now auto-chained (curation is near-fully
# automatic). The few genes still wanting a human look land in 04_output/curation/review_suggested.csv.
JOB_OUTPUT=$(sbatch --parsable \
    --dependency=afterok:${JOB_CURATE} \
    "${MAIL_OVERRIDES[@]}" \
    "${SCRIPT_DIR}/steps/04_output.sh")
ok "04_output.sh        JOBID=${JOB_OUTPUT}  (after ${JOB_CURATE})"

echo ""
echo -e "${YELLOW}[OPTIONAL, recommended]${NC} phynteny + synteny integration (one merged step):"
echo "    sbatch --dependency=afterok:${JOB_OUTPUT} steps/05_phynteny.sh"
echo "  (04_output/final_annotations_table.csv is complete after step 04; 05_phynteny"
echo "   adds phynteny categories + C1/Cro synteny naming into the integrated table & GBK)"
echo ""
step_banner "Submission complete — check  squeue -u \$USER"