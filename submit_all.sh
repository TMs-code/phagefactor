#!/usr/bin/env bash
# =============================================================================
# submit_all.sh — Orchestrated SLURM submission for phageFACTor
#
# AUTO-DETECTS input type from input/fasta/ and selects the correct pipeline:
#
#   GENOME MODE  (sequences > 5kb, nucleotide):
#     00c (pharokka, FASTA → GBK) → 01 (phold genome) → 01c merge (auto)
#     → 02 (FoldSeek) → 03 (compare) → 04 (curate + build output table)
#     → [optional, recommended] 05_phynteny
#
#   PROTEIN MODE (sequences are AA, < 5kb):
#     00p split (fixed PROTEIN_BATCH_SIZE-protein batches, inline)
#     → 01p phold-proteins ARRAY (one batch per task, CPU-default)
#     → 01p merge (afterok) → 02 → 03 → 04 (curate + build output table)
#     → [optional, recommended] 05_phynteny
#     All input proteins treated as hypothetical.
#
# Usage:
#   cd phagefactor/
#   bash submit_all.sh                          # auto-detect input type, NEW run
#   INPUT_MODE=genome bash submit_all.sh        # force genome mode
#   INPUT_MODE=protein bash submit_all.sh       # force protein mode
#   SKIP_PHAROKKA=1 bash submit_all.sh          # genome mode, skip pharokka only
#                                                # (still a NEW run -- no RESUME needed,
#                                                #  phold onward runs fresh)
#   REDETECT=1 bash submit_all.sh               # ignore cached .input_type
#
# Resuming an EXISTING run (SKIP_PHOLD / SKIP_FOLDSEEK / SKIP_COMPARE): each
# of these reuses output already sitting in some earlier run dir, so each
# REQUIRES you to say which run — with RESUME=1 (the most recent one) or
# PHAGEFACTOR_RUN_DIR=/path/to/run (a specific one). Sourcing an env file is
# NOT enough by itself; without RESUME=1 or PHAGEFACTOR_RUN_DIR this script
# will refuse to start rather than silently create a new, empty run:
#   RESUME=1 SKIP_PHOLD=1     bash submit_all.sh   # continue from step 02 (FoldSeek)
#   RESUME=1 SKIP_FOLDSEEK=1  bash submit_all.sh   # continue from step 03 (compare)
#   RESUME=1 SKIP_COMPARE=1   bash submit_all.sh   # continue from step 04 (curate +
#                                                    build output, one job) -- the one
#                                                    to use when you've only edited
#                                                    curation rules and want fresh
#                                                    output against an existing run's
#                                                    FoldSeek + comparison output
# (SKIP_FOLDSEEK / SKIP_COMPARE each imply the earlier skips too.)
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
export PHAGEFACTOR_HOME="${SCRIPT_DIR}"   # passed to every sbatch job so it finds config.sh regardless of launch dir
source "${SCRIPT_DIR}/config.sh"

SKIP_PHAROKKA="${SKIP_PHAROKKA:-0}"
SKIP_PHOLD="${SKIP_PHOLD:-0}"
SKIP_FOLDSEEK="${SKIP_FOLDSEEK:-0}"
SKIP_COMPARE="${SKIP_COMPARE:-0}"
REDETECT="${REDETECT:-0}"
# Each later skip implies the earlier ones -- SKIP_COMPARE means "everything
# through 03_compare is already done in this run dir", so pharokka/phold/
# FoldSeek must all be treated as skipped too.
[[ "${SKIP_COMPARE}" == "1" ]]  && SKIP_FOLDSEEK=1
[[ "${SKIP_FOLDSEEK}" == "1" ]] && SKIP_PHOLD=1
# SEARCH_MODE=local (default, step 02 local FoldSeek/afdb50) | webapi (no local DB,
# search via search.foldseek.com on the LOGIN node — see step-02 section below).
SEARCH_MODE="${SEARCH_MODE:-local}"

# --- Guard: any SKIP_* beyond SKIP_PHAROKKA needs an EXPLICIT existing run ---
# SKIP_PHOLD / SKIP_FOLDSEEK / SKIP_COMPARE all mean "earlier steps already ran
# in some run dir, reuse their output" — that's meaningless without telling
# this invocation WHICH run dir. Without this guard, forgetting to set
# PHAGEFACTOR_RUN_DIR or RESUME=1 silently falls through to "allocate a fresh,
# EMPTY run dir" below, and the skipped steps then submit against nothing.
if [[ "${SKIP_PHOLD}" == "1" && -z "${PHAGEFACTOR_RUN_DIR:-}" && "${RESUME:-0}" != "1" ]]; then
    die "SKIP_PHOLD / SKIP_FOLDSEEK / SKIP_COMPARE resume an EXISTING run, but
no run was specified. Set one of:
  RESUME=1                          # reuse the most recent run
  PHAGEFACTOR_RUN_DIR=/path/to/run  # reuse a specific run
alongside the SKIP_* flag, e.g.:
  RESUME=1 SKIP_COMPARE=1 bash submit_all.sh
Just sourcing an env file is not enough on its own — it does not set either of
these, so without one of them this would silently start a brand-new, empty run."
fi

# --- Resolve THIS run's output dir ONCE, then share it with every step -------
# Default: allocate a FRESH versioned dir (run, run_002, ...) so a new submission
# never overwrites a previous run. Override with:
#   PHAGEFACTOR_RUN_DIR=/path   → reuse an EXACT run dir (resume a specific run)
#   RESUME=1                    → reuse the MOST RECENT run (<run_base>/.current_run)
# Combine either with SKIP_PHAROKKA=1 / SKIP_PHOLD=1 / SKIP_FOLDSEEK=1 /
# SKIP_COMPARE=1 to skip already-done steps.
if [[ -z "${PHAGEFACTOR_RUN_DIR:-}" ]]; then
    if [[ "${RESUME:-0}" == "1" && -f "${RUN_BASE}/.current_run" ]]; then
        export PHAGEFACTOR_RUN_DIR="$(cat "${RUN_BASE}/.current_run")"
        echo "[resume] RESUME=1 — reusing most recent run: ${PHAGEFACTOR_RUN_DIR}"
    else
        export PHAGEFACTOR_RUN_DIR="$(next_run_dir "${RUN_BASE}" "${RUN_NAME}")"
    fi
    source "${SCRIPT_DIR}/config.sh"     # re-derive OUTPUT/LOGS dirs under the run dir
fi
mkdir -p "${RUN_DIR}" "${LOGS_DIR}" "${RUN_BASE}"
printf '%s\n' "${RUN_DIR}" > "${RUN_BASE}/.current_run"   # so later manual steps follow it

# --- Guard: the resolved run dir must actually HAVE the output a SKIP_* flag
# expects to reuse (catches a wrong/stale PHAGEFACTOR_RUN_DIR or RESUME=1
# pointing at a run that never got that far). Checks the earliest required
# file for the deepest skip requested.
if [[ "${SKIP_COMPARE}" == "1" && ! -f "${COMPARISON_DIR}/comparison_per_gene.csv" ]]; then
    die "SKIP_COMPARE=1 but ${COMPARISON_DIR}/comparison_per_gene.csv does not
exist — this run dir (${RUN_DIR}) never finished step 03 (compare). Point
PHAGEFACTOR_RUN_DIR / RESUME=1 at a run that has, or drop SKIP_COMPARE=1."
elif [[ "${SKIP_FOLDSEEK}" == "1" && ! -f "${FOLDSEEK_DIR}/3di_tokens/best_hit.csv" ]]; then
    die "SKIP_FOLDSEEK=1 but ${FOLDSEEK_DIR}/3di_tokens/best_hit.csv does not
exist — this run dir (${RUN_DIR}) never finished step 02 (FoldSeek). Point
PHAGEFACTOR_RUN_DIR / RESUME=1 at a run that has, or drop SKIP_FOLDSEEK=1."
elif [[ "${SKIP_PHOLD}" == "1" && ! -f "${PHOLD_AA_FASTA}" ]]; then
    die "SKIP_PHOLD=1 but ${PHOLD_AA_FASTA} does not exist — this run dir
(${RUN_DIR}) never finished phold + the merge. Point PHAGEFACTOR_RUN_DIR /
RESUME=1 at a run that has, or drop SKIP_PHOLD=1."
fi

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[SUBMITTED]${NC}  $*"; }
info() { echo -e "              $*"; }
mode() { echo -e "${CYAN}[MODE]${NC}       $*"; }

step_banner "phageFACTor — submitting SLURM jobs"
info "input : ${INPUT_DIR}"
info "run   : ${RUN_DIR}   (all outputs + logs land here)"

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
    # shared modules may live in scripts/ or scripts/lib/
    DETECT_PY="${SCRIPTS_DIR}/lib/detect_input_type.py"
    [[ -f "${DETECT_PY}" ]] || DETECT_PY="${SCRIPTS_DIR}/detect_input_type.py"
    DETECT_ARGS=()
    [[ "${REDETECT}" == "1" ]] && DETECT_ARGS+=(--redetect)

    # Call python3 explicitly (some submit nodes have no python->python3
    # symlink) and surface any failure with a clear message rather than
    # swallowing stderr.
    INPUT_MODE=$(python3 "${DETECT_PY}" \
        --fasta-dir "${FASTA_DIR}" \
        "${DETECT_ARGS[@]}") \
        || die "Auto-detection failed (see error above). Debug with:
  python3 ${SCRIPTS_DIR}/detect_input_type.py --fasta-dir ${FASTA_DIR} --verbose
or override explicitly: INPUT_MODE=protein bash submit_all.sh"

    mode "Auto-detected input type: ${INPUT_MODE}"
fi

# --- PANGENOME sub-mode: translate a bulk nucleotide CDS set to protein --------
# INPUT_MODE=pangenome → translate <nt fasta> into protein-mode input (clean ids
# + id map), then run exactly as protein mode. This is the one-command DNA path;
# detect_input_type.py otherwise deliberately REFUSES to guess on bulk nt CDS.
# Source nt FASTA: $PANGENOME_NT_FASTA, else the first file in $INPUT_DIR/pangenome_nt/.
if [[ "${INPUT_MODE}" == "pangenome" ]]; then
    # Source nt FASTA, in order: $PANGENOME_NT_FASTA, else input/pangenome_nt/*,
    # else a nucleotide FASTA sitting in input/fasta/ (the "just drop it in like
    # any job" workflow). In that last case the original nt file is MOVED to
    # input/pangenome_nt/ after translation so the protein split (00p) does not
    # re-process the DNA.
    NT_SRC="${PANGENOME_NT_FASTA:-}"
    [[ -z "${NT_SRC}" ]] && NT_SRC=$(ls "${INPUT_DIR}"/pangenome_nt/*.fa "${INPUT_DIR}"/pangenome_nt/*.fna "${INPUT_DIR}"/pangenome_nt/*.fasta 2>/dev/null | head -1 || true)
    [[ -z "${NT_SRC}" ]] && NT_SRC=$(ls "${FASTA_DIR}"/*.fa "${FASTA_DIR}"/*.fna "${FASTA_DIR}"/*.fasta 2>/dev/null | head -1 || true)
    [[ -n "${NT_SRC}" && -f "${NT_SRC}" ]] || die "INPUT_MODE=pangenome but no nucleotide FASTA found. Put it in ${FASTA_DIR}/ (or ${INPUT_DIR}/pangenome_nt/, or set PANGENOME_NT_FASTA=)."
    mode "PANGENOME: translating ${NT_SRC} -> protein-mode input (table ${PANGENOME_TABLE:-11})"
    mkdir -p "${FASTA_DIR}" "${INPUT_DIR}/pangenome_nt"
    _NT_BASE=$(basename "${NT_SRC}"); _NT_STEM="${_NT_BASE%.*}"
    python3 "${SCRIPTS_DIR}/00n_translate_pangenome.py" \
        --in "${NT_SRC}" \
        --out-faa "${FASTA_DIR}/${_NT_STEM}_proteins.faa" \
        --id-map "${INPUT_DIR}/pangenome_id_map.csv" \
        --table "${PANGENOME_TABLE:-11}" --prefix "${PANGENOME_PREFIX:-pang_}" \
        || die "pangenome translation failed (see error above)"
    if [[ "$(cd "$(dirname "${NT_SRC}")" && pwd)" == "$(cd "${FASTA_DIR}" && pwd)" ]]; then
        mv "${NT_SRC}" "${INPUT_DIR}/pangenome_nt/" \
            && info "Moved nucleotide source out of fasta/ -> ${INPUT_DIR}/pangenome_nt/${_NT_BASE}"
    fi
    export PANGENOME_ID_MAP="${INPUT_DIR}/pangenome_id_map.csv"
    INPUT_MODE="protein"
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
# --chdir sends each job's SLURM .out/.err (the '#SBATCH -o logs/...' paths, which
# are relative) into ${RUN_DIR}/logs alongside everything else for this run.
MAIL_OVERRIDES=( --mail-user="${EMAIL}" --mail-type="${MAIL_TYPE:-BEGIN,END,FAIL}" --chdir="${RUN_DIR}" )

# ═══════════════════════════════════════════════════════════════════════════════
# PROTEIN MODE: split → phold-proteins array (fixed-size batches) → merge →
#               then steps 02–04 chained
#
# Submits a fixed-size batch array (mirroring the genome-mode split →
# 01_phold_array → 01c_merge_phold chain) rather than one single large GPU
# job. Design choices: PROTEIN_BATCH_SIZE-protein batches (default 50,
# config.sh) keep batch composition reproducible across datasets regardless
# of total protein count, with the last batch simply left partial; CPU is the
# default for the initial submission (PHOLD_PROTEINS_USE_GPU=0, decoupled
# from genome-mode PHOLD_USE_GPU) to avoid the GPU-autotune failure mode
# described in steps/01_phold_array.sh.
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

        BATCH_LIST="${PROTEIN_BATCH_LIST}"
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

    # Steps 02–04 (same as genome mode, chained from the merge job, but
    # skip-aware: SKIP_COMPARE/SKIP_FOLDSEEK let you resume mid-chain when
    # only curation rules changed -- see the usage note at the top of this file).
    # steps/04_curate.sh runs BOTH curation and output-building in one job
    # (they always ran back-to-back, so no separate 04_output job/dependency).
    if [[ "${SKIP_COMPARE}" == "1" ]]; then
        info "SKIP_COMPARE=1 — skipping FoldSeek + compare, starting at 04_curate"
        JOB_CURATE=$(sbatch --parsable "${MAIL_OVERRIDES[@]}" "${SCRIPT_DIR}/steps/04_curate.sh")
        ok "04_curate.sh          JOBID=${JOB_CURATE}"
    elif [[ "${SKIP_FOLDSEEK}" == "1" ]]; then
        info "SKIP_FOLDSEEK=1 — skipping FoldSeek, starting at 03_compare"
        JOB_COMPARE=$(sbatch --parsable "${MAIL_OVERRIDES[@]}" "${SCRIPT_DIR}/steps/03_compare.sh")
        ok "03_compare.sh         JOBID=${JOB_COMPARE}"
        JOB_CURATE=$(sbatch --parsable \
            --dependency=afterok:${JOB_COMPARE} \
            "${MAIL_OVERRIDES[@]}" \
            "${SCRIPT_DIR}/steps/04_curate.sh")
        ok "04_curate.sh          JOBID=${JOB_CURATE}  (after ${JOB_COMPARE})"
    else
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
    fi

    # --- Pan-genome homology-reference: only when an id map exists (00n was run) --
    PANGENOME_ID_MAP="${PANGENOME_ID_MAP:-${INPUT_DIR}/pangenome_id_map.csv}"
    if [[ -f "${PANGENOME_ID_MAP}" ]]; then
        JOB_PANGREF=$(sbatch --parsable \
            --dependency=afterok:${JOB_CURATE} \
            "${MAIL_OVERRIDES[@]}" \
            "${SCRIPT_DIR}/steps/08_pangenome_reference.sh")
        ok "08_pangenome_reference.sh JOBID=${JOB_PANGREF}  (after ${JOB_CURATE}; pan-genome mode)"
    fi

    echo ""
    echo -e "${YELLOW}[OPTIONAL, recommended]${NC} phynteny + synteny integration (one merged step):"
    echo "    sbatch --chdir=${RUN_DIR} --dependency=afterok:${JOB_CURATE} steps/05_phynteny.sh"
    echo "  (this run's outputs are in ${RUN_DIR}; 05_phynteny auto-follows it via"
    echo "   ${RUN_BASE}/.current_run, and --chdir keeps its logs there too.)"
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
# Reached by falling through from the phold block (DEP_FS = the merge job),
# directly via SKIP_PHOLD=1 (DEP_FS empty -> step 02 starts immediately),
# or further in via SKIP_FOLDSEEK=1 (start at 03_compare) / SKIP_COMPARE=1
# (start at 04_curate) -- both require an EXISTING run's earlier-step outputs,
# see RESUME=1 / PHAGEFACTOR_RUN_DIR in the usage note at the top of this file.
# ---------------------------------------------------------------------------
echo ""

# --- WEBAPI search mode: login-node step, not part of the sbatch chain --------
# Only relevant when step 02 is actually about to run locally.
if [[ "${SEARCH_MODE}" == "webapi" && "${SKIP_FOLDSEEK}" != "1" ]]; then
    echo -e "${YELLOW}[WEBAPI MODE]${NC} Step 02 runs against search.foldseek.com (no local DB)."
    echo "  It must run on the LOGIN node (compute nodes have no internet) AFTER the"
    echo "  phold merge${JOB_MERGE:+ (job ${JOB_MERGE})} completes:"
    echo "    bash steps/02w_foldseek_webapi.sh        # writes best_hit.csv/top3.csv"
    echo "  then chain compare -> curate (this run = ${RUN_DIR}; steps auto-follow it"
    echo "  via ${RUN_BASE}/.current_run, and --chdir keeps their logs there):"
    echo "    J03=\$(sbatch --parsable --chdir=${RUN_DIR} steps/03_compare.sh)"
    echo "    J04=\$(sbatch --parsable --chdir=${RUN_DIR} --dependency=afterok:\$J03 steps/04_curate.sh)"
    echo "    sbatch --chdir=${RUN_DIR} --dependency=afterok:\$J04 steps/05_phynteny.sh   # optional"
    echo ""
    step_banner "Submitted through phold merge. WebAPI search 02w is a manual login-node step."
    exit 0
fi

info "Submitting steps 02–04 with dependency chain..."
echo ""

# steps/04_curate.sh runs BOTH curation and output-building in one job (they
# always ran back-to-back, so no separate 04_output job/dependency).
if [[ "${SKIP_COMPARE}" == "1" ]]; then
    info "SKIP_COMPARE=1 — skipping FoldSeek + compare, starting at 04_curate"
    JOB_CURATE=$(sbatch --parsable "${MAIL_OVERRIDES[@]}" "${SCRIPT_DIR}/steps/04_curate.sh")
    ok "04_curate.sh        JOBID=${JOB_CURATE}"
elif [[ "${SKIP_FOLDSEEK}" == "1" ]]; then
    info "SKIP_FOLDSEEK=1 — skipping FoldSeek, starting at 03_compare"
    JOB_COMPARE=$(sbatch --parsable "${MAIL_OVERRIDES[@]}" "${SCRIPT_DIR}/steps/03_compare.sh")
    ok "03_compare.sh       JOBID=${JOB_COMPARE}"
    JOB_CURATE=$(sbatch --parsable \
        --dependency=afterok:${JOB_COMPARE} \
        "${MAIL_OVERRIDES[@]}" \
        "${SCRIPT_DIR}/steps/04_curate.sh")
    ok "04_curate.sh        JOBID=${JOB_CURATE}  (after ${JOB_COMPARE})"
else
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
fi

echo ""
echo -e "${YELLOW}[OPTIONAL, recommended]${NC} phynteny + synteny integration (one merged step):"
echo "    sbatch --chdir=${RUN_DIR} --dependency=afterok:${JOB_CURATE} steps/05_phynteny.sh"
echo "  (this run = ${RUN_DIR}; 05_phynteny auto-follows it via ${RUN_BASE}/.current_run,"
echo "   and --chdir keeps its logs there. ${RUN_DIR}/04_output/final_annotations_table.csv"
echo "   is already complete after step 04; 05 adds phynteny + C1/Cro synteny naming.)"
echo ""
step_banner "Submission complete — check  squeue -u \$USER"