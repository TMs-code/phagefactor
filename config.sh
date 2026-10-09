#!/usr/bin/env bash
# =============================================================================
# config.sh — phageFACTor configuration (bash layer for the SLURM step scripts)
# =============================================================================
# The Python scripts read config/config.yaml (via scripts/lib/config.py). This bash
# file mirrors the few values the steps/*.sh SLURM wrappers need. Paths are
# resolved relative to the repo root, overridable with environment variables.
# Source it from a step:  source "$(dirname "$0")/../config.sh"
# =============================================================================

# =============================================================================
# SECTION 1 — PROJECT + SLURM  ★ edit / export to match your machine
# =============================================================================

# Repo root = where the CODE lives (this file's dir). Override with PHAGEFACTOR_ROOT.
PHAGEFACTOR_ROOT="${PHAGEFACTOR_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PIPELINE_DIR="${PHAGEFACTOR_ROOT}"

# INPUT — folder holding fasta/ + prophage_list.txt (data in). May be anywhere;
# override with PHAGEFACTOR_INPUT. Default: <repo>/input.
INPUT_DIR="${PHAGEFACTOR_INPUT:-${PIPELINE_DIR}/input}"

# RUN — ALL outputs + logs go under a versioned run dir. Override the PARENT with
# PHAGEFACTOR_RUN_BASE (default <repo>/runs). The first run is <run_base>/<run_name>;
# if it exists the pipeline makes <run_name>_002, _003, ... (never overwrites).
# submit_all.sh resolves the exact dir once and exports it as PHAGEFACTOR_RUN_DIR
# so every SLURM step writes to the SAME place.
RUN_BASE="${PHAGEFACTOR_RUN_BASE:-${PIPELINE_DIR}/runs}"
RUN_NAME="${PHAGEFACTOR_RUN_NAME:-run}"
if [[ -n "${PHAGEFACTOR_RUN_DIR:-}" ]]; then
    RUN_DIR="${PHAGEFACTOR_RUN_DIR}"                 # set by submit_all (authoritative)
elif [[ -f "${RUN_BASE}/.current_run" ]]; then
    RUN_DIR="$(cat "${RUN_BASE}/.current_run")"      # follow the most recent run
else
    RUN_DIR="${RUN_BASE}/${RUN_NAME}"                # first run / fallback
fi

# SLURM notification email (leave empty to disable). Override with SLURM_EMAIL.
EMAIL="${SLURM_EMAIL:-}"
MAIL_TYPE="${SLURM_MAIL_TYPE:-FAIL}"        # subset of BEGIN,END,FAIL,REQUEUE,ALL

# SLURM partition / QOS / account — set to your scheduler's names.
PARTITION="${SLURM_PARTITION:-common}"
QOS="${SLURM_QOS:-fast}"
ACCOUNT="${SLURM_ACCOUNT:-}"                 # empty if no project account needed

# =============================================================================
# SECTION 2 — DERIVED PATHS  (auto-computed; no edit needed)
# =============================================================================

# --- Input (from INPUT_DIR) ---
PROPHAGE_LIST="${INPUT_DIR}/prophage_list.txt"      # one prophage name per line
FASTA_DIR="${INPUT_DIR}/fasta"                      # genome mode: one nucleotide FASTA per prophage

# --- Outputs (ALL under the resolved RUN_DIR) ---
# Pharokka output (consumed by phold): <NAME>/<NAME>.gbk
PHAROKKA_OUT_DIR="${RUN_DIR}/00c_pharokka"
GBK_DIR="${PHAROKKA_OUT_DIR}"

PHOLD_OUT_DIR="${RUN_DIR}/01_phold"
FOLDSEEK_DIR="${RUN_DIR}/02_foldseek"
COMPARISON_DIR="${RUN_DIR}/03_comparison"
# 04_output/ = deliverables; 04_output/curation/ = curated + review_suggested;
# 05_phynteny/ = phynteny + integration; 05_phynteny/run/ = phynteny log + input FASTA.
OUTPUT_DIR="${RUN_DIR}/04_output"
CURATION_DIR="${OUTPUT_DIR}/curation"
PHYNTENY_DIR="${RUN_DIR}/05_phynteny"
LOGS_DIR="${RUN_DIR}/logs"
SPLIT_DIR="${RUN_DIR}/split"                        # generated split/batch/gene_metadata
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"               # code stays with the repo

# Combined phold outputs (produced by 01c_merge_phold.py after the array job)
PHOLD_COMBINED_TSV="${PHOLD_OUT_DIR}/combined/phold_all.tsv"
PHOLD_3DI_FASTA="${PHOLD_OUT_DIR}/combined/phold_3di.fasta"
PHOLD_AA_FASTA="${PHOLD_OUT_DIR}/combined/phold_aa.fasta"

# Protein-mode generated artifacts (kept OUT of INPUT_DIR — under the run dir)
PROTEIN_COMBINED_FAA="${SPLIT_DIR}/all_proteins_combined.faa"
PROTEIN_BATCH_DIR="${SPLIT_DIR}/protein_batches"
PROTEIN_BATCH_LIST="${SPLIT_DIR}/protein_batch_list.txt"

# =============================================================================
# SECTION 3 — TOOL DATABASES  (edit once after install, or export the vars)
# =============================================================================

# Database root. Override with PHAGEFACTOR_DB_ROOT. (For WebAPI search mode you
# do not need local FoldSeek DBs — see docs/databases.md.)
DB_ROOT="${PHAGEFACTOR_DB_ROOT:-${PIPELINE_DIR}/databases}"

PHAROKKA_DB="${PHAROKKA_DB:-${DB_ROOT}/pharokka_db}"

# FoldSeek databases (afdb50 + afdb-swissprot + pdb100; see steps/00b_setup_databases.sh)
FOLDSEEK_DB_ROOT="${FOLDSEEK_DB_ROOT:-${DB_ROOT}/foldseek_dbs}"
FOLDSEEK_DB_PDB100="${FOLDSEEK_DB_ROOT}/pdb100_db/pdb100"
FOLDSEEK_DB_SWISSPROT="${FOLDSEEK_DB_ROOT}/afdb_swissprot_db/afdb_swissprot"
FOLDSEEK_DB_AFDB50="${FOLDSEEK_DB_ROOT}/afdb50_db/afdb50"

# =============================================================================
# SECTION 4 — COMPUTE RESOURCES
# =============================================================================

THREADS="${THREADS:-16}"

# SLURM exports SLURM_CPUS_PER_TASK inside a batch job. Off-SLURM (run_local.sh,
# or running a single step by hand) it is UNSET, and every step reads it bare as
# ${SLURM_CPUS_PER_TASK} -- which aborts immediately under `set -u`. Default it
# here, once, rather than editing every step or faking the variable in the
# runner: one core is left free so the machine stays usable while phold runs.
if [[ -z "${SLURM_CPUS_PER_TASK:-}" ]]; then
    _ncpu=$(nproc 2>/dev/null || echo 4)
    SLURM_CPUS_PER_TASK=$(( _ncpu > 1 ? _ncpu - 1 : 1 ))
    export SLURM_CPUS_PER_TASK
fi
MEM_PHAROKKA="16G"
MEM_PHOLD="32G"
MEM_FOLDSEEK="64G"
MEM_PYTHON="16G"

TIME_PHAROKKA="00:30:00"    # per prophage (Phanotate gene calling)
TIME_PHOLD="02:00:00"       # per prophage / per batch
TIME_FOLDSEEK="08:00:00"
TIME_PHYNTENY="02:00:00"
TIME_PYTHON="01:00:00"      # steps 03–05

# phold is the only GPU-accelerable step. Default CPU (0): robust everywhere and
# fast enough (~15–30 min per 50-protein batch on 8 CPUs). Set to 1 only with a
# CUDA-enabled env (verify: python -c "import torch; print(torch.cuda.is_available())").
PHOLD_USE_GPU="${PHOLD_USE_GPU:-0}"
PHOLD_PROTEINS_USE_GPU="${PHOLD_PROTEINS_USE_GPU:-0}"
PHOLD_GPU_PARTITION="${SLURM_GPU_PARTITION:-gpu}"
PHOLD_GPU_QOS="${SLURM_GPU_QOS:-gpu}"
PHOLD_GPU_GRES="gpu:1"
PHOLD_AUTOTUNE=1

# Protein mode: fixed batch size for the phold array (last batch may be partial).
PROTEIN_BATCH_SIZE="${PROTEIN_BATCH_SIZE:-50}"

# Scratch dir for large temp files (cluster scratch if available, else /tmp).
SCRATCH_TMPDIR="${SCRATCH_TMPDIR:-${TMPDIR:-/tmp}}"

# =============================================================================
# SECTION 5 — CONDA / MICROMAMBA ENVIRONMENT
# =============================================================================

MAMBA_ENV="${MAMBA_ENV:-phagefactor}"            # env containing pharokka + phold + foldseek
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${HOME}/.mamba}"

# =============================================================================
# SECTION 6 — INTERNAL HELPERS  (do not edit)
# =============================================================================

log()         { echo "[$(date '+%H:%M:%S')] $*"; }
die()         { echo "[ERROR] $*" >&2; exit 1; }
warn()        { echo "[WARN]  $*" >&2; }

step_banner() {
    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "  $*"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "════════════════════════════════════════════════════════"
}

check_file() { [[ -f "$1" ]] || die "Required file not found: $1"; }
check_dir()  { [[ -d "$1" ]] || die "Required directory not found: $1"; }

# Return the next FREE versioned run dir under base $1 with base name $2:
#   <base>/<name>, then <base>/<name>_002, _003, ...  (never overwrites an existing one).
next_run_dir() {
    local base="$1" name="$2" i cand
    cand="${base}/${name}"
    [[ ! -e "${cand}" ]] && { printf '%s\n' "${cand}"; return; }
    for ((i = 2; i < 1000; i++)); do
        cand=$(printf '%s/%s_%03d' "${base}" "${name}" "${i}")
        [[ ! -e "${cand}" ]] && { printf '%s\n' "${cand}"; return; }
    done
    die "Too many run dirs under ${base} (>=999)"
}

activate_env() {
    # SLURM jobs start with a minimal shell (~/.bashrc is NOT sourced), so load
    # the micromamba shell hook inline before activating.
    if [[ "${MAMBA_DEFAULT_ENV:-}" == "${MAMBA_ENV}" ]]; then return; fi

    # ALREADY-ACTIVE ENV WINS. On a laptop the env is usually activated by hand
    # (conda activate ...), and there may be no micromamba at all -- this used to
    # die outright. If the tools we need are already on PATH, there is nothing to
    # do; PHAGEFACTOR_SKIP_ENV=1 forces this path.
    if [[ "${PHAGEFACTOR_SKIP_ENV:-0}" == "1" ]] \
       || { command -v phold &>/dev/null && command -v python3 &>/dev/null; }; then
        return
    fi

    local mamba_bin="${MAMBA_ROOT_PREFIX}/bin/micromamba"
    if [[ -x "${mamba_bin}" ]]; then
        export MAMBA_ROOT_PREFIX
        eval "$("${mamba_bin}" shell hook --shell bash)"
        micromamba activate "${MAMBA_ENV}" \
            || die "Cannot activate env '${MAMBA_ENV}'. Run steps/00_install_env.sh first."
        return
    fi
    # conda fallback, so a conda-only machine works without installing micromamba
    if command -v conda &>/dev/null; then
        # shellcheck disable=SC1091
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "${MAMBA_ENV}" \
            || die "Cannot activate conda env '${MAMBA_ENV}'. Create it (environment.yml),
       activate it yourself and re-run, or set MAMBA_ENV=<your env name>."
        return
    fi
    die "No micromamba at ${mamba_bin} and no conda on PATH.
       Activate an env with pharokka/phold/foldseek yourself and re-run, or run
       steps/00_install_env.sh."
}

require_tool() {
    command -v "$1" &>/dev/null \
        || die "Tool not found: $1 — is the '${MAMBA_ENV}' env active?"
}
