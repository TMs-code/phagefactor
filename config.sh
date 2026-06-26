#!/usr/bin/env bash
# =============================================================================
# config.sh — phageFACTor configuration (bash layer for the SLURM step scripts)
# =============================================================================
# The Python scripts read config/config.yaml (via config/config.py). This bash
# file mirrors the few values the steps/*.sh SLURM wrappers need. Paths are
# resolved relative to the repo root, overridable with environment variables.
# Source it from a step:  source "$(dirname "$0")/../config.sh"
# =============================================================================

# =============================================================================
# SECTION 1 — PROJECT + SLURM  ★ edit / export to match your machine
# =============================================================================

# Repo root (the directory that contains this file). Override with PHAGEFACTOR_ROOT.
PHAGEFACTOR_ROOT="${PHAGEFACTOR_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PIPELINE_DIR="${PHAGEFACTOR_ROOT}"

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

INPUT_DIR="${PIPELINE_DIR}/input"
PROPHAGE_LIST="${INPUT_DIR}/prophage_list.txt"      # one prophage name per line
FASTA_DIR="${INPUT_DIR}/fasta"                      # genome mode: one nucleotide FASTA per prophage

# Pharokka output (consumed by phold): <NAME>/<NAME>.gbk
PHAROKKA_OUT_DIR="${PIPELINE_DIR}/00c_pharokka"
GBK_DIR="${PHAROKKA_OUT_DIR}"

PHOLD_OUT_DIR="${PIPELINE_DIR}/01_phold"
FOLDSEEK_DIR="${PIPELINE_DIR}/02_foldseek"
COMPARISON_DIR="${PIPELINE_DIR}/03_comparison"
# 04_output/ = deliverables; 04_output/curation/ = curated + review_suggested;
# 05_phynteny/ = phynteny + integration; 05_phynteny/run/ = phynteny log + input FASTA.
OUTPUT_DIR="${PIPELINE_DIR}/04_output"
CURATION_DIR="${OUTPUT_DIR}/curation"
PHYNTENY_DIR="${PIPELINE_DIR}/05_phynteny"
LOGS_DIR="${PIPELINE_DIR}/logs"
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"

# Combined phold outputs (produced by 01c_merge_phold.py after the array job)
PHOLD_COMBINED_TSV="${PHOLD_OUT_DIR}/combined/phold_all.tsv"
PHOLD_3DI_FASTA="${PHOLD_OUT_DIR}/combined/phold_3di.fasta"
PHOLD_AA_FASTA="${PHOLD_OUT_DIR}/combined/phold_aa.fasta"

# =============================================================================
# SECTION 3 — TOOL DATABASES  (edit once after install, or export the vars)
# =============================================================================

# Database root. Override with PHAGEFACTOR_DB_ROOT. (For WebAPI search mode you
# do not need local FoldSeek DBs — see docs/databases.md.)
DB_ROOT="${PHAGEFACTOR_DB_ROOT:-${PIPELINE_DIR}/databases}"

PHAROKKA_DB="${PHAROKKA_DB:-${DB_ROOT}/pharokka_db}"

# FoldSeek databases (afdb50 + afdb-swissprot + pdb100; see steps/00b_download_foldseek_dbs.sh)
FOLDSEEK_DB_ROOT="${FOLDSEEK_DB_ROOT:-${DB_ROOT}/foldseek_dbs}"
FOLDSEEK_DB_PDB100="${FOLDSEEK_DB_ROOT}/pdb100_db/pdb100"
FOLDSEEK_DB_SWISSPROT="${FOLDSEEK_DB_ROOT}/afdb_swissprot_db/afdb_swissprot"
FOLDSEEK_DB_AFDB50="${FOLDSEEK_DB_ROOT}/afdb50_db/afdb50"

# =============================================================================
# SECTION 4 — COMPUTE RESOURCES
# =============================================================================

THREADS="${THREADS:-16}"
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

activate_env() {
    # SLURM jobs start with a minimal shell (~/.bashrc is NOT sourced), so load
    # the micromamba shell hook inline before activating.
    if [[ "${MAMBA_DEFAULT_ENV:-}" == "${MAMBA_ENV}" ]]; then return; fi
    local mamba_bin="${MAMBA_ROOT_PREFIX}/bin/micromamba"
    [[ -x "${mamba_bin}" ]] \
        || die "micromamba not found at ${mamba_bin}. Run steps/00_install_env.sh first."
    export MAMBA_ROOT_PREFIX
    eval "$("${mamba_bin}" shell hook --shell bash)"
    micromamba activate "${MAMBA_ENV}" \
        || die "Cannot activate env '${MAMBA_ENV}'. Run steps/00_install_env.sh first."
}

require_tool() {
    command -v "$1" &>/dev/null \
        || die "Tool not found: $1 — is the '${MAMBA_ENV}' env active?"
}
