#!/usr/bin/env bash
# =============================================================================
# 00b_install_pharokka_db.sh — One-time Pharokka database download
#
# Two modes:
#   1. INTERACTIVE on a login node (recommended):
#        bash steps/00b_install_pharokka_db.sh
#   2. SLURM:
#        sbatch steps/00b_install_pharokka_db.sh
#
# Why curl and not `install_databases.py -o ${PHAROKKA_DB}` ?
#   Behind some authenticated HTTP proxies, pharokka's Python `requests`
#   downloader mis-parses proxy URLs whose password contains `+`/`=` (which
#   must be URL-encoded); curl tolerates them. Once the tarball is on disk we
#   untar it and pharokka is happy. With no proxy, either route works.
# =============================================================================

#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --qos=fast
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -J pharokka_db_dl
#SBATCH -o logs/pharokka_db_dl.out
#SBATCH -e logs/pharokka_db_dl.err
#SBATCH --time=01:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email@example.com

set -euo pipefail

# Resolve project root (works interactively and under sbatch).
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
fi
source "${PROJECT_DIR}/config.sh"

PHAROKKA_TAR_URL="https://zenodo.org/record/17110353/files/pharokka_v1.8.0_databases.tar.gz"
TAR_PATH="${PHAROKKA_DB}/pharokka_v1.8.0_databases.tar.gz"

mkdir -p "${PHAROKKA_DB}"

log "Pharokka DB target dir : ${PHAROKKA_DB}"
log "Host                   : $(hostname)"
log "Proxy env (truncated)  : ${https_proxy:0:50}..."

# ---------------------------------------------------------------------------
# Download with curl (honours any HTTP(S) proxy env vars set by your shell init)
# ---------------------------------------------------------------------------
if [[ -f "${PHAROKKA_DB}/VERSION_1_8_0" ]] || [[ -f "${PHAROKKA_DB}/phrog_annot_v4.tsv" ]]; then
    log "Pharokka DB appears already installed -- skipping download."
else
    log "Downloading ${PHAROKKA_TAR_URL}..."
    curl --location \
         --output "${TAR_PATH}" \
         --retry 3 \
         --retry-delay 5 \
         "${PHAROKKA_TAR_URL}"

    log "Download complete.  Tarball size:"
    ls -lh "${TAR_PATH}"
fi

# ---------------------------------------------------------------------------
# Extract (idempotent: tar -k skips existing files)
# ---------------------------------------------------------------------------
if [[ -f "${TAR_PATH}" ]]; then
    log "Extracting ${TAR_PATH}..."
    cd "${PHAROKKA_DB}"
    # Pharokka's Zenodo tarball wraps the DB inside pharokka_v1.8.0_databases/.
    # Always strip the wrapping dir so files land directly in ${PHAROKKA_DB}.
    tar -xzf "${TAR_PATH}" --strip-components=1
    log "Extract complete.  Removing tarball."
    rm -f "${TAR_PATH}"
fi

# ---------------------------------------------------------------------------
# Verify
# ---------------------------------------------------------------------------
log "Files in ${PHAROKKA_DB}:"
ls -la "${PHAROKKA_DB}" | head -20

REQUIRED_FILES=( "VERSION_1_8_0" "phrog_annot_v4.tsv" )
for f in "${REQUIRED_FILES[@]}"; do
    if [[ ! -f "${PHAROKKA_DB}/${f}" ]]; then
        warn "Expected file missing: ${PHAROKKA_DB}/${f}"
    fi
done

log "Pharokka DB installation finished."
