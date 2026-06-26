#!/usr/bin/env bash
# =============================================================================
# 00b_download_foldseek_dbs.sh — One-time FoldSeek + Baktfold database setup
# =============================================================================
# Run INTERACTIVELY on a login node (downloads, not compute).
#
# What this does:
#   1. pdb100          — PDB experimental structures (~4 GB, via foldseek)
#   2. afdb-swissprot  — SwissProt predicted structures (~1.5 GB, via foldseek)
#   3. Baktfold DB     — bacterial+phage AFDB clusters, installed via baktfold tool
#
# NOTE: afdb-swissprot ≠ Baktfold DB.
#   afdb-swissprot = ~570k manually-curated SwissProt proteins → best functional names
#   Baktfold DB    = AFDB+PDB+SwissProt+CATH structurally clustered → bacterial depth
#   Keep both; they are complementary.
#
# Baktfold (https://github.com/gbouras13/baktfold) is a standalone annotation tool
# that manages its own databases.  It is NOT a raw FoldSeek tarball — its database
# is installed via  `baktfold install -d <dir>`.
# =============================================================================
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

activate_env

# Verify tools
require_tool foldseek
require_tool baktfold

TMPDIR="${SCRATCH_TMPDIR}"
export TMPDIR
mkdir -p "${TMPDIR}"

mkdir -p "${FOLDSEEK_DB_ROOT}"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[OK]${NC}  $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }

# ---------------------------------------------------------------------------
# 1. pdb100
# ---------------------------------------------------------------------------
step_banner "Downloading pdb100"

PDB_DIR="${FOLDSEEK_DB_ROOT}/pdb100_db"
if [[ -f "${PDB_DIR}/pdb100.dbtype" ]]; then
    ok "pdb100 already exists at ${PDB_DIR}"
else
    mkdir -p "${PDB_DIR}"
    foldseek databases PDB "${PDB_DIR}/pdb100" "${TMPDIR}/foldseek_tmp"
    ok "pdb100 downloaded → ${PDB_DIR}"
fi

# ---------------------------------------------------------------------------
# 2. afdb-swissprot
# ---------------------------------------------------------------------------
step_banner "Downloading afdb-swissprot"

SP_DIR="${FOLDSEEK_DB_ROOT}/afdb_swissprot_db"
if [[ -f "${SP_DIR}/afdb_swissprot.dbtype" ]]; then
    ok "afdb-swissprot already exists at ${SP_DIR}"
else
    mkdir -p "${SP_DIR}"
    foldseek databases Alphafold/Swiss-Prot "${SP_DIR}/afdb_swissprot" "${TMPDIR}/foldseek_tmp"
    ok "afdb-swissprot downloaded → ${SP_DIR}"
fi

# ---------------------------------------------------------------------------
# 3. Baktfold database (installed via baktfold tool)
# ---------------------------------------------------------------------------
# Baktfold is a full annotation tool (like phold but for bacteria/general genomes).
# It bundles SwissProt, AFDB clusters, PDB, and CATH into its own FoldSeek-format DB.
# Database is managed by `baktfold install`; do NOT try to download it manually.
#
# baktfold v0.3.0: https://github.com/gbouras13/baktfold/releases/tag/v0.3.0
#   Install:  micromamba install baktfold -c bioconda -c conda-forge  (in pharokka env)
#   Database: baktfold install -d <BAKTFOLD_DB_DIR> -t <threads>
# ---------------------------------------------------------------------------
step_banner "Installing Baktfold database"

BAKTFOLD_DB_DIR="${FOLDSEEK_DB_ROOT}/baktfold_db"
# baktfold install writes a marker file; check for it
if [[ -f "${BAKTFOLD_DB_DIR}/baktfold.dbtype" ]] || \
   [[ -n "$(ls "${BAKTFOLD_DB_DIR}"/*.dbtype 2>/dev/null | head -1)" ]]; then
    ok "Baktfold DB already installed at ${BAKTFOLD_DB_DIR}"
else
    mkdir -p "${BAKTFOLD_DB_DIR}"
    echo "  Running: baktfold install -d ${BAKTFOLD_DB_DIR} -t 8"
    echo "  (This downloads ~5–10 GB — takes 10–20 min depending on network)"
    baktfold install -d "${BAKTFOLD_DB_DIR}" -t 8
    ok "Baktfold DB installed → ${BAKTFOLD_DB_DIR}"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "════════════════════════════════════════════════════════"
echo "  FoldSeek + Baktfold databases summary"
echo "════════════════════════════════════════════════════════"
echo ""
echo "  DB root: ${FOLDSEEK_DB_ROOT}"
du -sh "${FOLDSEEK_DB_ROOT}"/* 2>/dev/null || true
echo ""
echo "Ensure databases.foldseek_db_root in config/config.yaml points here."
echo "Then run: sbatch steps/02_foldseek_3di.sh"
