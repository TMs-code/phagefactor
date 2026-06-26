#!/usr/bin/env bash
# =============================================================================
# 00d_download_afdb50.sh — Download AlphaFold DB (UniProt50) via FoldSeek
# =============================================================================
#
# Downloads the pre-built FoldSeek AFDBClusters/UniProt50 database (~80 GB).
# This DB:
#   - Covers ~2.3M representative bacterial/archaeal/viral/eukaryotic proteins
#     clustered at 50% sequence identity by MMseqs2 (sequence-based, not structural)
#   - Comes with FULL embedded taxonomy (taxid + taxname for every hit)
#   - Is the same source as the local pipeline's afdb50 (afdb-swissprot uses only
#     SwissProt ~570k; afdb50 is the broader TrEMBL+SwissProt representative set)
#
# Why afdb50 instead of (or in addition to) baktfold AFDBClusters:
#   - baktfold uses structural clustering (FoldSeek, 90% qcov + e<0.01):
#       • Cluster reps chosen by highest pLDDT — member annotations NOT searchable
#       • No taxonomy in installed DB (baktfold install doesn't run createtaxdb)
#   - afdb50 uses sequence clustering (MMseqs2, 50% identity):
#       • Better annotation diversity per cluster (sequence clusters are coarser)
#       • Taxonomy embedded — same_host correctly set for Campylobacter hits
#       • ~2.3M entries vs ~8M in structural clusters (faster search)
#
# Optional: restrict searches to bacteria+archaea+viruses only.
# After download, the search command in 02d_foldseek_3di.py supports a
# per-DB taxon filter via foldseek_taxon_filter in config/config.yaml.
# Set FOLDSEEK_TAXON_FILTER["afdb50"] = "2,2157,10239"  (Bacteria,Archaea,Viruses)
# to skip eukaryotic hits at search time (requires FoldSeek --taxon-list flag).
#
# Run INTERACTIVELY on a login node (large download, not compute):
#   bash steps/00d_download_afdb50.sh
# =============================================================================
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
fi
source "${PROJECT_DIR}/config.sh"

activate_env
require_tool foldseek

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[OK]${NC}  $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }

# ---------------------------------------------------------------------------
# Destination: sibling of baktfold_db/ in the same FoldSeek DB root
# ---------------------------------------------------------------------------
AFDB50_DIR="${FOLDSEEK_DB_ROOT}/afdb50_db"
AFDB50_PREFIX="${AFDB50_DIR}/afdb50"

if [[ -f "${AFDB50_PREFIX}.dbtype" ]]; then
    ok "afdb50 already exists at ${AFDB50_DIR}"
    echo "To force re-download, delete ${AFDB50_DIR} and re-run this script."
    exit 0
fi

mkdir -p "${AFDB50_DIR}"

# CRITICAL: TMPDIR must be on the SAME filesystem as AFDB50_DIR.
# foldseek downloads + extracts into TMPDIR, then renames files to the output
# prefix.  If TMPDIR is on a different filesystem (e.g. local node scratch),
# the final step becomes a full cross-filesystem copy of ~350 GB — very slow
# and likely to be killed by SLURM.  Using a subfolder of the destination
# makes the rename instantaneous.
TMPDIR="${AFDB50_DIR}/tmp_download"
export TMPDIR
mkdir -p "${TMPDIR}"

# Check available disk space (need ~400 GB free: ~350 GB tarball + extracted)
AVAIL_GB=$(df -BG "${AFDB50_DIR}" | awk 'NR==2 {gsub("G","",$4); print $4}')
echo "  Available space on destination filesystem: ${AVAIL_GB} GB"
if [[ "${AVAIL_GB}" -lt 400 ]]; then
    warn "Less than 400 GB available — download may fail."
    warn "Proceeding anyway; monitor with: watch -n 60 'du -sh ${AFDB50_DIR}'"
fi

step_banner "Downloading AlphaFold DB UniProt50 (~350 GB) — allow 2–4 hours"
echo ""
echo "  Destination : ${AFDB50_DIR}"
echo "  Temp dir    : ${TMPDIR}  (same filesystem — no cross-FS copy needed)"
echo "  FoldSeek cmd: foldseek databases Alphafold/UniProt50 ${AFDB50_PREFIX} ${TMPDIR}"
echo ""
echo "  Taxonomy is embedded — no createtaxdb step needed."
echo "  Monitor progress in another terminal:"
echo "    watch -n 60 'du -sh ${AFDB50_DIR}; ls -lh ${AFDB50_DIR}'"
echo ""

foldseek databases Alphafold/UniProt50 "${AFDB50_PREFIX}" "${TMPDIR}"

# Clean up the temp subfolder (foldseek may already remove it)
rm -rf "${TMPDIR}" 2>/dev/null || true

ok "afdb50 downloaded → ${AFDB50_DIR}"
echo ""

# ---------------------------------------------------------------------------
# Verify taxonomy companion exists
# ---------------------------------------------------------------------------
if [[ -f "${AFDB50_PREFIX}_taxonomy" ]]; then
    ok "Taxonomy file present: ${AFDB50_PREFIX}_taxonomy"
else
    warn "Taxonomy file NOT found at ${AFDB50_PREFIX}_taxonomy"
    warn "Hits will have empty taxname — run 00c_add_taxonomy_baktfold.py to fix."
fi

# ---------------------------------------------------------------------------
# Post-download: register the DB in config/config.yaml
# ---------------------------------------------------------------------------
echo ""
echo "════════════════════════════════════════════════════════"
echo "  Next steps:"
echo "════════════════════════════════════════════════════════"
echo ""
echo "1. Add the new DB under databases.foldseek_local_dbs in config/config.yaml"
echo "   (paths are relative to foldseek_db_root):"
echo "   foldseek_local_dbs:"
echo "     afdb-swissprot: afdb_swissprot_db/afdb_swissprot"
echo "     pdb100:         pdb100_db/pdb100"
echo "     afdb50:         afdb50_db/afdb50            # <-- NEW"
echo ""
echo "2. (Optional) Restrict to bacteria+archaea+viruses at search time."
echo "   Under databases.foldseek_taxon_filter in config/config.yaml:"
echo "     afdb50: \"2,2157,10239\"   # Bacteria,Archaea,Viruses"
echo "   02d_foldseek_3di.py passes --taxon-list when this is set."
echo ""
echo "3. Re-run step 02 with --force to regenerate hits including the new DB."
echo ""
du -sh "${AFDB50_DIR}"
