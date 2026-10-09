#!/usr/bin/env bash
# =============================================================================
# 00b_setup_databases.sh — One-time database setup (pick what you need)
#
# Replaces four separate scripts with one, controlled by flags, so you only
# install what your search mode actually needs:
#   --pharokka    Pharokka annotation DB (genome mode; ~1-5 GB, curl+tar)
#   --pdb100      FoldSeek pdb100 (experimental PDB structures, ~4 GB)
#   --swissprot   FoldSeek afdb-swissprot (curated SwissProt AFDB, ~1.5 GB)
#   --afdb50      FoldSeek afdb50 = AlphaFold DB UniProt50-minimal (~40-60 GB
#                 as shipped by FoldSeek's fetcher; ~350 GB during download
#                 before final layout — see disk-space check below)
#   --index       Precompute the afdb50 3Di k-mer index (needs a FAT node,
#                 ~180-200 GB RAM — see resource note below; run this as its
#                 own sbatch submission, not lumped with the others)
#   --all         Everything except --index (that one you opt into deliberately
#                 because of its RAM requirement)
#
# WebAPI search mode needs NONE of this — skip this script entirely and use
# SEARCH_MODE=webapi (see docs/modes.md).
#
# Usage (interactive, login node — downloads, not compute):
#   bash steps/00b_setup_databases.sh --pharokka --pdb100 --swissprot --afdb50
#
# The --index step needs a fat-memory node and should be its OWN submission,
# with resources overridden on the sbatch command line to match your
# cluster's large-memory partition (the file's own #SBATCH block is sized for
# the lightweight download operations, not this):
#   sbatch --partition=<your-fatnode-partition> --mem=200G --time=04:00:00 \
#       steps/00b_setup_databases.sh --index
# =============================================================================

#SBATCH -N 1
#SBATCH --partition=common
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -J pf_setup_db
#SBATCH -o logs/00b_setup_databases.out
#SBATCH -e logs/00b_setup_databases.err
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
# --mail-user omitted on purpose: submit_all.sh passes it (SLURM_EMAIL). For a
# manual run, add --mail-user=you@your.institute (some clusters reject others).

set -euo pipefail

# Resolve project root: SLURM_SUBMIT_DIR if sbatch-submitted (BASH_SOURCE points
# into the spool dir once queued), else this script's own parent dir.
PROJECT_DIR="${PHAGEFACTOR_HOME:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
source "${PROJECT_DIR}/config.sh"

DO_PHAROKKA=0; DO_PDB100=0; DO_SWISSPROT=0; DO_AFDB50=0; DO_INDEX=0
for arg in "$@"; do
    case "${arg}" in
        --pharokka)  DO_PHAROKKA=1 ;;
        --pdb100)    DO_PDB100=1 ;;
        --swissprot) DO_SWISSPROT=1 ;;
        --afdb50)    DO_AFDB50=1 ;;
        --index)     DO_INDEX=1 ;;
        --all)       DO_PHAROKKA=1; DO_PDB100=1; DO_SWISSPROT=1; DO_AFDB50=1 ;;
        *) die "Unknown flag: ${arg}. Valid: --pharokka --pdb100 --swissprot --afdb50 --index --all" ;;
    esac
done
if [[ $((DO_PHAROKKA + DO_PDB100 + DO_SWISSPROT + DO_AFDB50 + DO_INDEX)) -eq 0 ]]; then
    die "No flags given. See the usage note at the top of this script, or run with --all."
fi

# =============================================================================
# --pharokka — Pharokka annotation database
# =============================================================================
# curl (not `install_databases.py -o`) because behind some authenticated HTTP
# proxies, pharokka's Python `requests` downloader mis-parses proxy URLs whose
# password contains `+`/`=` (which must be URL-encoded); curl tolerates them.
if [[ "${DO_PHAROKKA}" == "1" ]]; then
    step_banner "Pharokka database"
    PHAROKKA_TAR_URL="https://zenodo.org/record/17110353/files/pharokka_v1.8.0_databases.tar.gz"
    TAR_PATH="${PHAROKKA_DB}/pharokka_v1.8.0_databases.tar.gz"
    mkdir -p "${PHAROKKA_DB}"

    if [[ -f "${PHAROKKA_DB}/VERSION_1_8_0" ]] || [[ -f "${PHAROKKA_DB}/phrog_annot_v4.tsv" ]]; then
        log "Pharokka DB appears already installed at ${PHAROKKA_DB} -- skipping."
    else
        log "Downloading ${PHAROKKA_TAR_URL}..."
        curl --location --output "${TAR_PATH}" --retry 3 --retry-delay 5 "${PHAROKKA_TAR_URL}"
        log "Download complete: $(ls -lh "${TAR_PATH}")"
        log "Extracting..."
        ( cd "${PHAROKKA_DB}" && tar -xzf "${TAR_PATH}" --strip-components=1 )
        rm -f "${TAR_PATH}"
        for f in VERSION_1_8_0 phrog_annot_v4.tsv; do
            [[ -f "${PHAROKKA_DB}/${f}" ]] || warn "Expected file missing: ${PHAROKKA_DB}/${f}"
        done
        log "Pharokka DB installed -> ${PHAROKKA_DB}"
    fi
fi

# =============================================================================
# --pdb100 / --swissprot / --afdb50 — FoldSeek databases
# =============================================================================
if [[ "${DO_PDB100}" == "1" || "${DO_SWISSPROT}" == "1" || "${DO_AFDB50}" == "1" ]]; then
    activate_env
    require_tool foldseek
    mkdir -p "${FOLDSEEK_DB_ROOT}"
fi

if [[ "${DO_PDB100}" == "1" ]]; then
    step_banner "FoldSeek pdb100"
    PDB_DIR="${FOLDSEEK_DB_ROOT}/pdb100_db"
    if [[ -f "${PDB_DIR}/pdb100.dbtype" ]]; then
        log "pdb100 already exists at ${PDB_DIR} -- skipping."
    else
        mkdir -p "${PDB_DIR}"
        TMPDIR="${SCRATCH_TMPDIR}" foldseek databases PDB "${PDB_DIR}/pdb100" "${SCRATCH_TMPDIR}/foldseek_tmp"
        log "pdb100 downloaded -> ${PDB_DIR}"
    fi
fi

if [[ "${DO_SWISSPROT}" == "1" ]]; then
    step_banner "FoldSeek afdb-swissprot"
    SP_DIR="${FOLDSEEK_DB_ROOT}/afdb_swissprot_db"
    if [[ -f "${SP_DIR}/afdb_swissprot.dbtype" ]]; then
        log "afdb-swissprot already exists at ${SP_DIR} -- skipping."
    else
        mkdir -p "${SP_DIR}"
        TMPDIR="${SCRATCH_TMPDIR}" foldseek databases Alphafold/Swiss-Prot "${SP_DIR}/afdb_swissprot" "${SCRATCH_TMPDIR}/foldseek_tmp"
        log "afdb-swissprot downloaded -> ${SP_DIR}"
    fi
fi

if [[ "${DO_AFDB50}" == "1" ]]; then
    step_banner "FoldSeek afdb50 (AlphaFold DB UniProt50-minimal, ~350 GB during download)"
    AFDB50_DIR="${FOLDSEEK_DB_ROOT}/afdb50_db"
    AFDB50_PREFIX="${AFDB50_DIR}/afdb50"
    if [[ -f "${AFDB50_PREFIX}.dbtype" ]]; then
        log "afdb50 already exists at ${AFDB50_DIR} -- skipping."
        log "(To force re-download, delete ${AFDB50_DIR} and re-run with --afdb50.)"
    else
        mkdir -p "${AFDB50_DIR}"
        # TMPDIR must be on the SAME filesystem as AFDB50_DIR: foldseek downloads
        # + extracts into TMPDIR then renames into the output prefix, and a
        # cross-filesystem rename becomes a full ~350 GB copy instead of an
        # instant move.
        AFDB_TMP="${AFDB50_DIR}/tmp_download"
        mkdir -p "${AFDB_TMP}"
        AVAIL_GB=$(df -BG "${AFDB50_DIR}" | awk 'NR==2 {gsub("G","",$4); print $4}')
        log "Available space on destination filesystem: ${AVAIL_GB} GB (need ~400 GB free)"
        [[ "${AVAIL_GB}" -lt 400 ]] && warn "Less than 400 GB available -- download may fail."
        log "Downloading (allow 2-4 hours)... destination=${AFDB50_DIR}, temp=${AFDB_TMP}"
        TMPDIR="${AFDB_TMP}" foldseek databases Alphafold/UniProt50 "${AFDB50_PREFIX}" "${AFDB_TMP}"
        rm -rf "${AFDB_TMP}" 2>/dev/null || true
        log "afdb50 downloaded -> ${AFDB50_DIR}"
        if [[ -f "${AFDB50_PREFIX}_taxonomy" ]]; then
            log "Taxonomy file present: ${AFDB50_PREFIX}_taxonomy"
        else
            warn "Taxonomy file NOT found at ${AFDB50_PREFIX}_taxonomy -- hits will have empty taxname."
        fi
        log "Add it under databases.foldseek_local_dbs in config/config.yaml:"
        log "    afdb50: afdb50_db/afdb50"
    fi
fi

# =============================================================================
# --index — Precompute FoldSeek's 3Di k-mer index for afdb50
# =============================================================================
# Needs a fat-memory node (~180-200 GB RAM). See the usage note at the top of
# this file for how to override --partition/--mem/--time on the sbatch
# command line -- the file's own #SBATCH block is sized for the lightweight
# download flags above, not this.
if [[ "${DO_INDEX}" == "1" ]]; then
    activate_env
    step_banner "FoldSeek createindex for afdb50"
    AFDB50_DB="${FOLDSEEK_DB_ROOT}/afdb50_db/afdb50"
    IDX_TMP="${SCRATCH_TMPDIR:-/local/scratch/tmp}/fsindex_${SLURM_JOB_ID:-manual}"
    mkdir -p "${IDX_TMP}"

    for f in "${AFDB50_DB}_ss.idx.0" "${AFDB50_DB}_ss.idx.index.0"; do
        if [[ -f "${f}" ]]; then
            log "Removing incomplete index fragment: ${f}"
            rm -f "${f}"
        fi
    done

    log "DB=${AFDB50_DB}  tmp=${IDX_TMP}  threads=${SLURM_CPUS_PER_TASK:-16}"
    # --index-exclude 2 skips Cα coordinates (only needed for TM-score/
    # --sort-by-structure-bits scoring; our ProstT5 query DB has no _ca, so
    # excluding it saves ~100 GB of index).
    foldseek createindex "${AFDB50_DB}" "${IDX_TMP}" \
        --index-exclude 2 --threads "${SLURM_CPUS_PER_TASK:-16}"
    rm -rf "${IDX_TMP}" 2>/dev/null || true

    log "Index files:"
    for sfx in ".idx" ".idx.dbtype" ".idx.index" "_ss.idx" "_ss.idx.dbtype" "_ss.idx.index"; do
        f="${AFDB50_DB}${sfx}"
        if [[ -f "${f}" ]]; then
            printf "  OK   %-25s %s\n" "${f##*/}" "$(du -sh "${f}" | cut -f1)"
        else
            printf "  MISS %s\n" "${f##*/}"
        fi
    done
fi

log ""
log "Done. Point config/config.yaml databases.foldseek_db_root (or FOLDSEEK_DB_ROOT)"
log "at ${FOLDSEEK_DB_ROOT}, confirm databases.foldseek_local_dbs, then:"
log "  sbatch steps/02_foldseek_3di.sh"
