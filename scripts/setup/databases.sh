#!/usr/bin/env bash
# =============================================================================
# databases.sh — download the databases phageFACTor needs (pick what you need).
# Called by `phagefactor databases`, which sets the locations below; it can
# also be run directly with them exported.
#
#   --pharokka   Pharokka annotation DB (genome mode; ~1-5 GB)
#   --phold      PHold DB (~2 GB; otherwise phold's own default location is used)
#   --pdb100     FoldSeek pdb100 (experimental PDB structures, ~4 GB)
#   --swissprot  FoldSeek afdb-swissprot (curated SwissProt AFDB, ~1.5 GB)
#   --afdb50     FoldSeek afdb50 = AlphaFold DB UniProt50-minimal
#                (~40-60 GB, ~350 GB free space needed during the download)
#   --taxdump    NCBI taxonomy (nodes.dmp + merged.dmp, ~70 MB download): the
#                domain of FoldSeek hits that carry no lineage (search_mode=webapi,
#                runs searched before v1.1). Local searches read the lineage from
#                the FoldSeek databases themselves.
#   --index      afdb50 3Di k-mer index: needs ~200 GB RAM, run it as its own
#                job, e.g.  sbatch --mem=200G --time=04:00:00 --wrap "phagefactor databases --index"
#   --all        everything except --index
# search_mode=webapi needs none of the FoldSeek databases (docs/modes.md).
#
# Locations (env): PHAROKKA_DB, PHOLD_DB, FOLDSEEK_DB_ROOT, TAXDUMP_DIR, SCRATCH_TMPDIR
# Tools: foldseek / phold on PATH, else run through $PF_TOOLS ("micromamba run -n <env>").
# =============================================================================
set -euo pipefail

: "${PHAROKKA_DB:?PHAROKKA_DB not set}" "${PHOLD_DB:?PHOLD_DB not set}" "${FOLDSEEK_DB_ROOT:?FOLDSEEK_DB_ROOT not set}"
SCRATCH_TMPDIR="${SCRATCH_TMPDIR:-${TMPDIR:-/tmp}}"
PF_TOOLS="${PF_TOOLS:-}"

log()  { echo "[$(date '+%H:%M:%S')] $*"; }
warn() { echo "[WARN]  $*" >&2; }
die()  { echo "[ERROR] $*" >&2; exit 1; }
banner() { echo; echo "== $* =="; }
tool() {  # tool <name>: the command to run it with
    if command -v "$1" &>/dev/null; then echo "$1"
    elif [[ -n "${PF_TOOLS}" ]]; then echo "${PF_TOOLS} $1"
    else die "$1 not found: activate the phageFACTor tools env (phagefactor install-envs)"; fi
}

DO_PHAROKKA=0; DO_PHOLD=0; DO_PDB100=0; DO_SWISSPROT=0; DO_AFDB50=0; DO_INDEX=0; DO_TAXDUMP=0
for arg in "$@"; do
    case "${arg}" in
        --pharokka)  DO_PHAROKKA=1 ;;
        --phold)     DO_PHOLD=1 ;;
        --pdb100)    DO_PDB100=1 ;;
        --swissprot) DO_SWISSPROT=1 ;;
        --afdb50)    DO_AFDB50=1 ;;
        --index)     DO_INDEX=1 ;;
        --taxdump)   DO_TAXDUMP=1 ;;
        --all)       DO_PHAROKKA=1; DO_PHOLD=1; DO_PDB100=1; DO_SWISSPROT=1; DO_AFDB50=1; DO_TAXDUMP=1 ;;
        *) die "unknown flag: ${arg}. Valid: --pharokka --phold --pdb100 --swissprot --afdb50 --taxdump --index --all" ;;
    esac
done
(( DO_PHAROKKA + DO_PHOLD + DO_PDB100 + DO_SWISSPROT + DO_AFDB50 + DO_INDEX + DO_TAXDUMP > 0 )) \
    || die "no database selected (e.g. --pharokka --phold --pdb100 --swissprot, or --all)"

# --- pharokka: curl, not install_databases.py: behind some authenticated proxies
# pharokka's `requests` downloader mis-parses proxy passwords with '+' or '='.
if (( DO_PHAROKKA )); then
    banner "Pharokka database -> ${PHAROKKA_DB}"
    if [[ -f "${PHAROKKA_DB}/VERSION_1_8_0" || -f "${PHAROKKA_DB}/phrog_annot_v4.tsv" ]]; then
        log "already installed"
    else
        mkdir -p "${PHAROKKA_DB}"
        TAR="${PHAROKKA_DB}/pharokka_v1.8.0_databases.tar.gz"
        curl --location --output "${TAR}" --retry 3 --retry-delay 5 \
            "https://zenodo.org/record/17110353/files/pharokka_v1.8.0_databases.tar.gz"
        ( cd "${PHAROKKA_DB}" && tar -xzf "${TAR}" --strip-components=1 )
        rm -f "${TAR}"
        for f in VERSION_1_8_0 phrog_annot_v4.tsv; do
            [[ -f "${PHAROKKA_DB}/${f}" ]] || warn "expected file missing: ${PHAROKKA_DB}/${f}"
        done
        log "installed"
    fi
fi

if (( DO_PHOLD )); then
    banner "PHold database -> ${PHOLD_DB}"
    mkdir -p "${PHOLD_DB}"
    PHOLD=$(tool phold)
    ${PHOLD} install -d "${PHOLD_DB}"
fi

fs_download() {  # fs_download <foldseek name> <dir> <prefix>
    local name="$1" dir="${FOLDSEEK_DB_ROOT}/$2" prefix="$3" tmp
    banner "FoldSeek ${name} -> ${dir}"
    if [[ -f "${dir}/${prefix}.dbtype" ]]; then log "already installed"; return; fi
    mkdir -p "${dir}"
    # tmp on the SAME filesystem: foldseek renames from tmp into place, a
    # cross-filesystem rename becomes a full copy (~350 GB for afdb50)
    tmp="${dir}/tmp_download"; mkdir -p "${tmp}"
    FOLDSEEK=$(tool foldseek)
    TMPDIR="${tmp}" ${FOLDSEEK} databases "${name}" "${dir}/${prefix}" "${tmp}"
    rm -rf "${tmp}"
    log "installed"
}
if (( DO_PDB100 )); then fs_download PDB pdb100_db pdb100; fi
if (( DO_SWISSPROT )); then fs_download Alphafold/Swiss-Prot afdb_swissprot_db afdb_swissprot; fi
if (( DO_AFDB50 )); then
    mkdir -p "${FOLDSEEK_DB_ROOT}"
    avail=$(df -BG "${FOLDSEEK_DB_ROOT}" | awk 'NR==2 {gsub("G","",$4); print $4}')
    if (( avail < 400 )); then warn "${avail} GB free on ${FOLDSEEK_DB_ROOT}; afdb50 needs ~400 GB during the download"; fi
    fs_download Alphafold/UniProt50 afdb50_db afdb50
    [[ -f "${FOLDSEEK_DB_ROOT}/afdb50_db/afdb50_taxonomy" ]] \
        || warn "no afdb50_taxonomy file: hits will have no taxon names"
fi

if (( DO_TAXDUMP )); then
    : "${TAXDUMP_DIR:?TAXDUMP_DIR not set}"
    banner "NCBI taxdump -> ${TAXDUMP_DIR}"
    if [[ -f "${TAXDUMP_DIR}/nodes.dmp" ]]; then
        log "already installed"
    else
        mkdir -p "${TAXDUMP_DIR}"
        TAR="${TAXDUMP_DIR}/taxdump.tar.gz"
        curl --location --output "${TAR}" --retry 3 --retry-delay 5 \
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        # names.dmp is not needed: only the parent tree is walked
        ( cd "${TAXDUMP_DIR}" && tar -xzf "${TAR}" nodes.dmp merged.dmp )
        rm -f "${TAR}"
        [[ -f "${TAXDUMP_DIR}/nodes.dmp" ]] || warn "nodes.dmp missing in ${TAXDUMP_DIR}"
        log "installed"
    fi
fi

if (( DO_INDEX )); then
    DB="${FOLDSEEK_DB_ROOT}/afdb50_db/afdb50"
    banner "FoldSeek createindex for ${DB}"
    [[ -f "${DB}.dbtype" ]] || die "afdb50 not found at ${DB} (download it with --afdb50 first)"
    rm -f "${DB}_ss.idx.0" "${DB}_ss.idx.index.0"          # fragments of an interrupted run
    tmp="${SCRATCH_TMPDIR}/fsindex_${SLURM_JOB_ID:-$$}"; mkdir -p "${tmp}"
    # --index-exclude 2: no C-alpha coordinates (the ProstT5 query DB has none;
    # saves ~100 GB of index)
    FOLDSEEK=$(tool foldseek)
    ${FOLDSEEK} createindex "${DB}" "${tmp}" --index-exclude 2 --threads "${SLURM_CPUS_PER_TASK:-16}"
    rm -rf "${tmp}"
    ls -lh "${DB}"*idx* 2>/dev/null || warn "no index files written"
fi

echo
log "Done. Check with: phagefactor check"
