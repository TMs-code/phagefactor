#!/usr/bin/env python3
"""
02d_foldseek_3di.py
===================
Local FoldSeek annotation pipeline using ProstT5 3di structural tokens.

Instead of folding proteins with ESMFold/AF3, this script uses the 3di
structural tokens computed by phold (stored in phold_3di.fasta) as direct
query input to a local FoldSeek installation.  This removes the 400aa
bottleneck and gives full coverage of all 313 hypothetical proteins,
including the 16 >400aa proteins still waiting for AF3 weights.

Scientific basis:
  ProstT5 3Di tokens give annotation quality close to actual 3D structures in
  FoldSeek searches, at a fraction of the cost. phold already runs ProstT5
  internally — we reuse its 3Di output rather than folding proteins ourselves.

FoldSeek query mode:
  Primary  : --prostt5-input    (FoldSeek >= 8, recommended)
  Fallback : --seq-type 1       (3di-only alignment, any FoldSeek version)
  Optional : --combined         (3Di + AA; requires manual DB construction,
                                 pass --combined flag to enable)

Steps:
  1. Filter phold_3di.fasta / phold_aa.fasta to hypothetical proteins only
  2. Write query FASTA(s) to 02_foldseek/3di_tokens/query_fastas/
  3. For each configured local database: run foldseek easy-search
     (skips if per-DB result file already exists — safe to re-run)
  4. Split combined m8 results into per-gene files
  5. Parse + score hits with the same logic as 02b_foldseek_pipeline.py
  6. Write summary CSVs to 02_foldseek/3di_tokens/

Output files:
  02_foldseek/3di_tokens/best_hit.csv   -- best informative hit per gene
  02_foldseek/3di_tokens/top3.csv       -- top-3 informative hits per gene
  02_foldseek/3di_tokens/all_hits.csv   -- all raw hits (manual inspection)
  02_foldseek/3di_tokens/raw_hits/      -- per-gene .m8 files

Usage:
  cd phagefactor/
  python scripts/02d_foldseek_3di.py

  # Parse only (no new foldseek calls):
  python scripts/02d_foldseek_3di.py --parse_only

  # Use combined 3Di+AA mode (slower DB construction, better sensitivity):
  python scripts/02d_foldseek_3di.py --combined

  # Override number of CPU threads:
  python scripts/02d_foldseek_3di.py --threads 16

  # Restrict to a single database for quick testing:
  python scripts/02d_foldseek_3di.py --dbs afdb-swissprot

Prerequisites:
  foldseek installed and on PATH (or set FOLDSEEK_CMD in config.py)
  At least one local FoldSeek database downloaded (see FOLDSEEK_LOCAL_DBS in config.py)
  For combined mode: both phold_3di.fasta AND phold_aa.fasta must exist
"""

import sys
import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).parent
_PROJECT_DIR = _SCRIPTS_DIR.parent
sys.path.insert(0, str(_PROJECT_DIR))

from config import (
    HYPO_GENE_LIST,
    PHOLD_OUT_DIR,
    PHOLD_3DI_FASTA, PHOLD_AA_FASTA,
    FOLDSEEK_3DI_DIR, FOLDSEEK_3DI_BEST, FOLDSEEK_3DI_TOP3, FOLDSEEK_3DI_ALL,
    FOLDSEEK_EVALUE_MAX, FOLDSEEK_SCORE_OVERRIDE,
    PROSTT5_MASK_THRESHOLD,
    FOLDSEEK_TAXON_FILTER,
    HOST_GENUS,
    FOLDSEEK_CMD, FOLDSEEK_THREADS, FOLDSEEK_LOCAL_DBS,
)
from utils import log, section

try:
    import pandas as pd
    from Bio import SeqIO
except ImportError as e:
    print(f"Missing dependency: {e}\n  pip install pandas biopython")
    sys.exit(1)


# =============================================================================
# IMPORT SCORING FUNCTIONS (local module)
# =============================================================================
# These helpers previously lived in an earlier 02b_foldseek_pipeline.py (alongside
# ESMFold + FoldSeek-Web-API logic). We keep the local-FoldSeek path
# self-contained by hosting them in foldseek_scoring.py.

from foldseek_scoring import (
    _is_informative_fs,
    _phage_boost_factor,
    _is_same_host_hit,
    _compute_fs_confidence,
    _upgrade_description,
    _DEFENSE_PATTERN,
    _is_promiscuous_fold_hit,
    _has_eukaryotic_description,
    build_best_and_top3,
)


# =============================================================================
# DIRECTORIES
# =============================================================================

QUERY_FASTA_DIR   = FOLDSEEK_3DI_DIR / "query_fastas"
QUERY_DB_DIR      = FOLDSEEK_3DI_DIR / "query_db"
PER_DB_HITS_DIR   = FOLDSEEK_3DI_DIR / "per_db_hits"    # one .m8 per DB
RAW_HITS_DIR      = FOLDSEEK_3DI_DIR / "raw_hits"       # one .m8 per gene (merged)

QUERY_3DI_FASTA   = QUERY_FASTA_DIR / "hypo_3di.fasta"
QUERY_AA_FASTA    = QUERY_FASTA_DIR / "hypo_aa.fasta"


# =============================================================================
# LOCAL FOLDSEEK OUTPUT FORMAT
# =============================================================================
# 15-column format for local FoldSeek searches.  We carry the same kind of
# information as the 21-column web API m8 format used by 02b, minus the
# fields that need structural inputs convertalis can't compute from our
# sequence-only queryDB:
#
#   - lddt   : per-alignment LDDT.  Needs queryDB_ca (alpha-carbon coords).
#              Our query side has no structures -- only ProstT5 3Di tokens.
#   - prob   : ProstT5/structure-derived homology probability.  Grouped in
#              foldseek's help under the same family as lddt/qca/tca/rmsd
#              (`foldseek convertalis --help` shows it next to those).
#              Requesting it also triggers the _ca file lookup -- THIS is
#              the column that was making our convertalis fail even after
#              we dropped lddt and added a _ca stub.
#   - taxid/taxname : OPTIONAL per-DB.  Require a <db>_taxonomy file or
#              a global NCBI taxdump.  afdb-swissprot and pdb100 have it;
#              baktfold-afdb does not.  We build the per-DB format on the
#              fly and pad missing rows with empty tax cells so all merged
#              m8 lines have a uniform column layout.
#
# Notes on remaining differences from web API m8 (M8_COLS in utils.py):
#   fident  : fraction 0-1 (web API column "pident" is percentage 0-100)
#   qcov    : query coverage fraction 0-1 (web API "qcov_aa" is absolute AA count)
#   theader : full target FASTA header "ACCESSION DESCRIPTION" -- same split as web API
#   taxname : empty if target DB has no taxonomy (treated as NA downstream)

LOCAL_FMT_CORE = (
    "query,theader,fident,alnlen,mismatch,gapopen,"
    "qstart,qend,tstart,tend,"
    "evalue,bits,qcov"
)
LOCAL_FMT_TAX = "taxid,taxname"
LOCAL_COLS_CORE = LOCAL_FMT_CORE.split(",")   # 13
LOCAL_COLS_TAX  = LOCAL_FMT_TAX.split(",")    # 2
LOCAL_COLS = LOCAL_COLS_CORE + LOCAL_COLS_TAX  # 15
# We append an extra column 'source_db' at write time recording which
# FOLDSEEK_LOCAL_DBS entry each hit came from (afdb-swissprot / pdb100 /
# baktfold-afdb).  parse_local_m8_file tolerates rows with or without this
# 16th column.
LOCAL_COLS_OUT = LOCAL_COLS + ["source_db"]    # 16
N_LOCAL_COLS = len(LOCAL_COLS_OUT)             # 16


def _db_has_taxonomy(db_path: Path) -> bool:
    """Return True if the FoldSeek target DB has a `_taxonomy` companion file."""
    return Path(str(db_path) + "_taxonomy").exists()


def _pad_m8_taxonomy(m8_path: Path, n_target_cols: int) -> None:
    """
    Pad each non-comment row in an m8 file to exactly `n_target_cols` columns
    by appending empty tab-separated cells.  Used after convertalis against a
    DB without taxonomy so the resulting file is column-compatible with the
    DBs that DO have taxonomy.
    """
    if not m8_path.exists() or m8_path.stat().st_size == 0:
        return
    out_lines = []
    changed = False
    for line in m8_path.read_text().splitlines():
        if not line or line.startswith("#"):
            out_lines.append(line)
            continue
        parts = line.split("\t")
        if len(parts) < n_target_cols:
            parts = parts + [""] * (n_target_cols - len(parts))
            changed = True
        out_lines.append("\t".join(parts))
    if changed:
        m8_path.write_text("\n".join(out_lines) + "\n")


# Eukaryote keyword set — mirrors step 03's _EUKARYOTE_KEYWORDS
_EUKA_KW = frozenset({
    "homo sapiens", "mus musculus", "saccharomyces", "arabidopsis",
    "drosophila", "caenorhabditis", "danio rerio", "xenopus",
    "gallus gallus", "bos taurus", "rattus norvegicus",
    "metazoa", "viridiplantae", "fungi", "eukaryot",
    "mammalia", "chordata", "synthetic construct",
})


def _is_eukaryotic_taxname(taxname: str) -> bool:
    """
    Return True if the taxname looks like a eukaryote (or synthetic construct).
    Used to apply FOLDSEEK_TAXON_FILTER at parse time for DBs like afdb50.
    """
    if not taxname or not isinstance(taxname, str):
        return False
    t = taxname.lower()
    return any(kw in t for kw in _EUKA_KW)


def _build_taxon_filter_set(taxon_filter_str: str) -> frozenset:
    """
    Parse comma-separated taxon IDs string (e.g. '2,2157,10239') into a frozenset
    of strings for fast membership testing.  Returns empty frozenset if blank.
    """
    if not taxon_filter_str:
        return frozenset()
    return frozenset(t.strip() for t in taxon_filter_str.split(",") if t.strip())


# -----------------------------------------------------------------------------
# Source-DB tag normalization
# -----------------------------------------------------------------------------
# Each entry in FOLDSEEK_LOCAL_DBS is a single-source FoldSeek database
# (afdb-swissprot, pdb100, baktfold-afdb).  We do not search any "bundled"
# DBs anymore, so the source_db name is itself the right tag for the
# foldseek_subdb column -- no accession-prefix parsing needed.
#
# Historical note:  an earlier design treated baktfold_db/ as a SINGLE bundled
# DB and tried to split hits into baktfold/{afdb50,swissprot,pdb,cath} by
# accession prefix.  We now know baktfold_db/ contains four INDEPENDENT
# FoldSeek databases.  Of those, we keep only AFDBClusters (= baktfold-afdb)
# because baktfold/swissprot and baktfold/pdb overlap with our standalone
# afdb_swissprot_db / pdb100_db, and baktfold/cath is out of scope.
# -----------------------------------------------------------------------------

def _classify_subdb(source_db: str, accession: str) -> str:
    """
    Return a clean source-DB tag for the `foldseek_subdb` column.

    Each FOLDSEEK_LOCAL_DBS entry is now a single-source DB, so we just
    normalize the source name (lowercase, underscores -> dashes).  The
    `accession` argument is kept in the signature for backward compatibility
    with previous bundled-DB logic but is ignored.
    """
    src = (source_db or "").strip().lower().replace("_", "-")
    return src or "unknown"


# =============================================================================
# STEP 1 — PREPARE QUERY FASTAS
# =============================================================================

def _parse_3di_id(header: str) -> str:
    """
    Extract locus_tag from a phold_3di.fasta header.

    phold_3di.fasta headers have the form:
        >CCIE1:THDMTYMJ_CDS_0001
    We want: THDMTYMJ_CDS_0001
    """
    bare = header.lstrip(">").strip()
    if ":" in bare:
        return bare.split(":", 1)[1]
    return bare


def _load_prostt5_probs(phold_out_dir: Path) -> dict:
    """
    Load per-residue ProstT5 3Di prediction probabilities from all per-prophage
    NDJSON files produced by phold (e.g. CCIE1_prostT5_3di_all_probabilities.json).

    Format (one JSON object per line):
        {"seq_id": "CCIE1:WRIRHRXM_CDS_0001", "probability": [97.75, 37.89, ...]}

    Returns dict: full_seq_id (with prophage prefix) -> list[float].
    Returns empty dict if no files found (masking gracefully disabled).
    """
    import json as _json

    probs: dict = {}
    # Recursive glob: genome mode writes 01_phold/<PROPHAGE>/..._probabilities.json
    # (1 level deep) but protein mode writes 01_phold/proteins/predictions/..._probabilities.json
    # (2 levels deep). The non-recursive "*/*..." pattern only matched the genome-mode
    # layout and silently disabled masking in protein mode.
    json_files = sorted(phold_out_dir.glob("**/*_prostT5_3di_all_probabilities.json"))
    if not json_files:
        log(f"  NOTE: No prostT5 per-residue probability files found under "
            f"{phold_out_dir} — per-residue masking disabled.")
        return probs

    for jf in json_files:
        try:
            with open(jf) as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        obj = _json.loads(line)
                        probs[obj["seq_id"]] = obj["probability"]
                    except (_json.JSONDecodeError, KeyError):
                        continue
        except OSError:
            log(f"  WARNING: could not read {jf}")

    log(f"  ProstT5 probs loaded: {len(probs):,} sequences "
        f"from {len(json_files)} files "
        f"(mask threshold: {PROSTT5_MASK_THRESHOLD}/100)")
    return probs


def _apply_3di_mask(seq: str, probs: list, threshold: float) -> str:
    """
    Replace 3Di characters whose ProstT5 confidence is below *threshold* with
    '*' (the FoldSeek wildcard/mask token).  Positions with no probability
    value (length mismatch) are left unchanged.

    Matches baktfold's --mask-threshold behaviour (default 25/100).
    """
    if not probs or len(probs) != len(seq):
        return seq
    return "".join("*" if p < threshold else c for c, p in zip(seq, probs))


def prepare_query_fastas(hypo_locus_tags: set, force: bool = False) -> tuple:
    """
    Filter phold_3di.fasta (and optionally phold_aa.fasta) to the hypothetical
    proteins, and write query FASTA files to QUERY_FASTA_DIR.

    Applies per-residue ProstT5 masking (PROSTT5_MASK_THRESHOLD) if probability
    files are available.  Masking replaces low-confidence 3Di tokens with '*'
    (FoldSeek mask), matching baktfold's --mask-threshold behaviour.

    Returns:
        (n_3di, n_aa)  — number of sequences written for each file.
        n_aa is 0 if phold_aa.fasta is not available.
    """
    QUERY_FASTA_DIR.mkdir(parents=True, exist_ok=True)

    if not PHOLD_3DI_FASTA.exists():
        log(f"ERROR: phold_3di.fasta not found at {PHOLD_3DI_FASTA}")
        log("  Run phold on the pharokka GBK first (step 01).")
        sys.exit(1)

    # Decide whether existing query FASTAs can be reused.  A stale empty file
    # from a previous failed run must NOT be reused -- we'd otherwise loop
    # back into the "0 hypothetical proteins" failure mode.
    def _seq_count(fa: Path) -> int:
        if not fa.exists() or fa.stat().st_size == 0:
            return 0
        try:
            return sum(1 for _ in SeqIO.parse(str(fa), "fasta"))
        except Exception:
            return 0

    # ---- Load per-residue ProstT5 probabilities for masking ----
    # Must happen before writing the 3di FASTA so masking is applied consistently.
    prostt5_probs: dict = {}
    if PROSTT5_MASK_THRESHOLD > 0:
        prostt5_probs = _load_prostt5_probs(PHOLD_OUT_DIR)

    # ---- 3di sequences ----
    existing_3di = _seq_count(QUERY_3DI_FASTA)
    if force or existing_3di == 0:
        written_3di  = 0
        masked_count = 0
        mask_stats_rows = []   # collected for per-gene CSV summary

        with open(QUERY_3DI_FASTA, "w") as fout:
            for rec in SeqIO.parse(str(PHOLD_3DI_FASTA), "fasta"):
                gene_id = _parse_3di_id(rec.id)
                if gene_id not in hypo_locus_tags:
                    continue
                seq = str(rec.seq)
                if PROSTT5_MASK_THRESHOLD > 0 and prostt5_probs:
                    # The prob dict is keyed by the FULL header "PROPHAGE:gene_id"
                    # (as stored in the phold JSON).  rec.id still has the prefix.
                    probs_for_gene = prostt5_probs.get(rec.id, [])
                    masked_seq = _apply_3di_mask(seq, probs_for_gene, PROSTT5_MASK_THRESHOLD)
                    n_masked  = sum(c == "*" for c in masked_seq)
                    n_total   = len(seq)
                    mean_prob = (sum(probs_for_gene) / len(probs_for_gene)
                                 if probs_for_gene else float("nan"))
                    if masked_seq != seq:
                        masked_count += 1
                    mask_stats_rows.append({
                        "gene_id":        gene_id,
                        "n_residues":     n_total,
                        "n_masked":       n_masked,
                        "frac_masked":    round(n_masked / n_total, 4) if n_total else 0,
                        "mean_prob":      round(mean_prob, 2),
                        "threshold":      PROSTT5_MASK_THRESHOLD,
                        "probs_available": bool(probs_for_gene),
                    })
                else:
                    masked_seq = seq
                    if PROSTT5_MASK_THRESHOLD > 0:
                        mask_stats_rows.append({
                            "gene_id": gene_id, "n_residues": len(seq),
                            "n_masked": 0, "frac_masked": 0.0,
                            "mean_prob": float("nan"),
                            "threshold": PROSTT5_MASK_THRESHOLD,
                            "probs_available": False,
                        })
                # Write with clean locus_tag as ID (no prophage prefix)
                fout.write(f">{gene_id}\n{masked_seq}\n")
                written_3di += 1

        # ---- Write per-gene masking stats CSV --------------------------------
        if mask_stats_rows and PROSTT5_MASK_THRESHOLD > 0:
            stats_csv = QUERY_FASTA_DIR / "masking_stats.csv"
            import csv as _csv
            with open(stats_csv, "w", newline="") as sfh:
                writer = _csv.DictWriter(sfh, fieldnames=list(mask_stats_rows[0]))
                writer.writeheader()
                writer.writerows(mask_stats_rows)

            # Summary statistics
            probs_rows = [r for r in mask_stats_rows if r["probs_available"]]
            total_res  = sum(r["n_residues"] for r in probs_rows)
            total_mask = sum(r["n_masked"]   for r in probs_rows)
            frac_genes = masked_count / len(probs_rows) if probs_rows else 0
            frac_res   = total_mask / total_res if total_res else 0
            mean_gene_prob = (sum(r["mean_prob"] for r in probs_rows if r["mean_prob"] == r["mean_prob"]) /
                              len([r for r in probs_rows if r["mean_prob"] == r["mean_prob"]])
                              if probs_rows else float("nan"))
            log(f"  Masking stats (threshold {PROSTT5_MASK_THRESHOLD}/100):")
            log(f"    Genes with prob data : {len(probs_rows)}/{len(mask_stats_rows)}")
            log(f"    Genes with ≥1 masked : {masked_count} ({frac_genes*100:.1f}%)")
            log(f"    Residues masked      : {total_mask}/{total_res} ({frac_res*100:.1f}%)")
            log(f"    Mean gene confidence : {mean_gene_prob:.1f}/100")
            log(f"    Full stats → {stats_csv}")

        found_ids = set()
        for rec in SeqIO.parse(str(QUERY_3DI_FASTA), "fasta"):
            found_ids.add(rec.id)
        missing_3di = sorted(hypo_locus_tags - found_ids)
        if missing_3di:
            log(f"  WARNING: {len(missing_3di)} hypothetical proteins not found "
                f"in phold_3di.fasta:")
            for g in missing_3di[:10]:
                log(f"    {g}")
            if len(missing_3di) > 10:
                log(f"    ... and {len(missing_3di) - 10} more")
        mask_note = (f", {masked_count} proteins had ≥1 residue masked"
                     if PROSTT5_MASK_THRESHOLD > 0 and prostt5_probs else "")
        log(f"  3di query FASTA: {written_3di}/{len(hypo_locus_tags)} genes written"
            f"{mask_note} → {QUERY_3DI_FASTA}")
    else:
        written_3di = existing_3di
        log(f"  3di query FASTA exists ({written_3di} seqs) — reuse")

    # ---- AA sequences (optional, for combined mode) ----
    written_aa = 0
    if PHOLD_AA_FASTA.exists():
        existing_aa = _seq_count(QUERY_AA_FASTA)
        if force or existing_aa == 0:
            with open(QUERY_AA_FASTA, "w") as fout:
                for rec in SeqIO.parse(str(PHOLD_AA_FASTA), "fasta"):
                    gene_id = _parse_3di_id(rec.id)
                    if gene_id not in hypo_locus_tags:
                        continue
                    fout.write(f">{gene_id}\n{str(rec.seq)}\n")
                    written_aa += 1
            log(f"  AA query FASTA:  {written_aa}/{len(hypo_locus_tags)} genes written"
                f" → {QUERY_AA_FASTA}")
        else:
            written_aa = existing_aa
            log(f"  AA query FASTA exists ({written_aa} seqs) — reuse")
    else:
        log(f"  phold_aa.fasta not found at {PHOLD_AA_FASTA} — AA query unavailable "
            f"(combined mode disabled)")

    return written_3di, written_aa


# =============================================================================
# STEP 2 — BUILD COMBINED QUERY DB (for --combined mode)
# =============================================================================

def _run_mmseqs(args: list, desc: str = "") -> subprocess.CompletedProcess:
    """Run an mmseqs2 sub-command and raise on non-zero exit."""
    cmd = ["mmseqs"] + [str(a) for a in args]
    log(f"  [{desc}] {' '.join(cmd[:6])}{'...' if len(cmd) > 6 else ''}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log(f"\n  STDERR:\n{result.stderr[-2000:]}")
        raise RuntimeError(f"mmseqs {desc} failed (exit {result.returncode})")
    return result


def _create_minimal_ca_stub(db_base: Path) -> None:
    """
    Foldseek 10.x's `convertalis` opens the queryDB's `_ca` (alpha-carbon
    coordinates) file at startup and exits with
        "No datafile could be found for ...queryDB_ca!"
    if it's missing -- there is NO flag to disable this check (verified by
    reading `foldseek convertalis --help`).

    `_ca` only ever gets written when createdb consumes real PDB/mmCIF
    input.  Our pipeline starts from phold's 3Di tokens (no structures),
    so neither `foldseek createdb --prostt5-model` nor `mmseqs createdb`
    produces a `_ca` file -- both leave the queryDB at the "Aminoacid"
    tier.  We work around the existence check by manufacturing a minimal
    stub:

      <db>_ca         : empty file (0 bytes -- no CA data)
      <db>_ca.dbtype  : copy of a target DB's _ca.dbtype, so convertalis
                        sees the correct foldseek "Calpha" magic value
      <db>_ca.index   : one row per query with offset=0 length=0, meaning
                        "no CA bytes for this key" -- foldseek reads zero
                        bytes if it ever looks one up

    This is safe iff we never request a CA-derived column in
    --format-output (lddt, qca, tca, rmsd, qtmscore, ttmscore, alntmscore,
    complex*tmscore).  Our LOCAL_FMT_CORE does not.
    """
    from config import FOLDSEEK_LOCAL_DBS

    ca_main   = Path(str(db_base) + "_ca")
    ca_dbtype = Path(str(db_base) + "_ca.dbtype")
    ca_index  = Path(str(db_base) + "_ca.index")
    aa_index  = Path(str(db_base) + ".index")

    # 1. Pick a reference _ca.dbtype from one of our configured target DBs.
    #    We just want the 4 magic bytes that say "this is a Calpha DB" --
    #    they should be identical across all foldseek structure DBs.
    ref_ca_dbtype = None
    for db_root in FOLDSEEK_LOCAL_DBS.values():
        candidate = Path(str(db_root) + "_ca.dbtype")
        if candidate.exists():
            ref_ca_dbtype = candidate
            break
    if ref_ca_dbtype is None:
        raise RuntimeError(
            "Cannot build queryDB_ca stub: no <db>_ca.dbtype found among "
            "FOLDSEEK_LOCAL_DBS targets.  Check FOLDSEEK_LOCAL_DBS paths "
            "in config.py."
        )

    # 2. Empty data file
    ca_main.touch()

    # 3. dbtype byte: byte-identical copy from the reference
    shutil.copy2(str(ref_ca_dbtype), str(ca_dbtype))
    try:
        b = ca_dbtype.read_bytes()
        dbtype_int = int.from_bytes(b, byteorder="little") if b else None
        log(f"    queryDB_ca.dbtype = {dbtype_int} (copied from "
            f"{ref_ca_dbtype.parent.name}/{ref_ca_dbtype.name})")
    except OSError:
        pass

    # 4. Index: one entry per query, all zero-length.  Keys must match
    #    queryDB.index so foldseek's hash lookups don't miss.
    if not aa_index.exists():
        raise RuntimeError(
            f"queryDB.index not found at {aa_index} -- cannot build _ca "
            f"index without knowing the query keys.  Did mmseqs createdb "
            f"fail silently?"
        )
    n_entries = 0
    with open(aa_index) as fin, open(ca_index, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            key = line.split("\t", 1)[0]
            fout.write(f"{key}\t0\t0\n")
            n_entries += 1
    log(f"    queryDB_ca.index = {n_entries} zero-length entries (one per query)")


def build_combined_query_db(foldseek_bin: str) -> Path:
    """
    Build a FoldSeek query database (queryDB / queryDB_ss / queryDB_h)
    using the **official ProstT5 -> foldseek recipe** documented by
    Michael Heinzinger (ProstT5 author) and endorsed by the foldseek
    maintainer in steineggerlab/foldseek#511:

        https://github.com/mheinzinger/ProstT5/blob/main/scripts/generate_foldseek_db.py

    The recipe -- 50 lines of canonical Python -- is:

      1. Parse the AA FASTA and 3Di FASTA into matching {id -> seq} dicts.
         IDs must align 1:1 (we drop any AA entry without a 3Di counterpart
         and vice versa, with warnings).
      2. Write three TSV files with integer keys 1..N (same key order in
         all three), formatted as:
                <key>\\t<sequence>\\n
      3. Run `foldseek tsv2db` three times:
                <db>     from aa.tsv     --output-dbtype 0  (AA, dbtype 0)
                <db>_ss  from 3di.tsv    --output-dbtype 0  (3Di as AA-format)
                <db>_h   from header.tsv --output-dbtype 12 (headers)
      4. Clean up the temp TSVs.

    No `_ca` file is created.  This works on FoldSeek 10.x as long as the
    --format-output requested in convertalis does NOT include any
    structure-derived column (lddt/lddtfull/qca/tca/t/u/qtmscore/ttmscore/
    alntmscore/rmsd/prob/complex*).  Our LOCAL_FMT_CORE deliberately omits
    all of these.

    Why this beats the mmseqs route we tried before:
      * It's the maintainer-endorsed path -- zero surprise factor.
      * tsv2db is foldseek-native, no cross-tool format risk.
      * Integer-key alignment is enforced by construction (we write key K
        for the same ID in all three TSVs).
      * No _ss-suffix manual move; tsv2db writes the right files directly.
      * No _ca stub needed (because we don't request CA-derived columns).
    """
    QUERY_DB_DIR.mkdir(parents=True, exist_ok=True)
    db_base = QUERY_DB_DIR / "queryDB"

    # Wipe any stale DB files from prior runs (incl. zero-byte stubs from
    # previous foldseek-createdb / mmseqs attempts).
    for f in list(QUERY_DB_DIR.glob("queryDB*")):
        try:
            if f.is_dir():
                shutil.rmtree(f, ignore_errors=True)
            else:
                f.unlink()
        except OSError:
            pass

    log(f"  Building queryDB via foldseek tsv2db (official ProstT5 recipe)...")

    # ---- 1. Parse AA + 3Di FASTAs into matching dicts ----
    if not QUERY_AA_FASTA.exists() or QUERY_AA_FASTA.stat().st_size == 0:
        raise RuntimeError(
            f"AA FASTA not found at {QUERY_AA_FASTA}.  The official "
            f"tsv2db recipe requires real AA sequences (not a placeholder)."
        )
    if not QUERY_3DI_FASTA.exists() or QUERY_3DI_FASTA.stat().st_size == 0:
        raise RuntimeError(
            f"3Di FASTA not found at {QUERY_3DI_FASTA}."
        )

    sequences_aa: dict = {}
    for rec in SeqIO.parse(str(QUERY_AA_FASTA), "fasta"):
        sequences_aa[rec.id] = str(rec.seq)

    sequences_3di: dict = {}
    extra_3di = 0
    for rec in SeqIO.parse(str(QUERY_3DI_FASTA), "fasta"):
        if rec.id not in sequences_aa:
            extra_3di += 1
            continue
        sequences_3di[rec.id] = str(rec.seq).upper()

    if extra_3di:
        log(f"    WARNING: {extra_3di} 3Di entries had no AA counterpart -- skipped")

    # Drop AA entries that have no matching 3Di -- keeps dicts strictly aligned.
    missing_3di = [k for k in list(sequences_aa.keys()) if k not in sequences_3di]
    if missing_3di:
        log(f"    WARNING: {len(missing_3di)} AA entries had no 3Di counterpart "
            f"-- dropping them from the query DB:")
        for g in missing_3di[:5]:
            log(f"      {g}")
        if len(missing_3di) > 5:
            log(f"      ... and {len(missing_3di) - 5} more")
        for k in missing_3di:
            del sequences_aa[k]

    n_matched = len(sequences_aa)
    log(f"    AA+3Di matched entries: {n_matched}")
    if n_matched == 0:
        raise RuntimeError(
            "No AA/3Di pairs to write to the query DB.  Check that "
            "phold_aa.fasta and phold_3di.fasta share the same IDs after "
            "prepare_query_fastas's _parse_3di_id step."
        )

    # ---- 2. Write the three TSV files with integer keys 1..N ----
    tsv_dir = QUERY_DB_DIR / "_tsv_tmp"
    tsv_dir.mkdir(parents=True, exist_ok=True)
    aa_tsv  = tsv_dir / "aa.tsv"
    ss_tsv  = tsv_dir / "3di.tsv"
    h_tsv   = tsv_dir / "header.tsv"

    with open(aa_tsv, "w") as f_aa, \
         open(ss_tsv, "w") as f_ss, \
         open(h_tsv,  "w") as f_h:
        for i, gene_id in enumerate(sequences_aa.keys()):
            key = i + 1   # 1-indexed, matches the official script
            f_aa.write(f"{key}\t{sequences_aa[gene_id]}\n")
            f_ss.write(f"{key}\t{sequences_3di[gene_id]}\n")
            f_h.write (f"{key}\t{gene_id}\n")

    log(f"    Wrote {n_matched}-row TSVs to {tsv_dir}/")

    # ---- 3. tsv2db x3 ----
    # dbtype 0  = AA (also used for 3Di stored as AA-format -- foldseek
    #             identifies 3Di by the _ss filename suffix)
    # dbtype 12 = generic DB (used for headers)
    _run(foldseek_bin,
         ["tsv2db", str(aa_tsv), str(db_base),
          "--output-dbtype", "0"],
         desc="tsv2db (AA)")
    _run(foldseek_bin,
         ["tsv2db", str(ss_tsv), str(db_base) + "_ss",
          "--output-dbtype", "0"],
         desc="tsv2db (3Di)")
    _run(foldseek_bin,
         ["tsv2db", str(h_tsv), str(db_base) + "_h",
          "--output-dbtype", "12"],
         desc="tsv2db (headers)")

    # ---- 4. Clean up temp TSVs ----
    shutil.rmtree(tsv_dir, ignore_errors=True)

    # ---- 5. Sanity checks ----
    expected = ["", ".dbtype", ".index",
                "_ss", "_ss.dbtype", "_ss.index",
                "_h",  "_h.dbtype",  "_h.index"]
    for suffix in expected:
        f = Path(str(db_base) + suffix)
        if not f.exists():
            raise RuntimeError(f"Missing expected DB file: {f}")

    # Key alignment cross-check: queryDB.index and queryDB_ss.index must
    # have the same keys in the same order (they do by construction here,
    # but verify in case tsv2db re-orders).
    try:
        aa_keys = [line.split("\t", 1)[0]
                   for line in (QUERY_DB_DIR / "queryDB.index").read_text().splitlines()
                   if line]
        ss_keys = [line.split("\t", 1)[0]
                   for line in (QUERY_DB_DIR / "queryDB_ss.index").read_text().splitlines()
                   if line]
        log(f"    queryDB entries: {len(aa_keys)} AA / {len(ss_keys)} 3Di")
        if aa_keys != ss_keys:
            raise RuntimeError(
                "queryDB and queryDB_ss have different keys -- tsv2db "
                "should have produced identical key sequences.  Aborting "
                "to avoid misaligned search results."
            )
        log(f"    queryDB / queryDB_ss key alignment: OK ({len(aa_keys)} entries)")
    except OSError as e:
        log(f"    Could not verify key alignment: {e}")

    log("  Query DB ready (no _ca by design -- format-output is sequence-only).")
    return db_base


# =============================================================================
# STEP 3 — RUN FOLDSEEK SEARCH
# =============================================================================

def _run(foldseek_bin: str, args: list, desc: str = "") -> subprocess.CompletedProcess:
    """Run a foldseek sub-command and raise on non-zero exit."""
    cmd = [foldseek_bin] + [str(a) for a in args]
    log(f"  [{desc}] {' '.join(cmd[:6])}{'...' if len(cmd) > 6 else ''}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log(f"\n  STDERR:\n{result.stderr[-2000:]}")
        raise RuntimeError(f"foldseek {desc} failed (exit {result.returncode})")
    return result


def _detect_prostt5_support(foldseek_bin: str) -> bool:
    """
    Return True if the installed foldseek binary supports --prostt5-input.
    Checks by running `foldseek easy-search --help` and grepping for the flag.
    """
    try:
        result = subprocess.run(
            [foldseek_bin, "easy-search", "--help"],
            capture_output=True, text=True, timeout=10
        )
        return "--prostt5-input" in (result.stdout + result.stderr)
    except Exception:
        return False


def run_foldseek_easy_search(
    foldseek_bin: str,
    db_name: str,
    db_path: Path,
    out_m8: Path,
    combined_query_db: Path,
    threads: int,
    evalue: float = 1.0,
    max_seqs: int = 300,
) -> bool:
    """
    Run a foldseek search against one target database.

    Query mode selection (auto-detected):
      1. If combined_query_db is provided (the normal path, built via
         `mmseqs createdb` in build_combined_query_db): use
         `foldseek search` + `convertalis` with alignment-type 2 (3Di+AA)
         when real AA is available, else alignment-type 1 (3Di-only).
      2. Fallback only: if combined_query_db is None AND foldseek supports
         `--prostt5-input` (older builds), use easy-search with it.
         On the pharokka-env foldseek 10.x build this fallback cannot fire
         -- the binary has neither --prostt5-input nor --seq-type, so the
         createdb route must succeed.

    Returns True on success, False on failure.
    """
    if out_m8.exists() and out_m8.stat().st_size > 0:
        n = sum(1 for l in out_m8.read_text().splitlines()
                if l and not l.startswith("#"))
        if n > 0:
            log(f"  {db_name}: result file exists ({n} hits) — skip")
            return True
        # zero-row cache from a previous failed/empty run: remove and re-search
        log(f"  {db_name}: cached m8 has 0 hits — re-running search")
        try:
            out_m8.unlink()
        except OSError:
            pass

    if not db_path.exists() and not Path(str(db_path) + ".dbtype").exists():
        log(f"  {db_name}: database not found at {db_path} — SKIP")
        log(f"    Download with: foldseek databases {db_name} {db_path} /tmp/foldseek_tmp")
        return False

    log(f"\n  Searching against {db_name} ({db_path})...")

    # Per-DB format-output: only include taxid,taxname if the target DB has
    # a _taxonomy file.  Otherwise convertalis falls back to a global NCBI
    # taxdump and errors "names.dmp, nodes.dmp, merged.dmp ... could not be
    # found!".  We pad the resulting m8 to a uniform column count afterwards.
    has_tax = _db_has_taxonomy(db_path)
    fmt = LOCAL_FMT_CORE + ("," + LOCAL_FMT_TAX if has_tax else "")
    n_target_cols = len(LOCAL_COLS_CORE) + len(LOCAL_COLS_TAX)

    with tempfile.TemporaryDirectory(prefix="foldseek_3di_") as tmpdir:
        tmp = Path(tmpdir)

        if combined_query_db is not None:
            # ---------- createdb route: search + convertalis ----------
            # alignment-type 2 = 3Di+AA (when real AA FASTA was used)
            # alignment-type 1 = 3di-only  (when AA placeholder was used)
            aln_type = "2" if QUERY_AA_FASTA.exists() else "1"
            mode_desc = "3Di+AA" if aln_type == "2" else "3di-only"
            tax_desc  = "+tax" if has_tax else "no-tax"
            result_db = tmp / "resultDB"

            # ---- Pre-search DB file audit --------------------------------
            # Log which DB files exist, their sizes, and whether the
            # precomputed index is available (index avoids in-RAM k-mer table).
            db_suffixes = [
                ("",               "AA sequences      "),
                ("_ss",            "3Di sequences     "),
                ("_h",             "Headers/descriptions"),
                ("_taxonomy",      "Taxonomy          "),
                ("_ca",            "Cα coords (delete!)"),
                ("_seq_ca.1",      "Member Cα (delete!)"),
                (".idx",           "AA precomp. index "),
                ("_ss.idx",        "3Di precomp. index"),  # ← FoldSeek needs this for 3Di+AA
                ("_ss.idx.0",      "3Di index chunk0  "),  # indicates incomplete index if present w/o _ss.idx
            ]
            log(f"  DB file audit for {db_name}:")
            # FoldSeek needs BOTH afdb50.idx (AA) AND afdb50_ss.idx (3Di) for 3Di+AA search
            has_aa_idx  = Path(str(db_path) + ".idx").exists()
            has_3di_idx = Path(str(db_path) + "_ss.idx").exists()
            has_incomplete_3di = (Path(str(db_path) + "_ss.idx.0").exists()
                                  and not has_3di_idx)
            for sfx, label in db_suffixes:
                fpath = Path(str(db_path) + sfx)
                if fpath.exists():
                    size_gb = fpath.stat().st_size / 1e9
                    log(f"    {label}: {fpath.name}  ({size_gb:.1f} GB)")
                else:
                    log(f"    {label}: {fpath.name}  [not found]")
            has_precomputed_idx = has_aa_idx and has_3di_idx
            if has_incomplete_3di:
                log(f"  ✗ INCOMPLETE 3Di index detected (afdb50_ss.idx.0 exists but not afdb50_ss.idx)")
                log(f"    Fix: rm {db_path}_ss.idx.0 {db_path}_ss.idx.index.0")
                log(f"         foldseek createindex {db_path} tmp --index-exclude 2 --threads 16")
            if has_precomputed_idx:
                log(f"  ✓ Precomputed index found (AA + 3Di) — low RAM search mode")
            elif has_aa_idx and not has_3di_idx:
                log(f"  ✗ AA index present but 3Di index missing — search will OOM on large DBs")
                log(f"    Fix: foldseek createindex {db_path} tmp --index-exclude 2 --threads 16")
            elif not has_aa_idx:
                log(f"  ✗ No precomputed index — k-mer table built in RAM "
                    f"(may OOM for large DBs; run: foldseek createindex {db_path} tmp --index-exclude 2)")

            try:
                # Check if --sort-by-structure-bits is supported by this foldseek build.
                # It reduces RAM from ~151 GB → ~35 GB for large DBs (afdb50).
                # Requires a precomputed index built with: createindex --index-exclude 2
                # Not supported in all builds — detect via help output.
                _sbsb_supported = "--sort-by-structure-bits" in subprocess.run(
                    [foldseek_bin, "search", "--help"],
                    capture_output=True, text=True
                ).stdout + subprocess.run(
                    [foldseek_bin, "search", "--help"],
                    capture_output=True, text=True
                ).stderr

                search_args = [
                    "search",
                    str(combined_query_db), str(db_path), str(result_db), str(tmp),
                    "--alignment-type", aln_type,
                    "-e", str(evalue),
                    "--max-seqs", str(max_seqs),
                    "--threads", str(threads),
                    # NOTE: do NOT add --db-load-mode 2 (mmap) here. With <index-size
                    # RAM the 196GB afdb50 index thrashes from disk (~6min -> >1h).
                    # Keep default load-to-RAM + --mem 256G in steps/02_foldseek_3di.sh.
                ]
                if _sbsb_supported:
                    search_args += ["--sort-by-structure-bits", "0"]
                    log(f"  --sort-by-structure-bits 0 supported → Cα not loaded into RAM")
                else:
                    log(f"  --sort-by-structure-bits not in this build; searching without it")
                    log(f"  (if OOM: ensure createindex was run with --index-exclude 2)")

                log(f"  Full search command: foldseek {' '.join(search_args)}")
                _run(foldseek_bin, search_args, desc=f"search ({db_name}, {mode_desc})")

                _run(foldseek_bin, [
                    "convertalis",
                    str(combined_query_db), str(db_path), str(result_db), str(out_m8),
                    "--format-output", fmt,
                ], desc=f"convertalis ({db_name}, {tax_desc})")

                # Normalize column count: pad with empty taxid/taxname if this
                # DB has no taxonomy, so merged per-gene m8s have uniform layout.
                if not has_tax:
                    _pad_m8_taxonomy(out_m8, n_target_cols)
                    log(f"    Padded {out_m8.name} with empty taxid/taxname "
                        f"({n_target_cols} cols total)")

                n_hits = sum(1 for l in out_m8.read_text().splitlines()
                             if l and not l.startswith("#"))
                log(f"  {db_name}: {mode_desc} {tax_desc} done -> {n_hits} hits -> {out_m8.name}")
                return True

            except RuntimeError as e:
                log(f"  {db_name}: createdb route FAILED: {e}")
                log(f"  Falling back to easy-search (may fail on foldseek 10.x)...")
                # Fall through to easy-search approach

        # ---------- easy-search fallback (only viable on ProstT5 builds) ----------
        # On FoldSeek 10.x without ProstT5 (e.g. the pharokka micromamba env),
        # easy-search has NO usable 3Di entry point: --prostt5-input is absent
        # and --seq-type was removed.  The createdb route above is the only
        # path that works, so reaching this block is itself a problem -- we
        # log it loudly.
        has_prostt5 = _detect_prostt5_support(foldseek_bin)
        query_fasta = str(QUERY_3DI_FASTA)

        if not has_prostt5:
            log(f"  {db_name}: createdb route failed AND this foldseek build "
                f"has no --prostt5-input.  No viable easy-search path on "
                f"FoldSeek 10.x without ProstT5.  Skipping {db_name}.")
            return False

        log(f"    Mode: --prostt5-input (easy-search fallback)")
        try:
            _run(foldseek_bin, [
                "easy-search",
                query_fasta, str(db_path), str(out_m8), str(tmp),
                "-e", str(evalue),
                "--max-seqs", str(max_seqs),
                "--threads", str(threads),
                "--format-output", fmt,
                "--prostt5-input",
            ], desc=f"easy-search ({db_name})")

            n = sum(1 for l in out_m8.read_text().splitlines()
                    if l and not l.startswith("#"))
            log(f"  {db_name}: done -> {n} hits -> {out_m8.name}")
            return True

        except RuntimeError as e:
            log(f"  {db_name}: easy-search --prostt5-input FAILED: {e}")
            return False


# =============================================================================
# STEP 4 — SPLIT COMBINED M8 INTO PER-GENE FILES
# =============================================================================

def split_and_merge_hits(per_db_m8_files: list, raw_hits_dir: Path,
                          gene_ids_in_scope: set, force: bool = False):
    """
    Read all per-DB m8 files, deduplicate across databases, and write one
    .m8 file per gene to raw_hits_dir/.

    Deduplication: if the same target accession appears in multiple DBs,
    keep the hit with the better (lower) evalue.
    """
    raw_hits_dir.mkdir(parents=True, exist_ok=True)

    if not per_db_m8_files:
        log("  No per-DB m8 files to merge.")
        return

    # Read all DB results into a single DataFrame
    dfs = []
    for m8_path in per_db_m8_files:
        if not m8_path.exists() or m8_path.stat().st_size == 0:
            continue
        content = m8_path.read_text()
        lines   = [l for l in content.splitlines() if l and not l.startswith("#")]
        if not lines:
            continue
        rows = []
        # Pad each row to the DATA column count (LOCAL_COLS = 16).  The
        # source_db tag is added as a separate column below, NOT padded in.
        # Using N_LOCAL_COLS here would mismatch the DataFrame columns and
        # raise ValueError.
        n_data_cols = len(LOCAL_COLS)
        for line in lines:
            parts = line.split("\t")
            parts = parts[:n_data_cols] + [""] * (n_data_cols - len(parts))
            rows.append(parts)
        df = pd.DataFrame(rows, columns=LOCAL_COLS)
        df["_source_db"] = m8_path.stem  # track which DB the hit came from
        dfs.append(df)

    if not dfs:
        log("  No hits found in any per-DB m8 file.")
        return

    combined = pd.concat(dfs, ignore_index=True)

    # Restrict to genes in scope
    combined = combined[combined["query"].isin(gene_ids_in_scope)].copy()
    log(f"  Total hits across all DBs: {len(combined)} "
        f"(genes: {combined['query'].nunique()})")

    # Deduplicate: same query + same accession → keep lowest evalue
    # Important: keep _source_db so the per-gene m8 file records which
    # FOLDSEEK_LOCAL_DBS entry each hit came from.
    combined["accession"] = combined["theader"].str.split(" ", n=1).str[0]
    combined["evalue_num"] = pd.to_numeric(combined["evalue"], errors="coerce")
    combined = (
        combined
        .sort_values(["query", "accession", "evalue_num"])
        .drop_duplicates(subset=["query", "accession"], keep="first")
        .drop(columns=["accession", "evalue_num"])
    )
    # Rename to the final 17th column name
    combined = combined.rename(columns={"_source_db": "source_db"})
    log(f"  After deduplication: {len(combined)} hits")

    # Write per-gene files (same filename convention as 02b: {gene_id}.m8)
    # 17-column TSV: 16 data m8 columns + source_db
    n_written = 0
    for gene_id, group in combined.groupby("query"):
        out_path = raw_hits_dir / f"{gene_id}.m8"
        if force or not out_path.exists():
            group[LOCAL_COLS_OUT].to_csv(
                str(out_path), sep="\t", index=False, header=False
            )
            n_written += 1

    log(f"  Per-gene m8 files written: {n_written} → {raw_hits_dir}/")


# =============================================================================
# STEP 5 — PARSE LOCAL M8 FILES
# =============================================================================

def parse_local_m8_file(m8_path: Path, gene_id: str) -> pd.DataFrame:
    """
    Parse a single local FoldSeek .m8 file into a DataFrame that matches
    the schema expected by build_best_and_top3() (from 02b).

    Local format differences vs web API (M8_COLS in utils.py):
      - 'fident' is fraction 0-1 (web API 'pident' is percentage 0-100)
      - 'qcov'   is fraction 0-1 (web API 'qcov_aa' is absolute residue count)
      - No 'lddt' column (dropped -- requires queryDB_ca which mmseqs2's
        createdb doesn't write; foldseek_scoring.py never reads it).
      - Only 16 data columns + source_db tag = 17 (web API has 21).
    """
    content = m8_path.read_text().strip()
    if not content:
        return pd.DataFrame()
    lines = [l for l in content.splitlines() if l and not l.startswith("#")]
    if not lines:
        return pd.DataFrame()

    rows = []
    for line in lines:
        parts = line.split("\t")
        parts = parts[:N_LOCAL_COLS] + [""] * (N_LOCAL_COLS - len(parts))
        rows.append(parts)

    df = pd.DataFrame(rows, columns=LOCAL_COLS_OUT)
    df["gene"] = gene_id

    # Derive accession + description from theader ("ACCESSION DESCRIPTION")
    df["accession"]   = df["theader"].str.split(" ", n=1).str[0]
    df["description"] = (
        df["theader"].str.split(" ", n=1).str[1].fillna("").str.strip()
    )

    # Source-DB tag: normalized name of the FOLDSEEK_LOCAL_DBS entry the hit
    # came from (afdb-swissprot / pdb100 / baktfold-afdb).
    df["foldseek_subdb"] = df.apply(
        lambda r: _classify_subdb(r.get("source_db", ""), r["accession"]),
        axis=1,
    )

    # Numeric conversions (lddt and prob no longer present in format-output)
    for col in ["fident", "evalue", "bits", "qcov", "alnlen"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Rename to match internal schema
    df = df.rename(columns={
        "fident": "pident",       # 0-1 fraction; display-only, not used in filtering
        "evalue": "evalue",
        "bits":   "score",
        "qcov":   "qcov_frac",    # 0-1 fraction; corrected to AA count in main()
    })
    # lddtfull no longer derivable from format-output; keep the column for
    # schema stability (downstream code doesn't read it but won't break if
    # it's there as NaN).
    df["lddtfull"]   = float("nan")
    df["evalue_raw"] = df["evalue"]
    df["score_raw"]  = df["score"]
    # qcov_aa placeholder — will be overwritten in main() once aa_length is known
    df["qcov_aa"]    = df["qcov_frac"]

    return df


def parse_all_local_hits(raw_hits_dir: Path, gene_ids_in_scope: set,
                          hypo_meta: pd.DataFrame = None) -> pd.DataFrame:
    """
    Parse all per-gene m8 files in raw_hits_dir and return a combined
    DataFrame, identical in schema to the output of parse_all_hits() in 02b.

    hypo_meta is used to convert qcov_frac → qcov_aa (AA count) when
    aa_length is available.
    """
    dfs = []
    for m8_path in sorted(raw_hits_dir.glob("*.m8")):
        gene_id = m8_path.stem
        if gene_ids_in_scope and gene_id not in gene_ids_in_scope:
            continue
        df = parse_local_m8_file(m8_path, gene_id)
        if not df.empty:
            dfs.append(df)

    if not dfs:
        log("  No local hits to parse.")
        return pd.DataFrame()

    hits = pd.concat(dfs, ignore_index=True)

    # Convert qcov_frac → qcov_aa (residue count) using aa_length
    if hypo_meta is not None and "aa_length" in hypo_meta.columns:
        aa_map = hypo_meta.set_index("locus_tag")["aa_length"].to_dict()
        hits["aa_length"] = hits["gene"].map(aa_map)
        hits["qcov_aa"] = (
            hits["qcov_frac"] * hits["aa_length"].fillna(0)
        ).round().astype("Int64")
    else:
        hits["qcov_aa"] = hits["qcov_frac"]  # keep as fraction if length unavailable

    # ---- Per-DB taxon filter (FOLDSEEK_TAXON_FILTER) -------------------------
    # Applied before quality scoring.  For DBs with embedded taxonomy (e.g.
    # afdb50) and a non-empty filter string, drop hits whose taxname indicates
    # a eukaryote (or synthetic construct).  This removes false positive
    # structural matches to eukaryotic folds early — before they influence
    # top3 selection or composite scoring.
    # Note: we use keyword-based eukaryote detection (fast, no taxdump needed)
    # rather than exact taxid matching, which would require a full lineage lookup.
    if FOLDSEEK_TAXON_FILTER:
        before = len(hits)
        keep_mask = pd.Series(True, index=hits.index)
        for db_name, filter_str in FOLDSEEK_TAXON_FILTER.items():
            if not filter_str:
                continue  # no filter configured for this DB
            is_this_db = hits.get("foldseek_subdb", pd.Series("", index=hits.index)) == db_name
            is_euka    = hits["taxname"].apply(_is_eukaryotic_taxname)
            # Drop hits from this DB that are eukaryotic
            keep_mask = keep_mask & ~(is_this_db & is_euka)
        hits = hits[keep_mask].copy()
        n_dropped = before - len(hits)
        if n_dropped:
            log(f"  Taxon filter: dropped {n_dropped} eukaryotic hits "
                f"({', '.join(k for k,v in FOLDSEEK_TAXON_FILTER.items() if v)} DBs)")

    # Quality filter (mirrors parse_all_hits in 02b)
    quality_ok = (
        (hits["evalue"].fillna(999.0) < FOLDSEEK_EVALUE_MAX) |
        (hits["score"].fillna(0.0)    >= FOLDSEEK_SCORE_OVERRIDE)
    )
    hits["informative"]           = hits["description"].apply(_is_informative_fs) & quality_ok
    hits["phage_boost"]           = hits["description"].apply(_phage_boost_factor)
    hits["composite_score"]       = hits["score"].fillna(0.0) * hits["phage_boost"]
    hits["same_host"]             = hits["taxname"].apply(_is_same_host_hit)
    hits["defense_flag"]          = hits["description"].apply(
        lambda d: bool(_DEFENSE_PATTERN.search(d)) if isinstance(d, str) else False
    )
    # Pre-flag promiscuous folds and eukaryotic descriptions so build_best_and_top3
    # can propagate these flags to the best-hit output (→ step 04 review routing).
    hits["promiscuous_fold_flag"] = hits["description"].apply(
        lambda d: _is_promiscuous_fold_hit(str(d) if d else "")
    )
    hits["eukaryotic_desc_flag"]  = hits["description"].apply(
        lambda d: _has_eukaryotic_description(str(d) if d else "")
    )

    hits = hits.sort_values(
        ["gene", "composite_score", "evalue"],
        ascending=[True, False, True]
    )
    return hits


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Local FoldSeek annotation via ProstT5 3di tokens"
    )
    parser.add_argument(
        "--parse_only", action="store_true",
        help="Skip FoldSeek search; only parse existing raw hits"
    )
    parser.add_argument(
        "--combined", action="store_true",
        help="(legacy / no-op) — combined 3Di+AA is now the default. Kept for "
             "backward compatibility with old submit_all.sh scripts."
    )
    parser.add_argument(
        "--easy-search", action="store_true",
        help="Bypass the createdb route and use easy-search with --seq-type 1 "
             "(only works on older foldseek builds)."
    )
    parser.add_argument(
        "--threads", type=int, default=FOLDSEEK_THREADS,
        help=f"CPU threads for foldseek (default: {FOLDSEEK_THREADS})"
    )
    parser.add_argument(
        "--evalue", type=float, default=FOLDSEEK_EVALUE_MAX,
        help=f"E-value cutoff for foldseek search (default: {FOLDSEEK_EVALUE_MAX} from config). "
             "Hits above this evalue are not returned at all, reducing noise and search RAM. "
             "The Python quality gate (FOLDSEEK_EVALUE_MAX) applies the same threshold "
             "post-hoc, so changing this only affects what foldseek itself fetches."
    )
    parser.add_argument(
        "--max_seqs", type=int, default=300,
        help="Max hits to keep per query (default: 300)"
    )
    parser.add_argument(
        "--dbs", nargs="+", default=None,
        help="Restrict to specific database names (default: all in FOLDSEEK_LOCAL_DBS)"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run all steps even if output files exist"
    )
    args = parser.parse_args()

    section("STEP 02d — LOCAL FOLDSEEK WITH PROSTT5 3di TOKENS")

    # Create output dirs
    for d in [FOLDSEEK_3DI_DIR, PER_DB_HITS_DIR, RAW_HITS_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    # ---- Check foldseek binary ----
    foldseek_bin = FOLDSEEK_CMD
    if not args.parse_only:
        if not shutil.which(foldseek_bin):
            log(f"ERROR: '{foldseek_bin}' not found on PATH.")
            log("  Install foldseek: https://github.com/steineggerlab/foldseek#installation")
            log("  Or set FOLDSEEK_CMD in config.py to the full path.")
            sys.exit(1)
        try:
            result = subprocess.run([foldseek_bin, "version"],
                                    capture_output=True, text=True, timeout=5)
            log(f"FoldSeek version: {result.stdout.strip() or result.stderr.strip()}")
        except Exception:
            pass

    # ---- Load hypothetical gene list ----
    # Auto-built from pharokka TSVs (pharokka-hypothetical rule, 313 targets)
    # if missing, stale, or inconsistent with gene_metadata.csv.
    #
    # Idempotency: regenerate if:
    #   (a) file is missing or unreadable / empty
    #   (b) gene_metadata.csv is absent (single source of truth)
    #   (c) count in hypothetical_genes.csv != hypothetical count in gene_metadata
    #       (catches old phold-hypothetical runs that produced 244 instead of 313)
    #   (d) 'prophage' column is absent (old schema, pre-cluster-enrichment)
    #   (e) --force flag
    from config import (PHOLD_COMBINED_TSV, PHOLD_AA_FASTA, UNINFORMATIVE_STRINGS,
                        PHAROKKA_CDS_TSV_GLOB, GENE_METADATA_CSV)

    _probe = None
    needs_regen = (not HYPO_GENE_LIST.exists()) or args.force
    if HYPO_GENE_LIST.exists() and not needs_regen:
        try:
            _probe = pd.read_csv(str(HYPO_GENE_LIST))
            if len(_probe) == 0:
                log(f"{HYPO_GENE_LIST.name} exists but is empty -- regenerating.")
                needs_regen = True
            elif "prophage" not in _probe.columns:
                log(f"{HYPO_GENE_LIST.name} missing 'prophage' column (old schema) -- regenerating.")
                needs_regen = True
        except Exception:
            log(f"{HYPO_GENE_LIST.name} unreadable -- regenerating.")
            needs_regen = True

    # Check (b): gene_metadata.csv must exist
    needs_regen = needs_regen or (not GENE_METADATA_CSV.exists())

    # Check (c): count consistency with gene_metadata
    if not needs_regen and _probe is not None and GENE_METADATA_CSV.exists():
        try:
            _meta = pd.read_csv(str(GENE_METADATA_CSV))
            _expected = int(_meta["is_hypothetical"].sum())
            if len(_probe) != _expected:
                log(f"{HYPO_GENE_LIST.name} has {len(_probe)} entries but "
                    f"gene_metadata has {_expected} hypotheticals -- regenerating.")
                needs_regen = True
        except Exception:
            pass

    if needs_regen:
        log(f"Auto-generating gene_metadata.csv + {HYPO_GENE_LIST.name} ...")
        HYPO_GENE_LIST.parent.mkdir(parents=True, exist_ok=True)
        GENE_METADATA_CSV.parent.mkdir(parents=True, exist_ok=True)

        if not PHOLD_AA_FASTA.exists():
            log(f"ERROR: {PHOLD_AA_FASTA} not found. Run steps/01d_merge_3di.sh first.")
            sys.exit(1)

        # phold --prefix <NAME> writes headers like '>CCIE1:LOCUS_CDS_0001';
        # strip the prophage prefix so the key is the bare locus_tag (matches
        # both phold's cds_id and pharokka's gene column).
        def _strip_prefix(h):
            return h.split(":", 1)[1] if ":" in h else h

        aa_lengths = {}
        cur, n = None, 0
        with open(str(PHOLD_AA_FASTA)) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if cur is not None:
                        aa_lengths[cur] = n
                    cur = _strip_prefix(line[1:].split()[0]); n = 0
                else:
                    n += len(line)
            if cur is not None:
                aa_lengths[cur] = n

        # ---- Build the full 507-row gene_metadata.csv from PHAROKKA output ----
        # IMPORTANT (design decision 2026-06-01): the hypothetical target set is
        # defined by the ORIGINAL PHAROKKA product (annot column), NOT phold's
        # post-rescue product.  This (a) matches the local pipeline's 313-gene
        # target set, (b) recovers the merged enrichment for genes phold
        # annotated, and (c) is the rule the future mode-1 (new prophages) and
        # mode-2 (single proteins) runs will use: always FoldSeek the
        # pharokka-hypotheticals.
        import glob as _glob
        def _is_hypo(p):
            s = str(p).strip().lower() if p is not None else ""
            return (s == "" or s in UNINFORMATIVE_STRINGS)

        pharokka_tsvs = sorted(_glob.glob(PHAROKKA_CDS_TSV_GLOB))
        meta_rows = []
        if pharokka_tsvs:
            for tsv in pharokka_tsvs:
                proph = Path(tsv).parent.name
                pdf = pd.read_csv(tsv, sep="\t", low_memory=False)
                gcol = "gene" if "gene" in pdf.columns else pdf.columns[0]
                acol = "annot" if "annot" in pdf.columns else None
                for _, r in pdf.iterrows():
                    lt = _strip_prefix(str(r[gcol]))
                    annot = str(r[acol]).strip() if acol else ""
                    hyp = _is_hypo(annot)
                    meta_rows.append(dict(
                        prophage=proph, locus_tag=lt,
                        function=("unknown function" if hyp else annot),
                        product=(annot if annot else "hypothetical protein"),
                        aa_length=int(aa_lengths.get(lt, 0)),
                        is_hypothetical=bool(hyp),
                        start=r.get("start"), end=r.get("stop", r.get("end")),
                        strand=r.get("strand")))
            meta_df = pd.DataFrame(meta_rows).sort_values(["prophage", "locus_tag"])
            src = "pharokka annot"
        else:
            # Fallback: no pharokka tsvs (e.g. mode-2 single proteins) -> define
            # hypothetical from phold product, still emitting the full schema.
            log("  WARNING: no pharokka *_cds_final_merged_output.tsv found "
                f"(glob: {PHAROKKA_CDS_TSV_GLOB}); falling back to phold product.")
            phold_df = pd.read_csv(str(PHOLD_COMBINED_TSV), sep="\t", low_memory=False)
            prod_col = "product" if "product" in phold_df.columns else "phold_product"
            id_col = "cds_id" if "cds_id" in phold_df.columns else "locus_tag"
            pcol = "prophage" if "prophage" in phold_df.columns else None
            rows = []
            for _, r in phold_df.iterrows():
                raw_id = str(r[id_col])
                lt = _strip_prefix(raw_id)
                # In protein mode, step 01p prefixes FASTA headers as
                # "phiNP:GENE_ID" so the prophage is recoverable here.
                # _strip_prefix already splits on ":" to give bare locus_tag.
                proph = raw_id.split(":", 1)[0] if ":" in raw_id else \
                        (str(r[pcol]) if pcol else "NA")
                prod = str(r[prod_col]).strip()
                rows.append(dict(
                    prophage=proph, locus_tag=lt,
                    function="unknown function",
                    product=(prod if prod else "hypothetical protein"),
                    aa_length=int(aa_lengths.get(lt, 0)),
                    is_hypothetical=True,  # protein mode: ALL genes are hypothetical
                    start=r.get("start"), end=r.get("end"), strand=r.get("strand")))
            meta_df = pd.DataFrame(rows).sort_values(["prophage", "locus_tag"])
            src = "phold product (fallback)"

        meta_df.to_csv(str(GENE_METADATA_CSV), index=False)
        log(f"  Wrote {GENE_METADATA_CSV.name}: {len(meta_df)} CDS "
            f"({int(meta_df['is_hypothetical'].sum())} hypothetical, source: {src})")

        # ---- Derive the FoldSeek target list = hypotheticals only ----
        out_df = (meta_df[meta_df["is_hypothetical"]]
                  [["prophage", "locus_tag", "aa_length"]].drop_duplicates())
        out_df.to_csv(str(HYPO_GENE_LIST), index=False)
        log(f"  Wrote {HYPO_GENE_LIST.name}: {len(out_df)} hypothetical target genes")

        if len(out_df) == 0:
            log("ERROR: 0 hypothetical proteins after auto-generation.")
            sys.exit(1)

    hypo_meta = pd.read_csv(str(HYPO_GENE_LIST))
    all_hypo_genes     = list(hypo_meta["locus_tag"])
    hypo_locus_tag_set = set(all_hypo_genes)
    log(f"Target genes: {len(all_hypo_genes)} hypothetical proteins "
        f"(pharokka-hypothetical rule)")

    if len(all_hypo_genes) == 0:
        log(f"ERROR: {HYPO_GENE_LIST} contains no entries.")
        log("  Re-run with --force to rebuild the gene list.")
        sys.exit(1)

    # ---- Select databases ----
    if args.dbs:
        dbs_to_search = {k: v for k, v in FOLDSEEK_LOCAL_DBS.items() if k in args.dbs}
        missing_keys = set(args.dbs) - set(FOLDSEEK_LOCAL_DBS.keys())
        if missing_keys:
            log(f"WARNING: --dbs specified unknown DB name(s): {missing_keys}")
            log(f"  Available: {list(FOLDSEEK_LOCAL_DBS.keys())}")
    else:
        dbs_to_search = FOLDSEEK_LOCAL_DBS

    if not dbs_to_search and not args.parse_only:
        log("ERROR: No local databases configured in FOLDSEEK_LOCAL_DBS (config.py).")
        log("  Download databases first, e.g.:")
        log("    foldseek databases afdb-swissprot ~/foldseek_dbs/afdb_swissprot_db /tmp/tmp")
        sys.exit(1)

    log(f"Databases to search: {list(dbs_to_search.keys())}")

    # ---- STEP 1: Prepare query FASTAs ----
    if not args.parse_only:
        section("STEP 1 — Prepare query FASTAs")
        n_3di, n_aa = prepare_query_fastas(hypo_locus_tag_set, force=args.force)
        log(f"  3di sequences available: {n_3di}")
        log(f"  AA  sequences available: {n_aa}")

        if n_3di == 0:
            log("ERROR: No 3di sequences found for hypothetical proteins.")
            sys.exit(1)

    # ---- STEP 2: Build query DB (createdb route) ----
    # On foldseek 10.x easy-search rejects --seq-type and --prostt5-input is
    # only available in GPU/ProstT5 builds, so we route everything through
    # createdb + search + convertalis.  When AA FASTA is available we use
    # alignment-type 2 (3Di+AA); otherwise alignment-type 1 (3di-only).
    combined_query_db = None
    if not args.parse_only and not args.easy_search:
        section("STEP 2 — Build query database (createdb route)")
        combined_query_db = build_combined_query_db(foldseek_bin)
        if combined_query_db is None:
            log("  Query DB build failed — falling back to easy-search "
                "(may fail on foldseek 10.x).")

    # ---- STEP 3: Run FoldSeek searches ----
    per_db_m8_files = []
    if not args.parse_only:
        section("STEP 3 — FoldSeek local searches")
        failed_dbs = []
        for db_name, db_path in dbs_to_search.items():
            out_m8 = PER_DB_HITS_DIR / f"{db_name}.m8"
            ok = run_foldseek_easy_search(
                foldseek_bin=foldseek_bin,
                db_name=db_name,
                db_path=db_path,
                out_m8=out_m8,
                combined_query_db=combined_query_db,
                threads=args.threads,
                evalue=args.evalue,
                max_seqs=args.max_seqs,
            )
            if ok:
                per_db_m8_files.append(out_m8)
            else:
                failed_dbs.append(db_name)

        if failed_dbs:
            log(f"\nFailed databases: {failed_dbs}")
            log("  Check database paths in config.py and re-run.")

        if not per_db_m8_files:
            log("No successful searches.  Exiting.")
            sys.exit(1)

    else:
        # Parse mode: collect all existing per-DB m8 files
        per_db_m8_files = sorted(PER_DB_HITS_DIR.glob("*.m8"))
        if not per_db_m8_files:
            log(f"No per-DB m8 files found in {PER_DB_HITS_DIR}.")
            log("  Run without --parse_only first to execute FoldSeek searches.")
            sys.exit(1)
        log(f"Parse-only mode: found {len(per_db_m8_files)} per-DB m8 files")

    # ---- STEP 4: Merge + split into per-gene files ----
    section("STEP 4 — Merge and split by gene")
    split_and_merge_hits(
        per_db_m8_files=per_db_m8_files,
        raw_hits_dir=RAW_HITS_DIR,
        gene_ids_in_scope=hypo_locus_tag_set,
        force=args.force,
    )

    # ---- STEP 5: Parse all hits ----
    section("STEP 5 — Parse and score hits")
    hits = parse_all_local_hits(
        raw_hits_dir=RAW_HITS_DIR,
        gene_ids_in_scope=hypo_locus_tag_set,
        hypo_meta=hypo_meta,
    )

    if hits.empty:
        log("No hits parsed.  Did FoldSeek find any matches?")
        sys.exit(0)

    log(f"Total hits: {len(hits)} across {hits['gene'].nunique()} genes")
    inf_count = hits["informative"].sum()
    log(f"  Informative (description OK + quality filter): {inf_count}")

    # ---- Build best + top3 ----
    # Attach aa_length so build_best_and_top3 can compute qcov_frac
    if "aa_length" in hypo_meta.columns:
        aa_map = hypo_meta.set_index("locus_tag")["aa_length"].to_dict()
        hits["aa_length"] = hits["gene"].map(aa_map)

    best, top3 = build_best_and_top3(hits, all_hypo_genes)

    # ---- Propagate foldseek_subdb (source-DB tag) ----
    # Map (gene, accession) -> foldseek_subdb from the parsed hits and attach
    # to best/top3 so downstream comparison + curation can see which DB the
    # hit came from (afdb-swissprot / pdb100 / baktfold-afdb).
    if "foldseek_subdb" in hits.columns:
        subdb_map = (
            hits[["gene", "accession", "foldseek_subdb"]]
            .drop_duplicates(["gene", "accession"])
        )
        if "accession" in best.columns:
            best = best.merge(subdb_map, on=["gene", "accession"], how="left")
        if "accession" in top3.columns:
            top3 = top3.merge(subdb_map, on=["gene", "accession"], how="left")

    # ---- Add prophage + foldseek_top3 columns (cluster enrichment) ----
    # prophage: from hypo_meta if present (built from gene_metadata.csv), else
    # derived from the per-prophage phold split. Carried so every output table
    # is traceable to its prophage (the cluster runs 8 contigs, not 1).
    prophage_map = {}
    if "prophage" in hypo_meta.columns:
        prophage_map = hypo_meta.set_index("locus_tag")["prophage"].to_dict()

    def _attach_prophage(df):
        if not prophage_map or "gene" not in df.columns:
            return df
        cols = list(df.columns)
        df = df.copy()
        df.insert(0, "prophage", df["gene"].map(prophage_map))
        return df

    # foldseek_top3: pipe-joined top-3 descriptions per gene, merged onto best
    # so the single-best table also shows the runner-up structural calls.
    if not top3.empty and "description" in top3.columns:
        top3_str = (
            top3.assign(_d=top3["description"].fillna("").astype(str))
                .groupby("gene")["_d"]
                .apply(lambda s: " | ".join([x for x in s if x][:3]))
                .rename("foldseek_top3")
                .reset_index()
        )
        best = best.merge(top3_str, on="gene", how="left")
        best["foldseek_top3"] = best["foldseek_top3"].fillna("NA")

    best = _attach_prophage(best)
    top3 = _attach_prophage(top3)

    # ---- Write outputs ----
    best.to_csv(str(FOLDSEEK_3DI_BEST), index=False)
    top3.to_csv(str(FOLDSEEK_3DI_TOP3), index=False)
    hits.drop(columns=["informative"], errors="ignore").to_csv(
        str(FOLDSEEK_3DI_ALL), index=False
    )

    log(f"\nOutputs written:")
    log(f"  Best hit : {FOLDSEEK_3DI_BEST}")
    log(f"  Top 3    : {FOLDSEEK_3DI_TOP3}")
    log(f"  All hits : {FOLDSEEK_3DI_ALL}")

    # ---- Summary report ----
    section("3Di FOLDSEEK SUMMARY REPORT")
    n_inf    = best["informative_hit_found"].sum() if "informative_hit_found" in best.columns else 0
    n_total  = len(all_hypo_genes)
    n_done   = best["score"].notna().sum() if "score" in best.columns else 0
    n_no_hit = n_total - n_done
    pct      = 100 * n_inf // max(n_done, 1)

    log(f"Target genes       : {n_total}  (all hypotheticals, incl. >400aa)")
    log(f"Genes with any hit : {n_done}")
    log(f"Informative hits   : {n_inf}  ({pct}%)")
    log(f"No hit at all      : {n_no_hit}")
    log(f"Query mode         : {'3Di+AA (combined)' if combined_query_db else '3di-only (prostt5/seq-type 1)'}")

    # Confidence distribution
    if "foldseek_confidence" in best.columns:
        conf_counts = best["foldseek_confidence"].value_counts()
        log("\nConfidence distribution:")
        for tier in ["CONFIDENT", "GOOD", "BORDERLINE", "WEAK", "NO_HIT"]:
            n = conf_counts.get(tier, 0)
            log(f"  {tier:<12}: {n}")

    # Same-host hits (potential AMG/moron)
    if "same_host" in best.columns:
        same_host = best[best["same_host"] == True]
        if not same_host.empty:
            log(f"\nSame-host hits ({HOST_GENUS}) — potential AMG/moron candidates: "
                f"{len(same_host)}")
            for _, r in same_host.iterrows():
                log("  " + str(r["gene"]) + ": "
                    + str(r.get("description", "?"))[:60]
                    + f" (evalue={r.get('evalue', '?')})")

    # Defense-system hits
    if "defense_flag" in best.columns:
        defense = best[best["defense_flag"] == True]
        if not defense.empty:
            log(f"\nDefense-system hits (check DefenseFinder/PADLOC): {len(defense)}")
            for _, r in defense.iterrows():
                log("  " + str(r["gene"]) + ": " + str(r.get("description", "?"))[:60])

    # Genes still unresolved after 3di search
    no_inf = best[best.get("informative_hit_found", pd.Series(False)) == False]["gene"].tolist() \
        if "informative_hit_found" in best.columns else []
    if no_inf:
        log(f"\nGenes without informative 3di hit ({len(no_inf)}) — candidates for manual curation:")
        log(f"  {no_inf[:20]}")
        if len(no_inf) > 20:
            log(f"  ... and {len(no_inf) - 20} more")

    log("\nDone.  Results are in 02_foldseek/3di_tokens/")
    log("Next step: run 03_compare_annotations.py with 3di results for comparison with phold.")


if __name__ == "__main__":
    main()
