#!/usr/bin/env python3
"""
01p_merge_phold_proteins.py
===========================
Merge per-batch phold proteins-{predict,compare} outputs (written by the
01p_phold_proteins_array.sh SLURM array — one batch of PROTEIN_BATCH_SIZE
proteins per task) into the single combined layout that steps 02-06 expect.

Mirrors 01c_merge_phold.py (genome mode), but additionally consolidates two
NEWLY load-bearing pieces of evidence (both surfaced in the 2026-06 fixes):

  - sub_db_tophits/*_cds_predictions.tsv   (ACR/VFDB/CARD/NetFlaX/DefenseFinder
    structured hits — 03_compare_annotations.py's _load_subdb_hits() globs
    PHOLD_OUT_DIR recursively, so leaving these scattered under
    01_phold/proteins/batches/*/compare/sub_db_tophits/ would technically
    still be found, but consolidating avoids globbing hundreds of small files
    and keeps a single predictable location for inspection)

  - *_prostT5_3di_all_probabilities.json   (per-residue masking probabilities —
    02d_foldseek_3di.py's _load_prostt5_probs() also globs recursively, so
    these are likewise discoverable in place; we leave them where phold wrote
    them rather than duplicating large per-batch JSON blobs, and just confirm
    that the recursive glob can see them)

Reads (per batch, from 01_phold/proteins/batches/<batch>/):
  compare/phold_per_cds_predictions.tsv
  compare/phold_3di.fasta
  compare/phold_aa.fasta
  compare/sub_db_tophits/{acr,vfdb,card,netflax,defensefinder}_cds_predictions.tsv
  predictions/*_prostT5_3di_all_probabilities.json   (left in place, just counted)

Writes:
  01_phold/phold_per_cds_predictions.tsv          (merged, single header)
  01_phold/combined/phold_all.tsv                 (= same content, expected by 03)
  01_phold/combined/phold_3di.fasta               (concatenated, expected by 02)
  01_phold/combined/phold_aa.fasta                (concatenated)
  01_phold/proteins/compare/sub_db_tophits/{source}_cds_predictions.tsv
                                                  (concatenated per source)
  input/hypothetical_gene_list.txt                (regenerated, all protein IDs)
  input/split/gene_metadata.csv                   (regenerated; columns match
                                                   what 03_compare_annotations.py
                                                   expects)

Usage:
  cd phagefactor/
  python scripts/01p_merge_phold_proteins.py
"""

import sys
import re
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).parent
_PROJECT_DIR = _SCRIPTS_DIR.parent
sys.path.insert(0, str(_PROJECT_DIR))

from config import PHOLD_OUT_DIR, PHOLD_COMB_DIR, INPUT_DIR
from utils import log, section

try:
    import pandas as pd
except ImportError:
    print("pandas required. Install with: pip install pandas")
    sys.exit(1)


BATCHES_DIR = PHOLD_OUT_DIR / "proteins" / "batches"
MERGED_COMPARE_DIR = PHOLD_OUT_DIR / "proteins" / "compare"
MERGED_SUBDB_DIR = MERGED_COMPARE_DIR / "sub_db_tophits"
HYPO_LIST = INPUT_DIR / "hypothetical_gene_list.txt"
GENE_META = INPUT_DIR / "split" / "gene_metadata.csv"

_SUBDB_SOURCES = ("acr", "vfdb", "card", "netflax", "defensefinder")


def _find_batches():
    """Return sorted list of (batch_name, batch_dir) for completed batches."""
    if not BATCHES_DIR.is_dir():
        return []
    out = []
    for d in sorted(BATCHES_DIR.iterdir()):
        if not d.is_dir():
            continue
        marker = d / ".batch_done"
        compare_tsv = d / "compare" / "phold_per_cds_predictions.tsv"
        if marker.exists() and compare_tsv.exists():
            out.append((d.name, d))
        else:
            log(f"  WARNING: {d.name} missing .batch_done or compare TSV — skipping "
                f"(re-submit array task if this batch should have completed)")
    return out


# ---------------------------------------------------------------------------
# 1) Merge phold_per_cds_predictions.tsv across batches
# ---------------------------------------------------------------------------
def _merge_predictions_tsv(batches):
    frames = []
    for name, bdir in batches:
        tsv = bdir / "compare" / "phold_per_cds_predictions.tsv"
        try:
            df = pd.read_csv(str(tsv), sep="\t", low_memory=False)
        except Exception as exc:
            log(f"  WARNING: could not read {tsv}: {exc}")
            continue
        df["_batch"] = name
        frames.append(df)
    if not frames:
        return None
    merged = pd.concat(frames, ignore_index=True)
    merged = merged.drop(columns=["_batch"])
    return merged


# ---------------------------------------------------------------------------
# 2) Concatenate FASTA files (3Di + AA) across batches
# ---------------------------------------------------------------------------
def _concat_fasta(batches, relpath, out_path):
    n_seqs = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out:
        for name, bdir in batches:
            src = bdir / "compare" / relpath
            if not src.exists():
                log(f"  WARNING: missing {src}")
                continue
            with open(src) as fh:
                for line in fh:
                    if line.startswith(">"):
                        n_seqs += 1
                    out.write(line)
    return n_seqs


# ---------------------------------------------------------------------------
# 3) Concatenate sub_db_tophits/*_cds_predictions.tsv per source
#
# Locus tags are globally unique (gene IDs are prophage-prefixed and batches
# partition the input disjointly), so straight concatenation is safe — no
# dedup needed here (03_compare_annotations._load_subdb_hits() also dedupes
# defensively by bitscore, so even residual overlap would be handled).
# ---------------------------------------------------------------------------
def _merge_subdb_tophits(batches):
    counts = {}
    MERGED_SUBDB_DIR.mkdir(parents=True, exist_ok=True)
    for source in _SUBDB_SOURCES:
        frames = []
        for name, bdir in batches:
            f = bdir / "compare" / "sub_db_tophits" / f"{source}_cds_predictions.tsv"
            if not f.exists() or f.stat().st_size == 0:
                continue
            try:
                df = pd.read_csv(str(f), sep="\t", low_memory=False)
            except Exception as exc:
                log(f"  WARNING: could not read {f}: {exc}")
                continue
            if df.empty:
                continue
            frames.append(df)
        out_path = MERGED_SUBDB_DIR / f"{source}_cds_predictions.tsv"
        if frames:
            merged = pd.concat(frames, ignore_index=True)
            merged.to_csv(str(out_path), sep="\t", index=False)
            counts[source] = len(merged)
        else:
            # Touch an empty file so downstream globbing/inspection sees a
            # consistent layout (matches the empty acr/netflax files phold
            # itself produces when a sub-DB has no hits).
            out_path.touch()
            counts[source] = 0
    return counts


# ---------------------------------------------------------------------------
# 4) Confirm per-residue probability JSONs are discoverable (left in place)
# ---------------------------------------------------------------------------
def _count_probability_jsons():
    return list(PHOLD_OUT_DIR.glob("**/*_prostT5_3di_all_probabilities.json"))


# ---------------------------------------------------------------------------
# 5) Regenerate hypothetical_gene_list.txt + gene_metadata.csv
#    (same logic previously inline in 01p_phold_proteins.sh, generalized to
#    read from the merged combined/ files instead of a single-job output)
# ---------------------------------------------------------------------------
def _write_hypothetical_list(aa_fasta):
    ids = []
    with open(aa_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    HYPO_LIST.parent.mkdir(parents=True, exist_ok=True)
    HYPO_LIST.write_text("\n".join(ids) + ("\n" if ids else ""))
    return len(ids)


def _write_gene_metadata(merged_tsv_path, aa_fasta):
    # AA lengths
    aa_lengths = {}
    cur, n = None, 0
    with open(aa_fasta) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur is not None:
                    aa_lengths[cur] = n
                hdr = line[1:].split()[0]
                cur = hdr.split(":", 1)[1] if ":" in hdr else hdr
                n = 0
            else:
                n += len(line.strip())
        if cur is not None:
            aa_lengths[cur] = n

    df = pd.read_csv(str(merged_tsv_path), sep="\t", low_memory=False)
    gcol = "gene" if "gene" in df.columns else df.columns[0]

    rows = []
    for _, r in df.iterrows():
        gene_raw = str(r[gcol])
        lt = gene_raw.split(":", 1)[-1]
        proph = gene_raw.split(":", 1)[0] if ":" in gene_raw else "input"
        aal = aa_lengths.get(lt, 0)
        rows.append(dict(prophage=proph, locus_tag=lt, function="unknown function",
                         aa_length=aal, is_hypothetical=True))

    GENE_META.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(str(GENE_META), index=False)
    return len(rows)


def main():
    section("MERGING PROTEIN-MODE PHOLD BATCHES")

    batches = _find_batches()
    if not batches:
        log(f"ERROR: no completed batches found under {BATCHES_DIR}")
        log("  (each batch needs a .batch_done marker + compare/phold_per_cds_predictions.tsv)")
        sys.exit(1)
    log(f"Found {len(batches)} completed batch(es): {', '.join(n for n, _ in batches)}")

    # 1) Merge predictions TSV -> 01_phold/ root + combined/
    section("Merging phold_per_cds_predictions.tsv")
    merged = _merge_predictions_tsv(batches)
    if merged is None or merged.empty:
        log("ERROR: failed to merge any per-batch prediction TSVs")
        sys.exit(1)
    PHOLD_COMB_DIR.mkdir(parents=True, exist_ok=True)
    root_tsv = PHOLD_OUT_DIR / "phold_per_cds_predictions.tsv"
    combined_tsv = PHOLD_COMB_DIR / "phold_all.tsv"
    merged.to_csv(str(root_tsv), sep="\t", index=False)
    merged.to_csv(str(combined_tsv), sep="\t", index=False)
    log(f"Merged {len(merged)} gene rows from {len(batches)} batches")
    log(f"  -> {root_tsv}")
    log(f"  -> {combined_tsv}")

    # 2) Concatenate 3Di + AA FASTAs -> combined/
    section("Concatenating 3Di / AA FASTAs")
    n_3di = _concat_fasta(batches, "phold_3di.fasta", PHOLD_COMB_DIR / "phold_3di.fasta")
    n_aa = _concat_fasta(batches, "phold_aa.fasta", PHOLD_COMB_DIR / "phold_aa.fasta")
    log(f"  phold_3di.fasta -> {PHOLD_COMB_DIR / 'phold_3di.fasta'}  ({n_3di} sequences)")
    log(f"  phold_aa.fasta  -> {PHOLD_COMB_DIR / 'phold_aa.fasta'}  ({n_aa} sequences)")
    if n_3di != len(merged) or n_aa != len(merged):
        log(f"  WARNING: sequence counts (3di={n_3di}, aa={n_aa}) don't match "
            f"merged gene rows ({len(merged)}) — check for missing/failed batches")

    # 3) Merge sub_db_tophits per source (newly load-bearing — see Task #15)
    section("Merging sub_db_tophits (ACR/VFDB/CARD/NetFlaX/DefenseFinder)")
    subdb_counts = _merge_subdb_tophits(batches)
    for source, n in subdb_counts.items():
        log(f"  {source:<14} {n} hit(s) -> {MERGED_SUBDB_DIR / (source + '_cds_predictions.tsv')}")
    log("(03_compare_annotations._load_subdb_hits() globs PHOLD_OUT_DIR recursively, "
        "so this consolidated location is found automatically — no config change needed)")

    # 4) Confirm per-residue probability JSONs remain discoverable (newly load-bearing — Task #16)
    section("Checking per-residue masking-probability JSONs")
    json_files = _count_probability_jsons()
    log(f"  Found {len(json_files)} *_prostT5_3di_all_probabilities.json file(s) "
        f"under {PHOLD_OUT_DIR} (left in place under proteins/batches/*/predictions/)")
    log("(02d_foldseek_3di._load_prostt5_probs() now globs '**/...' recursively — "
        "the Task #16 fix — so these are discovered regardless of batch-array depth)")

    # 5) Regenerate hypothetical_gene_list.txt + gene_metadata.csv
    section("Regenerating hypothetical_gene_list.txt and gene_metadata.csv")
    aa_fasta = PHOLD_COMB_DIR / "phold_aa.fasta"
    n_hypo = _write_hypothetical_list(aa_fasta)
    log(f"  hypothetical_gene_list.txt -> {HYPO_LIST}  ({n_hypo} protein IDs)")
    n_meta = _write_gene_metadata(combined_tsv, aa_fasta)
    log(f"  gene_metadata.csv -> {GENE_META}  ({n_meta} rows)")

    log("")
    log("=== Protein-mode batch merge complete ===")
    log(f"  Batches merged   : {len(batches)}")
    log(f"  Genes (rows)     : {len(merged)}")
    log(f"  Annotation TSV   : {root_tsv}  (combined copy: {combined_tsv})")
    log("")
    log("Next: sbatch steps/02_foldseek_3di.sh")
    log("Done.")


if __name__ == "__main__":
    main()
