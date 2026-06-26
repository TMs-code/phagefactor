#!/usr/bin/env python3
"""
01c_merge_phold.py
==================
Post-process PHold output and write a single merged TSV.

Supports two input layouts automatically:

  Layout A -- single combined run (01_run_phold_local.sh / 01b_run_phold_hpc.sh):
    01_phold/combined/phold_per_cds_predictions.tsv

  Layout B -- per-sample run (01_run_pharokka_phold.sh):
    01_pharokka/<sample>/phold/phold_per_cds_predictions.tsv
    (one TSV per prophage, merged here)

Writes:
  01_phold/combined/phold_all.tsv        (tab-separated, with 'prophage' column)
  01_phold/combined/phold_all_simple.csv (key columns only, for quick inspection)

Usage:
  cd phagefactor/
  python scripts/01c_merge_phold.py
"""

import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).parent
_PROJECT_DIR = _SCRIPTS_DIR.parent
sys.path.insert(0, str(_PROJECT_DIR))

from config import (PHOLD_COMB_DIR, PHOLD_COMBINED_TSV,
                    PHAROKKA_OUT_DIR, PROPHAGE_NAMES, PHOLD_TSV_FILENAME,
                    PHOLD_OUT_DIR)
from utils import log, section, load_phold_tsv, is_informative

try:
    import pandas as pd
except ImportError:
    print("pandas required. Install with: pip install pandas")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Input discovery
# ---------------------------------------------------------------------------

def _find_combined_tsv():
    """Return path to the combined PHold TSV if it exists (layout A)."""
    p = PHOLD_COMB_DIR / PHOLD_TSV_FILENAME
    if p.exists() and p.stat().st_size > 0:
        return p
    return None


def _find_per_sample_tsvs():
    """
    Return list of (sample_name, Path) pairs from the pharokka per-sample
    output directory (layout B: 01_pharokka/<sample>/phold/phold_per_cds_predictions.tsv).
    """
    hits = []
    if not PHAROKKA_OUT_DIR.exists():
        return hits
    for sample_dir in sorted(PHAROKKA_OUT_DIR.iterdir()):
        if not sample_dir.is_dir():
            continue
        tsv = sample_dir / "phold" / PHOLD_TSV_FILENAME
        if tsv.exists() and tsv.stat().st_size > 0:
            hits.append((sample_dir.name, tsv))
    return hits


def _find_cluster_phold_tsvs():
    """
    Cluster layout (Layout C): phold writes per-prophage outputs to
    01_phold/<NAME>/<NAME>_per_cds_predictions.tsv (phold --prefix <NAME>).
    Return list of (sample_name, Path) pairs.
    """
    hits = []
    if not PHOLD_OUT_DIR.exists():
        return hits
    for sample_dir in sorted(PHOLD_OUT_DIR.iterdir()):
        if not sample_dir.is_dir() or sample_dir.name == "combined":
            continue
        sample_name = sample_dir.name
        # The phold prefix is the prophage name
        tsv = sample_dir / f"{sample_name}_per_cds_predictions.tsv"
        if tsv.exists() and tsv.stat().st_size > 0:
            hits.append((sample_name, tsv))
            continue
        # Fall back: any *_per_cds_predictions.tsv in the dir
        matches = list(sample_dir.glob("*_per_cds_predictions.tsv"))
        if matches:
            hits.append((sample_name, matches[0]))
    return hits


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    section("STEP 01c -- PROCESS PHOLD OUTPUT")

    PHOLD_COMB_DIR.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 1. Discover input                                                    #
    # ------------------------------------------------------------------ #
    combined_tsv   = _find_combined_tsv()
    per_sample     = _find_per_sample_tsvs()       # Layout B (legacy)
    cluster_sample = _find_cluster_phold_tsvs()    # Layout C (this pipeline)

    if combined_tsv:
        log(f"Layout A detected: single combined TSV")
        log(f"  {combined_tsv}")
        merged = load_phold_tsv(combined_tsv)
        log(f"  {len(merged)} CDS rows, {len(merged.columns)} columns")
        # prophage name comes from contig_id (= LOCUS from GenBank record)
        if "contig_id" not in merged.columns:
            log("ERROR: 'contig_id' column missing in PHold output.")
            log(f"  Available: {list(merged.columns)}")
            sys.exit(1)
        merged["prophage"] = merged["contig_id"]

    elif per_sample:
        log(f"Layout B detected: {len(per_sample)} per-sample TSV(s) in {PHAROKKA_OUT_DIR}/")
        dfs = []
        for sample_name, tsv_path in per_sample:
            df = load_phold_tsv(tsv_path)
            # Use contig_id if available, else fall back to sample dir name
            if "contig_id" in df.columns:
                df["prophage"] = df["contig_id"]
            else:
                df["prophage"] = sample_name
            log(f"  {sample_name}: {len(df)} CDS")
            dfs.append(df)
        merged = pd.concat(dfs, ignore_index=True)
        log(f"  Total: {len(merged)} CDS rows")

    elif cluster_sample:
        log(f"Layout C detected: {len(cluster_sample)} per-prophage TSV(s) in {PHOLD_OUT_DIR}/")
        dfs = []
        for sample_name, tsv_path in cluster_sample:
            df = load_phold_tsv(tsv_path)
            if "contig_id" in df.columns:
                df["prophage"] = df["contig_id"]
            else:
                df["prophage"] = sample_name
            log(f"  {sample_name}: {len(df)} CDS  ({tsv_path.name})")
            dfs.append(df)
        merged = pd.concat(dfs, ignore_index=True)
        log(f"  Total: {len(merged)} CDS rows")

    else:
        log("ERROR: No PHold output found.")
        log(f"  Looked for:")
        log(f"    Layout A: {PHOLD_COMB_DIR / PHOLD_TSV_FILENAME}")
        log(f"    Layout B: {PHAROKKA_OUT_DIR}/<sample>/phold/{PHOLD_TSV_FILENAME}")
        log(f"    Layout C: {PHOLD_OUT_DIR}/<NAME>/<NAME>_per_cds_predictions.tsv")
        log("")
        log("  Run one of the following first:")
        log("    bash scripts/01_run_phold_local.sh       (local, combined GB)")
        log("    sbatch scripts/01b_run_phold_hpc.sh      (SLURM, combined GB)")
        log("    bash scripts/01_run_pharokka_phold.sh    (Pharokka + PHold per FASTA)")
        log("    sbatch steps/01_phold_array.sh           (cluster array job)")
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # 2. Validate prophage names                                           #
    # ------------------------------------------------------------------ #
    found_prophages = sorted(merged["prophage"].unique())
    log("\nProphages found in PHold output: " + str(found_prophages))
    unexpected = [p for p in found_prophages if p not in PROPHAGE_NAMES]
    missing    = [p for p in PROPHAGE_NAMES if p not in found_prophages]
    if unexpected:
        log("  WARNING: Unexpected prophage IDs: " + str(unexpected))
        log("  (contig_id values may differ from config.PROPHAGE_NAMES -- check LOCUS names in GB)")
    if missing:
        log("  WARNING: No PHold output for: " + str(missing))

    # ------------------------------------------------------------------ #
    # 3. Per-prophage statistics                                           #
    # ------------------------------------------------------------------ #
    log("\n-- Per-prophage annotation summary --")
    log("{:<10} {:>6} {:>12} {:>6} {:>10} {:>9} {:>8}".format(
        "Prophage", "Total", "Informative", "%Inf", "Conf=high", "Conf=med", "Conf=low"))
    log("-" * 65)

    for name in found_prophages:
        grp = merged[merged["prophage"] == name]
        n = len(grp)
        n_inf = grp["product"].apply(
            lambda x: is_informative(str(x) if x is not None else "")
        ).sum()
        n_high = (grp["annotation_confidence"] == "high").sum()
        n_med  = (grp["annotation_confidence"] == "medium").sum()
        n_low  = (grp["annotation_confidence"] == "low").sum()
        pct    = 100 * n_inf / n if n > 0 else 0
        log("{:<10} {:>6} {:>12} {:>5.0f}% {:>10} {:>9} {:>8}".format(
            name, n, n_inf, pct, n_high, n_med, n_low))

    n_total     = len(merged)
    n_inf_total = merged["product"].apply(
        lambda x: is_informative(str(x) if x is not None else "")
    ).sum()
    log("-" * 65)
    pct_total = 100 * n_inf_total / n_total if n_total > 0 else 0
    log("{:<10} {:>6} {:>12} {:>5.0f}%".format("TOTAL", n_total, n_inf_total, pct_total))

    log("\n-- Annotation method breakdown --")
    if "annotation_method" in merged.columns:
        log(merged["annotation_method"].value_counts().fillna("none").to_string())

    # ------------------------------------------------------------------ #
    # 4. Write outputs                                                     #
    # ------------------------------------------------------------------ #
    merged.to_csv(str(PHOLD_COMBINED_TSV), sep="\t", index=False)
    log("\nPHold combined TSV -> " + str(PHOLD_COMBINED_TSV))
    log("  " + str(len(merged)) + " rows, " + str(len(merged.columns)) + " columns")

    simple_cols = ["prophage", "cds_id", "contig_id", "phrog",
                   "function", "product", "annotation_confidence",
                   "bitscore", "evalue", "fident", "qCov", "prostt5_confidence"]
    simple_cols = [c for c in simple_cols if c in merged.columns]
    simple_out = PHOLD_COMB_DIR / "phold_all_simple.csv"
    merged[simple_cols].to_csv(str(simple_out), index=False)
    log("Simplified CSV    -> " + str(simple_out))

    log("\n-> Next steps (cluster):")
    log("     bash steps/01d_merge_3di.sh           # concat 3di FASTAs")
    log("     SKIP_PHOLD=1 bash submit_all.sh       # chains 02 (FoldSeek 3di) -> 03 -> 04")


if __name__ == "__main__":
    main()
