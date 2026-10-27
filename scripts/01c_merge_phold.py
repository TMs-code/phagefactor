#!/usr/bin/env python3
"""
01c_merge_phold.py
==================
Genome mode: merge the per-prophage PHold outputs (01_phold/<NAME>/) into
01_phold/combined/:

  phold_all.tsv          all CDS predictions (tab-separated, with 'prophage' column)
  phold_all_simple.csv   key columns only, for quick inspection
  phold_3di.fasta        3Di sequences of every prophage (input of step 02)
  phold_aa.fasta         amino-acid sequences of every prophage

Prophages are taken in the order of <input>/prophage_list.txt.
"""

import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).parent
_PROJECT_DIR = _SCRIPTS_DIR.parent
sys.path.insert(0, str(_PROJECT_DIR))
# scripts/lib/ holds the shared modules (config, utils, lexicon, ...).
# Harmless if they still sit next to this file.
_LIB_DIR = Path(__file__).resolve().parent / "lib"
if _LIB_DIR.is_dir():
    sys.path.insert(0, str(_LIB_DIR))

from config import (PHOLD_COMB_DIR, PHOLD_COMBINED_TSV, PROPHAGE_NAMES,
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

def _find_cluster_phold_tsvs():
    """
    phold --prefix <NAME> writes 01_phold/<NAME>/<NAME>_per_cds_predictions.tsv.
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


def _merge_fasta(suffix: str, out: Path) -> None:
    """Concatenate 01_phold/<NAME>/<NAME>_<suffix>.fasta in prophage_list order."""
    found, n_seq = 0, 0
    with open(out, "w") as fh:
        for name in PROPHAGE_NAMES:
            src = next((c for c in (PHOLD_OUT_DIR / name / f"{name}_{suffix}.fasta",
                                    PHOLD_OUT_DIR / name / f"{name}_phold_{suffix}.fasta",
                                    PHOLD_OUT_DIR / name / f"phold_{suffix}.fasta")
                        if c.is_file()), None)
            if src is None:
                log(f"  WARNING: no {suffix} FASTA for {name}")
                continue
            text = src.read_text()
            fh.write(text)
            found += 1
            n_seq += text.count(">")
    log(f"{suffix} FASTA: {found} prophages, {n_seq} sequences -> {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    section("STEP 01c -- PROCESS PHOLD OUTPUT")

    PHOLD_COMB_DIR.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 1. Discover input                                                    #
    # ------------------------------------------------------------------ #
    cluster_sample = _find_cluster_phold_tsvs()
    if not cluster_sample:
        log(f"ERROR: no PHold output found ({PHOLD_OUT_DIR}/<NAME>/<NAME>_per_cds_predictions.tsv)")
        sys.exit(1)
    log(f"{len(cluster_sample)} per-prophage TSV(s) in {PHOLD_OUT_DIR}/")
    dfs = []
    for sample_name, tsv_path in cluster_sample:
        df = load_phold_tsv(tsv_path)
        df["prophage"] = df["contig_id"] if "contig_id" in df.columns else sample_name
        log(f"  {sample_name}: {len(df)} CDS  ({tsv_path.name})")
        dfs.append(df)
    merged = pd.concat(dfs, ignore_index=True)
    log(f"  Total: {len(merged)} CDS rows")

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

    # ------------------------------------------------------------------ #
    # 5. Concatenate the per-prophage 3Di / AA FASTAs                      #
    # ------------------------------------------------------------------ #
    for suffix in ("3di", "aa"):
        _merge_fasta(suffix, PHOLD_COMB_DIR / f"phold_{suffix}.fasta")


if __name__ == "__main__":
    main()
