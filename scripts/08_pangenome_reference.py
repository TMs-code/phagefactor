#!/usr/bin/env python3
# =============================================================================
# 08_pangenome_reference.py — pan-genome mode post-processing.
#
# Rejoin phageFACTor's final_annotations_table (keyed on the clean pang_NNNNN
# ids from 00n_translate_pangenome.py) back to the biological group names +
# protein sequences, producing a compact HOMOLOGY-REFERENCE set:
#
#   <prefix>.csv / .xlsx   one row per group:
#     group_id, representative_locus, clean_id, aa_length, aa_sequence,
#     final_product, short_name, annotation_source, agreement,
#     (evidence_category / needs_review / fs_* if present)
#
#   <prefix>.faa           the SAME representatives as a FASTA, header =
#                          ">group_id  final_product", ready for
#                          `diamond makedb` / `mmseqs createdb`. A homology hit's
#                          sseqid == group_id joins straight back to the CSV row,
#                          so you transfer final_product/short_name onto whatever
#                          query protein hit it.
#
# WHY A FASTA TOO: DIAMOND / MMseqs2 / BLAST index sequences, not CSVs. The CSV
# carries the annotation; the .faa is what you actually search against. Keeping
# the header id == group_id == CSV key makes annotation transfer a plain join.
#
# USAGE
#   python scripts/08_pangenome_reference.py \
#       --final    <run_dir>/04_output/final_annotations_table.csv \
#       --id-map   pangenome_id_map.csv \
#       --faa      input/fasta/<name>_proteins.faa \
#       --out-prefix <run_dir>/04_output/pangenome_reference
# =============================================================================

import argparse
import re
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas required (pip install pandas).", file=sys.stderr)
    sys.exit(1)


def log(msg):
    print(f"[08_pangenome_ref] {msg}", flush=True)


def read_faa(path: Path) -> dict:
    """clean_id -> aa sequence. Indexes each record under BOTH its full header
    token AND the part after a trailing ':' (00p prefixes headers as
    'PROPH:clean_id'), so the id map's clean_id matches either layout."""
    seqs, cur, buf = {}, None, []
    if not path or not Path(path).exists():
        return seqs

    def _store(h, s):
        seqs[h] = s
        if ":" in h:
            seqs.setdefault(h.split(":", 1)[1], s)

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    _store(cur, "".join(buf))
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        _store(cur, "".join(buf))
    return seqs


# annotation columns we surface, in order, IF present in the final table
_ANNOT_COLS = ["final_product", "short_name", "final_function",
               "annotation_source", "agreement", "evidence_category",
               "needs_review", "best_hit_kingdom",
               "fs_description", "fs_confidence", "phold_product"]


def _fasta_safe(text: str) -> str:
    """collapse whitespace so a product string is a legal one-line FASTA desc."""
    return re.sub(r"\s+", " ", str(text)).strip()


def main():
    ap = argparse.ArgumentParser(description="Build the pan-genome homology-"
                                 "reference set (group + sequence + annotation).")
    ap.add_argument("--final", required=True, type=Path,
                    help="final_annotations_table.csv from the run")
    ap.add_argument("--id-map", required=True, type=Path,
                    help="pangenome_id_map.csv from 00n_translate_pangenome.py")
    ap.add_argument("--faa", type=Path, default=None,
                    help="protein FASTA (attaches aa_sequence + enables the .faa output)")
    ap.add_argument("--out-prefix", required=True, type=Path,
                    help="output path prefix (.csv/.xlsx/.faa appended)")
    ap.add_argument("--key", default="locus_tag",
                    help="column in --final holding the clean id (default locus_tag)")
    ap.add_argument("--no-faa", action="store_true",
                    help="do not emit the DIAMOND/MMseqs-ready .faa")
    args = ap.parse_args()

    for p in (args.final, args.id_map):
        if not p.exists():
            log(f"ERROR: not found: {p}")
            sys.exit(1)

    idmap = pd.read_csv(str(args.id_map), dtype=str).fillna("")
    if "clean_id" not in idmap.columns:
        log(f"ERROR: {args.id_map} lacks a 'clean_id' column "
            f"(has: {list(idmap.columns)})")
        sys.exit(1)
    final = pd.read_csv(str(args.final), dtype=str).fillna("")
    if args.key not in final.columns:
        log(f"ERROR: final table lacks key column '{args.key}' "
            f"(has: {list(final.columns)})")
        sys.exit(1)

    seqs = read_faa(args.faa) if args.faa else {}

    # Left-join from the id map so EVERY input gene appears, even any the
    # pipeline may have dropped (those get blank annotations, flagged below).
    merged = idmap.merge(final, how="left", left_on="clean_id", right_on=args.key)

    if "final_product" in merged:
        n_missing = int((merged["final_product"] == "").sum()
                        + merged["final_product"].isna().sum())
        if n_missing:
            log(f"WARNING: {n_missing} genes from the id map have NO/blank row in the "
                "final table (dropped/failed batch?) — kept with blank annotation.")

    present = [c for c in _ANNOT_COLS if c in merged.columns]
    out = pd.DataFrame({
        "group_id":             merged["group"].where(merged["group"].astype(bool),
                                                       merged["clean_id"]),
        "representative_locus": merged["orig_id"],
        "clean_id":             merged["clean_id"],
        "aa_length":            merged.get("aa_len", ""),
    })
    for c in present:
        out[c] = merged[c].fillna("")
    if seqs:
        out["aa_sequence"] = merged["clean_id"].map(seqs).fillna("")

    # one row per group; representative-per-group is 1:1 for a Roary pan-genome,
    # but if a group repeats keep the most-informative (non-hypothetical) first.
    if "final_product" in out.columns:
        out["_hypo"] = (out["final_product"].isin(["", "hypothetical protein"])).astype(int)
        out = out.sort_values(["group_id", "_hypo"]).drop(columns="_hypo")
    n_before = len(out)
    out = out.drop_duplicates(subset="group_id", keep="first").reset_index(drop=True)
    if len(out) != n_before:
        log(f"Collapsed {n_before} rows -> {len(out)} unique groups.")

    csv_path = args.out_prefix.with_suffix(".csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(str(csv_path), index=False)
    log(f"CSV  -> {csv_path}  ({len(out)} groups, {len(out.columns)} cols)")

    try:
        xlsx_path = args.out_prefix.with_suffix(".xlsx")
        out.to_excel(str(xlsx_path), index=False)
        log(f"XLSX -> {xlsx_path}")
    except Exception as exc:
        log(f"(xlsx skipped: {exc} — install openpyxl for Excel output)")

    # DIAMOND/MMseqs-ready FASTA (header id == group_id == CSV key)
    if not args.no_faa and "aa_sequence" in out.columns:
        faa_path = args.out_prefix.with_suffix(".faa")
        n = 0
        with open(faa_path, "w") as fh:
            for _, r in out.iterrows():
                seq = r.get("aa_sequence", "")
                if not seq:
                    continue
                desc = _fasta_safe(r.get("final_product", "")) if "final_product" in out.columns else ""
                fh.write(f">{r['group_id']} {desc}\n{seq}\n")
                n += 1
        log(f"FAA  -> {faa_path}  ({n} sequences; `diamond makedb --in {faa_path.name}`)")
    elif not args.no_faa:
        log("(.faa skipped: no --faa given, so no sequences to write)")

    if "final_product" in out.columns:
        annotated = int((~out["final_product"].isin(["", "hypothetical protein"])).sum())
        log(f"Annotated (non-hypothetical): {annotated}/{len(out)} "
            f"({100*annotated//max(len(out),1)}%)")


if __name__ == "__main__":
    main()
