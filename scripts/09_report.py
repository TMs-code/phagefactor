#!/usr/bin/env python3
"""
09_report.py -- interactive, single-file HTML report of a phageFACTor run.

One self-contained file (data embedded, no internet needed): a genome map per
prophage coloured by function, a searchable/filterable table and, per gene, the
evidence behind the call (PHold, FoldSeek top-3, curation reason, phynteny).

  python scripts/09_report.py --run_dir <run>            # uses the standard layout
  python scripts/09_report.py --final <table.csv> [--gbk <updated.gb>] --out report.html

With --run_dir it reads
  04_output/final_annotations_table.csv
  05_phynteny/final_annotations_integrated.csv   (preferred when present: + phynteny)
  04_output/updated_prophages.gb                 (gene coordinates, when present)
and writes 04_output/report.html (the workflow runs it as its last step).
"""
import argparse
import datetime as _dt
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "lib" / "report_template.html"


def _coords_from_gbk(gbk: Path):
    """locus_tag -> (contig, start, end, strand); contigs -> length."""
    try:
        from Bio import SeqIO
    except ImportError:
        print("  biopython not available: map drawn without coordinates", file=sys.stderr)
        return {}, {}
    pos, contigs = {}, {}
    for rec in SeqIO.parse(str(gbk), "genbank"):
        contigs[rec.id] = len(rec.seq)
        for f in rec.features:
            if f.type != "CDS":
                continue
            for lt in f.qualifiers.get("locus_tag", []):
                pos[lt] = (rec.id, int(f.location.start) + 1, int(f.location.end),
                           1 if (f.location.strand or 1) >= 0 else -1)
    return pos, contigs


def build(final: Path, out: Path, gbk: Path = None, integrated: Path = None,
          title: str = "", run_label: str = ""):
    table = final
    has_phynteny = False
    if integrated and integrated.is_file():
        if integrated.stat().st_mtime >= final.stat().st_mtime:
            table, has_phynteny = integrated, True
        else:   # curation re-run after phynteny: the integrated table is stale
            print(f"  NOTE: {integrated.name} is older than {final.name}; phynteny columns "
                  f"not shown (re-run phynteny to include them)", file=sys.stderr)
    df = pd.read_csv(table, dtype=str).fillna("")
    if df.empty:
        raise SystemExit(f"empty table: {table}")

    pos, contigs = ({}, {})
    if gbk and gbk.is_file():
        pos, contigs = _coords_from_gbk(gbk)
    has_coords = bool(pos) and bool(df["locus_tag"].isin(pos.keys()).mean() > 0.5)

    genes = []
    offsets = {}
    for i, r in enumerate(df.to_dict("records")):
        g = {k: v for k, v in r.items() if v != ""}
        lt = r.get("locus_tag", "")
        if has_coords and lt in pos:
            c, s, e, st = pos[lt]
            g.update(_contig=c, _start=s, _end=e, _strand=st)
        else:
            # no GenBank: genes in input order, length from the protein (x3 nt)
            p = r.get("prophage", "") or "all"
            o = offsets.get(p, 0)
            try:
                L = max(int(float(r.get("aa_length") or 100)) * 3, 60)
            except ValueError:
                L = 300
            g.update(_contig=p, _start=o + 1, _end=o + L, _strand=0)
            offsets[p] = o + L + 40
        g["_i"] = i
        genes.append(g)
    if not has_coords:
        contigs = {k: v for k, v in offsets.items()}

    mode = "genome" if "pharokka_product" in df.columns else "protein"
    meta = {
        "title": title or "phageFACTor report",
        "run": run_label or str(final.parent.parent),
        "table": str(table),
        "mode": mode,
        "generated": _dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "n": len(genes),
        "has_coords": has_coords,
        "has_phynteny": has_phynteny,
        "contigs": [{"id": k, "length": int(v)} for k, v in contigs.items()],
    }
    data = json.dumps({"meta": meta, "genes": genes}, separators=(",", ":"))
    data = data.replace("</", "<\\/")          # never close the <script> early
    html = TEMPLATE.read_text(encoding="utf-8")
    html = html.replace("/*__PF_DATA__*/null", data, 1)
    html = html.replace("__PF_TITLE__", meta["title"].replace("<", "&lt;"), 1)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    print(f"report: {out}  ({len(genes)} genes, coordinates: {'GenBank' if has_coords else 'input order'}, "
          f"phynteny: {'yes' if has_phynteny else 'no'})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run_dir", type=Path, help="run dir with the standard layout")
    ap.add_argument("--final", type=Path, help="final_annotations_table.csv")
    ap.add_argument("--integrated", type=Path, help="05_phynteny/final_annotations_integrated.csv")
    ap.add_argument("--gbk", type=Path, help="updated_prophages.gb (gene coordinates)")
    ap.add_argument("--out", type=Path, help="output HTML (default <run>/04_output/report.html)")
    ap.add_argument("--title", default="", help="report title")
    a = ap.parse_args()

    if a.run_dir:
        r = a.run_dir.resolve()
        a.final = a.final or r / "04_output" / "final_annotations_table.csv"
        a.integrated = a.integrated or r / "05_phynteny" / "final_annotations_integrated.csv"
        a.gbk = a.gbk or r / "04_output" / "updated_prophages.gb"
        a.out = a.out or r / "04_output" / "report.html"
        label = str(r)
    else:
        if not a.final:
            ap.error("give --run_dir or --final")
        a.out = a.out or a.final.parent / "report.html"
        label = str(a.final.resolve().parent.parent)
    if not a.final.is_file():
        raise SystemExit(f"not found: {a.final}")
    build(a.final, a.out, a.gbk, a.integrated, a.title, label)


if __name__ == "__main__":
    main()
