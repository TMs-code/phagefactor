#!/usr/bin/env python3
"""
check_databases.py -- inventory which databases are installed, and what each
one unlocks. Run it before a DB download so you only fetch what is missing.

  phagefactor check            (or: python scripts/tools/check_databases.py)

Paths follow scripts/lib/config.py, so PHAGEFACTOR_DB_ROOT, PHAROKKA_DB,
FOLDSEEK_DB_ROOT and the location recorded by `phagefactor databases` are honoured.
"""
import os
import sys
from pathlib import Path



def human(n):
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.0f}{u}"
        n /= 1024
    return f"{n:.0f}PB"


def dir_size(p: Path, cap_files=4000):
    """Approximate: stop after cap_files entries so a 60 GB DB does not take
    a minute just to be measured."""
    tot = n = 0
    try:
        for f in p.rglob("*"):
            if f.is_file():
                tot += f.stat().st_size
                n += 1
                if n >= cap_files:
                    return tot, True
    except OSError:
        pass
    return tot, False


def main():
    # Import config the way the pipeline does, so this reports the paths the
    # pipeline will ACTUALLY use rather than a second guess at them.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
    import config as C
    # FOLDSEEK_LOCAL_DBS values are <root>/<subdir>/<db>, so the root is 2 up
    fs_root = Path(next(iter(C.FOLDSEEK_LOCAL_DBS.values()))).parents[1]
    pharokka = Path(C.PHAROKKA_DB)
    db_root = fs_root.parent

    print(f"  DB root        : {db_root}")
    print(f"  FoldSeek DBs   : {fs_root}")
    print()

    checks = [
        ("pharokka_db", pharokka, pharokka / "phrog_annot_v4.tsv", "~1-5 GB",
         "genome mode (gene calling + PHROG/CARD/VFDB)", "--pharokka"),
        ("pdb100", fs_root / "pdb100_db", fs_root / "pdb100_db" / "pdb100", "~4 GB",
         "experimental structures", "--pdb100"),
        ("afdb-swissprot", fs_root / "afdb_swissprot_db",
         fs_root / "afdb_swissprot_db" / "afdb_swissprot", "~1.5 GB",
         "curated SwissProt AFDB", "--swissprot"),
        ("afdb50", fs_root / "afdb50_db", fs_root / "afdb50_db" / "afdb50", "~40-60 GB",
         "THE important one: ~91% of FoldSeek-derived calls", "--afdb50"),
        ("taxdump", Path(C.TAXDUMP_DIR), Path(C.TAXDUMP_DIR) / "nodes.dmp", "~250 MB",
         "hit domain without a FoldSeek lineage (webapi mode, older runs)", "--taxdump"),
    ]

    print(f"  {'database':<16}{'status':<12}{'size':>10}   unlocks")
    print(f"  {'-'*16}{'-'*12}{'-'*10}   {'-'*40}")
    missing = []
    for name, d, marker, approx, what, flag in checks:
        present = d.is_dir() and (marker.exists() or any(d.glob(f"{marker.name}*")))
        if present:
            sz, capped = dir_size(d)
            size = f"{'>' if capped else ''}{human(sz)}"
            status = "INSTALLED"
        else:
            size, status = approx, "missing"
            missing.append((name, flag, approx))
        print(f"  {name:<16}{status:<12}{size:>10}   {what}")

    print()
    if not missing:
        print("  All databases present — local search can use every DB.")
        return 0

    print("  MISSING:")
    for name, flag, approx in missing:
        print(f"    {name:<16}{flag:<14}{approx}")

    small = [f for n, f, _ in missing if n in ("pdb100", "afdb-swissprot", "pharokka_db")]
    if small:
        print("\n  Laptop-friendly set (~5.5 GB for the two FoldSeek DBs) — enough for a")
        print("  full reproducible run, though pdb100+swissprot supplied only ~9% of")
        print("  hits in our 17-phage benchmark:")
        print(f"    phagefactor databases {' '.join(small)}")
    if any(n == "afdb50" for n, _, _ in missing):
        print("\n  afdb50 is where the annotations actually come from, but it is 40-60 GB")
        print("  and its index wants a fat-memory node. Cluster only:")
        print("    phagefactor databases --afdb50")
        print('    sbatch --mem=200G --time=04:00:00 --wrap "phagefactor databases --index"')
    print("\n  Restrict a run to the DBs you have:")
    print('    phagefactor run --config foldseek_dbs="afdb-swissprot pdb100"')
    return 0


if __name__ == "__main__":
    sys.exit(main())
