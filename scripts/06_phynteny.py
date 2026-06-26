#!/usr/bin/env python3
"""
06_phynteny.py — Run Phynteny on the final curated GenBank
===========================================================
Runs as the final annotation step (step 06) AFTER the full curation pipeline
(steps 01–05) to predict PHROG functional categories from synteny context.

Why here and not at step 01:
  - Avoids circularity: phold_function_cat is NOT used as input here; instead
    the fully-curated final_product annotations are the phynteny input.
  - Runs once on the combined curated GBK, not per-prophage inside phold array.
  - Preserves phynteny independence from phold's category calls.

What phynteny predicts:
  - PHROG functional category (9 categories: head/packaging, tail, integration,
    transcription regulation, DNA metabolism, lysis, connector, moron/AMG, other)
  - Per-CDS probability (0–1)
  - Threshold: --threshold 0.8 (configurable, recommended 0.8 for low FP rate)

Output columns added to per-gene table:
  phynteny_function_cat   — predicted category (str) or 'unknown' if below threshold
  phynteny_probability    — raw model probability (float 0–1)

Writes (into <output_dir> = 05_phynteny/):
  phynteny_predictions.csv     — per-CDS phynteny output
  phynteny.tsv / phynteny.gbk  — phynteny_transformer per-CDS table + annotated GBK
  run/phynteny_run.log         — phynteny_transformer's verbose run log
  run/phynteny_input.fasta     — the FASTA phynteny_transformer wrote/used

Usage:
  python scripts/06_phynteny.py \\
      --input  04_output/updated_prophages.gb \\
      --output 05_phynteny/ \\
      --threshold 0.8
"""

import sys
import argparse
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).parent
_PROJECT_DIR = _SCRIPTS_DIR.parent
sys.path.insert(0, str(_PROJECT_DIR))

from utils import log, section

try:
    import pandas as pd
except ImportError:
    print("pandas required: pip install pandas")
    sys.exit(1)


def parse_args():
    p = argparse.ArgumentParser(description="Run Phynteny on curated GenBank")
    p.add_argument("--input",     required=True,  help="Curated GenBank (.gb/.gbk)")
    p.add_argument("--output",    required=True,  help="Output directory")
    p.add_argument("--threshold", type=float, default=0.8,
                   help="Phynteny acceptance threshold (default 0.8)")
    return p.parse_args()


def run_phynteny(input_gb: Path, output_dir: Path, threshold: float):
    """
    Run phynteny on the input GenBank.
    Phynteny predicts PHROG functional categories from ±5 gene neighborhood context.
    """
    import shutil, subprocess
    # 2026-06: switched to Phynteny Transformer (susiegriggo/Phynteny_transformer).
    # Console script `phynteny_transformer`. Verified CLI (README):
    #   phynteny_transformer <input.gbk> -o <outdir> [-m <models_dir>]
    # Needs PHROG-annotated GenBank (our step-05 GBK qualifies post-phold) and
    # pre-downloaded models (`install_models`, see 00_install_env.sh Step 5).
    # No 120-gene cap (unlike classic LSTM phynteny). No --min_confidence flag —
    # it emits a phynteny score + confidence per CDS, filtered in
    # parse_phynteny_output(threshold). Outputs:
    #   <outdir>/phynteny_transformer.gbk
    #   <outdir>/phynteny_per_cds_funcions.tsv   (sic — upstream filename typo)
    cli = "phynteny_transformer"
    if shutil.which(cli) is None:
        log(f"ERROR: '{cli}' CLI not found on PATH.")
        log("  Activate the phynteny env and install (note: setuptools is required):")
        log("    micromamba activate phynteny && pip install -U setuptools phynteny_transformer && install_models")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    log(f"Input GenBank: {input_gb}")
    log(f"Output dir:    {output_dir}")
    log(f"Post-hoc confidence filter: {threshold}")

    # -f/--force: overwrite the output dir if it already exists (phynteny_transformer
    # refuses otherwise and exits 1 — the usual cause of a re-run failure).
    cmd = [cli, str(input_gb), "-o", str(output_dir), "-f"]
    # install_models drops weights into <site-packages>/phynteny_utils/models and
    # tells you to pass them with -m. Auto-detect and add it so we don't depend on
    # phynteny finding them implicitly.
    try:
        import importlib.util as _ilu
        _spec = _ilu.find_spec("phynteny_utils")
        if _spec and _spec.submodule_search_locations:
            _models = Path(list(_spec.submodule_search_locations)[0]) / "models"
            if _models.is_dir() and any(_models.iterdir()):
                cmd += ["-m", str(_models)]
                log(f"Models:        {_models}")
    except Exception as _e:
        log(f"  (could not auto-locate phynteny models: {_e}; relying on default)")
    log(f"Running: {' '.join(cmd)}")
    # Capture output so we can echo phynteny's OWN error verbatim on failure
    # (avoids the previous misleading generic hints).
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        log("---- phynteny stdout ----\n" + result.stdout.rstrip())
    if result.stderr:
        log("---- phynteny stderr ----\n" + result.stderr.rstrip())
    if result.returncode != 0:
        log(f"ERROR: {cli} exited with code {result.returncode}. "
            "Read phynteny's stdout/stderr above for the actual cause. Common ones: "
            "output dir exists (now handled by -f); missing models (run install_models); "
            "'No module named pkg_resources' (pip install 'setuptools<82'); "
            "numpy/sklearn ABI clash 'dtype size changed' (pip install 'scikit-learn>=1.4').")
        sys.exit(result.returncode)


def parse_phynteny_output(output_dir: Path, threshold: float) -> pd.DataFrame:
    """
    Parse phynteny output CSV/TSV and filter by threshold.
    Returns DataFrame with columns: locus_tag, phynteny_function_cat, phynteny_probability
    """
    # Phynteny Transformer writes `phynteny_per_cds_funcions.tsv` (sic). Prefer it
    # explicitly, then fall back to any phynteny table for forward/backward compat.
    candidates = list(output_dir.glob("phynteny_per_cds_funcions.tsv")) + \
                 list(output_dir.glob("*per_cds*.tsv")) + \
                 list(output_dir.glob("*phynteny*.tsv")) + \
                 list(output_dir.glob("*phynteny*.csv"))
    # de-dup while preserving order
    seen = set(); candidates = [c for c in candidates if not (c in seen or seen.add(c))]
    if not candidates:
        log("WARNING: No phynteny output CSV found — check phynteny output format")
        return pd.DataFrame(columns=["locus_tag", "phynteny_function_cat", "phynteny_probability"])

    _f = candidates[0]
    df = pd.read_csv(str(_f), sep="\t" if _f.suffix == ".tsv" else ",")
    log(f"Phynteny output ({_f.name}): {len(df)} predictions (before threshold filter)")

    # Phynteny Transformer schema (verified 2026-06):
    #   ID,start,end,strand,phrog_id,phrog_category,phynteny_category,
    #   phynteny_score,phynteny_confidence,phage
    # NOTE: `ID` is POSITIONAL ("phiNP_0"), NOT the GBK locus_tag — so downstream
    # joins to the final table must use (phage,start) coordinates, not locus_tag.
    # Map by EXACT known names (the old substring heuristic collided: "ID"/"phrog_id"
    # both matched "id", "phrog_category"/"phynteny_category" both matched "category").
    def pick(*names):
        for n in names:
            if n in df.columns:
                return n
        return None
    c_id   = pick("ID", "id", "locus_tag", "cds_id")
    c_cat  = pick("phynteny_category", "predicted_category", "category")
    c_conf = pick("phynteny_confidence", "confidence", "phynteny_probability", "prob")
    c_phage = pick("phage", "contig", "record", "prophage")
    if c_cat is None or c_conf is None:
        log(f"WARNING: could not find category/confidence columns. Actual columns: {list(df.columns)}")
        return pd.DataFrame(columns=["phynteny_id", "phynteny_function_cat", "phynteny_probability"])

    out = pd.DataFrame({
        "phynteny_id":           df[c_id] if c_id else range(len(df)),
        "phynteny_function_cat": df[c_cat].astype(str),
        "phynteny_probability":  pd.to_numeric(df[c_conf], errors="coerce").fillna(0.0),
    })
    for extra in ("start", "end", "strand"):
        if extra in df.columns:
            out[extra] = df[extra]
    if c_phage:
        out["phage"] = df[c_phage]

    n_before = len(out)
    out.loc[out["phynteny_probability"] < threshold, "phynteny_function_cat"] = "unknown"
    n_accepted = (out["phynteny_function_cat"] != "unknown").sum()
    log(f"Threshold {threshold}: {n_accepted}/{n_before} predictions accepted "
        "(join to final table on (phage,start) — phynteny IDs are positional, not locus_tags)")
    return out


def _relocate_run_artifacts(output_dir: Path):
    """Move phynteny_transformer's verbose log and the FASTA it writes into a
    05_phynteny/run/ subfolder with clearer, content-describing names, keeping the
    deliverables (phynteny.tsv/.gbk, phynteny_predictions.csv) at the top level.
    The tool writes phynteny.log and phynteny.fasta into the output dir."""
    import shutil
    run_dir = output_dir / "run"
    for pattern, newname in (("phynteny*.log",   "phynteny_run.log"),
                             ("phynteny*.fasta", "phynteny_input.fasta")):
        for src in sorted(output_dir.glob(pattern)):
            run_dir.mkdir(parents=True, exist_ok=True)
            try:
                shutil.move(str(src), str(run_dir / newname))
                log(f"  relocated {src.name} -> run/{newname}")
            except OSError as e:
                log(f"  (could not relocate {src.name}: {e})")


def main():
    args = parse_args()
    input_gb  = Path(args.input)
    output_dir = Path(args.output)

    if not input_gb.exists():
        log(f"ERROR: Input GenBank not found: {input_gb}")
        sys.exit(1)

    section("STEP 06 -- PHYNTENY")

    run_phynteny(input_gb, output_dir, args.threshold)

    section("PARSING PHYNTENY OUTPUT")
    phynteny_df = parse_phynteny_output(output_dir, args.threshold)

    if phynteny_df.empty:
        log("No phynteny predictions to process.")
        return

    pred_csv = output_dir / "phynteny_predictions.csv"
    phynteny_df.to_csv(str(pred_csv), index=False)
    log(f"Predictions saved → {pred_csv}")

    # Summary by category
    section("PHYNTENY SUMMARY")
    cat_counts = phynteny_df["phynteny_function_cat"].value_counts()
    log(f"\nCategory distribution (threshold={args.threshold}):")
    log(cat_counts.to_string())

    n_assigned = (phynteny_df["phynteny_function_cat"] != "unknown").sum()
    log(f"\nTotal assigned (>= {args.threshold}): {n_assigned} / {len(phynteny_df)}")

    # Tidy: move phynteny's verbose log + input FASTA into run/ with clearer names.
    _relocate_run_artifacts(output_dir)

    section("DONE")
    log(f"Output: {output_dir}/  (run/ holds phynteny_run.log + phynteny_input.fasta)")


if __name__ == "__main__":
    main()
