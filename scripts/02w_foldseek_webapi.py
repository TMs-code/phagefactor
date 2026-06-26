#!/usr/bin/env python3
"""
02w_foldseek_webapi.py — FoldSeek search via the public web server (NO local DB)
================================================================================
Optional alternative to the local 02d step for users who don't have the disk/RAM
for the afdb50 DB + index. Submits the hypothetical-protein AMINO-ACID FASTA to
the FoldSeek web server (https://search.foldseek.com), which runs ProstT5 itself
to get 3Di and searches PDB/AFDB — so no local database is required.

Adapted from an earlier 02b_foldseek_pipeline.py (same ticket/
poll/backoff machinery), with the needed fix: on a per-job TIMEOUT we wait and
RESUBMIT rather than failing, since the public server queue is the bottleneck.

Output: one concatenated `webapi_hits.m8` (FoldSeek BTAB/m8). Feed it into the
existing parse+score path (the scoring/confidence/phage-boost logic in
foldseek_scoring.py / 02d STEP 5) so downstream 03/04/05 are unchanged.

  python 02w_foldseek_webapi.py --fasta input/all_proteins_combined.faa \
      --out 02_foldseek/3di_tokens/webapi_hits.m8

⚠️ The public FoldSeek server API params can change — verify `mode`/`database[]`
   against https://search.foldseek.com if submissions start failing.
"""

import io
import sys
import time
import tarfile
import argparse
from pathlib import Path

try:
    import requests
except ImportError:
    raise SystemExit("pip install requests")

API_URL = "https://search.foldseek.com/api"
DATABASES = ["pdb100", "afdb50", "afdb-swissprot"]
POLL_SECS = 5
JOB_TIMEOUT = 300          # seconds before we treat a single job as stuck
MAX_RESUBMITS = 4          # resubmit a timed-out/queued job this many times
RATE_BACKOFF = [30, 60, 120, 240, 300, 300]   # 429 backoff schedule


def log(msg, end="\n"):
    print(msg, end=end, flush=True)


def _read_fasta(path):
    name, seq, out = None, [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        out.append((name, "".join(seq)))
    return out


def _submit(name, seq):
    fasta = f">{name}\n{seq}\n".encode()
    for attempt, wait in enumerate(RATE_BACKOFF, 1):
        r = requests.post(
            f"{API_URL}/ticket",
            files={"q": (f"{name}.fasta", fasta, "application/octet-stream")},
            data=[("mode", "3diaa")] + [("database[]", db) for db in DATABASES],
            timeout=60,
        )
        if r.status_code == 429:
            log(f"  429 rate-limited; wait {wait}s ({attempt}/{len(RATE_BACKOFF)})")
            time.sleep(wait)
            continue
        r.raise_for_status()
        return r.json()["id"]
    raise RuntimeError(f"rate-limited too long for {name}")


def _poll_and_download(name, ticket):
    deadline = time.time() + JOB_TIMEOUT
    while time.time() < deadline:
        st = requests.get(f"{API_URL}/ticket/{ticket}", timeout=30).json().get("status", "")
        if st == "COMPLETE":
            r = requests.get(f"{API_URL}/result/download/{ticket}", timeout=120)
            with tarfile.open(fileobj=io.BytesIO(r.content)) as tar:
                m8 = []
                for mem in tar.getmembers():
                    if mem.name.endswith(".m8"):
                        f = tar.extractfile(mem)
                        if f:
                            m8.append(f.read().decode("utf-8", "replace"))
                return "".join(m8)
        if st == "ERROR":
            raise RuntimeError(f"server ERROR for {name}")
        time.sleep(POLL_SECS)
    raise TimeoutError(name)


def search_one(name, seq):
    """Submit + poll, resubmitting on timeout (the public-queue pain point)."""
    for attempt in range(1, MAX_RESUBMITS + 1):
        try:
            ticket = _submit(name, seq)
            return _poll_and_download(name, ticket)
        except TimeoutError:
            wait = 60 * attempt
            log(f"  [{name}] timeout; resubmit {attempt}/{MAX_RESUBMITS} after {wait}s")
            time.sleep(wait)
    log(f"  [{name}] FAILED after {MAX_RESUBMITS} resubmits — skipping")
    return ""


# =============================================================================
# Bridge: web m8 -> best_hit.csv / top3.csv  (closes the loop to step 03)
# =============================================================================
# Web-server m8 column order (M8_COLS in utils.py, 21 cols):
_WEB_M8_COLS = ["query", "target_raw", "pident", "alnlen", "mismatch", "gapopen",
                "qstart", "qend", "tstart", "tend", "prob", "evalue_raw",
                "score_raw", "lddtfull", "qcov_aa", "qaln", "taln", "tcoords",
                "tseq", "taxid", "taxname"]


def webm8_to_best_hit(m8_path, best_csv, top3_csv, all_csv=None):
    """Parse the web-API m8 and reuse the SAME scoring path as local mode
    (foldseek_scoring.build_best_and_top3) so step 03 is identical for both
    search modes. ⚠ Verify the server's m8 column order against _WEB_M8_COLS if
    descriptions/taxa look wrong — the public server's format can drift."""
    import sys as _sys
    from pathlib import Path as _P
    _sys.path.insert(0, str(_P(__file__).parent))
    import pandas as pd
    from foldseek_scoring import (
        _is_informative_fs, _phage_boost_factor, _is_same_host_hit,
        _is_promiscuous_fold_hit, _has_eukaryotic_description, build_best_and_top3,
    )
    try:
        from config import FOLDSEEK_EVALUE_MAX, FOLDSEEK_SCORE_OVERRIDE
    except Exception:
        FOLDSEEK_EVALUE_MAX, FOLDSEEK_SCORE_OVERRIDE = 0.1, 200

    rows = [l.split("\t") for l in _P(m8_path).read_text().splitlines()
            if l and not l.startswith("#")]
    if not rows:
        log("No rows in m8 — nothing to bridge."); return
    df = pd.DataFrame([r[:len(_WEB_M8_COLS)] + [""] * (len(_WEB_M8_COLS) - len(r))
                       for r in rows], columns=_WEB_M8_COLS)
    # map to the schema build_best_and_top3 expects
    df["gene"]        = df["query"]
    df["accession"]   = df["target_raw"].str.split(" ", n=1).str[0]
    df["description"] = df["target_raw"].str.split(" ", n=1).str[1].fillna("").str.strip()
    for c in ("pident", "evalue_raw", "score_raw", "qcov_aa"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.rename(columns={"evalue_raw": "evalue", "score_raw": "score"})
    quality_ok = (df["evalue"].fillna(999) < FOLDSEEK_EVALUE_MAX) | \
                 (df["score"].fillna(0) >= FOLDSEEK_SCORE_OVERRIDE)
    df["informative"]           = df["description"].apply(_is_informative_fs) & quality_ok
    df["phage_boost"]           = df["description"].apply(_phage_boost_factor)
    df["composite_score"]       = df["score"].fillna(0) * df["phage_boost"]
    df["same_host"]             = df["taxname"].apply(_is_same_host_hit)
    df["defense_flag"]          = False
    df["promiscuous_fold_flag"] = df["description"].apply(lambda d: _is_promiscuous_fold_hit(str(d)))
    df["eukaryotic_desc_flag"]  = df["description"].apply(lambda d: _has_eukaryotic_description(str(d)))
    df = df.sort_values(["gene", "composite_score", "evalue"], ascending=[True, False, True])
    best, top3 = build_best_and_top3(df, sorted(df["gene"].unique()))
    best.to_csv(best_csv, index=False)
    top3.to_csv(top3_csv, index=False)
    if all_csv:
        df.to_csv(all_csv, index=False)
    log(f"Bridged web m8 -> {best_csv} ({len(best)} genes) + {top3_csv}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", help="hypothetical-protein AA FASTA (search step)")
    ap.add_argument("--out", required=True, help="concatenated web m8 path")
    ap.add_argument("--from-m8", action="store_true",
                    help="skip the search; only bridge an existing m8 at --out")
    ap.add_argument("--best-hit", help="write best_hit.csv here (enables the bridge)")
    ap.add_argument("--top3", help="write top3.csv here (with --best-hit)")
    ap.add_argument("--all-hits", help="optional all_hits.csv")
    a = ap.parse_args()
    out = Path(a.out)

    if not a.from_m8:
        if not a.fasta:
            ap.error("--fasta is required unless --from-m8 is given")
        seqs = _read_fasta(a.fasta)
        log(f"Submitting {len(seqs)} proteins to {API_URL} (DBs: {', '.join(DATABASES)})")
        out.parent.mkdir(parents=True, exist_ok=True)
        n_ok = 0
        with open(out, "w") as fh:
            for i, (name, seq) in enumerate(seqs, 1):
                log(f"[{i}/{len(seqs)}] {name}", end=" ")
                m8 = search_one(name, seq)
                if m8:
                    fh.write(m8); n_ok += 1; log("ok")
        log(f"Search done: {n_ok}/{len(seqs)} proteins returned hits -> {out}")

    # Bridge to best_hit.csv/top3.csv so step 03 sees the usual inputs.
    if a.best_hit:
        webm8_to_best_hit(out, a.best_hit, a.top3 or str(out.with_name("top3.csv")),
                          a.all_hits)
    else:
        log("No --best-hit given; m8 written but not bridged. Re-run with "
            "--from-m8 --best-hit <path> to produce best_hit.csv for step 03.")


if __name__ == "__main__":
    main()
