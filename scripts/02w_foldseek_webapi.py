#!/usr/bin/env python3
"""
02w_foldseek_webapi.py — FoldSeek search via the public web server (no local DB)

STATUS: on hold. search.foldseek.com accepts a STRUCTURE, so queries must be
folded first, and the ESMFold API (api.esmatlas.com) is currently unavailable.
Use SEARCH_MODE=local instead.

Folds with scripts/lib/esmfold_cache.py (sequence-keyed cache), searches each
structure, and writes one concatenated m8 that feeds the same scoring path as
local mode, so 03/04/05 are unchanged.

ESMFold caps at ~400 aa; longer proteins are reported by name and get no
structural evidence.

  python 02w_foldseek_webapi.py --fasta <aa.faa> --out <webapi_hits.m8>
"""

import io
import re
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
MAX_RESUBMITS = 4          # resubmit a TIMED-OUT/queued job this many times
RATE_BACKOFF = [30, 60, 120, 240, 300, 300]   # 429 backoff schedule


class ServerJobError(RuntimeError):
    """Server returned status=ERROR. This is DETERMINISTIC (the server caches by
    sequence, so resubmitting the identical query just replays the cached error)
    -> do NOT hammer with retries; clean the sequence or skip the protein."""


_STD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _clean_seq(seq: str):
    """Return (clean_seq, note). ProstT5 on the FoldSeek server rejects stop '*',
    gaps and some non-standard residues -> a single such char makes the WHOLE job
    ERROR. Strip '*'/'-'/whitespace, map rare residues (UZOBJ) to X, uppercase.
    The note names EXACTLY what was changed, so the cause is pinpointed in the log."""
    s = re.sub(r"\s", "", seq).upper()
    from collections import Counter
    stripped = Counter(ch for ch in s if ch in "*-")
    s = s.replace("*", "").replace("-", "")
    rare = Counter(ch for ch in s if ch in "UZOBJ")
    s = re.sub(r"[UZOBJ]", "X", s)
    bad = Counter(ch for ch in s if ch not in _STD_AA and ch != "X")
    if bad:
        s = re.sub("[" + re.escape("".join(bad)) + "]", "X", s)
    parts = []
    if stripped:
        parts.append("removed " + ", ".join(f"{n}×'{c}'" for c, n in stripped.items()))
    for label, cnt in (("mapped→X", rare), ("mapped→X", bad)):
        if cnt:
            parts.append(label + " " + ", ".join(f"{n}×{c}" for c, n in cnt.items()))
    return s, "; ".join(parts)


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


def _submit(name, seq, dbs=None):
    dbs = dbs or DATABASES
    fasta = f">{name}\n{seq}\n".encode()
    for attempt, wait in enumerate(RATE_BACKOFF, 1):
        r = requests.post(
            f"{API_URL}/ticket",
            files={"q": (f"{name}.fasta", fasta, "application/octet-stream")},
            data=[("mode", "3diaa")] + [("database[]", db) for db in dbs],
            timeout=60,
        )
        if r.status_code == 429:
            log(f"  429 rate-limited; wait {wait}s ({attempt}/{len(RATE_BACKOFF)})")
            time.sleep(wait)
            continue
        r.raise_for_status()
        return r.json()["id"]
    raise RuntimeError(f"rate-limited too long for {name}")


def _submit_fasta(fasta_text, dbs, label):
    """Submit a (possibly multi-sequence) FASTA as ONE FoldSeek job. The public
    server accepts many query sequences per ticket, so 37 proteins = 1 job, not 37
    — vastly faster and kinder to the queue than one job per protein."""
    dbs = dbs or DATABASES
    data = fasta_text.encode()
    for attempt, wait in enumerate(RATE_BACKOFF, 1):
        r = requests.post(
            f"{API_URL}/ticket",
            files={"q": ("query.fasta", data, "application/octet-stream")},
            data=[("mode", "3diaa")] + [("database[]", db) for db in dbs],
            timeout=120,
        )
        if r.status_code == 429:
            log(f"  429 rate-limited; wait {wait}s ({attempt}/{len(RATE_BACKOFF)})")
            time.sleep(wait); continue
        r.raise_for_status()
        return r.json()["id"]
    raise RuntimeError(f"rate-limited too long for {label}")


def _poll_and_download(label, ticket, timeout=None):
    timeout = timeout or JOB_TIMEOUT
    start = time.time(); deadline = start + timeout; nextlog = start + 30; last = None
    while time.time() < deadline:
        last = requests.get(f"{API_URL}/ticket/{ticket}", timeout=30).json()
        st = last.get("status", "")
        if st == "COMPLETE":
            r = requests.get(f"{API_URL}/result/download/{ticket}", timeout=180)
            with tarfile.open(fileobj=io.BytesIO(r.content)) as tar:
                m8 = []
                for mem in tar.getmembers():
                    if mem.name.endswith(".m8"):
                        f = tar.extractfile(mem)
                        if f:
                            m8.append(f.read().decode("utf-8", "replace"))
                return "".join(m8)
        if st == "ERROR":
            raise ServerJobError(f"server ERROR for {label} (ticket={ticket}); response={last}")
        # heartbeat so a slow public-queue job doesn't look hung
        if time.time() >= nextlog:
            log(f"  [{label}] still running… {int(time.time()-start)}s elapsed (status={st or 'PENDING'})")
            nextlog += 30
        time.sleep(POLL_SECS)
    raise TimeoutError(f"{label} (ticket={ticket}); last status={last}")


def search_batch(pairs, dbs=None):
    """Search a list of (name, seq) as ONE job. Returns (m8_text, [(name,err),...]).
    On a deterministic server ERROR, BISECTS the batch to isolate the offending
    sequence so one bad protein can't sink the other 36. Transient errors
    (timeout/HTTP/conn) retry the whole batch with backoff, then bisect."""
    if not pairs:
        return "", []
    fasta = "".join(f">{n}\n{s}\n" for n, s in pairs)
    label = (pairs[0][0] if len(pairs) == 1
             else f"{len(pairs)} seqs [{pairs[0][0]} … {pairs[-1][0]}]")
    timeout = min(JOB_TIMEOUT + 15 * len(pairs), 2400)
    last_err = None
    for attempt in range(1, MAX_RESUBMITS + 1):
        try:
            ticket = _submit_fasta(fasta, dbs, label)
            return _poll_and_download(label, ticket, timeout), []
        except ServerJobError as e:
            if len(pairs) == 1:
                log(f"  [{label}] server ERROR — skipping this protein. {e}")
                return "", [(pairs[0][0], f"ServerJobError: {e}")]
            mid = len(pairs) // 2
            log(f"  [{label}] batch ERROR — bisecting to isolate the bad sequence "
                f"({len(pairs)} -> {mid} + {len(pairs)-mid})")
            m8a, fa = search_batch(pairs[:mid], dbs)
            m8b, fb = search_batch(pairs[mid:], dbs)
            return m8a + m8b, fa + fb
        except TimeoutError as e:
            last_err = f"timeout ({e})"
        except Exception as e:                      # HTTP / conn / JSON / tar
            last_err = f"{type(e).__name__}: {e}"
        wait = 30 * attempt
        log(f"  [{label}] attempt {attempt}/{MAX_RESUBMITS} failed — {last_err}; retry in {wait}s")
        time.sleep(wait)
    # transient errors exhausted -> bisect (or give up on a single protein)
    if len(pairs) == 1:
        return "", [(pairs[0][0], last_err or "failed after retries")]
    mid = len(pairs) // 2
    log(f"  [{label}] still failing — bisecting {len(pairs)} -> {mid} + {len(pairs)-mid}")
    m8a, fa = search_batch(pairs[:mid], dbs)
    m8b, fb = search_batch(pairs[mid:], dbs)
    return m8a + m8b, fa + fb


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
        _is_informative_fs, _phage_boost_factor, _is_same_host_hit, _host_boost_factor,
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
    df["same_host"]             = df["taxname"].apply(_is_same_host_hit)
    # host-aware ranking tie-breaker (optional; foldseek.host_boost, 1.0 disables)
    df["host_boost"]            = df.apply(
        lambda r: _host_boost_factor(r["description"], r["same_host"]), axis=1)
    df["composite_score"]       = (df["score"].fillna(0)
                                   * df["phage_boost"] * df["host_boost"])
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
    ap.add_argument("--strict", action="store_true",
                    help="abort on the first failed protein (old behaviour). Default: "
                         "skip failures and keep going.")
    ap.add_argument("--dbs", default=",".join(DATABASES),
                    help="comma-separated FoldSeek server DBs (default: %(default)s). "
                         "Use to test a single DB, e.g. --dbs afdb50")
    ap.add_argument("--chunk-size", type=int, default=50,
                    help="proteins submitted per FoldSeek job (default 50). One job "
                         "for a whole prophage is far faster than one job per protein.")
    ap.add_argument("--no-canary", action="store_true",
                    help="skip the known-good-protein probe that aborts early on a "
                         "systematic server failure.")
    a = ap.parse_args()
    out = Path(a.out)
    dbs = [d.strip() for d in a.dbs.split(",") if d.strip()]

    if not a.from_m8:
        if not a.fasta:
            ap.error("--fasta is required unless --from-m8 is given")
        seqs = _read_fasta(a.fasta)
        out.parent.mkdir(parents=True, exist_ok=True)
        fail_log = out.with_name("webapi_failures.log")

        log(f"Submitting {len(seqs)} proteins to {API_URL} (DBs: {', '.join(dbs)})")

        # --- preflight: reachability + VALIDATE the requested DB names ----------
        # A single invalid database[] value makes the server ERROR the whole job
        # (the symptom you hit). So we fetch the server's DB list and check.
        try:
            pr = requests.get(f"{API_URL}/databases/all", timeout=30)
            log(f"  [preflight] API reachable (GET /databases/all -> HTTP {pr.status_code})")
            try:
                data = pr.json()
                items = data.get("databases", data) if isinstance(data, dict) else data
                avail = {d.get("path") for d in items if isinstance(d, dict)}
                avail.discard(None)
                if avail:
                    log(f"  [preflight] server DBs: {sorted(avail)}")
                    bad = [d for d in dbs if d not in avail]
                    if bad:
                        log(f"  [preflight][WARN] requested DB(s) NOT offered by server: {bad} "
                            f"-> this alone makes every job ERROR. Fix --dbs to use the names above.")
            except Exception as e:
                log(f"  [preflight] could not parse DB list ({e}); continuing.")
        except Exception as e:
            log(f"  [preflight][WARN] could not reach {API_URL} ({type(e).__name__}: {e}). "
                f"On Maestro this MUST run on the LOGIN node (compute nodes have no internet).")

        # --- clean sequences up front; drop too-short ones with a reason -------
        pairs, failures = [], []
        for name, seq in seqs:
            cseq, note = _clean_seq(seq)
            if len(cseq) < 10:
                failures.append((name, f"too short after cleaning ({len(cseq)} aa)"))
                log(f"  [{name}] SKIP — too short ({len(cseq)} aa)")
            else:
                pairs.append((name, cseq))
                if note:
                    log(f"  [{name}] {note}")

        # --- submit in CHUNKS (one job per chunk), not one job per protein -----
        fl = open(fail_log, "w"); fl.write("protein\terror\n")
        for nm, er in failures:
            fl.write(f"{nm}\t{er}\n")
        fl.flush()
        # --- fold, then search ------------------------------------------
        # The server takes a structure, so fold first. Both steps are cached by
        # sequence in esmfold_cache, so re-runs and repeats cost nothing.
        from esmfold_cache import fold_and_search, report as _web_report

        m8_all, rep = fold_and_search(pairs, dbs)
        with open(out, "w") as fh:
            fh.write(m8_all)
        _web_report(rep)

        for nm in rep["fold_failed"]:
            failures.append((nm, "ESMFold failed")); fl.write(f"{nm}\tESMFold failed\n")
        for nm in rep["search_failed"]:
            failures.append((nm, "FoldSeek search failed")); fl.write(f"{nm}\tFoldSeek search failed\n")
        for nm, ln in rep["too_long"]:
            failures.append((nm, f"too long for ESMFold ({ln} aa)"))
            fl.write(f"{nm}\ttoo long for ESMFold ({ln} aa)\n")
        fl.flush()
        if a.strict and failures:
            fl.close(); raise SystemExit(f"--strict: {len(failures)} failure(s)")
        fl.close()

        # count distinct queries that returned hits
        n_q = len({l.split("\t", 1)[0] for l in Path(out).read_text().splitlines() if l.strip()})
        log(f"\nSearch done: {n_q}/{len(seqs)} proteins returned hits -> {out}")
        if failures:
            log(f"  {len(failures)} protein(s) failed/skipped — details -> {fail_log}")
            log("  first few: " + "; ".join(f"{n} [{e[:70]}]" for n, e in failures[:3]))
            if n_q == 0:
                log("  ALL proteins failed. If the preflight flagged a bad DB name, fix --dbs. "
                    "Else test one DB on one protein: --dbs afdb50 --chunk-size 1")

    # Bridge to best_hit.csv/top3.csv so step 03 sees the usual inputs.
    if a.best_hit:
        webm8_to_best_hit(out, a.best_hit, a.top3 or str(out.with_name("top3.csv")),
                          a.all_hits)
    else:
        log("No --best-hit given; m8 written but not bridged. Re-run with "
            "--from-m8 --best-hit <path> to produce best_hit.csv for step 03.")


if __name__ == "__main__":
    main()
