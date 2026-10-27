#!/usr/bin/env python3
"""
esmfold_cache.py -- fold sequences with ESMFold and search them against the
public Foldseek server, caching both.

search.foldseek.com accepts a STRUCTURE (not a sequence), so web-mode queries
must be folded first.

Caching is keyed by sha256 of the SEQUENCE, so an identical protein is folded
and searched once regardless of gene name, phage or run.

A search that completes with no rows is cached as an empty .m8 -- a real result,
never retried. Failures (ERROR/timeout/429) cache nothing and retry next run.

ESMFold caps at ~400 aa, which excludes large structural proteins (tape measure,
tail fibre, portal, terminase large). Those are reported by name, never dropped
silently.
"""
import gzip
import hashlib
import io
import os
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    import requests
except ImportError:
    raise SystemExit("pip install requests")

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
FOLDSEEK_API = "https://search.foldseek.com/api"

ESMFOLD_MAX_LEN = int(os.environ.get("ESMFOLD_MAX_LEN", "400"))
ESMFOLD_DELAY = float(os.environ.get("ESMFOLD_DELAY", "1.5"))   # politeness
WEB_WORKERS = int(os.environ.get("PHAGEFACTOR_WEB_WORKERS", "4"))  # proteins in flight
SUBMIT_DELAY = float(os.environ.get("FOLDSEEK_SUBMIT_DELAY", "4"))
POLL_SECS = 5
MAX_WAIT = int(os.environ.get("FOLDSEEK_MAX_WAIT", "300"))
RATE_BACKOFF = [30, 60, 120, 240, 300, 300]     # 429 schedule
ESMFOLD_BACKOFF = [5, 15, 30]                   # 5xx/timeout schedule
ERROR_RETRIES = int(os.environ.get("FOLDSEEK_ERROR_RETRIES", "2"))


def log(msg, end="\n"):
    print(msg, end=end, flush=True)


def _key(seq: str) -> str:
    return hashlib.sha256(seq.strip().upper().encode()).hexdigest()[:24]


def _write_atomic(path: Path, text: str) -> None:
    """Write via a temporary file + rename: jobs sharing the cache never read a
    partial file. A .gz path is written gzip-compressed."""
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    data = text.encode()
    tmp.write_bytes(gzip.compress(data) if path.suffix == ".gz" else data)
    os.replace(tmp, path)


def cache_root() -> Path:
    """Persistent across runs: structures do not depend on a run."""
    env = os.environ.get("PHAGEFACTOR_WEB_CACHE")
    if env:
        return Path(env)
    db_root = os.environ.get("PHAGEFACTOR_DB_ROOT")
    if db_root:
        return Path(db_root) / "web_cache"
    # the resolved database root (project config > ~/.config/phagefactor
    # settings > <project>/databases), never the code directory
    from config import DB_ROOT
    return Path(DB_ROOT) / "web_cache"


# =============================================================================
# ESMFold
# =============================================================================

def fold(seq: str, root: Optional[Path] = None) -> Tuple[Optional[str], str]:
    """Return (pdb_text, status). status in {cached, folded, too_long, failed}."""
    root = root or cache_root()
    sdir = root / "structures"
    sdir.mkdir(parents=True, exist_ok=True)
    pdb_path = sdir / f"{_key(seq)}.pdb"

    if pdb_path.exists() and pdb_path.stat().st_size > 100:
        return pdb_path.read_text(), "cached"
    if len(seq) > ESMFOLD_MAX_LEN:
        return None, "too_long"

    # 5xx are transient. HTTP 504 is the API gateway giving up at ~29s while
    # ESMFold is still folding, so it correlates with sequence LENGTH. Retry
    # with backoff; cache nothing on failure so the next run tries again.
    for attempt, wait in enumerate(ESMFOLD_BACKOFF, 1):
        try:
            r = requests.post(
                ESMFOLD_URL, data=seq,
                headers={"Content-Type": "application/x-www-form-urlencoded"},
                timeout=180,
            )
        except Exception as e:
            log(f"      ESMFold EXC {type(e).__name__}: {e} (len={len(seq)}) "
                f"[{attempt}/{len(ESMFOLD_BACKOFF)}]")
            time.sleep(wait)
            continue
        finally:
            time.sleep(ESMFOLD_DELAY)

        if r.status_code == 200 and len(r.text) > 100 and "ATOM" in r.text:
            _write_atomic(pdb_path, r.text)
            return r.text, "folded"
        if r.status_code in (429, 500, 502, 503, 504):
            log(f"      ESMFold HTTP {r.status_code} (len={len(seq)} aa) — "
                f"retry in {wait}s [{attempt}/{len(ESMFOLD_BACKOFF)}]")
            time.sleep(wait)
            continue
        log(f"      ESMFold HTTP {r.status_code}: {r.text[:100]}")
        return None, "failed"
    return None, "failed"


# =============================================================================
# Foldseek web search
# =============================================================================

def _submit(pdb_text: str, dbs: List[str], label: str) -> str:
    for attempt, wait in enumerate(RATE_BACKOFF, 1):
        r = requests.post(
            f"{FOLDSEEK_API}/ticket",
            files={"q": (f"{label}.pdb", pdb_text.encode(), "application/octet-stream")},
            data=[("mode", "3diaa")] + [("database[]", d) for d in dbs],
            timeout=90,
        )
        if r.status_code == 429:
            log(f"      429 rate-limited; waiting {wait}s ({attempt}/{len(RATE_BACKOFF)})")
            time.sleep(wait)
            continue
        r.raise_for_status()
        return r.json()["id"]
    raise RuntimeError(f"rate-limited too long for {label}")


def _poll(ticket: str) -> str:
    deadline = time.time() + MAX_WAIT
    while time.time() < deadline:
        st = requests.get(f"{FOLDSEEK_API}/ticket/{ticket}", timeout=30).json().get("status", "")
        if st in ("COMPLETE", "ERROR"):
            return st
        time.sleep(POLL_SECS)
    return "TIMEOUT"


# qaln, taln, tcoords, tseq: unused downstream and most of the m8 size; kept as
# empty columns so the layout stays the server's
_DROP_COLS = (15, 16, 17, 18)


def _slim_m8(text: str) -> str:
    out = []
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) <= max(_DROP_COLS):
            continue        # the server's taxonomy report rows, not hits
        for i in _DROP_COLS:
            f[i] = ""
        out.append("\t".join(f))
    return "\n".join(out) + ("\n" if out else "")


def _download_m8(ticket: str) -> str:
    r = requests.get(f"{FOLDSEEK_API}/result/download/{ticket}", timeout=180, stream=True)
    r.raise_for_status()
    out = []
    with tarfile.open(fileobj=io.BytesIO(r.content)) as tar:
        for m in tar.getmembers():
            if m.name.endswith(".m8") or "alis" in m.name:
                f = tar.extractfile(m)
                if f:
                    out.append(f.read().decode("utf-8", errors="replace"))
    return "".join(out)


def search(seq: str, pdb_text: str, dbs: List[str], label: str,
           root: Optional[Path] = None) -> Tuple[Optional[str], str]:
    """Return (m8_text, status). status in {cached, searched, empty, failed}.

    An EMPTY cached file means "searched, no homologs" -- a real answer, never
    retried. Nothing is written on ERROR/TIMEOUT so those retry next run.
    """
    root = root or cache_root()
    # results depend on WHICH DBs were searched, so the DB set is part of the key
    dbkey = hashlib.sha256(",".join(sorted(dbs)).encode()).hexdigest()[:8]
    mdir = root / "searches" / dbkey
    mdir.mkdir(parents=True, exist_ok=True)
    m8_path = mdir / f"{_key(seq)}.m8.gz"
    plain = mdir / f"{_key(seq)}.m8"           # uncompressed entries are still read

    if m8_path.exists() or plain.exists():
        txt = (gzip.decompress(m8_path.read_bytes()).decode() if m8_path.exists()
               else plain.read_text())
        return txt, ("empty" if not txt.strip() else "cached")

    for attempt in range(1, ERROR_RETRIES + 2):
        try:
            ticket = _submit(pdb_text, dbs, label)
            st = _poll(ticket)
        except Exception as e:
            log(f"      search EXC {type(e).__name__}: {e}")
            st = "EXC"
        if st == "COMPLETE":
            try:
                m8 = _slim_m8(_download_m8(ticket))
            except Exception as e:
                log(f"      download EXC {type(e).__name__}: {e}")
                return None, "failed"
            _write_atomic(m8_path, m8)      # empty file == "no hits", cached
            time.sleep(SUBMIT_DELAY)
            return m8, ("empty" if not m8.strip() else "searched")
        # ERROR / TIMEOUT / EXC are transient -> retry, cache NOTHING
        if attempt <= ERROR_RETRIES:
            log(f"      {st} for {label} — retry {attempt}/{ERROR_RETRIES}")
            time.sleep(10 * attempt)
    return None, "failed"


# =============================================================================
# Driver
# =============================================================================

# A 26-aa peptide: short enough that the backend cannot plausibly time out, so
# a failure here means the SERVICE is down rather than the sequence being hard.
_PROBE_SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERL"


def preflight(root: Optional[Path] = None) -> bool:
    """Fold a tiny peptide to establish the endpoint is alive before spending
    hours on a real dataset. Returns True if ESMFold answered."""
    log("  [preflight] testing ESMFold with a 27 aa peptide…")
    pdb, st = fold(_PROBE_SEQ, root)
    if pdb:
        log(f"  [preflight] ESMFold OK ({st})")
        return True
    log("  [preflight] ESMFold did NOT answer for a trivially short peptide.")
    return False


def fold_and_search(pairs: Iterable[Tuple[str, str]], dbs: List[str],
                    root: Optional[Path] = None) -> Tuple[str, Dict[str, list]]:
    """Fold + search every (name, seq). Returns (combined_m8, report).

    The m8 query column is rewritten to the GENE NAME: each protein is folded
    and submitted separately, so the server's query id is our cache key, not
    something downstream can join on.
    """
    root = root or cache_root()
    pairs = list(pairs)
    rep: Dict[str, list] = {"too_long": [], "fold_failed": [], "search_failed": [],
                            "no_hits": [], "ok": [], "cached": []}
    chunks: List[str] = []
    log(f"  cache: {root}")
    log(f"  {len(pairs)} protein(s); ESMFold cap {ESMFOLD_MAX_LEN} aa; DBs {dbs}")
    if not preflight(root):
        raise SystemExit(
            "\n[02w] The ESMFold API (api.esmatlas.com) is not answering.\n"
            "  It could not fold a 27 aa peptide, so this is the SERVICE, not your\n"
            "  Web search needs structures, so without ESMFold it cannot run.\n"
            "\n  Use instead:\n"
            "   1. LOCAL search, small DBs (~5.5 GB, no folding -- phold already\n"
            "      wrote the 3Di):\n"
            "        phagefactor databases --pdb100 --swissprot\n"
            "        phagefactor run --config search_mode=local \\\n"
            "            foldseek_dbs=\"afdb-swissprot pdb100\"   (pharokka/phold are kept)\n"
            "   2. the cluster, where afdb50 is on disk\n"
            "   3. search_mode=none -- curation on pharokka + phold only\n"
            "\n  Retry web mode later; the structure cache means nothing is lost.\n")

    def _one(seq):
        pdb_text, fstat = fold(seq, root)
        if pdb_text is None:
            return fstat, None, None
        m8, sstat = search(seq, pdb_text, dbs, _key(seq), root)
        return fstat, m8, sstat

    # a few proteins in flight at once; results are handled in input order
    with ThreadPoolExecutor(max_workers=WEB_WORKERS) as pool:
        results = pool.map(_one, [seq for _, seq in pairs])
        for i, ((name, seq), (fstat, m8, sstat)) in enumerate(zip(pairs, results), 1):
            if fstat == "too_long":
                rep["too_long"].append((name, len(seq)))
                log(f"  [{i}/{len(pairs)}] {name}: {len(seq)} aa > {ESMFOLD_MAX_LEN} — SKIPPED")
                continue
            if fstat == "failed":
                rep["fold_failed"].append(name)
                log(f"  [{i}/{len(pairs)}] {name}: fold FAILED")
                continue
            if m8 is None:
                rep["search_failed"].append(name)
                log(f"  [{i}/{len(pairs)}] {name}: search FAILED (will retry next run)")
                continue
            if sstat == "empty":
                rep["no_hits"].append(name)
            elif sstat == "cached":
                rep["cached"].append(name)
            else:
                rep["ok"].append(name)

            # rewrite query id -> gene name
            for line in m8.splitlines():
                if line.strip():
                    parts = line.split("\t")
                    parts[0] = name
                    chunks.append("\t".join(parts))
            log(f"  [{i}/{len(pairs)}] {name}: fold={fstat} search={sstat} "
                f"({len([l for l in m8.splitlines() if l.strip()])} hits)")

    return "\n".join(chunks) + ("\n" if chunks else ""), rep


def report(rep: Dict[str, list]) -> None:
    """Print the summary. Long proteins are named, not just counted -- they are
    a systematic gap in the evidence, not noise."""
    log("\n  --- web search summary ---")
    log(f"    searched now      : {len(rep['ok'])}")
    log(f"    reused from cache : {len(rep['cached'])}")
    log(f"    no homologs (real): {len(rep['no_hits'])}")
    if rep["fold_failed"]:
        log(f"    fold FAILED       : {len(rep['fold_failed'])} (retried next run)")
    if rep["search_failed"]:
        log(f"    search FAILED     : {len(rep['search_failed'])} (retried next run)")
    if rep["too_long"]:
        log(f"\n    !! {len(rep['too_long'])} protein(s) EXCEED the {ESMFOLD_MAX_LEN} aa "
            f"ESMFold limit and have NO structural evidence.")
        log("       These are typically the large structural proteins (tape measure,")
        log("       tail fibre, portal, terminase large) — i.e. a biased gap, not a")
        log("       random one. Fold them elsewhere (ColabFold/AF3) and drop the .pdb")
        log(f"       into {cache_root()/'structures'}/<sha256-24>.pdb to include them.")
        for nm, ln in sorted(rep["too_long"], key=lambda x: -x[1]):
            log(f"         {nm}  {ln} aa")
