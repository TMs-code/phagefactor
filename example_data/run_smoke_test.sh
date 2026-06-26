#!/usr/bin/env bash
# =============================================================================
# run_smoke_test.sh — end-to-end check on the genome_2campy example
#
# Runs the pipeline on the 2 Campylobacter prophages and verifies the output
# table is produced with the expected schema and a plausible annotated fraction.
# Uses WebAPI search mode by default so NO local FoldSeek database is needed.
#
# Requires the 'pharokka' env (phold + foldseek + pharokka) — see steps/00_install_env.sh.
# Usage:  bash example_data/run_smoke_test.sh
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "${HERE}")"
cd "${ROOT}"

echo "[smoke] staging genome_2campy into input/"
mkdir -p input/fasta
cp -f example_data/genome_2campy/input/fasta/*.fasta input/fasta/
cp -f example_data/genome_2campy/input/prophage_list.txt input/ 2>/dev/null || true
rm -f input/.input_type

echo "[smoke] running pipeline (WebAPI search mode — no local DB)"
SEARCH_MODE=webapi bash submit_all.sh || {
  echo "[smoke] submit_all returned non-zero (cluster/SLURM not available?)."
  echo "        On a workstation, run the steps directly; see README quick-start."
}

TABLE="04_output/final_annotations_table.csv"
echo "[smoke] checking ${TABLE}"
test -f "${TABLE}" || { echo "FAIL: ${TABLE} not produced"; exit 1; }

python3 - "${TABLE}" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
need = {"prophage","locus_tag","final_product","final_function",
        "annotation_source","needs_review","fs_top3_kingdoms"}
missing = need - set(rows[0].keys()) if rows else need
assert not missing, f"FAIL: missing columns {missing}"
n = len(rows)
annotated = sum(1 for r in rows if r["final_product"] not in ("", "hypothetical protein"))
print(f"PASS: {n} CDS, {annotated} annotated ({100*annotated//max(n,1)}%), schema OK")
assert n > 0, "FAIL: empty table"
PY
echo "[smoke] done."
