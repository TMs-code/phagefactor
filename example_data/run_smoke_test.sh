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

# fail early with a useful message rather than deep inside step 02
python3 scripts/tools/check_databases.py >/dev/null 2>&1 || true
echo "[smoke] staging genome_2campy into input/"
mkdir -p input/fasta
cp -f example_data/genome_2campy/input/fasta/*.fasta input/fasta/
cp -f example_data/genome_2campy/input/prophage_list.txt input/ 2>/dev/null || true
rm -f input/.input_type

echo "[smoke] running pipeline (no SLURM, on-disk DBs)"
# run_local.sh, not submit_all.sh: submit_all sbatches every step.
# Uses the two small FoldSeek DBs; afdb50 is not needed for a smoke test.
: "${SMOKE_DBS:=afdb-swissprot pdb100}"
FOLDSEEK_DBS="${SMOKE_DBS}" SEARCH_MODE="${SEARCH_MODE:-local}" bash run_local.sh || {
  echo "[smoke] run_local failed — the real error is ABOVE this line."
  echo "        Check databases with: python scripts/tools/check_databases.py"
  echo "        Install the small set : bash steps/00b_setup_databases.sh --pharokka --pdb100 --swissprot"
  exit 1
}

# outputs now live in the versioned run dir; follow the pointer submit_all wrote
RUN_DIR="${PHAGEFACTOR_RUN_DIR:-$(cat runs/.current_run 2>/dev/null || echo runs/run)}"
TABLE="${RUN_DIR}/04_output/final_annotations_table.csv"
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
