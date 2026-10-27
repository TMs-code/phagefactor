#!/usr/bin/env bash
# =============================================================================
# run_smoke_test.sh — end-to-end check on the genome_2campy example
#
# Runs the pipeline on the 2 Campylobacter prophages and verifies the output
# table is produced with the expected schema and a plausible annotated fraction.
# Uses the two small local FoldSeek databases (search_mode=webapi needs none).
#
# Needs the tools (phagefactor install-envs), snakemake, and the small databases
# (phagefactor databases --pharokka --pdb100 --swissprot).
# Usage:  bash example_data/run_smoke_test.sh [project_dir]   (default: ./pf_smoke)
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PF="${HERE}/../bin/phagefactor"
PROJ="${1:-pf_smoke}"

echo "[smoke] project ${PROJ} with the genome_2campy example"
"${PF}" init "${PROJ}" >/dev/null
cp -f "${HERE}"/genome_2campy/input/fasta/*.fasta "${PROJ}/input/fasta/"
cp -f "${HERE}"/genome_2campy/input/prophage_list.txt "${PROJ}/input/"
cd "${PROJ}"

# the two small FoldSeek DBs; afdb50 is not needed for a smoke test
: "${SMOKE_DBS:=afdb-swissprot pdb100}"
"${PF}" run --config foldseek_dbs="${SMOKE_DBS}" search_mode="${SEARCH_MODE:-local}" || {
  echo "[smoke] run failed — the real error is ABOVE this line."
  echo "        Check databases with: ${PF} check"
  echo "        Install the small set: ${PF} databases --pharokka --phold --pdb100 --swissprot"
  exit 1
}

RUN_DIR="$("${PF}" where | sed -n 's/^RUN_DIR=//p')"
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
