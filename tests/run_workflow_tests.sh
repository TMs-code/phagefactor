#!/usr/bin/env bash
# =============================================================================
# tests/run_workflow_tests.sh — end-to-end tests of the Snakemake workflow
# =============================================================================
# Runs workflow/Snakefile with MOCK tools (tests/mock_bin: pharokka.py, phold,
# foldseek, phynteny_transformer, micromamba) on the shipped example data, so
# it needs no database, no GPU and no SLURM, and takes ~2 minutes.
#
# What it checks
#   * every input mode (genome, protein, pangenome) and search mode (local, none)
#     runs to the final table;
#   * the `phagefactor` launcher: init + run from a project dir, with the code
#     dir read-only (nothing may be written into it);
#   * a run dir without Snakemake metadata (older version) is adopted as up to date;
#   * INCREMENTAL re-runs: each kind of edit re-runs exactly the expected steps
#     (lexicon -> parse+compare+curate; curate threshold -> curate only;
#     evalue -> search onward; one input FASTA -> that prophage onward; one
#     protein -> one phold batch);
#   * a failed job resumes without redoing finished work;
#   * stale protein batches from an earlier, larger split are not merged.
# The mocks prove plumbing, not annotation quality.
#
# Usage:   bash tests/run_workflow_tests.sh [work_dir]
# Needs:   snakemake >= 8.20 on PATH, python3 with pandas, biopython, openpyxl, pyyaml
#          (the driver env envs/pf_snakemake.yaml + pandas/biopython/openpyxl is enough)
# The repo is COPIED into work_dir first; your checkout is never modified.
# =============================================================================
set -uo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE="${1:-$(mktemp -d -t pf_wf_tests.XXXX)}"
# everything goes in a dedicated sub-folder that is wiped first, so a previous
# (possibly interrupted) test run can never leak edited inputs into this one
W="$(mkdir -p "$BASE" && cd "$BASE" && pwd)/pf_wf_tests"
rm -rf "$W"; mkdir -p "$W"
echo "[tests] work dir: $W"

# --- isolated copy of the repo + fake databases --------------------------------
rm -rf "$W/repo"; mkdir -p "$W/repo"
( cd "$SRC" && tar --exclude=.git --exclude=runs --exclude=.snakemake -cf - . ) | ( cd "$W/repo" && tar -xf - )
R="$W/repo"
chmod +x "$R"/tests/mock_bin/*   # exec bits are lost through Windows / core.filemode=false
mkdir -p "$W/db/pharokka_db" "$W/db/foldseek_dbs/afdb_swissprot_db" "$W/db/foldseek_dbs/pdb100_db" "$W/db/foldseek_dbs/afdb50_db" "$W/tmp"
touch "$W/db/foldseek_dbs/afdb_swissprot_db/afdb_swissprot.dbtype" \
      "$W/db/foldseek_dbs/afdb_swissprot_db/afdb_swissprot_taxonomy" \
      "$W/db/foldseek_dbs/pdb100_db/pdb100.dbtype" "$W/db/foldseek_dbs/afdb50_db/afdb50.dbtype"

# never inherit a real project's settings from the caller's shell
for v in PHAGEFACTOR_INPUT PHAGEFACTOR_RUN_DIR PHAGEFACTOR_PROJECT PHAGEFACTOR_CONFIG \
         PHAGEFACTOR_HOME PHOLD_DB INPUT_MODE SEARCH_MODE FOLDSEEK_DBS MAMBA_ENV \
         FOLDSEEK_DB_ROOT PHAROKKA_DB PROTEIN_BATCH_SIZE PROTEIN_FASTA_DIR PHOLD_MODE \
         MAMBA_EXE MAMBA_ROOT_PREFIX PHYNTENY_ENV PHYNTENY_CACHE SLURM_EMAIL RUN_PHYNTENY; do
    unset "$v"
done
export PATH="$R/tests/mock_bin:$PATH"
# phynteny is ON by default; sections [1]-[6] test the rest without it
export RUN_PHYNTENY=0
export PHYNTENY_CACHE="$W/phycache"     # shared ESM2 cache (mock), never the run dir
export PHAGEFACTOR_DB_ROOT="$W/db" SCRATCH_TMPDIR="$W/tmp" PHAGEFACTOR_HOST_GENUS=""
export XDG_CONFIG_HOME="$W/xdg"         # never read the caller's ~/.config/phagefactor
cd "$R"

N_PASS=0; N_FAIL=0
pass() { echo "  PASS  $*"; N_PASS=$((N_PASS + 1)); }
fail() { echo "  FAIL  $*"; N_FAIL=$((N_FAIL + 1)); }

smk() {  # smk <input_dir> <run_dir> [snakemake args...]
    local in="$1" run="$2"; shift 2
    snakemake --profile profiles/local --cores 4 --config input_dir="$in" run_dir="$run" "$@" \
        > "$W/last.log" 2>&1
}
plan() { # plan <input_dir> <run_dir> [args] -> "rule=n ..." of a dry run, or "NOTHING"
    local in="$1" run="$2"; shift 2
    snakemake --profile profiles/local -n --config input_dir="$in" run_dir="$run" "$@" 2>&1 | python3 -c '
import sys, re
txt = sys.stdin.read()
if "Nothing to be done" in txt:
    print("NOTHING"); sys.exit()
m = re.search(r"Job stats:\n.*?\n-+\s+-+\n(.*?)\ntotal", txt, re.S)
if not m:
    print("ERROR"); sys.exit()
jobs = dict(l.split() for l in m.group(1).strip().splitlines())
print(" ".join(f"{k}={v}" for k, v in sorted(jobs.items()) if k != "all"))'
}
expect_plan() { # expect_plan <label> <expected> <input> <run> [args]
    local label="$1" want="$2"; shift 2
    local got; got="$(plan "$@")"
    [[ "$got" == "$want" ]] && pass "$label: $got" || fail "$label: expected [$want] got [$got]"
}
same_tree() { # same_tree <label> <dirA> <dirB> : all data files byte-identical
    local label="$1" a="$2" b="$3" bad=0 f
    while IFS= read -r f; do
        if [[ ! -f "$b/$f" ]]; then echo "        only in A: $f"; bad=1
        elif ! cmp -s "$a/$f" "$b/$f"; then echo "        differs: $f"; bad=1; fi
    done < <(cd "$a" && find . -type f \( -name '*.csv' -o -name '*.tsv' -o -name '*.fasta' \
             -o -name '*.faa' -o -name '*.gb' -o -name '*.gbk' -o -name '*.m8' -o -name '*.txt' \) \
             -not -path './logs/*' -not -path './.snakemake/*' | sort)
    [[ $bad == 0 ]] && pass "$label" || fail "$label"
}
table_rows() { python3 -c "import pandas as pd,sys; print(len(pd.read_csv(sys.argv[1])))" "$1" 2>/dev/null || echo 0; }

cp -r example_data/genome_2campy/input  "$W/in_genome"
cp -r example_data/protein_2citro/input  "$W/in_protein"

# ------------------------------------------------------------------------------
echo "[1] genome mode"
smk "$W/in_genome" "$W/g_smk" && pass "snakemake genome run" || { fail "snakemake genome run"; tail -30 "$W/last.log"; }
n=$(table_rows "$W/g_smk/04_output/final_annotations_table.csv"); [[ $n -gt 0 ]] && pass "final table: $n rows" || fail "final table empty"
expect_plan "idempotent re-run" "NOTHING" "$W/in_genome" "$W/g_smk"
cp -a "$W/g_smk" "$W/g_nometa"; rm -rf "$W/g_nometa/.snakemake/metadata"
expect_plan "run dir without Snakemake metadata (older version) adopted" "NOTHING" "$W/in_genome" "$W/g_nometa"

echo "[2] incremental re-runs (genome)"
echo "# test edit" >> scripts/lib/lexicon.py
expect_plan "lexicon edited" "compare=1 curate=1 foldseek_parse=1 report=1" "$W/in_genome" "$W/g_smk"
smk "$W/in_genome" "$W/g_smk"
sed -i 's/subdb_close_min_fident: 0.4$/subdb_close_min_fident: 0.45/' config/config.yaml
expect_plan "curate threshold" "curate=1 report=1" "$W/in_genome" "$W/g_smk"
smk "$W/in_genome" "$W/g_smk"
sed -i 's/  fuzzy_strong_threshold: 0.35/  fuzzy_strong_threshold: 0.36/' config/config.yaml
expect_plan "compare threshold" "compare=1 curate=1 report=1" "$W/in_genome" "$W/g_smk"
smk "$W/in_genome" "$W/g_smk"
sed -i 's/low_info_malus: 0.60/low_info_malus: 0.55/' config/config.yaml
expect_plan "foldseek scoring setting" "compare=1 curate=1 foldseek_parse=1 report=1" "$W/in_genome" "$W/g_smk"
smk "$W/in_genome" "$W/g_smk"
sed -i 's/  evalue_max: 0.1 /  evalue_max: 0.2 /' config/config.yaml
expect_plan "search e-value" "compare=1 curate=1 foldseek_parse=1 foldseek_search=1 report=1" "$W/in_genome" "$W/g_smk"
smk "$W/in_genome" "$W/g_smk"
printf 'ACGT\n' >> "$W/in_genome/fasta/CJIE2.fasta"
expect_plan "one genome FASTA edited" \
  "compare=1 curate=1 foldseek_parse=1 foldseek_search=1 foldseek_targets=1 pharokka=1 phold=1 phold_merge=1 report=1" \
  "$W/in_genome" "$W/g_smk"
smk "$W/in_genome" "$W/g_smk"
expect_plan "settled again" "NOTHING" "$W/in_genome" "$W/g_smk"

echo "[3] failure + resume"
cp -r example_data/genome_2campy/input "$W/in_fail"
MOCK_FAIL=CJIE2 smk "$W/in_fail" "$W/g_fail" --keep-going && fail "forced failure did not fail" || pass "forced phold failure stops the run"
expect_plan "resume redoes only the failed prophage onward" \
  "compare=1 curate=1 foldseek_parse=1 foldseek_search=1 foldseek_targets=1 phold=1 phold_merge=1 report=1" \
  "$W/in_fail" "$W/g_fail"
smk "$W/in_fail" "$W/g_fail" && pass "resume completes" || fail "resume"

echo "[4] protein mode (precomputed sub-mode: rich metadata + host genome)"
smk "$W/in_protein" "$W/p_smk" && pass "snakemake protein run" || { fail "snakemake protein run"; tail -30 "$W/last.log"; }
[[ -f "$W/p_smk/04_output/updated_prophages.gb" ]] && pass "GenBank built" || fail "no GenBank"
python3 - "$W/in_protein/fasta/phiSM_proteins.faa" <<'PY'
import sys; p = sys.argv[1]; L = open(p).read().splitlines()
L[-1] = L[-1][:-1] + ("W" if L[-1][-1] != "W" else "A"); open(p, "w").write("\n".join(L) + "\n")
PY
smk "$W/in_protein" "$W/p_smk"
nb=$(grep -cE "^(local)?rule phold_proteins:" "$W/last.log" || true)
[[ "$nb" == "1" ]] && pass "one protein edited -> 1 phold batch re-run" || fail "one protein edited -> $nb phold batches re-run (expected 1)"
mkdir -p "$W/aside"; mv "$W/in_protein/fasta/phiSM_proteins.faa" "$W/aside/"
smk "$W/in_protein" "$W/p_smk" --config protein_batch_size=40 input_dir="$W/in_protein" run_dir="$W/p_smk"
n=$(table_rows "$W/p_smk/04_output/final_annotations_table.csv")
[[ "$n" == "105" ]] && pass "input shrunk: stale batches not merged (105 rows)" || fail "input shrunk: $n rows (expected 105)"
mv "$W/aside/phiSM_proteins.faa" "$W/in_protein/fasta/"
smk "$W/in_protein" "$W/p_smk" >/dev/null
G="$W/p_smk/01_phold/combined/all_prophages_combined.gbk"
echo "STALE-MARKER" >> "$G"                                   # a GenBank left by an earlier run
echo "" >> "$W/in_protein/fasta/prophage_windows.csv"           # window file edited
smk "$W/in_protein" "$W/p_smk"
grep -qE "^(local)?rule curate:" "$W/last.log" && ! grep -q "STALE-MARKER" "$G" \
    && pass "protein mode: window file edited -> curate re-run, GenBank rebuilt (not reused)" \
    || fail "protein mode: stale GenBank reused"

echo "[5] pangenome mode"
mkdir -p "$W/in_pan/fasta"
python3 - "$W/g_smk/04_output/updated_prophages.gb" "$W/in_pan/fasta/pan_reps.fa" <<'PY'
import sys
from Bio import SeqIO
n = 0
with open(sys.argv[2], "w") as out:
    for r in SeqIO.parse(sys.argv[1], "genbank"):
        for f in r.features:
            if f.type == "CDS":
                n += 1
                out.write(f">{f.qualifiers['locus_tag'][0]}\n{f.extract(r.seq)}\n")
PY
smk "$W/in_pan" "$W/pan_smk" && pass "pangenome run (auto-detected)" || { fail "pangenome run"; tail -30 "$W/last.log"; }
[[ -s "$W/pan_smk/04_output/pangenome_reference.csv" ]] && pass "pangenome_reference.csv" || fail "no pangenome_reference.csv"
[[ -f "$W/in_pan/fasta/pan_reps.fa" ]] && pass "nucleotide input left in place" || fail "input FASTA was moved"

echo "[6] search modes"
smk "$W/in_genome" "$W/none_smk" --config search_mode=none input_dir="$W/in_genome" run_dir="$W/none_smk" \
    && pass "search_mode=none run" || fail "search_mode=none run"
got=$(plan "$W/in_genome" "$W/web_smk" search_mode=webapi)
[[ "$got" == *"foldseek_webapi=1"* && "$got" != *"foldseek_search"* ]] && pass "search_mode=webapi DAG: $got" || fail "webapi DAG: $got"

echo "[7] phynteny"
got=$(RUN_PHYNTENY= plan "$W/in_genome" "$W/g_smk")
[[ "$got" == "phynteny=1 report=1" ]] && pass "phynteny ON by default (only it + the report left to run): $got" || fail "phynteny default: $got"
smk "$W/in_genome" "$W/g_smk" --config phynteny=true input_dir="$W/in_genome" run_dir="$W/g_smk" \
    && [[ -s "$W/g_smk/05_phynteny/final_annotations_integrated.csv" ]] \
    && grep -q "PHYNTENY RAN" "$W/last.log" && grep -q "PHYNTENY: ON" "$W/last.log" \
    && pass "phynteny + integration, flagged at start and end" || fail "phynteny"
[[ ! -e "$W/g_smk/cache" && ! -e "$W/g_smk/skipped_genomes.txt" \
   && -f "$W/g_smk/05_phynteny/run/skipped_genomes.txt" && -L "$W/g_smk/05_phynteny/run/cache" \
   && -f "$W/phycache/models--facebook--esm2_t33_650M_UR50D/weights.bin" ]] \
    && pass "phynteny: ESM2 cache shared, nothing written to the run root" \
    || fail "phynteny: cache/skipped_genomes.txt location ($(ls "$W/g_smk" | tr '\n' ' '))"
R_HTML="$W/g_smk/04_output/report.html"
python3 - "$R_HTML" "$W/g_smk/04_output/final_annotations_table.csv" <<'PY' && pass "report.html: all genes embedded, phynteny + GenBank coordinates included" || fail "report.html content"
import json, re, sys, pandas as pd
html = open(sys.argv[1]).read()
data = json.loads(re.search(r"const PF = (\{.*?\});\n", html, re.S).group(1))
n = len(pd.read_csv(sys.argv[2]))
assert data["meta"]["n"] == n == len(data["genes"]), (data["meta"]["n"], n)
assert data["meta"]["has_phynteny"] and data["meta"]["has_coords"], data["meta"]
assert "</script" not in json.dumps(data), "unescaped </script>"
PY
echo "<!-- edit -->" >> scripts/lib/report_template.html
expect_plan "report template edited -> report only" "report=1" "$W/in_genome" "$W/g_smk" phynteny=true
echo "# synteny" >> scripts/lib/synteny_hint.py
expect_plan "lib/synteny_hint.py edited -> phynteny only" "phynteny=1 report=1" "$W/in_genome" "$W/g_smk" phynteny=true
echo "# rules" >> scripts/lib/rules.py
expect_plan "lib/rules.py edited -> FoldSeek post-processing onward + phynteny" "compare=1 curate=1 foldseek_parse=1 phynteny=1 report=1" "$W/in_genome" "$W/g_smk" phynteny=true
smk "$W/in_protein" "$W/p_smk" --config phynteny=true input_dir="$W/in_protein" run_dir="$W/p_smk" \
    ; grep -qE "PHYNTENY (RAN|WAS SKIPPED)" "$W/last.log" \
    && pass "protein mode phynteny flagged: $(grep -oE 'PHYNTENY (RAN|WAS SKIPPED)' "$W/last.log" | head -1)" || fail "protein phynteny flag"

echo "[8] software environments"
# tools on PATH (conda install): run in the current env, no micromamba needed
mkdir -p "$W/nomamba_bin" "$W/mroot/bin"
for f in "$R"/tests/mock_bin/*; do [[ "$(basename "$f")" == micromamba ]] || ln -sf "$f" "$W/nomamba_bin/"; done
cp "$R/tests/mock_bin/micromamba" "$W/mroot/bin/micromamba"
NOMAMBA_PATH="$W/nomamba_bin:$(echo "$PATH" | tr ':' '\n' | grep -v "tests/mock_bin" | paste -sd:)"
PATH="$NOMAMBA_PATH" smk "$W/in_genome" "$W/env_auto" --config search_mode=none input_dir="$W/in_genome" run_dir="$W/env_auto" \
    && grep -q "env: (current)" "$W/last.log" && [[ -s "$W/env_auto/04_output/final_annotations_table.csv" ]] \
    && pass "tools on PATH: runs in the current env, without micromamba" || fail "tools_env auto"
# a named env: micromamba found as a binary ($MAMBA_EXE, $MAMBA_ROOT_PREFIX), never a shell function
PATH="$NOMAMBA_PATH" MAMBA_EXE="$W/mroot/bin/micromamba" \
    smk "$W/in_genome" "$W/mamba_exe" --config search_mode=none tools_env=phagefactor input_dir="$W/in_genome" run_dir="$W/mamba_exe" \
    && grep -q "env: phagefactor" "$W/last.log" && [[ -s "$W/mamba_exe/04_output/final_annotations_table.csv" ]] \
    && pass "tools_env=phagefactor via \$MAMBA_EXE (full run)" || fail "MAMBA_EXE lookup"
got=$(PATH="$NOMAMBA_PATH" MAMBA_ROOT_PREFIX="$W/mroot" plan "$W/in_genome" "$W/mamba_root" tools_env=phagefactor)
[[ "$got" == *"pharokka="* ]] && pass "found via \$MAMBA_ROOT_PREFIX/bin" || fail "MAMBA_ROOT_PREFIX lookup: $got"
NOTOOLS_PATH="$(echo "$PATH" | tr ':' '\n' | grep -v "tests/mock_bin" | paste -sd:)"
PATH="$NOTOOLS_PATH" MAMBA_ROOT_PREFIX="$W/nowhere" CONDA_EXE="" \
    snakemake --profile profiles/local -n --config input_dir="$W/in_genome" run_dir="$W/mamba_none" > "$W/last.log" 2>&1
grep -q "no micromamba or conda binary was found" "$W/last.log" && pass "clear error when the tools and micromamba are both missing" \
    || fail "missing micromamba: $(tail -3 "$W/last.log")"

echo "[9] phagefactor launcher, read-only code dir"
RO="$W/ro_code"; mkdir -p "$RO"; ( cd "$R" && tar --exclude=.snakemake -cf - . ) | ( cd "$RO" && tar -xf - )
find "$RO" -name __pycache__ -prune -exec rm -rf {} +
chmod -R a-w "$RO"
# __pycache__ excluded: Python writes it only where it can (root ignores a-w)
snap() { find "$RO" -name __pycache__ -prune -o -print | sort | md5sum; }
before=$(snap)
P="$W/proj"
"$RO/bin/phagefactor" init "$P" > "$W/last.log" 2>&1 && [[ -f "$P/config/config.yaml" && -f "$P/config/workflow.yaml" ]] \
    && pass "phagefactor init" || fail "phagefactor init: $(tail -3 "$W/last.log")"
cp -r example_data/genome_2campy/input/. "$P/input/"
( cd "$P" && "$RO/bin/phagefactor" run --cores 4 --config search_mode=none > "$W/last.log" 2>&1 ) \
    && [[ -s "$P/runs/run/04_output/final_annotations_table.csv" ]] \
    && pass "phagefactor run in a project dir -> <project>/runs/run" || { fail "phagefactor run"; tail -20 "$W/last.log"; }
[[ "$(snap)" == "$before" ]] && pass "nothing written into the code dir" || fail "code dir modified"
sed -i 's|^  run_dir: .*|  run_dir:        runs/second|' "$P/config/config.yaml"
[[ "$(cd "$P" && "$RO/bin/phagefactor" where | sed -n 's/^RUN_DIR=//p')" == "$P/runs/second" ]] \
    && pass "project config.yaml is the one read (run_dir)" || fail "project config not used"
( cd "$W" && "$RO/bin/phagefactor" run -n > "$W/last.log" 2>&1 ) && fail "run outside a project dir accepted" \
    || { grep -q "no config/" "$W/last.log" && pass "clear error outside a project dir" || fail "outside project: $(tail -2 "$W/last.log")"; }
"$RO/bin/phagefactor" --version | grep -q "phageFACTor $(cat "$R/VERSION")" && pass "--version" || fail "--version"
mkdir -p "$W/shared_db"
( cd "$P" && PHAGEFACTOR_DB_ROOT= "$RO/bin/phagefactor" databases "$W/shared_db" --bogus > /dev/null 2>&1 ); \
grep -q "db_root: $W/shared_db" "$XDG_CONFIG_HOME/phagefactor/settings.yaml" \
    && [[ "$(cd "$P" && PHAGEFACTOR_DB_ROOT= "$RO/bin/phagefactor" where | sed -n 's/^DB_ROOT=//p')" == "$W/shared_db" ]] \
    && pass "phagefactor databases <dir>: location recorded and used by every project" || fail "shared database location"
# per-user machine settings: below the project, above the code defaults
mkdir -p "$XDG_CONFIG_HOME/phagefactor"
cat >> "$XDG_CONFIG_HOME/phagefactor/settings.yaml" <<'YAML'
slurm_qos: userqos
resources:
  foldseek_search: {constraint: avx512}
YAML
( cd "$P" && "$RO/bin/phagefactor" run -n -p --forceall > "$W/last.log" 2>&1 )
grep -q "constraint=avx512" "$W/last.log" && grep -q "mem_mb=262144" "$W/last.log" && grep -q "qos=userqos" "$W/last.log" \
    && pass "user settings: resource added to one step, defaults kept, slurm_qos applied" || fail "user settings layering"
echo "slurm_qos: projectqos" >> "$P/config/workflow.yaml"
( cd "$P" && "$RO/bin/phagefactor" run -n -p --forceall > "$W/last.log" 2>&1 )
grep -q "qos=projectqos" "$W/last.log" && ! grep -q "qos=userqos" "$W/last.log" \
    && pass "project workflow.yaml overrides the user settings" || fail "project over user settings"
grep -q "^db_root: $W/shared_db" "$XDG_CONFIG_HOME/phagefactor/settings.yaml" \
    && pass "user settings file keeps db_root next to the workflow settings" || fail "settings.yaml db_root"
chmod -R u+w "$RO"

# ------------------------------------------------------------------------------
echo "[10] hit taxonomy (lineage > taxid > organism name)"
grep -rqs "Hit domain source: lineage=[1-9]" "$W/g_smk/logs" \
    && pass "local search: hit domain read from the FoldSeek lineage" || fail "no lineage-based hit domain in the 02 log"
mkdir -p "$W/taxdump"
printf '1\t|\t1\t|\n2\t|\t1\t|\n2759\t|\t1\t|\n10239\t|\t1\t|\n4751\t|\t2759\t|\n5000\t|\t4751\t|\n561\t|\t2\t|\n' > "$W/taxdump/nodes.dmp"
printf '9999\t|\t5000\t|\n' > "$W/taxdump/merged.dmp"
PHAGEFACTOR_TAXDUMP="$W/taxdump" python3 - "$R/scripts/lib" <<'PY' && pass "taxonomy: lineage prefixes, taxdump walk + merged ids, name fallback" || fail "taxonomy module"
import sys; sys.path.insert(0, sys.argv[1])
from taxonomy import hit_kingdom, kingdom_from_lineage
assert kingdom_from_lineage("-_cellular organisms;d_Eukaryota;k_Fungi") == "Eukaryote"
assert kingdom_from_lineage("-_cellular organisms;-_Eukaryota;k_Metazoa") == "Eukaryote"
assert kingdom_from_lineage("-_Viruses;c_Caudoviricetes") == "Virus"
assert hit_kingdom("5000", "Kipferlia bialata", "") == ("Eukaryote", "taxid")
assert hit_kingdom("9999", "", "") == ("Eukaryote", "taxid")             # merged id
assert hit_kingdom("561", "Homo sapiens", "d_Bacteria") == ("Bacteria", "lineage")
assert hit_kingdom("256318", "metagenome", "")[0] == "Unknown"          # no domain
assert hit_kingdom("", "Escherichia coli", "") == ("Bacteria", "name")
PY

echo
echo "[tests] ${N_PASS} passed, ${N_FAIL} failed   (work dir: $W)"
[[ $N_FAIL == 0 ]]
