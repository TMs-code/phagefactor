#!/usr/bin/env python3
"""
config.py — phageFACTor configuration loader
=========================================================================
Reads config/config.yaml ONCE and exposes the same constant names the analysis
scripts already import (`from config import FOLDSEEK_EVALUE_MAX, ...`), so the
step scripts need ZERO changes when moving between machines.

Path resolution:
  ROOT = ${PHAGEFACTOR_ROOT}  if set, else the repo root (two levels up from here)
  databases.foldseek_db_root = ${FOLDSEEK_DB_ROOT} env wins over the yaml value.

The large semantic sets (UNINFORMATIVE_STRINGS, GENERIC_WORDS,
COMPLEMENTARY_CATEGORY_MAP) stay in this file as code — they are vocabulary, not
configuration, and live next to the logic that uses them.
"""

import os
from pathlib import Path

try:
    import yaml
except ImportError as e:  # pragma: no cover
    raise SystemExit("phageFACTor config needs PyYAML:  pip install pyyaml") from e

# -----------------------------------------------------------------------------
# Load YAML
# -----------------------------------------------------------------------------
_CONF_DIR = Path(__file__).resolve().parent
_YAML = _CONF_DIR / "config.yaml"
with open(_YAML) as _fh:
    _C = yaml.safe_load(_fh)

# -----------------------------------------------------------------------------
# ROOT
# -----------------------------------------------------------------------------
_env_root = os.environ.get("PHAGEFACTOR_ROOT")
PROJECT_ROOT = Path(_env_root).resolve() if _env_root else _CONF_DIR.parent

def _p(rel: str) -> Path:
    return PROJECT_ROOT / rel

# -----------------------------------------------------------------------------
# PROJECT
# -----------------------------------------------------------------------------
HOST_GENUS = _C["project"]["host_genus"]
_NOTE_TAG  = _C["project"]["note_tag"]

# -----------------------------------------------------------------------------
# DIRECTORIES
# -----------------------------------------------------------------------------
_P = _C["paths"]
INPUT_DIR        = _p(_P["input_dir"])
FASTA_DIR        = _p(_P["fasta_dir"])
SPLIT_DIR        = _p(_P["split_dir"])
GBK_DIR          = INPUT_DIR / "gbk"
PHAROKKA_OUT_DIR = _p(_P["pharokka_out"])
PHOLD_OUT_DIR    = _p(_P["phold_out"])
FOLDSEEK_DIR     = _p(_P["foldseek_out"])
COMPARISON_DIR   = _p(_P["comparison_out"])
OUTPUT_DIR       = _p(_P["output_out"])          # 04_output (deliverables)
CURATION_DIR     = _p(_P["curation_out"])        # 04_output/curation
PHYNTENY_DIR     = _p(_P["phynteny_out"])        # 05_phynteny
PHYNTENY_RUN_DIR = _p(_P.get("phynteny_run", "05_phynteny/run"))
GOPHAGE_DIR      = _p(_P["gophage_out"])
SCRIPTS_DIR      = PROJECT_ROOT / "phagefactor"

# Combined phold outputs
PHOLD_COMB_DIR     = PHOLD_OUT_DIR / "combined"
PHOLD_COMBINED_TSV = PHOLD_COMB_DIR / "phold_all.tsv"
PHOLD_3DI_FASTA    = PHOLD_COMB_DIR / "phold_3di.fasta"
PHOLD_AA_FASTA     = PHOLD_COMB_DIR / "phold_aa.fasta"
RAW_GB             = PHOLD_COMB_DIR / "all_prophages_combined.gbk"

# Gene metadata + hypothetical target list
GENE_METADATA_CSV = SPLIT_DIR / "gene_metadata.csv"
PHAROKKA_CDS_TSV_GLOB = str(PHAROKKA_OUT_DIR / "*" / "*_cds_final_merged_output.tsv")
HYPO_TARGETS_DIR = FOLDSEEK_DIR / "targets"
HYPO_GENE_LIST   = HYPO_TARGETS_DIR / "hypothetical_genes.csv"

# Pre-CDS-computed sub-mode
_PC = _C.get("precomputed", {})
SOURCE_GENOME_FASTA    = _p(_PC.get("source_genome_fasta", "input/genome/source_genome.fasta"))
PROPHAGE_WINDOWS_CSV   = _p(_PC.get("prophage_windows_csv", "input/fasta/prophage_windows.csv"))
RICH_GENE_METADATA_CSV = _p(_PC.get("rich_gene_metadata_csv", "input/gene_metadata_rich.csv"))

# -----------------------------------------------------------------------------
# PROPHAGE NAMES  (read from input/prophage_list.txt at runtime if present)
# -----------------------------------------------------------------------------
_plist = INPUT_DIR / "prophage_list.txt"
PROPHAGE_NAMES = (
    [l.strip() for l in _plist.read_text().splitlines() if l.strip()]
    if _plist.exists() else []
)

# -----------------------------------------------------------------------------
# STEP 01 — PHOLD
# -----------------------------------------------------------------------------
PHOLD_THREADS      = 8
PHOLD_TSV_FILENAME = "phold_per_cds_predictions.tsv"

# -----------------------------------------------------------------------------
# STEP 02 — FOLDSEEK
# -----------------------------------------------------------------------------
_F = _C["foldseek"]
FOLDSEEK_3DI_DIR  = FOLDSEEK_DIR / "3di_tokens"
FOLDSEEK_3DI_BEST = FOLDSEEK_3DI_DIR / "best_hit.csv"
FOLDSEEK_3DI_TOP3 = FOLDSEEK_3DI_DIR / "top3.csv"
FOLDSEEK_3DI_ALL  = FOLDSEEK_3DI_DIR / "all_hits.csv"
FOLDSEEK_BEST_HIT = FOLDSEEK_3DI_BEST   # compat aliases for downstream scripts
FOLDSEEK_TOP3     = FOLDSEEK_3DI_TOP3
FOLDSEEK_ALL_HITS = FOLDSEEK_3DI_ALL

FOLDSEEK_CMD     = "foldseek"
FOLDSEEK_THREADS = _F["threads"]
FOLDSEEK_API_URL = _F["api_url"]
FOLDSEEK_EVALUE_MAX     = _F["evalue_max"]
FOLDSEEK_SCORE_OVERRIDE = _F["score_override"]
PROSTT5_MASK_THRESHOLD  = _F["prostt5_mask_threshold"]
FS_CONFIDENT_EVALUE  = _F["confident_evalue"]
FS_GOOD_EVALUE       = _F["good_evalue"]
FS_BORDERLINE_SCORE  = _F["borderline_score"]

# FoldSeek DB roots: env FOLDSEEK_DB_ROOT > yaml > blank
_db_root_str = os.environ.get("FOLDSEEK_DB_ROOT") or _C["databases"].get("foldseek_db_root") or ""
_DB_ROOT = Path(_db_root_str) if _db_root_str else None
FOLDSEEK_LOCAL_DBS = {
    name: (_DB_ROOT / rel if _DB_ROOT else Path(rel))
    for name, rel in _C["databases"]["foldseek_local_dbs"].items()
}
FOLDSEEK_TAXON_FILTER = dict(_C["databases"].get("foldseek_taxon_filter", {}))
PHAROKKA_DB = Path(_C["databases"].get("pharokka_db") or "")

# -----------------------------------------------------------------------------
# STEP 03 — COMPARISON THRESHOLDS
# -----------------------------------------------------------------------------
_CM = _C["compare"]
FUZZY_STRONG_THRESHOLD  = _CM["fuzzy_strong_threshold"]
FUZZY_PARTIAL_THRESHOLD = _CM["fuzzy_partial_threshold"]
PHOLD_TRUSTED_CONF = set(_CM["phold_trusted_conf"])
PHOLD_WEAK_CONF    = set(_CM["phold_weak_conf"])

# -----------------------------------------------------------------------------
# STEP 06 — PHYNTENY
# -----------------------------------------------------------------------------
PHYNTENY_THRESHOLD = _C["phynteny"]["threshold"]

# -----------------------------------------------------------------------------
# STEP 05 — OUTPUT
# -----------------------------------------------------------------------------
FINAL_ANNOTATIONS_TABLE = OUTPUT_DIR / "final_annotations_table.csv"
FINAL_ANNOTATIONS_XLSX  = OUTPUT_DIR / "final_annotations_table.xlsx"
UPDATED_GB              = OUTPUT_DIR / "updated_prophages.gb"
NOTE_TEMPLATE = (f"pipeline={_NOTE_TAG}; source={{source}}; {{evidence}}; "
                 "original_pharokka=hypothetical protein")

# =============================================================================
# SEMANTIC VOCABULARY  (code-level; unchanged from the validated pipeline)
# =============================================================================
UNINFORMATIVE_STRINGS = frozenset({
    "unknown function", "hypothetical protein", "no phold match",
    "na", "nan", "", "none", "n/a",
    "phage protein (unknown function)", "phage protein",
    "uncharacterized protein", "uncharacterized protein (fragment)",
    "predicted protein", "phage protein (fragment)", "defensefinder protein",
    "four helix bundle protein", "four-helix bundle protein",
    "three helix bundle protein", "helix bundle protein", "beta-barrel protein",
    "putative uncharacterized protein",
})

GENERIC_WORDS = frozenset({
    "protein", "domain", "containing", "related", "like", "putative",
    "probable", "hypothetical", "unknown", "function", "phage", "predicted",
    "subunit", "family", "superfamily", "type", "class", "component",
    "associated", "binding", "factor", "homolog", "gene", "conserved",
    "the", "and", "of", "or", "a", "an", "is", "in", "to", "na",
    "cds", "phrog", "phrogs",
})

# NOTE (2026-06-26): consumed two ways. (1) is_complementary() keys on the PHold
# *category* (left) and looks for any keyword (right) in the FoldSeek description.
# (2) 04_curate's _shares_function also treats each value-list as a SYNONYM GROUP
# (symmetric, category-agnostic): two descriptions that each contain a keyword from
# the SAME list are "related" (no review flag). Extend known biological equivalences
# here rather than adding a parallel constant.
COMPLEMENTARY_CATEGORY_MAP = {
    # cell-wall hydrolases (lysis): muramidase/amidase/glycosidase are all peptidoglycan
    # hydrolases -> phosphodiester glycosidase ~ amidase.
    "lysis":                   ["lysin", "endolysin", "holin", "spanin", "amidase",
                                "muramidase", "lysozyme", "glycosidase", "glucosaminidase",
                                "phosphodiester"],
    "tail":                    ["tail", "baseplate", "fiber", "spike", "needle", "tape measure"],
    # protease/peptidase/metallopeptidase as a SYNONYM GROUP (prohead/maturation
    # proteases are head&packaging); used by _shares_function for metallo-protease ~
    # metallopeptidase pairs.
    "head and packaging":      ["portal", "capsid", "terminase", "head", "major capsid",
                                "scaffold", "prohead", "morphogenesis",
                                "protease", "peptidase", "metallopeptidase", "metalloprotease"],
    "integration and excision":["integrase", "excisionase", "recombinase", "transposase",
                                "resolvase", "dna repair", "site-specific"],
    # a transcriptional regulator is a DNA-binding/HTH protein. Fur/ferric-uptake and
    # LexA are CAP-like HTH metal/SOS regulators (Holm 1994 PMID 7708205).
    "transcription regulation":["repressor", "activator", "anti-repressor", "cro",
                                "dna binding", "dna-binding", "helix-turn-helix", "hth",
                                "winged helix", "arc", "regulator", "sigma", "abrb", "lexa",
                                "fur", "ferric", "regulation"],
    # replisome~replication, parA/parB/Soj~chromosome partitioning~ribonuclease, RecT/RecA.
    # Soj is a ParA homolog (ParA-like ~ "Sporulation initiation inhibitor Soj").
    "DNA, RNA and nucleotide metabolism": ["replication", "polymerase", "ligase", "helicase",
                                "primase", "nuclease", "exonuclease", "endonuclease",
                                "recombinase", "recombination", "ribonuclease", "replisome",
                                "reca", "rect", "sbc", "methyltransferase",
                                "deoxyribosyltransferase", "partition", "chromosome partitioning",
                                "parb", "para", "soj", "korb", "atpase", "aaa", "atp-binding",
                                "primosome", "terminase"],
    # Fe-S redox enzymes carried as morons: Fe-S-cluster redox ~ Fe-S-oxidoreductase.
    "moron, auxiliary metabolic gene and host takeover": ["toxin", "antitoxin", "amg",
                                "secreted", "effector", "oxidoreductase", "redox",
                                "fe-s", "fes", "ferredoxin"],
    "connector":               ["connector", "adaptor", "portal", "neck", "head-tail"],
}
