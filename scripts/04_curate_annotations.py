#!/usr/bin/env python3
"""
04_curate_annotations.py
=========================
Automated curation and flagging of the PHold vs FoldSeek comparison.

Reads:
  03_comparison/comparison_per_gene.csv

Decision tree per gene (for the 262 hypothetical proteins):
  ---------------------------------------------------------------------
  Case                                     Auto-action       Flag?
  ---------------------------------------------------------------------
  strong/partial match, both inf           Auto-merge        No
  complementary (cat matches FS desc)      Auto-merge        No
  phold_only (PHold inf, FS uninf)         Use PHold         No
  foldseek_only (FS inf, PHold uninf)      Use FoldSeek      No
  different (both inf, diverge)            Flag              Yes
  both_uninformative                       "hypothetical"    No
  ---------------------------------------------------------------------

For flagged genes, an AI_suggestion column is pre-filled with a proposed
merged annotation, and an AI_explanation column explains the reasoning.
You review these in Excel, adjust if needed, then fill final_annotation.

Writes:
  04_curation/auto_curated.csv    -> confidently decided genes (no review needed)
  04_curation/needs_review.csv    -> flagged genes for manual review
      Columns include: ...(comparison cols)..., AI_suggestion, AI_explanation,
                       final_annotation  [<- fill this in Excel]

After manual review:
  -> Save needs_review.csv with final_annotation column filled
  -> Copy it back as 04_curation/needs_review.csv
  -> Run scripts/05_build_output.py

Usage:
  cd phagefactor/
  python scripts/04_curate_annotations.py
"""

import sys
import re
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).parent
_PROJECT_DIR = _SCRIPTS_DIR.parent
sys.path.insert(0, str(_PROJECT_DIR))

from config import (
    COMPARISON_DIR, CURATION_DIR,
    FUZZY_STRONG_THRESHOLD, FUZZY_PARTIAL_THRESHOLD,
    PHOLD_TRUSTED_CONF, PHOLD_WEAK_CONF,
    COMPLEMENTARY_CATEGORY_MAP,
    HOST_GENUS,
    GENERIC_WORDS,
)
from utils import log, section, clean_str, is_informative, tokenize, fuzzy_score

try:
    import pandas as pd
except ImportError:
    print("pandas required.")
    sys.exit(1)


# -----------------------------------------------------------------------------
# COMPLEMENTARY DETECTION
# -----------------------------------------------------------------------------

def is_complementary(phold_cat: str, fs_desc: str) -> bool:
    """
    Return True if PHold's functional category and the FoldSeek description
    are complementary (same biology, different vocabulary).

    E.g. PHold category='lysis', FS desc='endolysin' -> complementary
    """
    if not phold_cat or not fs_desc:
        return False
    cat_lower = phold_cat.lower()
    desc_lower = fs_desc.lower()
    for category, keywords in COMPLEMENTARY_CATEGORY_MAP.items():
        if category.lower() in cat_lower:
            for kw in keywords:
                if kw.lower() in desc_lower:
                    return True
    return False


# -----------------------------------------------------------------------------
# DEFENSEFINDER ANNOTATION HELPER
# -----------------------------------------------------------------------------

def _parse_defensefinder_name(tophit: str) -> str:
    """
    Convert a DefenseFinder gene ID to a readable system name.
    E.g. 'Gao_Mza_753'  -> 'Gao Mza defense system protein'
         'Avs_1_001'    -> 'Avs1 defense system protein'
    Returns empty string if the pattern does not match, the tophit is a raw
    number, or the name contains phage-gene-number segments (gp##, orf##) that
    indicate a phage protein rather than a defense-system entry.
    """
    if not tophit or not isinstance(tophit, str):
        return ""
    t = tophit.strip()
    # Reject sentinel / missing-value strings (e.g. pandas NaN read as "NA")
    if t.lower() in ("nan", "na", "none", "n/a", ""):
        return ""
    # Skip raw numeric scores like "55.96875"
    if re.match(r'^[\d.]+$', t):
        return ""
    # Reject bare phage-gene patterns (gp29, orf12 alone) but allow genuine
    # DefenseFinder multi-component system IDs like MMB_gp29_gp30_1094.
    # Only reject when the ENTIRE name starts with a phage-gene prefix.
    if re.match(r'^(?:gp|orf|phi)\d+', t, re.IGNORECASE):
        return ""
    # Strip trailing numeric serial (e.g. _1094, _753, _001)
    base = re.sub(r'_\d+$', '', t)
    if not base:
        return ""
    # Split on underscores: all-lowercase/digit tokens with gp/orf prefix = component IDs
    parts = base.split('_')
    # Separate the system acronym(s) from the component protein IDs
    system_parts, comp_parts = [], []
    for p in parts:
        if re.match(r'^(?:gp|orf)\d+', p, re.IGNORECASE):
            comp_parts.append(p)
        else:
            system_parts.append(p)
    if not system_parts:
        return ""
    system_name = ' '.join(system_parts)
    if comp_parts:
        return f"{system_name} ({'/'.join(comp_parts)}) defense system protein"
    else:
        return f"{system_name} defense system protein"


# -----------------------------------------------------------------------------
# POST-SELECTION ANNOTATION HELPERS
# -----------------------------------------------------------------------------

_DEFENSE_DESC_RE = re.compile(
    r'\b(anti)?toxin\b|toxin-antitoxin|abortive\s+infection|anti-?crispr|anti-?phage|'
    r'defen[cs]e|restriction|\bcbass\b|\bparis\b|\bRM\b\s*system|immunity', re.I)


def _desc_is_defense(desc: str) -> bool:
    """True if the description itself names a defense/TA function. Used to decide
    whether a defense flag (which 03_compare may set from ANY top-3 hit) should
    force review: only when the CHOSEN phold/FS call is itself defense -- an
    incidental TA entry in a non-chosen top-3 slot must not flag an otherwise
    agreeing gene (e.g. CEDLNMDB_00945, a transcriptional regulator)."""
    return bool(desc) and isinstance(desc, str) and bool(_DEFENSE_DESC_RE.search(desc))


def _strip_fragment(desc: str) -> str:
    """
    Remove trailing '(Fragment)' from an annotation string.
    The fragment suffix refers to the PDB structure being incomplete, not the
    gene -- it is irrelevant for functional annotation.
    """
    if not desc or not isinstance(desc, str):
        return desc
    return re.sub(r'\s*\(fragment\)\s*$', '', desc.strip(), flags=re.IGNORECASE).strip()


# Rules: (pattern_on_desc, pattern_to_find_in_top3, replacement_template)
# replacement_template: None -> use the matched top3 entry verbatim;
#                       str  -> format string with {match} = matched top3 entry
_UPGRADE_RULES = [
    # Prefer specific phage repressor over generic "Transcriptional regulator"
    (re.compile(r'^transcriptional regulator$', re.I),
     re.compile(r'bacteriophage\s+[Cc][Ii]\s+repressor|ci\s+repressor\s+protein', re.I),
     None),
    # Prefer GTPase Era (a bacterial ribosome assembly factor) over structural
    # "Four helix bundle protein" -- EST3 telomere replication protein is a
    # false structural homology from yeast; GTPase Era is the correct bacterial hit.
    (re.compile(r'^four helix bundle protein$', re.I),
     re.compile(r'\bGTPase Era\b', re.I),
     None),
    # Prefer Endonuclease I over the broader "Deoxyribonuclease"
    (re.compile(r'^deoxyribonuclease$', re.I),
     re.compile(r'\bendonuclease\s+I\b', re.I),
     None),
    # Prefer Phosphocholine transferase AnkX over generic ankyrin repeat
    (re.compile(r'^ankyrin repeat', re.I),
     re.compile(r'phosphocholine transferase|AnkX', re.I),
     None),
    # Biotin operon repressor + phage associated context -> compound name
    (re.compile(r'^winged helix.turn.helix domain', re.I),
     re.compile(r'biotin operon repressor', re.I),
     "{match} (phage associated)"),
    # Prefer antitoxin RelB over CopG when both are in top3;
    # otherwise expand CopG to full descriptive name
    (re.compile(r'^CopG family', re.I),
     re.compile(r'antitoxin RelB|RelB\b', re.I),
     "CopG family transcriptional regulator (putative antitoxin, RHH family)"),
    # Prefer HicB antitoxin qualifier for anti-repressor
    (re.compile(r'^anti.repressor$', re.I),
     re.compile(r'HicB family antitoxin|hicB', re.I),
     "anti-repressor (putative HicB antitoxin)"),
    # Resolvase beats generic DNA binding protein
    (re.compile(r'^DNA binding protein$', re.I),
     re.compile(r'resolvase', re.I),
     None),
    # Com family DNA-binding -> conservative generic name (top-3 don't agree on specifics)
    (re.compile(r'^Com family DNA-binding', re.I),
     None,   # unconditional
     "transcription regulator"),
    # Strip "domain-containing" from tail fiber descriptions
    (re.compile(r'^tail fiber domain.containing protein$', re.I),
     None,   # no top3 lookup needed -- unconditional rename
     "Tail fiber protein"),
    # YopX: strip verbose "domain-containing protein" suffix
    (re.compile(r'^YopX protein domain.containing protein', re.I),
     None,
     "YopX protein"),
    # "minor tail protein / TIGR04255" -> "minor tail protein"
    (re.compile(r'^minor tail protein\s*/\s*TIGR', re.I),
     None,
     "minor tail protein"),
    # "DNA helicase / AAA family ATPase" -> "DNA helicase" (helicases are ATPases)
    (re.compile(r'^DNA helicase\s*/\s*AAA family ATPase', re.I),
     None,
     "DNA helicase"),
    # gp16 (Mu-like) -> shorter display name with GemA note (user-curated form)
    (re.compile(r'Mu-like prophage protein gp16.*GemA DNA gyrase', re.I),
     None,
     "Mu-like gp16 (putative GemA DNA gyrase)"),
    # gp16 (Mu-like) without GemA yet -> add GemA note first, then shorten above
    (re.compile(r'^Mu-like prophage protein gp16$', re.I),
     None,
     "Mu-like gp16 (putative GemA DNA gyrase)"),
    # Bacteriophage CI repressor: fix lowercase "ci" -> "CI" (standard nomenclature)
    (re.compile(r'(?i)bacteriophage\s+ci\s+repressor', re.I),
     None,
     "Bacteriophage CI repressor protein"),
    # HTH cro/C1-type domain-containing protein -> concise form with domain in parens
    (re.compile(r'^HTH cro/C1-type domain.containing protein$', re.I),
     None,
     "Transcription regulator (HTH cro/C1-type domain)"),
    # Transcriptional repressor NrdR-like N-terminal domain-containing protein -> NrdR
    (re.compile(r'^Transcriptional repressor NrdR.*$', re.I),
     None,
     "Transcriptional repressor NrdR"),
    # Resolvase HTH domain-containing protein -> concise form
    (re.compile(r'^Resolvase HTH domain.containing protein$', re.I),
     None,
     "Resolvase (HTH domain)"),
]


def _apply_final_upgrades(desc: str, top3_str: str) -> str:
    """
    Apply post-selection annotation upgrades:
      1. Strip (Fragment) suffix.
      2. Apply _UPGRADE_RULES: prefer more informative top-3 description when
         the primary description matches a less-specific pattern.
    """
    if not desc or not isinstance(desc, str):
        return desc

    # Strip (Fragment) first
    desc = _strip_fragment(desc)

    # Parse top3 entries
    top3 = [t.strip() for t in top3_str.split("|") if t.strip() and t.strip() != "NA"] \
           if top3_str and top3_str not in ("NA", "") else []

    for desc_pat, top3_pat, template in _UPGRADE_RULES:
        if not desc_pat.search(desc):
            continue
        if top3_pat is None:
            # Unconditional rename (no top3 lookup)
            if template:
                desc = template
            break
        # Search top3 for a matching entry
        for t3 in top3:
            if top3_pat.search(t3):
                matched = _strip_fragment(t3)
                if template:
                    desc = template.format(match=matched)
                else:
                    desc = matched
                break
    return desc


# Map description text -> function category for proteins that have no PHold call
_FUNC_CAT_RULES = [
    (re.compile(r'transcription.{0,4}regulat|ci.repressor|sigma.factor|copG|NrdR|'
                r'anti.repressor|HTH.cro|bacteriophage.ci|Com.family.DNA.binding|'
                r'repressor.protein|winged.helix', re.I),
     "transcription regulation"),
    (re.compile(r'DNA.helicase|RNA.polymerase|DNA.gyrase|topoisomerase|'
                r'GemA.DNA.gyrase|gp16.*gyrase|replication.initiation|'
                r'DNA.transposition|primase|DNA.repair', re.I),
     "DNA, RNA and nucleotide metabolism"),
    (re.compile(r'tail.fiber|tail.protein|tail.assembly|baseplate|'
                r'head.decoration|major.tail|minor.tail|tape.measure|'
                r'Phage.tail', re.I),
     "tail"),
    (re.compile(r'portal.protein|head.morphogenesis|capsid|head.and.packaging', re.I),
     "head and packaging"),
    (re.compile(r'head.tail.connector|head.to.tail.connector|connector.protein|'
                r'gp6.like.head.tail|neck.protein', re.I),
     "connector"),
    (re.compile(r'RloG|defense.system|abortive.infection|anti.crispr|anti.phage', re.I),
     "defense"),
    (re.compile(r'integrase|excisionase|transposase|resolvase|invertase', re.I),
     "integration and excision"),
    (re.compile(r'beta.lactamase|antibiotic.resistance', re.I),
     "moron, auxiliary metabolic gene and host takeover"),
    (re.compile(r'deoxyribonuclease|endonuclease|restriction.endonuclease|'
                r'HNH.endonuclease', re.I),
     "DNA, RNA and nucleotide metabolism"),
    # 2026-06: host-derived metabolic enzymes, ribosomal proteins and metal-uptake
    # regulators carried by the prophage are auxiliary metabolic / host-takeover
    # morons. (Only applied when phold left the category 'unknown' -- see
    # _infer_function_cat -- so phold's own category calls are never overwritten.)
    (re.compile(r'ribosom|zinc.uptake|\bzur\b|fur.family|\bkinase\b|'
                r'oxidoreductase|dehydrogenase|transaldolase|aldolase|racemase|'
                r'reductase|permease|dehydratase|\bsynthase\b|isomerase|epimerase|'
                r'\bmutase\b|\bphosphatase\b|sialyltransferase', re.I),
     "moron, auxiliary metabolic gene and host takeover"),
    (re.compile(r'antitoxin|RelB|HicB|MazE', re.I),
     "other"),   # toxin-antitoxin systems
]

# FoldSeek descriptions carrying no real function beyond "some protein". Used to
# decide when a specialized sub-DB hit (T6SS/T4SS/VFDB/CARD) should win over a
# CONFIDENT-but-generic FoldSeek call even at LOW phold confidence (TagO, IcmN).
# Deliberately narrow so SPECIFIC FS calls (e.g. vanT -> "Alanine racemase") do
# NOT match and keep winning.
_FS_GENERIC = re.compile(
    r'\b(putative|hypothetical|uncharacteri[sz]ed|unknown)\b.*\bprotein\b|'
    r'exported protein|lipoprotein|^duf\d+|^upf\d+|'
    r'^(putative |conserved )?(outer |inner )?membrane protein$', re.I)


def _infer_function_cat(desc: str, current_cat: str) -> str:
    """
    Infer a function category from the annotation description when the current
    category is unknown.  Trusts existing non-unknown PHold categories.
    """
    if current_cat and current_cat.lower() not in ("unknown function", "na", "nan"):
        return current_cat
    if not desc or not isinstance(desc, str):
        return current_cat or "unknown function"
    for pattern, cat in _FUNC_CAT_RULES:
        if pattern.search(desc):
            return cat
    return current_cat or "unknown function"


# -----------------------------------------------------------------------------
# CURATION RULES A, B, C, E (user-approved 2026-06-08, reconstructed from the
# original proposal text + the user's terse refinements -- see
# project_curation_subdb_gate_gap.md / memory for the lost-proposal-recovery
# note). Rule D (sub-DB -> phold_product) was implemented separately, upstream
# in 03_compare_annotations.py, per the user's clarification that it belongs
# "during phold handling, before any curation/comparison with FS".
#
# All four rules below act as PRE-CHECKS at the top of the "different"
# (divergent) branch of merge_annotations(): if a rule matches, the gene is
# auto-merged/resolved and returned immediately, never reaching the
# needs_review_divergent flag. Order: B (suppress) -> A (FS upgrade) ->
# C (top-3 corroboration) -> E (generic-vs-specific). Each docstring carries
# the user's verbatim refinement so future-you can judge edge cases.
# -----------------------------------------------------------------------------

# ---- Rule A: "FS specificity upgrade" -------------------------------------
# User: "ok lets try."
# Original proposal: several "different" pairs are cases where phold gives a
# generic PHROG-category-level structural term (e.g. "tail protein",
# "baseplate wedge subunit protein", "endolysin") and FoldSeek's top-3 hits
# independently CONVERGE on the same specific, named phage gene product
# (e.g. "Phage protein D", "Phage protein GP46", "...protein gp45"). These
# aren't disagreements -- FS is just naming the same thing more precisely.
# Detect: phold call matches a curated generic-structural-term pattern AND
# >=2 FS top-3 entries mutually share a specific (non-generic) token. Merge as
# "<phold generic term> (<FS specific name>)", attribute phold+foldseek,
# classify concordant (no review flag).
_GENERIC_STRUCTURAL_TERMS = re.compile(
    r'^(tail protein|tail fiber protein|'
    r'baseplate (wedge|spike|hub)?\s*(subunit\s*)?protein|'
    r'endolysin|holin|spanin|portal protein|'
    r'major capsid protein|minor capsid protein|head protein|'
    r'tape measure protein|terminase(\s+(large|small))?\s*subunit|'
    r'connector protein|neck protein)$', re.I)

# Words to additionally ignore when looking for FS top-3 internal concordance:
# top-3 "specific names" are frequently short bare IDs ("D", "GP46", "gp45")
# that tokenize()'s len>2 filter would drop as noise -- so we build our own
# lighter token set here rather than reusing tokenize() directly.
_RULE_A_DROP_WORDS = GENERIC_WORDS | {"phage", "prophage", "homolog", "to", "of"}


def _rule_a_top3_internal_id(desc: str) -> set:
    """Tokenise an FS top-3 entry for Rule A concordance checking, keeping
    short specific identifiers (e.g. 'D', 'GP46') that tokenize() would drop."""
    if not desc or not isinstance(desc, str):
        return set()
    raw = set(re.findall(r"[a-zA-Z0-9]+", desc.lower()))
    return {t for t in raw if t not in _RULE_A_DROP_WORDS and len(t) > 1}


def _fs_top3_concordant_specific(top3_str: str):
    """Rule A helper: return the most specific FS top-3 description that is
    corroborated by >=1 OTHER top-3 entry sharing a specific identifier token
    (e.g. both mention 'gp46' / 'd' / 'gp45'), or None if no such internal
    concordance exists. This is the signal that FS's top-3 calls are mutually
    consistent on a precise name (not just noise)."""
    entries = [_strip_fragment(t.strip()) for t in (top3_str or "").split("|")
               if t.strip() and t.strip() != "NA"]
    if len(entries) < 2:
        return None
    tok = [_rule_a_top3_internal_id(e) for e in entries]
    for i in range(len(entries)):
        if not tok[i]:
            continue
        for j in range(len(entries)):
            if j != i and tok[j] and (tok[i] & tok[j]):
                return entries[i]
    return None


# ---- Rule B: suppress near-uninformative "phage protein" wrappers ---------
# User: "yes, but only if alone, if 'baseplate phage protein' or something
#        informative then need to keep the informative part"
# Original proposal: "Prophage protein" / "Phage protein" / "Putative
# (pro)phage protein" carry no functional content beyond "this is a phage
# gene" -- when phold has a specific call, don't flag divergent, just keep
# phold's call. REFINEMENT (critical): this suppression applies ONLY when the
# FS description is JUST the wrapper phrase. If a qualifier survives in front
# of it (e.g. "baseplate phage protein" -- "baseplate" IS informative), the
# rule must NOT fire; the qualifier has to remain visible, so we fall through
# to normal divergent handling instead of suppressing it.
_PHAGE_WRAPPER_ONLY_RE = re.compile(
    # A bare "(pro)phage [associated/related/derived] protein", any "gpNN" / "Phage
    # protein GPNN" / "Homolog to ... gpNN" -- carries no function beyond "phage
    # gene". Must START with it (anchored), so "baseplate phage protein", "Phage
    # tail protein", "Phage transcriptional activator" are NOT matched.
    r'^(putative\s+|conserved\s+)?'
    r'((bacterio)?(pro)?phage[\s\-]*(associated|related|derived)?\s*protein\b.*'
    r'|gp\d+\w*'
    r'|phage\s+protein\s+(gp)?\d+.*'
    r'|homolog\s+to\s+.*\bgp\d+.*)$', re.I)


# Broader "FS hit carries no usable function" detector (superset of the wrapper):
# gp/HK97/Mu-like prophage proteins, generic "phage <descriptor> protein",
# bare ATP-binding/ATPase, "Prophage X protein NN", uncharacterised. Used to keep
# phold (no review) when FoldSeek only returns one of these. Carefully excludes
# specific phage names (tail/portal/baseplate/holin/terminase/capsid/connector/fiber).
_FS_GENERIC_BROAD = re.compile(
    r'^(putative\s+|conserved\s+)?('
    r'(bacterio)?(pro)?phage[\s\-]*(associated|related|derived)?\s*protein\b.*'
    r'|gp\d+\w*|\S*\s*gp\d+\b.*|mu-?like\s+prophage.*|.*\bhk97\b.*'
    r'|phage\s+(nucleotide-?binding|minor\s+structural|structural|virion\s+morphogenesis|derived|associated)\s+protein.*'
    r'|atp-?binding protein|atpase'
    r'|dna-?damage-?inducible protein.*'
    r'|prophage\s+\w+\s+protein\s+\d+'
    r'|.*\b(ig-?like|asch|atp-grasp|duf\d+|zinc[\s\-]?(ribbon|finger))\s+domain.?containing\s+protein.*'
    r'|pentapeptide\s+repeat.*'        # structural repeat fold, no function
    r'|tigr\d+\s+family\s+protein.*'   # uncharacterised TIGRFAM family
    r'|chp\d+.*'                        # "conserved hypothetical protein" id
    r'|(uncharacteri[sz]ed|hypothetical)\s+protein.*'
    r')$', re.I)


# PDB/structure-title detector: pdb100 descriptions are sometimes the full paper/
# structure title (a sentence) rather than a gene name -- carries no protein
# identity. E.g. LDGKBLMO_01121's GOOD hit "Biophysical and cellular
# characterisation of a junctional epitope antibody that locks IL-6 and gp80...".
_FS_PDB_TITLE_RE = re.compile(
    r'\b(structure|characteri[sz]ation|biophysical|cryo-?em|in complex with|'
    r'bound to|implications for|crystal structure of)\b', re.I)


def _is_pdb_title(desc: str) -> bool:
    """True if the FS description reads like a structure/paper title, not a gene
    name (>=8 words, or a tell-tale structural-biology phrase)."""
    if not desc or not isinstance(desc, str):
        return False
    return len(desc.split()) >= 8 or bool(_FS_PDB_TITLE_RE.search(desc))


def _fs_uninformative(desc: str) -> bool:
    """FoldSeek hit carries no usable function (wrapper, broad-generic, or a
    structure/paper title that is not a gene name)."""
    if not desc or not isinstance(desc, str):
        return False
    return bool(_FS_GENERIC_BROAD.match(desc.strip())) or _is_pdb_title(desc)


def _is_phage_wrapper_only(desc: str) -> bool:
    """True only if `desc` is JUST a generic phage-protein wrapper phrase with
    no qualifier in front (anchored match -- "baseplate phage protein" does
    NOT match, satisfying the user's "only if alone" refinement)."""
    if not desc or not isinstance(desc, str):
        return False
    return bool(_PHAGE_WRAPPER_ONLY_RE.match(desc.strip()))


# ---- Rule C: scan FS top-3 (not just top-1) for corroboration -------------
# User: "yes. but for anti-term vs antiterm there is also just the hyphen that
#        maybe blocked the agreement scan? This rule C can only work for
#        Confident sub FS hits (we dont want phold comparing to non equivalent
#        (in terms of confidence) hits"
# Original proposal: "Anti-termination protein Q-like" (phold) was flagged
# divergent against FS's generic top-1 "Antitermination protein", but FS's
# THIRD-ranked hit ("Antiterminator Q protein of prophage CP-933K") directly
# corroborates phold's "Q-like" call -- merge_annotations() only ever compares
# against fs_desc (top-1). REFINEMENTS: (1) normalise hyphens before matching
# -- "anti-term" vs "antiterm" / "anti-termination" vs "antitermination" would
# otherwise silently fail token overlap; (2) gate to FS confidence == CONFIDENT
# only, since comparing phold against a low-confidence FS hit is comparing
# across non-equivalent evidence tiers.
def _normalize_hyphens(text: str) -> str:
    """Strip hyphens so 'anti-termination'/'antitermination' and
    'anti-term'/'antiterm' tokenise identically (Rule C hyphen fix)."""
    if not text or not isinstance(text, str):
        return text
    return text.replace("-", "")


def _top3_corroborates(p_desc: str, top3_str: str):
    """Rule C helper: scan ALL FS top-3 entries (hyphen-normalised) for
    keyword/fuzzy overlap with PHold's description. Returns the corroborating
    top-3 description, or None."""
    p_tok = tokenize(_normalize_hyphens(p_desc))
    if not p_tok:
        return None
    for entry in (t.strip() for t in (top3_str or "").split("|")):
        if not entry or entry == "NA":
            continue
        e_clean = _strip_fragment(entry)
        e_norm = _normalize_hyphens(e_clean)
        e_tok = tokenize(e_norm)
        if not e_tok:
            continue
        if (p_tok & e_tok) and fuzzy_score(_normalize_hyphens(p_desc), e_norm) >= FUZZY_PARTIAL_THRESHOLD:
            return e_clean
    return None


# ---- Rule E: generic-vs-specific preference, regardless of source ---------
# User: "OK. be careful this blends in well with other gates, notably the
#        curated merges. Like helicase preferred over more generic ATPase"
# Original proposal: "CII-like transcriptional activator" (phold, a real named
# phage regulator family) vs. "DNA-binding protein (Fragment)" (FS, generic) --
# the existing strong/partial substring-preference logic never reaches these
# because they land in "different". Proposal: if one side IS a pure generic
# descriptor and the other names a specific family, auto-resolve to the
# specific one. REFINEMENT (critical): must not collide with the EXISTING
# specific-over-generic merges that already live in _UPGRADE_RULES (e.g.
# "DNA helicase / AAA family ATPase" -> "DNA helicase", "minor tail protein /
# TIGR04255" -> "minor tail protein") -- those are COMBINED single-source
# strings, resolved by pattern upgrade, not cross-source disagreements. Rule E
# therefore only fires on an EXACT membership match against a small, pure
# generic-descriptor set -- it can never match a combined "X / generic-Y"
# string (which contains "/" and extra text), so it cannot double-fire with
# _UPGRADE_RULES.
_GENERIC_DESCRIPTORS = frozenset({
    "dna-binding protein", "dna binding protein",
    "membrane protein", "transmembrane protein", "integral membrane protein",
    "metal-binding protein", "nucleotide-binding protein",
    "atp-binding protein", "atpase", "binding protein",
    "hydrolase", "transferase", "lipoprotein",
    # 2026-06: phold terms that are too generic to keep when FS names a specific
    # protein (validated on genome+citro: 9 good flips, 0 regressions). The
    # winner-must-be-specific guard in _generic_vs_specific stops a generic FS
    # wrapper from ever beating a specific phold call.
    "virion structural protein", "structural protein", "minor structural protein",
    "phage protein", "prophage protein", "putative phage protein",
    "phage-related protein", "phage derived protein", "conserved protein",
    "membrane-flanked domain", "transcriptional regulator",
})


def _generic_vs_specific(p_desc: str, fs_desc: str):
    """Rule E helper: if exactly one of (p_desc, fs_desc) is a PURE member of
    _GENERIC_DESCRIPTORS (after stripping a trailing "(Fragment)") and the
    other is informative/specific, return (specific_desc, specific_source).
    Else None. Exact-membership keeps this from ever matching combined
    "X / generic-Y" strings already handled by _UPGRADE_RULES (see docstring)."""
    p_key  = _strip_fragment(clean_str(p_desc)).strip().lower()
    fs_key = _strip_fragment(clean_str(fs_desc)).strip().lower()
    p_generic  = p_key in _GENERIC_DESCRIPTORS
    fs_generic = fs_key in _GENERIC_DESCRIPTORS

    def _specific(desc, key):
        # winner must be genuinely specific: informative, NOT itself in the
        # generic set, and NOT a bare (pro)phage-protein wrapper. This stops a
        # generic FS wrapper ("Putative phage protein") from beating a specific
        # phold call when phold happens to be in the generic set.
        return (is_informative(desc) and key not in _GENERIC_DESCRIPTORS
                and not _is_phage_wrapper_only(desc))

    if p_generic and not fs_generic and _specific(fs_desc, fs_key):
        return (fs_desc, "foldseek")
    if fs_generic and not p_generic and _specific(p_desc, p_key):
        return (p_desc, "phold")
    return None


# -----------------------------------------------------------------------------
# ANNOTATION MERGING
# -----------------------------------------------------------------------------

def merge_annotations(row: pd.Series) -> tuple:
    """
    Apply the curation decision tree for one gene.

    Returns:
      (final_desc: str,
       final_function_cat: str,
       source: str,
       curation_action: str,
       needs_review: bool,
       ai_suggestion: str,
       ai_explanation: str)
    """
    agreement    = str(row.get("agreement", "both_uninformative")).lower()
    p_desc       = clean_str(row.get("phold_product", "NA"))
    p_cat        = clean_str(row.get("phold_function_cat", "NA"))
    p_conf       = str(row.get("phold_confidence", "none")).lower()
    p_evalue     = row.get("phold_evalue")
    p_inf        = bool(row.get("phold_inf", False))
    p_tophit     = clean_str(row.get("phold_tophit", "NA"))
    p_method     = clean_str(row.get("phold_method", "NA"))
    sub_source   = clean_str(row.get("subdb_source", "NA")).lower()
    sub_name     = clean_str(row.get("subdb_name", "NA"))
    p_phrog_l    = str(row.get("phold_phrog", "")).lower()

    # Sub-DB hit (ACR/VFDB/CARD/NetFlaX/DefenseFinder): phold's global top1 search
    # often only records a GENERIC placeholder in `product` when one of these wins
    # ("VFDB virulence factor protein", "CARD resistance protein", or leaves it
    # blank for DefenseFinder) -- the real structured identity was joined in
    # 03_compare from sub_db_tophits/*.tsv (see _load_subdb_hits there) and is
    # surfaced here as subdb_source/subdb_name. We trust phold_phrog (not the
    # placeholder string or a regex on the raw tophit ID) to tell us which sub-DB
    # phold actually flagged this gene against -- it survives the NaN-cleaning
    # that can blank out phold_product for DefenseFinder hits.
    #
    # NOTE: this supersedes the old _parse_defensefinder_name()-based approach,
    # which regex-parsed raw DefenseFinder tophit ID strings (fragile, DF-only).
    # Reading the structured sub_db_tophits table directly is more robust and
    # generalises to all five sub-DBs uniformly.
    subdb_hit = is_informative(sub_name) and sub_source == p_phrog_l
    if subdb_hit:
        p_desc = sub_name
        p_inf  = True   # it IS informative: a named, structurally-identified protein
        if sub_source in ("defensefinder", "netflax", "acr"):
            # PHold maps these to the generic "moron, auxiliary metabolic gene and
            # host takeover" PHROG category; override to "defense", the
            # biologically correct category for these toxin/antitoxin & immune-
            # system protein families.
            p_cat = "defense"
    elif p_phrog_l == "defensefinder" and (
        p_desc.lower() in ("defensefinder protein", "defense protein")
        or "defensefinder" in p_method.lower()
    ):
        # Backward-compat fallback for comparison_per_gene.csv files generated
        # before the subdb_name join existed (no subdb_source/subdb_name columns,
        # or sub_db_tophits/*.tsv was missing/empty for this run). Keeps the
        # older regex-based DefenseFinder ID parsing as a safety net.
        df_name = _parse_defensefinder_name(p_tophit)
        if df_name:
            p_desc = df_name
            p_inf  = True
            p_cat  = "defense"

    fs_desc      = clean_str(row.get("foldseek_description", "NA"))
    fs_conf      = str(row.get("foldseek_confidence", "NO_HIT") or "NO_HIT")
    # foldseek_inf from 03_compare already has the WEAK filter applied (03_compare
    # sets fs_inf=False for WEAK/NO_HIT confidence hits before writing the CSV).
    # We keep the raw flag for re-routing logic below.
    _fs_inf_raw  = bool(row.get("foldseek_inf", False))
    # Safety net: WEAK/NO_HIT confidence hits must not drive annotation decisions
    # even if 03_compare classified them as informative (e.g. if CSV was generated
    # by an older version without the WEAK filter).
    fs_inf       = _fs_inf_raw and (fs_conf not in ("WEAK", "NO_HIT"))
    fs_evalue    = row.get("foldseek_evalue")
    fs_score     = row.get("foldseek_score")
    fs_pident    = row.get("foldseek_pident")
    fz_score     = row.get("fuzzy_score")
    top3         = clean_str(row.get("foldseek_top3", "NA"))
    defense_flag   = bool(row.get("foldseek_defense", False))
    fs_taxname     = clean_str(row.get("foldseek_taxname", "NA"))
    fs_db          = clean_str(row.get("foldseek_db", "NA"))
    fs_promiscuous = bool(row.get("foldseek_promiscuous", False))
    fs_euka_desc   = bool(row.get("foldseek_euka_desc", False))

    # Effective defense flag for the REVIEW decision: only force review when the
    # chosen phold or FS description is itself a defense/TA call. A defense_flag
    # set purely by an incidental TA hit in a non-chosen top-3 slot must not flag
    # an otherwise-agreeing gene (CEDLNMDB_00945). The verbose _defense_note still
    # appears in the explanation regardless.
    _eff_defense = defense_flag and (_desc_is_defense(p_desc) or _desc_is_defense(fs_desc))

    # ----- Agreement re-routing gates (applied in order) ----------------------

    # Gate 14: DF override made p_inf=True but 03_compare still saw
    # "defensefinder protein"/NaN as uninformative -> reclassify as phold_only.
    if p_inf and agreement == "both_uninformative" and not fs_inf:
        agreement = "phold_only"

    # Gate 15: WEAK/NO_HIT confidence means 03_compare may have set an agreement
    # that relied on an unreliable FS hit.  Re-route to prevent Case 4 from
    # auto-accepting poor structural evidence.
    # (Primary filtering is in 03_compare; this is a safety net for old CSVs.)
    if _fs_inf_raw and not fs_inf:
        # FS was informative at parse time but confidence filter demoted it
        if agreement == "foldseek_only":
            agreement = "both_uninformative"
        elif agreement == "different" and p_inf:
            agreement = "phold_only"

    # Gate 16: DF upgrade (Fix 1) set p_inf=True AFTER 03_compare classified this
    # as "foldseek_only" (PHold appeared uninformative then).  Now both sources
    # are informative -> the gene has two competing annotations and needs review.
    if p_inf and agreement == "foldseek_only" and fs_inf:
        agreement = "different"

    # Build suffix notes that appear in any explanation
    _defense_note = (
        " [DEFENSE SYSTEM: check DefenseFinder/PADLOC for exact protein name]"
        if defense_flag else ""
    )

    # -- Case 1: strong or partial match -------------------------------------
    if agreement in ("strong", "partial"):
        # Prefer PHold description: PHold searches phage-specific PHROG/PHAGE-DB
        # databases and is generally more phage-appropriate than generic PDB hits.
        # Exception: if FS description contains p_desc as a substring (i.e. FS is
        # a more qualified version of the same annotation), use the FS description.
        # E.g. PHold "tail protein" + FS "GpE family phage tail protein" -> use FS.
        if is_informative(p_desc) and is_informative(fs_desc):
            if p_desc.lower() in fs_desc.lower():
                best_desc = fs_desc   # FS is more specific
            else:
                best_desc = p_desc    # default: trust PHold
        else:
            best_desc = p_desc if is_informative(p_desc) else fs_desc
        best_desc = _apply_final_upgrades(best_desc, top3)
        best_cat  = _infer_function_cat(best_desc, p_cat if is_informative(p_cat) else "unknown function")
        ev_str = f"phold_evalue={p_evalue:.2e}" if _is_valid_float(p_evalue) else ""
        fs_ev  = f"foldseek_evalue={fs_evalue:.2e}" if _is_valid_float(fs_evalue) else ""
        ev_full = "; ".join(x for x in [ev_str, fs_ev] if x)
        flag_review = _eff_defense
        return (
            best_desc, best_cat,
            ("both agree" if agreement == "strong" else "merged"),
            f"auto_merge_{agreement}",
            flag_review,
            best_desc,
            f"{agreement.capitalize()} match (fuzzy={fz_score:.2f}). "
            f"PHold='{p_desc}', FoldSeek='{fs_desc}' (cat='{p_cat}'). "
            f"{ev_full}{_defense_note}",
        )

    # -- Case 2: complementary -----------------------------------------------
    # Triggers if 03_compare classified this as complementary (top3-aware), OR
    # if this script's own complementary check agrees.
    # Use PHold product (most phage-specific) unless FS best-hit is the direct
    # functional match.
    if agreement == "complementary" or (p_inf and fs_inf and is_complementary(p_cat, fs_desc)):
        # Prefer the specific FoldSeek call ONLY when PHold is a pure generic
        # descriptor (Rule E, e.g. PHold "transcriptional regulator" -> FS
        # "RinA"/"Excisionase"); otherwise keep the specific PHold name even when
        # the categories are complementary (e.g. PHold "ParA-like partition
        # protein" vs FS "Sporulation initiation inhibitor Soj" -> keep ParA).
        _ge = _generic_vs_specific(p_desc, fs_desc)
        if _ge and _ge[1] == "foldseek":
            best_desc = _apply_final_upgrades(_ge[0], top3)
            note_src  = (f"PHold product='{p_desc}' is generic; preferring the "
                         f"specific FoldSeek call '{_ge[0]}' (cat='{p_cat}').")
        else:
            best_desc = p_desc if is_informative(p_desc) else fs_desc
            best_desc = _apply_final_upgrades(best_desc, top3)
            note_src  = (f"PHold product='{p_desc}' (cat='{p_cat}') is consistent with "
                         f"FoldSeek context '{fs_desc}'. Using PHold name.")
        best_cat    = _infer_function_cat(best_desc, p_cat if is_informative(p_cat) else "unknown function")
        flag_review = _eff_defense
        conf_note   = f" [FS confidence: {fs_conf}]" if fs_conf not in ("CONFIDENT","GOOD") else ""
        return (
            best_desc, best_cat,
            ("both agree" if agreement == "strong" else "merged"),
            "auto_merge_complementary",
            flag_review,
            best_desc,
            f"Complementary annotations. {note_src}{conf_note}"
            f"{_defense_note}",
        )

    # -- Case 3: phold_only ---------------------------------------------------
    if agreement == "phold_only":
        # PHold "low" confidence can still have good evalue/bitscore -- don't
        # demote it; just label the confidence level in the explanation.
        conf_label = p_conf if p_conf in PHOLD_TRUSTED_CONF | PHOLD_WEAK_CONF else "unvalidated"
        ev_str = f"; evalue={p_evalue:.2e}" if _is_valid_float(p_evalue) else ""
        best_desc = _apply_final_upgrades(p_desc, top3)
        best_cat  = _infer_function_cat(best_desc, p_cat if is_informative(p_cat) else "unknown function")
        return (
            best_desc, best_cat,
            "phold",
            "auto_phold_only",
            False,
            best_desc,
            f"PHold annotation only (confidence={conf_label}{ev_str}). "
            f"No informative FoldSeek hit. PHROG={row.get('phold_phrog', 'NA')}.",
        )

    # -- Case 4: foldseek_only ------------------------------------------------
    if agreement == "foldseek_only":
        ev_str  = f"; evalue={fs_evalue:.2e}" if _is_valid_float(fs_evalue) else ""
        sc_str  = f"; score={fs_score:.0f}" if _is_valid_float(fs_score) else ""
        pid_str = f"; pident={fs_pident:.1f}%" if _is_valid_float(fs_pident) else ""
        best_desc = _apply_final_upgrades(fs_desc, top3)
        best_cat  = _infer_function_cat(best_desc, "unknown function")

        # 2026-06: a promiscuous-fold / eukaryotic FoldSeek hit that is BORDERLINE
        # and the only evidence is too weak to even send to review -- it's a fold-level
        # false positive (beta-lactamase/TLD, BCL-6/ankyrin). Treat as no_informative_hit
        # (no flag), unless it's a defense hit. A CONFIDENT promiscuous/euka hit still
        # goes to review below (worth a look). Plain (non-promiscuous) borderline hits are
        # left untouched -> still auto-annotated, so this does NOT over-demote.
        if (fs_promiscuous or fs_euka_desc) and fs_conf == "BORDERLINE" and not defense_flag:
            return (
                "hypothetical protein", "unknown function",
                "no_informative_hit",
                "both_uninformative",
                False,
                "hypothetical protein",
                f"Only a BORDERLINE promiscuous/eukaryotic FoldSeek hit ('{fs_desc}', "
                f"{fs_taxname}) and no phold call -- fold-level false positive; treated "
                f"as no informative hit.",
            )

        # Promiscuous-fold gate: known structurally promiscuous folds (TLD/MBL,
        # motor proteins, eukaryotic ubiquitination machinery) should NOT be
        # auto-annotated even when FoldSeek confidence is high.  The structural
        # similarity is real but the FUNCTION cannot be inferred without top-3
        # agreement (Aravind 1999 PMID 11471255; Daiyasu 2001 PMID 11513844).
        # Eukaryotic-description gate: description keywords indicate a eukaryote-
        # specific protein (inferred from text when DB lacks embedded taxonomy).
        # In both cases, route to manual review with a clear explanation.
        if fs_promiscuous or fs_euka_desc:
            prom_note = ""
            if fs_promiscuous:
                prom_note += (
                    " [PROMISCUOUS FOLD: this description matches a known "
                    "structurally promiscuous fold (e.g. TLD/beta-lactamase, "
                    "motor protein). Structural similarity does not imply "
                    "functional identity — verify top-3 agreement before "
                    "accepting this annotation. Ref: Aravind 1999 PMID 11471255.]"
                )
            if fs_euka_desc:
                prom_note += (
                    " [EUKARYOTIC DESCRIPTION: description suggests a eukaryote-"
                    "specific protein despite phage context — likely a fold-level "
                    "false positive from afdb-swissprot/pdb100 (no embedded taxonomy).]"
                )
            return (
                best_desc, best_cat,
                "foldseek",
                "needs_review_promiscuous_fold",
                True,   # flag for review
                best_desc,
                f"FoldSeek annotation only (PHold no hit){ev_str}{sc_str}{pid_str}. "
                f"Organism: {fs_taxname}. DB: {fs_db}. [FS confidence: {fs_conf}]"
                f"{_defense_note}{prom_note}",
            )

        flag_review = defense_flag
        action = "needs_review_defense" if defense_flag else "auto_foldseek_only"
        return (
            best_desc, best_cat,
            "foldseek",
            action,
            flag_review,
            best_desc,
            f"FoldSeek annotation only (PHold no hit){ev_str}{sc_str}{pid_str}. "
            f"Organism: {fs_taxname}. DB: {fs_db}. [FS confidence: {fs_conf}]{_defense_note}",
        )

    # -- Case 5: different (both informative, diverge) -------------------------
    if agreement == "different":
        # Rule DEF (2026-06): a confident DEFENSE-system phold call (DefenseFinder /
        # PHROG defense) is HMM-validated and trusted -- keep it, no review, even when
        # FoldSeek names a different fold (transposase, RloG, ...). Defense modules
        # routinely reuse mobile-element folds, so a divergent FS is expected here.
        if p_cat == "defense" and is_informative(p_desc) and p_conf in ("high", "medium"):
            best_desc = _apply_final_upgrades(p_desc, top3)
            return (
                best_desc, "defense", "phold", "auto_phold_defense", False, best_desc,
                f"Confident defense-system call '{p_desc}' (DefenseFinder) kept over "
                f"FoldSeek '{fs_desc}'.{_defense_note}",
            )

        # Rule ATPase-family (2026-06: "DNA helicase still in review! unflag"):
        # phold helicase / ATP-dependent enzyme + FoldSeek AAA(+)-ATPase are the SAME
        # superfamily (helicases ARE AAA+ ATPases) -> not divergent, keep phold, no flag.
        # Belt-and-suspenders: the complementary map already covers this, but this rule
        # guarantees it even if the complementary path isn't reached.
        # ABC transporters and SMC/condensin proteins are also AAA(+)-family ATPases,
        # so phold "ABC transporter" vs FS "Chromosome segregation protein SMC" /
        # "AAA family ATPase" is the same nucleotide-binding superfamily, not a
        # disagreement ("transporter can be ATPases" -> CEDLNMDB_00938).
        if (re.search(r'helicase|atp-?dependent|\babc\b|transporter', str(p_desc), re.I)
                and re.search(r'\baaa\b|atpase|atp-?binding|\bsmc\b|chromosome segregation',
                              str(fs_desc), re.I)
                and not re.search(r'transpos', str(p_desc), re.I)):
            best_desc = _apply_final_upgrades(p_desc, top3)
            best_cat  = _infer_function_cat(
                best_desc, p_cat if is_informative(p_cat) else "DNA, RNA and nucleotide metabolism")
            return (
                best_desc, best_cat, "phold", "auto_phold_atpase_family", False, best_desc,
                f"phold '{p_desc}' and FoldSeek '{fs_desc}' are the same ATPase superfamily "
                f"(helicases are AAA+ ATPases); kept phold, no review.{_defense_note}",
            )

        # Rule F (2026-06): a BORDERLINE/WEAK FoldSeek hit must NOT override or be
        # flagged-divergent against an informative phold call. phold (even low conf)
        # is the more phage-appropriate evidence; only CONFIDENT/GOOD FS competes.
        # Big reducer of false "needs_review" entries. (sub-DB Rule 0 runs after, but
        # its p_desc is already the sub-DB name, so keeping phold here is consistent.)
        _fs_weak = fs_conf not in ("CONFIDENT", "GOOD")
        _fs_generic = _fs_uninformative(fs_desc)
        if (_fs_weak or _fs_generic) and is_informative(p_desc) \
                and not _is_phage_wrapper_only(p_desc) and not subdb_hit:
            best_desc = _apply_final_upgrades(p_desc, top3)
            best_cat  = _infer_function_cat(best_desc, p_cat if is_informative(p_cat) else "unknown function")
            _why = "below GOOD confidence" if _fs_weak else "a generic/gp/Mu-like wrapper"
            return (
                best_desc, best_cat, "phold", "auto_phold_over_weak_fs", False, best_desc,
                f"FoldSeek hit '{fs_desc}' is {_why}; kept phold '{p_desc}' without review.{_defense_note}",
            )

        # Rule TA (2026-06): a phold NetFlaX / toxin-antitoxin hit gives the ROLE
        # (toxin/antitoxin) but FoldSeek often names the actual protein -- merge them
        # ("AbrB family transcriptional regulator (antitoxin protein)") and don't
        # flag. Recurring pattern: FS functionalises phold TA/NetFlaX calls.
        _is_ta_phold = (p_phrog_l in ("netflax",) or
                        re.search(r'\b(anti)?toxin\b|netflax|abrb', str(p_desc), re.I))
        if _is_ta_phold and is_informative(fs_desc) and not _fs_uninformative(fs_desc):
            _role = re.search(r'\b(antitoxin|toxin)\b', str(p_desc), re.I)
            role_str = (_role.group(1).lower() + " protein") if _role else clean_str(p_desc)
            merged = f"{_apply_final_upgrades(fs_desc, top3)} ({role_str})"
            return (
                merged, "moron, auxiliary metabolic gene and host takeover",
                "merged", "auto_merge_TA_functionalised", False, merged,
                f"phold TA/NetFlaX role '{p_desc}' functionalised by FoldSeek "
                f"'{fs_desc}'.{_defense_note}",
            )

        # Rule TR (2026-06): phold "DNA transposition"/transposase + FoldSeek AAA-ATPase
        # / NTPase-KAP / P-loop NTPase -> the FS hit IS the transposition motor ATPase;
        # merge "DNA transposition (<motor>)" instead of flagging divergent.
        _is_transpos = bool(re.search(r'transpos', str(p_desc), re.I))
        _fs_ntpase = bool(re.search(r'\baaa\b.*atpase|aaa\s+family\s+atpase|ntpase\s+kap'
                                    r'|kap\s+family|p-?loop', str(fs_desc), re.I))
        if _is_transpos and _fs_ntpase:
            motor = ("AAA ATPase" if re.search(r'aaa', str(fs_desc), re.I)
                     else "NTPase KAP" if re.search(r'kap', str(fs_desc), re.I)
                     else "P-loop NTPase")
            merged = f"DNA transposition ({motor})"
            return (
                merged, "integration and excision",
                "merged", "auto_merge_transposition_motor", False, merged,
                f"phold transposition call '{p_desc}' + FoldSeek motor '{fs_desc}' "
                f"merged as a transposition ATPase.{_defense_note}",
            )

        # Rule MORPH (2026-06): two virion structural/assembly calls (e.g. tail
        # assembly chaperone + portal) are parts of the same morphogenesis module --
        # merge "phold (FS)" rather than flag divergent.
        _MORPHO = re.compile(
            r'\b(tail|head|portal|capsid|baseplate|tape|measure|major|minor|'
            r'fib(er|re)|sheath|collar|neck|scaffold|chaperone|prohead|assembly)\b', re.I)
        if (is_informative(p_desc) and is_informative(fs_desc)
                and _MORPHO.search(str(p_desc)) and _MORPHO.search(str(fs_desc))):
            fs_short = re.sub(r'^phage\s+', '', clean_str(fs_desc), flags=re.I)
            merged = f"{clean_str(p_desc)} ({fs_short})"
            best_cat = _infer_function_cat(p_desc, p_cat if is_informative(p_cat) else "tail")
            return (
                merged, best_cat, "merged", "auto_merge_morphogenesis", False, merged,
                f"Structural morphogenesis calls merged: phold '{p_desc}' + "
                f"FoldSeek '{fs_desc}'.{_defense_note}",
            )

        # Rule TAL (2026-06): PHold 'transaldolase' vs FS panB / 3-methyl-2-
        # oxobutanoate hydroxymethyltransferase -> keep transaldolase, no flag.
        if re.search(r'transaldolase', str(p_desc), re.I) and \
           re.search(r'3-methyl-2-oxobutanoate hydroxymethyltransferase|panB', str(fs_desc), re.I):
            return (
                "transaldolase",
                _infer_function_cat("transaldolase", p_cat if is_informative(p_cat) else "other"),
                "phold", "auto_phold_transaldolase", False, "transaldolase",
                f"Kept PHold 'transaldolase' over FoldSeek fold-level '{fs_desc}'.{_defense_note}",
            )

        # ---- Rules B / A / C / E pre-checks (user-approved 2026-06-08) -----
        # If any of these fire, the gene is auto-resolved here and never
        # reaches the needs_review_divergent flag below. Order matters:
        # B (cheap pattern suppression) -> A (FS specificity upgrade via
        # top-3 concordance) -> C (CONFIDENT-gated top-3 corroboration) ->
        # E (generic-vs-specific, narrowest/most conservative). See the rule
        # docstrings above merge_annotations() for the user's verbatim
        # refinements and the reasoning behind each gate.

        # Rule 0 (2026-06, low-conf sub-DB): a specialized sub-DB hit
        # (T6SS/T4SS/VFDB/CARD/...) was adopted into p_desc, but phold confidence
        # was low so 03 classified this "different" and a generic FoldSeek call
        # ("Putative ... protein", "lipoprotein") would otherwise win. When the FS
        # call is GENERIC, trust the structured sub-DB name even at low phold
        # confidence (TagO, IcmN). The _FS_GENERIC guard keeps SPECIFIC FS calls
        # (vanT -> "Alanine racemase") winning, so this never over-fires.
        if subdb_hit and (_fs_uninformative(fs_desc) or _FS_GENERIC.search(fs_desc or "")):
            best_desc = _apply_final_upgrades(p_desc, top3)
            best_cat  = _infer_function_cat(best_desc, p_cat if is_informative(p_cat) else "unknown function")
            return (
                best_desc, best_cat,
                "phold",
                "auto_merge_subdb_over_generic_fs",
                False,
                best_desc,
                f"Specialized sub-DB hit '{p_desc}' ({sub_source}) preferred over "
                f"generic FoldSeek call '{fs_desc}' despite low phold confidence "
                f"[sub-DB low-conf rule].{_defense_note}",
            )

        # Rule B: FS hit is JUST a "(pro)phage protein" wrapper -- no function
        # beyond "this is a phage gene" -- and PHold has a specific call. Keep
        # PHold's call; do NOT flag divergent. ("only if alone": an
        # informative qualifier in front, e.g. "baseplate phage protein",
        # fails the anchored regex and falls through untouched.)
        if _is_phage_wrapper_only(fs_desc) and is_informative(p_desc):
            best_desc = _apply_final_upgrades(p_desc, top3)
            best_cat  = _infer_function_cat(best_desc, p_cat if is_informative(p_cat) else "unknown function")
            return (
                best_desc, best_cat,
                "phold",
                "auto_merge_ruleB_suppress_fs_wrapper",
                False,
                best_desc,
                f"FoldSeek hit '{fs_desc}' is a near-uninformative phage-protein "
                f"wrapper (no function beyond 'this is a phage gene'); kept "
                f"PHold's specific call '{p_desc}' (cat='{p_cat}'). [Rule B]"
                f"{_defense_note}",
            )

        # Rule A: PHold gives a generic PHROG-category-level structural term
        # (e.g. "tail protein", "endolysin") and >=2 of FS's top-3 hits
        # mutually converge on the same specific named gene product (e.g.
        # both top-3 entries mention "gp46"). Not a disagreement -- FS is
        # just more precise. Merge as "<phold generic> (<FS specific>)".
        if is_informative(p_desc) and _GENERIC_STRUCTURAL_TERMS.match(p_desc):
            fs_specific = _fs_top3_concordant_specific(top3)
            if fs_specific:
                merged = f"{p_desc} ({fs_specific})"
                best_cat = _infer_function_cat(p_desc, p_cat if is_informative(p_cat) else "unknown function")
                return (
                    merged, best_cat,
                    ("both agree" if agreement == "strong" else "merged"),
                    "auto_merge_ruleA_fs_specificity_upgrade",
                    False,
                    merged,
                    f"PHold's generic structural call '{p_desc}' and FoldSeek's "
                    f"internally-concordant top-3 (independently converging on "
                    f"'{fs_specific}') describe the same protein at different "
                    f"levels of specificity -- not a disagreement. "
                    f"Top-3 FS: {top3}. [Rule A]{_defense_note}",
                )

        # Rule C: scan ALL of FS's top-3 (hyphen-normalised), not just top-1,
        # for corroboration of PHold's call -- e.g. PHold "Anti-termination
        # protein Q-like" vs FS top-1 "Antitermination protein" (flagged
        # divergent) but FS top-3 #3 "Antiterminator Q protein of prophage
        # CP-933K" actually confirms it; the literal hyphen difference
        # ("anti-term" vs "antiterm") was silently blocking the match.
        # Gated to FS confidence == CONFIDENT only (user: "we dont want phold
        # comparing to non equivalent (in terms of confidence) hits").
        # Rule C scans ALL of FS's top-3 (not just top-1) for corroboration of
        # PHold's call, INCLUDING low-confidence phold (user: "i did say i wanted
        # to include phold low"). Cross-method agreement (phold + any FS top-3 hit
        # naming the same thing) is treated as concordance. AraC case: phold low
        # "kinase" + FS top-3 "histidine kinase" -> keep "kinase" (two independent
        # methods agree on kinase). Gated to fs_conf==CONFIDENT so we never compare
        # against an overall-unreliable FS hit. (Per-hit FS confidence isn't in the
        # top-3 string, so we can't require per-entry confidence parity.)
        if fs_conf == "CONFIDENT" and is_informative(p_desc):
            corroborating = _top3_corroborates(p_desc, top3)
            if corroborating:
                best_desc = _apply_final_upgrades(p_desc, top3)
                best_cat  = _infer_function_cat(best_desc, p_cat if is_informative(p_cat) else "unknown function")
                return (
                    best_desc, best_cat,
                    ("both agree" if agreement == "strong" else "merged"),
                    "auto_merge_ruleC_top3_corroboration",
                    False,
                    best_desc,
                    f"PHold='{p_desc}' (cat='{p_cat}') is corroborated by "
                    f"FoldSeek's top-3 hit '{corroborating}' once hyphen "
                    f"variants are normalised (FS confidence=CONFIDENT). "
                    f"Treating as concordant rather than divergent. "
                    f"Top-3 FS: {top3}. [Rule C]{_defense_note}",
                )

        # Rule E: one side is a PURE generic descriptor ("DNA-binding
        # protein", "membrane protein", ...) and the other names a specific
        # protein/family (e.g. PHold "CII-like transcriptional activator" vs
        # FS "DNA-binding protein (Fragment)") -- prefer the specific call,
        # with correct source attribution. Exact-membership check only, so
        # this can never collide with _UPGRADE_RULES' combined-string
        # resolutions (e.g. "DNA helicase / AAA family ATPase" -> "DNA
        # helicase" -- that string contains "/" + extra text and will never
        # match _GENERIC_DESCRIPTORS by exact membership).
        rule_e = _generic_vs_specific(p_desc, fs_desc)
        if rule_e:
            specific_desc, specific_src = rule_e
            generic_desc = fs_desc if specific_src == "phold" else p_desc
            best_desc = _apply_final_upgrades(specific_desc, top3)
            best_cat  = _infer_function_cat(best_desc, p_cat if is_informative(p_cat) else "unknown function")
            return (
                best_desc, best_cat,
                specific_src,
                "auto_merge_ruleE_specific_over_generic",
                False,
                best_desc,
                f"One side is a pure generic descriptor ('{generic_desc}'), the "
                f"other names a specific protein ('{specific_desc}', "
                f"source={specific_src}); preferring the specific call. [Rule E]"
                f"{_defense_note}",
            )
        # ---- end Rules B/A/C/E pre-checks -----------------------------------

        fs_strong = (_is_valid_float(fs_score) and float(fs_score) >= 300 and
                     _is_valid_float(fs_evalue) and float(fs_evalue) <= 1e-5)
        p_trusted = p_conf in PHOLD_TRUSTED_CONF

        # Apply fragment/upgrade rules to both candidates
        p_final = _apply_final_upgrades(p_desc, top3)
        fs_final = _apply_final_upgrades(fs_desc, top3)

        # Attribute source based on which evidence the suggestion actually uses.
        # This is the key fix: rows in needs_review still carry proper source
        # attribution so the final table's annotation_source is phold / foldseek /
        # phold+foldseek (not the generic "flagged" / "AI_suggestion").
        if fs_strong and p_trusted:
            # Dedup: exact equality OR one description contains the other
            # (e.g. "panB (3-methyl-...)" contains "3-methyl-...") -> use the longer
            if p_final == fs_final:
                suggestion = p_final
                divergent_source = ("both agree" if agreement == "strong" else "merged")   # same call from both
            elif fs_final.lower() in p_final.lower():
                suggestion = p_final   # p_final is more specific
                divergent_source = ("both agree" if agreement == "strong" else "merged")
            elif p_final.lower() in fs_final.lower():
                suggestion = fs_final  # fs_final is more specific
                divergent_source = ("both agree" if agreement == "strong" else "merged")
            else:
                suggestion = f"{p_final} / {fs_final}"
                divergent_source = ("both agree" if agreement == "strong" else "merged")   # combined annotation
            explanation = (
                f"DIVERGENT: both PHold (conf={p_conf}) and FoldSeek "
                f"(score={fs_score:.0f}, evalue={fs_evalue:.2e}) are informative but disagree. "
                f"PHold='{p_desc}' (cat={p_cat}), FS='{fs_desc}' ({fs_taxname}). "
                f"Top-3 FS: {top3}. Suggestion: '{suggestion}' "
                f"(review if they describe the same function)."
                f"{_defense_note}"
            )
        elif fs_strong and not p_trusted:
            suggestion = fs_final
            divergent_source = "foldseek"
            explanation = (
                f"DIVERGENT: FoldSeek (score={fs_score:.0f}, evalue={fs_evalue:.2e}) "
                f"more reliable (PHold conf='{p_conf}'). "
                f"PHold='{p_desc}', FS='{fs_desc}'. "
                f"Suggestion: use FoldSeek '{fs_final}'."
                f"{_defense_note}"
            )
        elif p_trusted and not fs_strong:
            suggestion = p_final
            divergent_source = "phold"
            explanation = (
                f"DIVERGENT: PHold conf='{p_conf}' more reliable "
                f"(FS score={fs_score}, evalue={fs_evalue}). "
                f"PHold='{p_desc}' (cat={p_cat}), FS='{fs_desc}' ({fs_taxname}). "
                f"Top-3 FS: {top3}. Suggestion: use PHold '{p_final}'."
                f"{_defense_note}"
            )
        else:
            # Neither source is clearly stronger -- prefer PHold as more
            # phage-specific, but use FS if it is the only informative call.
            if is_informative(fs_desc) and not is_informative(p_desc):
                suggestion = fs_final
                divergent_source = "foldseek"
            else:
                suggestion = p_final
                divergent_source = "phold"
            explanation = (
                f"DIVERGENT -- low confidence both sources. "
                f"PHold='{p_desc}' (conf={p_conf}), "
                f"FS='{fs_desc}' (score={fs_score}, evalue={fs_evalue}, "
                f"taxname={fs_taxname}). "
                f"Top-3 FS: {top3}. Manual review strongly recommended."
                f"{_defense_note}"
            )
        # Apply upgrade rules to the assembled suggestion (catches combined strings
        # like "DNA helicase / AAA family ATPase" and "minor tail / TIGR04255")
        suggestion = _apply_final_upgrades(suggestion, top3)
        final_cat = _infer_function_cat(suggestion, p_cat if is_informative(p_cat) else "unknown function")

        # 2026-06: only REVIEW genuine disagreements. If phold & FS (top-1 or any
        # top-3) share a specific functional token, OR phold's category matches FS,
        # they aren't really divergent -> mark "both agree" and DON'T flag review
        # (big reducer of the review file). Genuine unrelated disagreements stay flagged.
        _STOP = {"protein", "domain", "family", "like", "dna", "rna", "binding",
                 "putative", "phage", "prophage", "type", "system", "containing",
                 "subunit", "terminal", "associated", "related", "homolog"}
        def _split_tokens(s):
            # split on every non-alphanumeric WITHOUT hyphen-merging, so compound
            # enzyme names keep their parts ('metallo','protease') for the
            # substring/prefix tests below.
            return {t for t in re.split(r"[^a-z0-9]+", str(s).lower())
                    if t and t not in _STOP}
        def _common_prefix_len(x, y):
            n = 0
            for cx, cy in zip(x, y):
                if cx != cy:
                    break
                n += 1
            return n
        def _shares_function(a, b):
            # 1) hyphen-normalised token overlap (anti-term ~ antitermination).
            ta = tokenize(_normalize_hyphens(a)) - _STOP
            tb = tokenize(_normalize_hyphens(b)) - _STOP
            if ta & tb:
                return True
            # 2) split-token substring OR shared-prefix>=6 -- catches compound and
            # typo'd enzyme words that hyphen-merge would otherwise hide:
            # metallo-protease ~ metallopeptidase ('metallo' prefix of both),
            # deoxyribosyltransferase ~ "Deoxyribosyltansferase" (typo, 13-char prefix).
            sa, sb = _split_tokens(a), _split_tokens(b)
            for x in sa:
                if len(x) < 5:
                    continue
                for y in sb:
                    if len(y) < 5:
                        continue
                    if x in y or y in x or _common_prefix_len(x, y) >= 6:
                        return True
            # 3) synonym groups sourced from COMPLEMENTARY_CATEGORY_MAP value-lists
            # (reuse the curated map symmetrically; category label ignored here):
            # Fe-S~oxidoreductase, glycosidase~amidase, ParA~Soj, Fur~LexA, etc.
            la, lb = str(a).lower(), str(b).lower()
            for _kws in COMPLEMENTARY_CATEGORY_MAP.values():
                if any(k in la for k in _kws) and any(k in lb for k in _kws):
                    return True
            return False
        _t3_share = any(_shares_function(p_desc, t) for t in str(top3).split("|"))
        related = (is_informative(p_desc) and is_informative(fs_desc) and
                   (_shares_function(p_desc, fs_desc) or _t3_share
                    or is_complementary(p_cat, fs_desc)))
        if related:
            # both methods point at the same biology -> keep the single cleaner
            # phold name instead of a verbose "phold / foldseek" string
            # (e.g. "metallo-protease", not "metallo-protease / Phage metallopeptidase...").
            if is_informative(p_final):
                suggestion = p_final
                final_cat  = _infer_function_cat(
                    suggestion, p_cat if is_informative(p_cat) else final_cat)
            divergent_source = "both agree" if suggestion in (p_final, fs_final) else "merged"
        flag_review = (not related) or _eff_defense
        action = "auto_merge_related" if related and not _eff_defense else "needs_review_divergent"
        return (
            suggestion, final_cat,
            divergent_source,
            action,
            flag_review,
            suggestion,
            explanation,
        )

    # -- Case 6: both_uninformative -------------------------------------------
    # Distinguish TRUE dark matter (no FoldSeek hit at all) from a hit that exists
    # but is uninformative (DUF/empty/below-threshold) -> separate annotation_source
    # so diagnostics can tell "nothing found" from "found a fold, no functional name".
    fs_had_hit = (str(fs_conf) not in ("NO_HIT", "", "nan", "none", "None")) or bool(_fs_inf_raw)
    src = "no_informative_hit" if fs_had_hit else "no_hit"
    weak_note = (f" FoldSeek returned a WEAK/unreliable hit ({fs_desc}) -- filtered out."
                 if _fs_inf_raw and fs_conf in ("WEAK",) and is_informative(fs_desc) else "")
    return (
        "hypothetical protein", "unknown function",
        src,
        "both_uninformative",
        False,
        "hypothetical protein",
        f"Neither PHold nor FoldSeek found an informative annotation.{weak_note}",
    )


def _is_valid_float(v) -> bool:
    import math
    try:
        f = float(v)
        return not math.isnan(f) and not math.isinf(f)
    except (TypeError, ValueError):
        return False


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():
    section("STEP 04 -- CURATE ANNOTATIONS")
    CURATION_DIR.mkdir(parents=True, exist_ok=True)

    # -- Load comparison table -------------------------------------------------
    comp_path = COMPARISON_DIR / "comparison_per_gene.csv"
    if not comp_path.exists():
        log(f"ERROR: {comp_path} not found. Run 03_compare_annotations.py first.")
        sys.exit(1)

    comp = pd.read_csv(str(comp_path))
    log(f"Loaded comparison table: {len(comp)} genes")
    log(f"Columns: {list(comp.columns)}")

    # -- Apply curation logic --------------------------------------------------
    section("APPLYING CURATION DECISION TREE")
    results = []
    for _, row in comp.iterrows():
        (final_desc, final_cat, source, action,
         flag, suggestion, explanation) = merge_annotations(row)

        results.append({
            **row.to_dict(),   # keep all comparison columns
            "final_product":    final_desc,
            "final_function":   final_cat,
            "best_source":      source,
            "curation_action":  action,
            "needs_review":     flag,
            "AI_suggestion":    suggestion,
            "AI_explanation":   explanation,
            "final_annotation": "",  # <- user fills this in for needs_review rows
        })

    curated = pd.DataFrame(results)

    # -- Split into auto-curated and needs-review ------------------------------
    auto_curated = curated[curated["needs_review"] == False].copy()
    needs_review = curated[curated["needs_review"] == True].copy()

    # For needs_review, pre-fill final_annotation with AI_suggestion so user
    # can simply accept or edit inline
    needs_review["final_annotation"] = needs_review["AI_suggestion"]

    # -- Statistics ------------------------------------------------------------
    section("CURATION STATISTICS")
    n_total = len(curated)
    n_auto  = len(auto_curated)
    n_flag  = len(needs_review)

    log(f"Total hypothetical genes  : {n_total}")
    log(f"Auto-curated (no review)  : {n_auto}  ({100*n_auto//n_total}%)")
    log(f"Flagged for review        : {n_flag}  ({100*n_flag//n_total}%)")

    log("\nAuto-curated action breakdown:")
    log(auto_curated["curation_action"].value_counts().to_string())

    if n_flag > 0:
        log(f"\nFlagged genes ({n_flag}):")
        for _, row in needs_review.iterrows():
            log(f"  {row['locus_tag']} ({row['prophage']}): "
                f"phold='{str(row['phold_product'])[:40]}' vs "
                f"fs='{str(row['foldseek_description'])[:40]}'")

    # Informativeness after curation
    n_now_inf = curated["final_product"].apply(is_informative).sum()
    n_still_hypo = (curated["final_product"] == "hypothetical protein").sum()
    log(f"\nAfter curation:")
    log(f"  Genes with informative annotation : {n_now_inf} / {n_total}  "
        f"({100*n_now_inf//n_total}%)")
    log(f"  Still 'hypothetical protein'      : {n_still_hypo} / {n_total}")

    # Per-prophage annotation rate
    log("\nPer-prophage informativeness after curation:")
    pp = curated.groupby("prophage").apply(
        lambda g: pd.Series({
            "n": len(g),
            "annotated": g["final_product"].apply(is_informative).sum(),
        })
    )
    if not pp.empty and "annotated" in pp.columns:
        pp["pct"] = (100 * pp["annotated"] / pp["n"]).round(0).astype(int)
        log(pp.to_string())
    else:
        log("  (no per-prophage data — prophage column may be unset)")

    # -- Save outputs ---------------------------------------------------------
    # 2026-06: emit ONE combined curated file (all genes + needs_review flag)
    # into the shared general-output folder (04_output). step05 reads this and emits
    # the single review_suggested.csv + final table there -- no duplicate review file.
    CURATION_DIR.mkdir(parents=True, exist_ok=True)
    curated_out = CURATION_DIR / "curated_annotations.csv"

    # final_annotation: review rows keep AI_suggestion (user edits), auto rows = final_product
    curated = curated.copy()
    curated["final_annotation"] = curated.apply(
        lambda r: (r["AI_suggestion"] if r["needs_review"] else r["final_product"]), axis=1)

    # Reorder columns so the curated file reads like the final table (key identity +
    # decision columns first, then phold evidence, then foldseek evidence, then rest).
    _lead = ["prophage", "locus_tag", "aa_length", "final_product", "final_function",
             "final_annotation", "best_source", "needs_review", "curation_action",
             "agreement", "pharokka_function"]
    _phold = [c for c in curated.columns if c.startswith("phold_") or c.startswith("subdb_")]
    _fs    = [c for c in curated.columns if c.startswith("foldseek_") or c.startswith("fuzzy")]
    _ai    = [c for c in curated.columns if c.startswith("AI_")]
    _seen  = set(_lead + _phold + _fs + _ai)
    _rest  = [c for c in curated.columns if c not in _seen]
    _order = [c for c in (_lead + _phold + _fs + _rest + _ai) if c in curated.columns]
    curated = curated[_order]

    curated.to_csv(str(curated_out), index=False)
    log(f"\ncurated_annotations.csv -> {curated_out}  ({len(curated)} rows, "
        f"{len(needs_review)} flagged for review; review_suggested.csv written by step05)")

    # -- Instructions for manual review ----------------------------------------
    if n_flag > 0:
        section("MANUAL REVIEW INSTRUCTIONS")
        log(f"{n_flag} gene(s) flagged. step05 writes the review subset to:")
        log(f"  {CURATION_DIR / 'review_suggested.csv'}")
        log(f"")
        log(f"To adjust a flagged call: edit the 'final_annotation' column in")
        log(f"  {curated_out}  (needs_review==True rows), then re-run step05.")
        log(f"  AI_suggestion is pre-filled; AI_explanation gives the reasoning.")
        log(f"  Key evidence: phold_product/confidence/evalue, foldseek_description/")
        log(f"  score/evalue/pident, foldseek_top3.")
        log(f"  -> Run: python scripts/05_build_output.py")
    else:
        log("\nNo manual review needed!")
        log("-> Run: python scripts/05_build_output.py")

if __name__ == "__main__":
    main()
