#!/usr/bin/env python3
"""
naming_rules.py — derive a concise `short_name` from `final_product`
====================================================================
Used by 05_build_output.py to add a `short_name` column placed BEFORE
`final_product`. It does NOT touch the curation decision tree (03/04) — it is a
pure, deterministic post-pass on the already-chosen final_product string, so it
is easy to audit and tweak and cannot change which annotation was selected.

Rule order (first match wins), distilled from per-gene review:
  1. Sub-DB "Name (Category): description"  -> Name        (RecN, Cya, Ibes,
       "Hemolysin HlyA", LapB)   [VFDB/CARD/DefenseFinder formatting]
  2. Slash alternatives "A / B / C"         -> clean(A)     (primary call;
       "transaldolase / Dihydrodipicolinate..." -> transaldolase;
       "baseplate wedge subunit / Phage protein GP46" -> baseplate wedge subunit)
  3. Decoration stripping (always): drop (Fragment), (Modular protein), (EC ...),
       ", Tyr-sensitive", "from bacteriophage X", ", catalytic domain",
       "Lambda family ", ", lambda family", "HK97 gpNN family ", "domain-
       containing protein"->core, leading "Phage/Prophage/Putative".
  4. Trailing gene symbol "<long descriptor> XxxN" -> symbol (DinI, TagO, BamE,
       YfiB, DapA) when >=3 descriptor words precede it.
  5. Head enzyme noun: collapse long enzyme names to their class word
       (aldolase, dehydrogenase, oxidoreductase, hydrolase, ...).
  6. Otherwise: the cleaned string, trimmed to <=5 words.

These are heuristics — a first pass. Edge cases the user wants worded a specific
way (e.g. "Adhesion protein LapB", "T4SSB protein IcmN") may need a small manual
override map later; rule 1/4 give the gene symbol which is a sensible default.
"""

import re

# CamelCase-style gene symbols: DinI, TagO, BamE, YfiB, DapA, HlyA, RecN, Zur,
# IcmN. Require an initial capital, some lowercase, then a capital or digit, so
# plain English words (e.g. "Repressor") don't match.
_GENE_SYMBOL = re.compile(r'^[A-Z][a-z]{1,4}[A-Z0-9][A-Za-z0-9]*$')

# enzyme/function class "head nouns" — collapse long names to these.
_HEAD_NOUNS = [
    "transaldolase", "aldolase", "dehydrogenase", "oxidoreductase", "reductase",
    "dehydratase", "racemase", "epimerase", "isomerase", "mutase", "synthetase",
    "synthase", "transferase", "kinase", "phosphatase", "hydrolase", "nuclease",
    "ligase", "permease", "peptidase", "protease", "phosphodiesterase",
    "cyclase", "deaminase", "hydroxylase", "decarboxylase", "carboxylase",
]
_HEAD_NOUN_RE = re.compile(r'\b(' + "|".join(_HEAD_NOUNS) + r')\b', re.I)

_DECORATION = [
    re.compile(r'\s*\((fragment|modular protein|partial)\)', re.I),
    re.compile(r'\s*\(EC[ :][\d.\-]+\)', re.I),
    re.compile(r'\s*,\s*(tyr|phe|trp|his)-?sensitive', re.I),
    re.compile(r'\s+from bacteriophage\s+\w+', re.I),
    re.compile(r'\s*,?\s*catalytic domain', re.I),
    re.compile(r'\s*,\s*contains\b.*$', re.I),
    re.compile(r'\s*,?\s*lambda family', re.I),
    re.compile(r'\bHK97\s+gp\d+\s+family\s+', re.I),
    re.compile(r'\s*\(ACLAME[^)]*\)', re.I),
    re.compile(r'\s*\(Modular protein\)', re.I),
    re.compile(r'\s+\(IPT/TIG domain\)', re.I),
]
_LEADING = re.compile(r'^(putative|prophage|phage|bacteriophage|probable)\s+', re.I)
_TRAIL_NUM = re.compile(r'\s+\d+$')


def _strip_decoration(s: str) -> str:
    for pat in _DECORATION:
        s = pat.sub("", s)
    s = s.replace("domain-containing protein", "protein")
    s = _TRAIL_NUM.sub("", s).strip(" /,;")
    # collapse repeated leading Phage/Putative (e.g. "Phage holin" -> "holin")
    prev = None
    while prev != s:
        prev = s
        s = _LEADING.sub("", s).strip()
    return s.strip() or prev


def make_short_name(final_product, phold_phrog=None) -> str:
    """Return a concise short_name for a final_product string."""
    if final_product is None:
        return ""
    desc = str(final_product).strip()
    if not desc or desc.lower() in ("hypothetical protein", "nan", "na", ""):
        return "hypothetical protein" if "hypoth" in desc.lower() else ""

    # Rule 1: sub-DB "Name (Category): description". Use the Name as short_name
    # ONLY when it is a real gene symbol (RecN, Cya, HlyA, LapB, Ibes) — not a
    # phenotype/system label like "Lateral flagella" or "Biofilm". Otherwise fall
    # through and shorten the functional description instead.
    m = re.match(r'^(.{2,40}?)\s*\([^)]+\):\s*(.+)$', desc)
    if m:
        name, descpart = m.group(1).strip(), m.group(2).strip()
        name_toks = name.split()
        looks_like_symbol = (
            any(_GENE_SYMBOL.match(t) or t.isupper() for t in name_toks)
            or (len(name_toks) == 1 and len(name) <= 5)   # Cya, Ibes, Lap
        )
        if looks_like_symbol:
            return name
        desc = descpart   # phenotype label -> shorten the description instead

    # Rule 2: slash alternatives -> take the first (the primary/chosen call)
    if " / " in desc:
        desc = desc.split(" / ")[0].strip()

    # Rule 3: strip decoration
    desc = _strip_decoration(desc)

    # Rule 4: trailing gene symbol after a long descriptor -> symbol
    toks = desc.split()
    if len(toks) >= 4 and _GENE_SYMBOL.match(toks[-1]):
        return toks[-1]

    # Rule 5: collapse long enzyme names to their head noun. Triggers when the
    # name is wordy (>3 tokens) OR carries a long chemical-substrate prefix
    # (any token >12 chars, e.g. "Phospho-2-dehydro-3-deoxyheptonate aldolase").
    if len(toks) > 3 or any(len(t) > 12 for t in toks):
        hits = _HEAD_NOUN_RE.findall(desc)
        if hits:
            return hits[-1].lower()

    # Rule 6: trim to <=5 words
    out = " ".join(toks[:5]) if len(toks) > 5 else desc
    # the word-trim can cut INSIDE a parenthetical ("tail fiber protein (Large
    # polyvalent") -> drop the dangling "(...": never leave an unbalanced bracket.
    if out.count("(") > out.count(")"):
        out = out[:out.rfind("(")].strip()
    out = out.strip(" -/,;")
    # never return a useless bare "protein" (e.g. from stripping "Prophage protein")
    if out.strip().lower() in ("", "protein"):
        return "phage protein"
    return out
