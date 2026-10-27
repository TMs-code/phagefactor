"""Shared helpers for the mock tools (tests only).

The mocks write files in the SAME layout and formats the real tools write, with
deterministic synthetic content, so the workflow and
all python steps can be exercised end-to-end without databases or GPUs.
They prove plumbing and parity, NOT annotation quality.
"""
import hashlib
import random

AA = "ACDEFGHIKLMNPQRSTVWY"
DI = "ACDEFGHIKLMNPQRSTVWY"

PRODUCTS = [  # phold products, PHROG function
    ("hypothetical protein", "unknown function", "No_PHROG", "low"),
    ("hypothetical protein", "unknown function", "No_PHROG", "low"),
    ("hypothetical protein", "unknown function", "No_PHROG", "low"),
    ("terminase large subunit", "head and packaging", "phrog_2", "high"),
    ("major capsid protein", "head and packaging", "phrog_3", "high"),
    ("tail fiber protein", "tail", "phrog_4", "medium"),
    ("integrase", "integration and excision", "phrog_1", "high"),
    ("endolysin", "lysis", "phrog_5", "medium"),
    ("Cro/C1-type HTH DNA-binding domain", "transcription regulation", "phrog_6", "medium"),
]

FS_DESCS = [  # FoldSeek target descriptions: informative, low-information, titles
    "Site-specific integrase",
    "Major capsid protein",
    "Phage tail tube protein",
    "SHOCT domain-containing protein",
    "XkdX family protein",
    "YodL-like",
    "Uncharacterized protein",
    "Crystal structure of a phage protein in complex with DNA",
    "Holin",
    "Type II toxin-antitoxin system HicA family toxin",
    "Lipoprotein",
    "DNA-binding protein",
    "Tetratricopeptide repeat protein",
    "Terminase small subunit",
]


def rng(*keys):
    h = hashlib.md5("|".join(map(str, keys)).encode()).hexdigest()
    return random.Random(int(h[:12], 16))


def fake_3di(seq, key):
    r = rng("3di", key)
    return "".join(r.choice(DI) for _ in seq)


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


PHOLD_COLS = ["cds_id", "phrog", "function", "product", "annotation_method",
              "annotation_confidence", "bitscore", "fident", "evalue", "qStart", "qEnd",
              "qLen", "qCov", "tStart", "tEnd", "tLen", "tCov", "annotation_source",
              "tophit_protein", "function_with_highest_bitscore_proportion",
              "top_bitscore_proportion_not_unknown", "prostt5_confidence"]


# genome-mode `phold run` column order (37 columns; 01c renames by POSITION)
GENOME_COLS = [
    "contig_id", "cds_id", "start", "end", "strand", "phrog", "function", "product",
    "annotation_method", "annotation_confidence", "transl_table", "bitscore", "fident",
    "evalue", "qStart", "qEnd", "qLen", "qCov", "tStart", "tEnd", "tLen", "tCov",
    "annotation_source", "tophit_protein", "function_with_highest_bitscore_proportion",
    "top_bitscore_proportion_not_unknown", "head_and_packaging_bitscore_proportion",
    "integration_and_excision_bitscore_proportion", "tail_bitscore_proportion",
    "moron_auxiliary_metabolic_gene_and_host_takeover_bitscore_proportion",
    "DNA_RNA_and_nucleotide_metabolism_bitscore_proportion", "connector_bitscore_proportion",
    "transcription_regulation_bitscore_proportion", "lysis_bitscore_proportion",
    "other_bitscore_proportion", "unknown_function_bitscore_proportion", "prostt5_confidence",
]


def phold_row(cds_id, seq, contig=None):
    r = rng("phold", cds_id)
    prod, func, phrog, conf = r.choice(PRODUCTS)
    hyp = prod == "hypothetical protein"
    row = {
        "cds_id": cds_id, "phrog": phrog, "function": func, "product": prod,
        "annotation_method": "" if hyp else "foldseek",
        "annotation_confidence": conf,
        "bitscore": "" if hyp else f"{r.uniform(50, 900):.1f}",
        "fident": "" if hyp else f"{r.uniform(0.15, 0.9):.3f}",
        "evalue": "" if hyp else f"{10 ** -r.uniform(3, 40):.3e}",
        "qStart": "", "qEnd": "", "qLen": len(seq), "qCov": "", "tStart": "", "tEnd": "",
        "tLen": "", "tCov": "", "annotation_source": "" if hyp else "phrog",
        "tophit_protein": "" if hyp else f"tp_{r.randint(1, 10**6)}",
        "function_with_highest_bitscore_proportion": func,
        "top_bitscore_proportion_not_unknown": "",
        "prostt5_confidence": f"{r.uniform(40, 80):.2f}",
    }
    if contig is not None:
        row["contig_id"] = contig
    return row


def write_tsv(path, rows, cols):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
