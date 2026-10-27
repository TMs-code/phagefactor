# Architecture

Two diagrams (Mermaid, rendered by GitHub; PNG copies in [`diagrams/`](diagrams/)). Prose version of the scoring and curation
logic: [`curation_gates.md`](curation_gates.md); where each word-list lives:
[`code_lists.md`](code_lists.md).

## 1. Pipeline

Three input modes converge on one shared spine. Step numbers are the script
prefixes in `scripts/`.

```mermaid
flowchart TD
    RUN["phagefactor run<br/>launch check: tools + databases of pending steps"] --> DET{"input mode<br/>(auto-detected)"}

    DET -->|genome| PK["00c Pharokka · per prophage<br/>Phanotate gene calls"]
    PK --> PH["01 Phold · per prophage<br/>ProstT5 3Di + phage DBs"]
    PH --> M1["01c merge"]

    DET -->|pangenome| TR["00n translate<br/>nt CDS → protein, clean ids"]
    TR --> SP
    DET -->|protein| SP["00p split<br/>50 proteins / batch"]
    SP --> PP["01p Phold proteins · per batch"]
    PP --> M2["01p merge"]

    M1 --> TG["02 targets<br/>hypothetical genes"]
    M2 --> TG
    TG --> SM{"search_mode"}
    SM -->|local| LS["02d FoldSeek search<br/>pdb100 · afdb-swissprot · afdb50"]
    LS --> LP["02d parse<br/>hit domain · scoring → best hit + top-3"]
    SM -->|webapi| WS["02w ESMFold → search.foldseek.com<br/>≤ 400 aa · cached · same scoring"]
    SM -->|none| NS["no FoldSeek evidence"]

    LP --> CO["03 compare<br/>Phold ⟷ FoldSeek · sub-DBs · agreement"]
    WS --> CO
    NS --> CO
    CO --> CU["04 curate + 05 output<br/>final table · GenBank · review list"]
    CU --> PY["06 Phynteny + 07 synteny<br/>(needs genome coordinates)"]
    CU --> PG["08 pangenome reference<br/>(pangenome mode)"]
    CU --> RE["09 report.html"]
    PY --> RE
```

Local and web-API search produce identically shaped `best_hit.csv` / `top3.csv`, so
everything downstream is the same whichever ran. The web API takes structures only,
hence the ESMFold step and the 400 aa limit.

## 2. Scoring and curation

How a set of FoldSeek hits becomes one call per gene.

```mermaid
flowchart TD
    subgraph S02["02 · score hits"]
        H["FoldSeek hits"] --> Q["quality + informative gate"]
        Q --> K["hit domain<br/>lineage · taxid · name"]
        K --> B["rank: score × phage boost"]
        B --> G["boost and eukaryote gates"]
        G --> T["best hit + top-3<br/>name upgrade"]
        T --> CF["confidence tier"]
    end
    subgraph S03["03 · compare"]
        SD["sub-DB substitution"] --> EU["eukaryote demote"]
        EU --> AG["agreement label"]
    end
    subgraph S04["04 · curate"]
        R1["sub-DB gate"] --> R2["agreement cases"]
        R2 --> R3["rule tree<br/>first match wins"]
        R3 --> R4["post-checks<br/>domain · head vs tail"]
        R4 --> OV["overrides.tsv"]
    end
    CF --> SD
    AG --> R1
    OV --> OUT["final table · review list · report"]
```

Two points the diagram is meant to make obvious:

1. **The informativeness gate removes most hits**, before any scoring happens.
2. **Trusted calls are never re-annotated.** Only Pharokka-hypothetical genes enter
   curation; a sequence-homology call is passed through untouched.
