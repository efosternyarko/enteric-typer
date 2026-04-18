# Fig 4 — Plasmid Replicon Overview

The plasmid overview figure summarises plasmid replicon diversity across all samples of a given species. Because the meaningfulness of drug-class attribution depends critically on assembly quality, the figure **automatically selects one of two layouts** based on the data. This page explains both layouts, when each is chosen, and how to interpret what you see.

---

## Background: replicons and plasmid detection

[PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/) screens each assembly for known plasmid replicon sequences. A replicon is the origin-of-replication region that defines a plasmid's incompatibility group (Inc type). The Inc type determines:

- Which other plasmids the plasmid is incompatible with (cannot stably coexist in the same cell)
- Broadly, the host range (which bacterial species can carry it)
- Broadly, the transfer machinery available (conjugation frequency, host range of transfer)

Common replicons in enteric pathogens:

| Replicon family | Examples | Common associations |
|---|---|---|
| IncF | IncFII(S), IncFIB(S), IncFII(29) | Dominant in *E. coli* and *Salmonella*; often carry resistance |
| IncI | IncI1 | Broad host range; common in *Salmonella* |
| IncN | IncN | Broad host range; associated with ESBL genes |
| IncX | IncX1, IncX3 | Carry carbapenemase and ESBL genes |
| IncL/M | IncL/M(pOXA-48) | OXA-48 carbapenemase plasmids |
| Col | Col156, ColRNAI | Small cryptic plasmids; Col156 is specific to *S. sonnei* pandemic clone |

---

## Linking replicons to AMR: the contig co-location approach

The `plasmid_amr_map.py` script links replicon types to AMR genes by **contig co-location**: if a replicon sequence and an AMR gene are detected on the same contig, they are counted as co-occurring. This is the most rigorous approach because it confirms physical proximity on the same DNA molecule.

### Why this matters for short-read vs long-read assemblies

| Assembly type | Plasmid contiguity | Co-location reliability |
|---|---|---|
| Long-read (Nanopore / PacBio) | Plasmids typically assemble as single closed contigs | High — replicon and AMR genes on same contig = same plasmid |
| Hybrid (long + short reads) | Near-complete plasmid reconstruction | High |
| Short-read (Illumina) only | Plasmids often fragmented across 5–50 small contigs | Low — replicon and AMR genes routinely on different contigs of the same plasmid |

For short-read assemblies, the contig carrying the IncFII replicon sequence may be only a few kb long, while the AMR gene cassettes sit on separate small contigs from the same plasmid. Same-contig co-location misses this relationship entirely, even when the biology is well-established (e.g. MDR IncFII(S) plasmids in *Salmonella* Typhimurium).

---

## Auto-layout decision

After computing same-contig co-location data, the pipeline calculates:

```
fraction = (isolates with ≥ 1 replicon-carrying contig that also has AMR)
           ÷ (isolates with any replicon detected)
```

If this fraction is **< 10%**, the figure switches to the **simplified layout**. If ≥ 10%, the **full layout** is used.

### Decision tree

```
Compute same-contig AMR fraction
            │
            ├── < 10%  →  Simplified layout (2 panels)
            │               Likely: short-read Illumina assembly
            │               Action: show replicon prevalence by Inc family
            │
            └── ≥ 10%  →  Full layout (3 panels)
                            Likely: long-read or reference-quality assembly
                            Action: show drug-class breakdown + bubble matrix
```

---

## Simplified layout (2 panels) — short-read assemblies

Used when same-contig AMR co-occurrence is too sparse to be meaningful (< 10% threshold).

### Panel A — Replicon prevalence by Inc/Col family

A horizontal bar chart showing the percentage of isolates that carry each replicon, with bars coloured by broad Inc/Col family:

| Colour | Family | Examples |
|---|---|---|
| Blue | IncF | IncFII, IncFIB, IncFII(S), IncFIB(S) |
| Orange | IncI | IncI1, IncI2 |
| Red | IncN | IncN, IncN2 |
| Teal | IncA/C | IncA/C, IncA/C2 |
| Green | IncH | IncHI1, IncHI2 |
| Yellow | IncL/M | IncL/M, IncL/M(pOXA-48) |
| Purple | IncX | IncX1, IncX3, IncX4 |
| Pink | IncP | IncP-1α, IncP-1β |
| Brown | IncQ | IncQ1, IncQ2 |
| Grey | Col | Col156, ColRNAI, Col(BS512) |
| Light grey | Other | Unclassified or rare types |

**What the bar length tells you:** prevalence (%). A replicon at 50% is carried by half of your isolates. This does not mean the AMR genes on that replicon are at the same prevalence — the replicon prevalence is a property of the plasmid backbone, not of any specific resistance cargo.

**What the colour tells you:** the Inc family grouping, which correlates with plasmid biology (transfer range, incompatibility, typical gene cargo). It does *not* tell you about AMR in this layout — for that, cross-reference Fig 3 and Fig 6.

### Panel B — SNP phylogeny + plasmid heatmap

The midpoint-rooted SNP tree (from IQ-TREE 2) with:
- **ST strip** — MLST sequence type colour per tip
- **PG strip** — Clermont phylogroup (for *E. coli*) or left blank for Salmonella/Shigella
- **Plasmid heatmap** — one column per replicon; blue = present, white = absent

This panel lets you see whether plasmid carriage is structured by phylogeny. If replicon X is confined to a single clade, the plasmid was likely acquired once and co-transmitted with the bacterial lineage. If it is scattered across unrelated clades, it has been acquired multiple times (horizontal gene transfer).

---

## Full layout (3 panels) — long-read or reference-quality assemblies

Used when same-contig AMR co-occurrence is ≥ 10% — meaning the assemblies are contiguous enough that co-location data is biologically meaningful.

### Panel A — Replicon prevalence by drug class

As for the simplified layout but now each bar is stacked by the AMR drug class found on that replicon type. Each isolate carrying a given replicon is assigned to exactly one drug class segment (the highest-priority class found on the replicon-carrying contig), so segments sum to total replicon prevalence without double-counting.

**Priority order** (highest to lowest): BETA-LACTAM → QUINOLONE → COLISTIN → AMINOGLYCOSIDE → TETRACYCLINE → MACROLIDE → PHENICOL → SULFONAMIDE → TRIMETHOPRIM → FOSFOMYCIN → STREPTOTHRICIN → NITROFURAN → No AMR

Grey ("No AMR") means the replicon was detected but no acquired AMR gene was on the same contig.

### Panel B — Bubble matrix of replicon–AMR co-occurrence

A grid with replicons on the y-axis and drug classes on the x-axis. Each bubble represents a (replicon, drug class) combination:
- **Bubble size** — proportional to the square root of the percentage of isolates carrying both the replicon and that drug class (area ∝ percentage)
- **Bubble colour** — drug class
- **Number inside** — percentage of isolates (shown when ≥ 5%)
- **No bubble** — no isolate has both that replicon and that drug class on the same contig

### Panel C — SNP phylogeny + plasmid heatmap

Same as Panel B in the simplified layout.

---

## Common interpretation questions

**"Why does a well-known MDR plasmid type (e.g. IncFII) show No AMR in the bars?"**

Almost certainly an assembly fragmentation issue. Switch to the simplified layout interpretation: the replicon prevalence bar tells you the plasmid backbone is there, but the AMR cargo assessment requires long-read data. See Fig 3 for AMR gene prevalence independently of plasmid attribution.

**"What does it mean when a replicon is at 100% prevalence?"**

The replicon is present in every isolate. This can mean:
- A highly conserved plasmid that co-evolved with this lineage (common for Col-type replicons)
- A large virulence plasmid carried by all isolates (e.g. pINV in *Shigella*)
- A chromosome-integrated plasmid sequence (rare but possible)

**"Why is Col156 at 100% in *S. sonnei*?"**

Col156 is the replicon of a small cryptic plasmid that is ubiquitous in the global pandemic *S. sonnei* ST152 lineage. Its presence is a lineage marker as much as a plasmid marker.

**"IncL/M(pOXA-48) appeared in my Salmonella dataset — should I be alarmed?"**

Yes, this warrants attention. IncL/M(pOXA-48) is the plasmid backbone associated with OXA-48-type carbapenemase genes. However, confirm by checking the AMRFinder results for `blaOXA-48` or related genes — the replicon alone does not prove carbapenemase carriage.
