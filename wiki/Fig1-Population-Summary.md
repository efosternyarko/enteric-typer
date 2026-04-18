# Fig 1 — Population Summary

Fig 1 is a four-panel overview of the population structure and resistance landscape for a single species run. It is the first figure produced and is designed to give a rapid, high-level picture of who is in your dataset before diving into the gene-level and phylogenetic detail.

---

## Layout

```
┌─────────────────────┬─────────────────────┐
│  A. MLST STs        │  B. Serotypes / IS  │
├─────────────────────┼─────────────────────┤
│  C. AMR prevalence  │  D. MDR burden      │
└─────────────────────┴─────────────────────┘
```

---

## Panel A — MLST sequence types

A horizontal bar chart of the top 15 most common MLST sequence types (STs) by isolate count. STs beyond the top 15 are pooled into "Other"; isolates with no ST prediction are labelled "Unknown".

### How bars are filled (species-specific)

**E. coli**: Each bar is stacked by Clermont phylogroup. If `clermont_phylogroup` is present in the results TSV, those assignments are used directly. If not, phylogroup is inferred from a curated ST → phylogroup lookup table. The title notes which approach was used ("(Clermont)" vs "(ST-inferred)").

| Phylogroup colour | Phylogroup |
|---|---|
| Blue | A |
| Green | B1 |
| Orange | B2 |
| Teal | C |
| Red | D |
| Purple | E |
| Pink | F |
| Grey | Unknown |

**Salmonella**: Bars are filled with a single-hue gradient (darker = more common ST). No stacking — the equivalent stratification is by serovar, shown in Panel B.

**Shigella**: Bars are stacked by *Shigella* species (*S. sonnei*, *S. flexneri*, *S. boydii*, *S. dysenteriae*).

### Interpreting Panel A

- **Dominant ST**: A single very tall bar indicates a clonal outbreak or collection enriched for one lineage. Multiple STs of roughly equal height suggests a diverse, community-level sample.
- **Phylogroup composition** (E. coli): B2 dominates extra-intestinal pathogenic *E. coli* (ExPEC); B1 is more common in commensal and livestock isolates; A is ubiquitous across niches.
- **"Other" is large**: Your dataset contains many rare STs. This is common in community surveillance; not a quality issue.
- **"Unknown" is large**: MLST failed for many samples. Common causes: highly fragmented assemblies (< 4 MLST alleles found), novel alleles, or contaminated assemblies.

---

## Panel B — Serotypes / IS element landscape

The content of Panel B differs by species.

### E. coli: O:H serotypes (ECTyper)

Top 15 O:H serotype combinations (e.g. O25:H4), stacked by K-locus group (G1–G4) when Kaptive K-locus data are available in the results TSV.

| K-group colour | Group |
|---|---|
| Blue | G1 |
| Orange | G2 |
| Green | G3 |
| Red | G4 |
| Purple | G1/G4 overlap |
| Grey | Unknown |

**K-locus groups** are broad structural groupings of capsule loci relevant to virulence and vaccine design. G1 and G4 groups are associated with ExPEC and neonatal meningitis isolates respectively.

### Salmonella: Serovars (SISTR)

Top 15 serovars by isolate count, stacked by MLST ST when available. This lets you see whether a given serovar (e.g. Typhimurium) is composed of a single dominant ST or multiple lineages.

### Shigella: IS element landscape (copy number heatmap)

A compact matrix with samples as rows and IS element families as columns. Cell values are copy numbers (number of IS element insertions per genome); darker = more copies. Samples are sorted by *Shigella* species and serotype.

Tracked IS element families include IS*1*, IS*2*, IS*4*, IS*66*, IS*100*, IS*200/605*, IS*256*, IS*630*, IS*3*, IS*21*, IS*5*, IS*110*, IS*1595*, IS*481*, and IS*91*.

IS elements are central to *Shigella* evolution: IS*600* and IS*1* contribute to pathoadaptation by disrupting genes (e.g. *cadA*, *ompT*) that are incompatible with intracellular lifestyle. IS*1* and IS*4* expansions can mark specific epidemic lineages.

---

## Panel C — AMR drug class prevalence

A stacked horizontal bar chart where each bar represents one clinical drug class and the bar length shows what percentage of isolates carry at least one acquired resistance gene in that class.

### Bar filling

Each bar is broken into segments, one per gene (up to the top 6 genes for that class). Segment width is proportional to the percentage of isolates carrying that gene. A light grey "Other" segment captures additional genes beyond the top 6.

**All genes shown are acquired genes only** (AMRrules-filtered). Intrinsic genes — such as chromosomal *blaEC* AmpC in *E. coli* — are excluded.

### Colour coding

Genes are assigned colours by global prevalence rank across all drug classes (most prevalent gene = vivid red, second = steel blue, etc., cycling through 12 distinct hues). The colour assignment is consistent within a run, so the same gene always gets the same colour across all panels that show individual genes.

### Interpreting Panel C

- **Bar length**: the percentage of isolates resistant to this drug class. A 50% beta-lactam bar means half your isolates carry an acquired beta-lactam resistance gene.
- **Segment widths**: which specific genes are driving resistance in each class. A broad single-gene segment (e.g. *bla*TEM-1 filling the entire beta-lactam bar) means one gene dominates. Multiple narrow segments mean gene diversity within the class.
- **Classes absent from the chart**: no isolate had a detected acquired gene in that class — not necessarily true susceptibility (phenotypic testing may still detect resistance via mechanisms not in the AMRFinder database).

---

## Panel D — MDR burden

A histogram where the x-axis is the number of drug classes with acquired resistance per isolate and the y-axis is the number of isolates. A dashed vertical line marks the mean number of resistant classes.

### Colour coding

- **Red bars** (x ≥ 3): MDR isolates
- **Blue bars** (x < 3): non-MDR isolates

The MDR threshold of ≥ 3 drug classes follows WHO/ECDC definitions. See [AMR Classification](AMR-Classification.md) for details.

### Interpreting Panel D

- A bimodal distribution (spike at 0, another at 4–6) indicates a population with a distinct MDR subgroup, often a single successful clone.
- A unimodal peak at 0–1 is typical for commensal collections or low-resistance settings.
- A long tail to the right (isolates with 8+ classes) suggests extensively drug-resistant (XDR) isolates are present.
- The mean is a single-number summary of the resistance burden of your collection as a whole.
