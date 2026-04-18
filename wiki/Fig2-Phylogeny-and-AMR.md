# Fig 2 — Phylogeny and AMR

Fig 2 is the annotated phylogenetic tree: it places every isolate in its evolutionary context and overlays virulence and AMR gene presence/absence as coloured heatmap panels. It is the figure that most directly answers the question "which lineages carry resistance and virulence genes, and do they cluster together?"

---

## Layout (left to right)

```
Phylogenetic tree  |  ST strip  |  PG strip*  |  Virulence heatmap  |  AMR genes heatmap
```

\* Phylogroup (PG) strip is present for *E. coli* only.

---

## The phylogenetic tree

The tree is drawn from the IQ-TREE 2 midpoint-rooted Newick file (see [SNP Phylogenetics](SNP-Phylogenetics.md) for how it is constructed). Branch lengths represent substitutions per site in the SKA2 SNP alignment.

- Branches are drawn in dark grey (#333333) with uniform line weight
- A scale bar at the top left gives the branch length in substitutions per site
- If ≤ 60 isolates are present, tip labels (sample names) are shown in monospace font
- If > 60 isolates, tip labels are omitted to avoid overlap; use the heatmap panels to identify samples

### Reading the tree

- **Branch length**: longer branches = more SNP divergence from the ancestor. Very long branches on one or two isolates often indicate contamination or mislabelling.
- **Topology**: closely clustered tips share recent common ancestry (same lineage or outbreak cluster). Dispersed identical-coloured ST strips indicate polyphyletic origins.
- **Scale bar**: helps interpret whether apparent clustering reflects true phylogenetic proximity. For context, within-outbreak isolates typically have branch lengths of 0–20 SNPs (see [SNP Phylogenetics](SNP-Phylogenetics.md)).

---

## ST strip

A single column of coloured rectangles, one per tip, showing the MLST sequence type. Each ST gets a unique colour drawn from a 12-colour qualitative palette; when > 12 STs are present the palette cycles. An ST legend is placed at the bottom of the figure.

**What to look for**: ST colour blocks that coincide with subtrees (clades) indicate clonal expansion. ST colour blocks scattered across unrelated clades indicate the same ST has been assigned to genomically distant isolates — uncommon but possible when the MLST scheme has low resolution for that lineage.

---

## Phylogroup strip (*E. coli* only)

A second thin column showing the Clermont phylogroup, using the canonical colour scheme:

| Colour | Phylogroup |
|---|---|
| Blue (#4e79a7) | A |
| Green (#59a14f) | B1 |
| Orange (#f28e2b) | B2 |
| Teal (#76b7b2) | C |
| Red (#e15759) | D |
| Purple (#b07aa1) | E |
| Pink (#ff9da7) | F |
| Grey (#bab0ac) | Unknown |

Phylogroup is taken from `clermont_phylogroup` in the results TSV (EzClermont output) when available; otherwise from a curated ST → phylogroup lookup.

**What to look for**: phylogroups should be largely consistent within major clades (because phylogroup and MLST reflect the same underlying phylogeny). Phylogroup inconsistencies within a clade can indicate recombination events or errors in phylogroup calling.

---

## Virulence heatmap

A binary heatmap (teal = present, light grey = absent) showing virulence determinants. The exact columns shown depend on the species.

### *E. coli* virulence columns

Sources used (in order of priority):

1. **Kleborate virulence determinants** — columns include: Stx1, Stx2, eae (intimin), ipaH, LT (heat-labile enterotoxin), ST-toxin (heat-stable enterotoxin). Only columns where at least one isolate is positive are shown.
2. **AMRFinder virulence genes** — up to 12 of the most prevalent virulence-associated genes reported by AMRFinder Plus under the `--plus` flag (e.g. *hlyA*, *cnf1*, *papA*, *sat*, *iutA*). Shown in addition to Kleborate columns.

### *Salmonella* virulence columns

AMRFinder virulence genes only (top 12 by prevalence across the dataset).

### *Shigella* virulence columns

1. ShigEiFinder outputs: ipaH status, virulence plasmid classification
2. pINV screen result
3. AMRFinder virulence genes (top 12)

### Column labels

Each column label shows the gene or marker name and, in parentheses, the percentage of isolates that are positive. Labels are in italic and rotated 40°.

---

## AMR genes heatmap

A binary heatmap (red = present, light grey = absent) showing the top 20 acquired AMR genes by prevalence across the dataset. Genes are ordered by drug class (highest-burden class leftmost), with drug class brackets and labels below the column names.

**Only acquired genes are shown** — intrinsic genes excluded by AMRrules (see [AMR Classification](AMR-Classification.md)).

### Reading the AMR heatmap

- A solid red vertical stripe means that gene is nearly universal in the dataset.
- A gene confined to a tight clade on the tree means it was acquired once and transmitted clonally; the same tree topology will be visible as a red block coinciding with a monophyletic group.
- A gene scattered across unrelated lineages means it has been acquired multiple times, suggesting a mobile element is circulating. Confirm by checking Fig 4 for a plasmid replicon co-occurring at similar prevalence.
- Comparing the virulence and AMR panels side by side shows whether virulence and resistance are co-selected in the same isolates.

### Drug class brackets

Below the gene tick labels, drug class names are shown with a bracket spanning all genes in that class. Classes are ordered by total prevalence (highest leftmost). This groups related resistance determinants (e.g. all beta-lactamases together) for easier visual parsing.

---

## When Fig 2 is not produced

Fig 2 requires a phylogenetic tree, which in turn requires ≥ 3 samples. If the species has fewer than 3 samples, or if `--skip_local_phylo` was passed, Fig 2 is skipped and a note is written to the pipeline log.

---

## Legend

The full figure has two legend areas:

- **Lower-left of tree panel**: phylogroup colours (E. coli) + virulence present/absent swatches + AMR present/absent swatches
- **Bottom centre of figure**: ST colour legend (multi-column, all observed STs)
