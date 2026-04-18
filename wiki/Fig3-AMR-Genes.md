# Fig 3 — AMR Genes

Fig 3 is a horizontal bar chart of the top 25 most prevalent acquired AMR genes and point mutations across all isolates of the species. It is the gene-level complement to the drug-class view in Fig 1 Panel C, and provides the most granular view of resistance determinants in your dataset.

---

## What is shown

Each bar represents one AMR gene or chromosomal point mutation. Bar length is the percentage of isolates that carry that gene. Bars are sorted from most to least prevalent (top to bottom).

The figure shows up to the top 25 genes by prevalence. If fewer than 25 distinct acquired genes are detected across the whole dataset, all are shown.

---

## Colour coding

| Colour | Meaning |
|---|---|
| Blue (#4e79a7) | Acquired gene (horizontally transferred resistance gene) |
| Light blue (#aec6cf) | Chromosomal point mutation |

Point mutations are distinguished from acquired genes by their name format: a point mutation is reported by AMRFinder as `geneName_X99Y` — a base gene name, underscore, then an amino acid substitution in the form (reference)(position)(alternate). For example, `gyrA_S83L` is a serine-to-leucine substitution at position 83 of GyrA, conferring fluoroquinolone resistance.

---

## Data source

Genes are taken from the `amrfinder_acquired_genes` column in the results TSV — the AMRrules-filtered list of acquired resistance determinants only. Intrinsic chromosomal genes (e.g. *blaEC*, efflux pump subunits) are excluded.

---

## Interpreting the figure

### Gene name conventions

AMRFinder uses the NCBI AMR gene nomenclature:

- `bla` prefix = beta-lactamase (e.g. *bla*TEM-1, *bla*CTX-M-15, *bla*OXA-48)
- `aac`, `aph`, `ant` = aminoglycoside-modifying enzymes
- `tet` = tetracycline resistance gene
- `mcr` = mobile colistin resistance
- `qnr` = quinolone resistance (plasmid-mediated)
- `sul` = sulfonamide resistance
- `dfrA` = trimethoprim resistance (dihydrofolate reductase)
- `cat`, `cml` = chloramphenicol/phenicol resistance
- `mph`, `erm` = macrolide resistance
- Point mutations: `gyrA_*`, `parC_*` = fluoroquinolone; `pmrA_*`, `mgrB_*` = colistin

### Common patterns

**High *bla*TEM-1 prevalence**: TEM-1 is the most common plasmid-mediated beta-lactamase globally. It confers ampicillin resistance but is typically not an ESBL (does not hydrolyse cephalosporins). Very common in *E. coli* and *Salmonella*.

**CTX-M genes**: Extended-spectrum beta-lactamases (ESBLs) that hydrolyse third-generation cephalosporins. CTX-M-15 is the dominant ESBL globally; CTX-M-27 is common in *E. coli* ST131. If a CTX-M gene is present, check Fig 4 for IncF plasmids — these are the most common CTX-M carriers.

**gyrA/parC mutations**: Quinolone resistance via chromosomal mutations. Often co-selected with other resistance — isolates with multiple quinolone mutations (e.g. *gyrA* S83L + D87N + *parC* S80I) typically have high-level fluoroquinolone resistance.

**mcr genes**: Plasmid-mediated colistin resistance. Even at low prevalence, warrants attention as colistin is a last-resort agent. Cross-reference Fig 4 for IncX4 or IncI2 plasmids (the main *mcr-1* and *mcr-2* carriers).

**aac(6')-Ib-cr**: Aminoglycoside-acetyltransferase variant that also inactivates fluoroquinolones. A gene detected under the AMG class that also contributes to QNL resistance — if you see this alongside *gyrA* mutations it suggests co-selection pressure.

---

## Relationship to other figures

| Question | Figure to check |
|---|---|
| Which drug classes are these genes in? | Fig 1 Panel C |
| Are these genes phylogenetically clustered? | Fig 2 AMR heatmap |
| Which plasmids carry these genes? | Fig 4 (long-read data) |
| How does gene carriage vary by ST? | Fig 6 |
| How does gene carriage vary by serovar/phylogroup? | Fig 7 |

---

## Output file

`{prefix}_fig3_amr_genes.pdf` and `.png`
