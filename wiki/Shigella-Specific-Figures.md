# Shigella-Specific Figures

In addition to the standard Figs 1–7, *Shigella* runs produce two additional figures (Figs 8 and 9) that capture aspects of *Shigella* biology not addressed by the cross-species framework: serotype composition and the invasion/IS element feature landscape.

---

## Background: why Shigella needs dedicated figures

*Shigella* species (*S. sonnei*, *S. flexneri*, *S. boydii*, *S. dysenteriae*) are phylogenetically embedded within *E. coli*, but their clinical presentation, typing scheme, and genomic features differ substantially:

- Typing is by **serotype** (O-antigen serotype for *S. flexneri*; biotype for *S. sonnei*) rather than by sequence type alone
- The **virulence plasmid (pINV)** is the primary determinant of pathogenicity — its presence and invasion gene complement are clinically important
- **IS element content** differs markedly between species and can distinguish epidemic lineages (e.g. IS*600* expansion in *S. sonnei* pandemic clade)
- MDR in *S. sonnei* is increasingly driven by integrative conjugative elements (ICEs) rather than plasmids

---

## Fig 8 — Shigella Serotype Composition

A horizontal stacked bar chart showing the number of isolates per *Shigella* species (one bar per species), stacked by serotype (ShigEiFinder).

### Species identification

Species is inferred in this priority order:

1. ShigEiFinder serotype prefix: `SF` = *S. flexneri*, `SS`/SONNEI = *S. sonnei*, `SB` = *S. boydii*, `SD` = *S. dysenteriae*
2. ShigEiFinder cluster / serotype field (keyword matching)
3. Mykrobe lineage: lineage3 = *S. sonnei* pandemic

### Serotype colours

Each serotype is assigned a colour from the matplotlib `tab20` palette. The total n is shown to the right of each bar.

### Interpreting Fig 8

- **Dominant *S. sonnei*** with one serotype bar: typical of a single-strain or clonal epidemic; common in many LMICs where the pandemic lineage ST152 dominates.
- **Diverse *S. flexneri* serotypes**: multiple serotypes within *S. flexneri* indicate diverse circulating strains; can complicate vaccine efficacy assessment since candidate vaccines (e.g. Shigella flexneri 2a-based) may not cover all circulating serotypes.
- **Presence of *S. dysenteriae* type 1**: historically the cause of epidemic dysentery; carriage of Stx genes; any Sd1 isolates warrant immediate attention.

---

## Fig 9 — Shigella Virulence and Invasion Feature Panel

A binary presence/absence heatmap with samples as rows and virulence/invasion features as columns. Samples are sorted by species, then serotype, then sample name.

A dark red species colour strip to the left of the heatmap marks which species each sample belongs to.

### Column groups

The columns are divided into three bracketed groups:

#### 1. ShigEiFinder markers

| Column | Description |
|---|---|
| ipaH | Invasion plasmid antigen H — the gold standard PCR target for Shigella/EIEC; present on pINV |
| Vir. plasmid | ShigEiFinder's classification of whether a virulence plasmid is detected |

**ipaH** is the single most important diagnostic marker: it is present on the virulence plasmid of essentially all pathogenic *Shigella* and is absent from commensal *E. coli*. An isolate typed as *Shigella* but lacking ipaH may have lost its virulence plasmid (plasmid-cured strains exist) or be a misclassification.

#### 2. pINV invasion genes

| Gene | Protein | Function |
|---|---|---|
| icsA / virG | IcsA/VirG | Surface-exposed autotransporter; enables actin-based motility |
| virF | VirF | AraC-family transcription factor; master virulence regulator |
| virB | VirB | Transcriptional activator; activated by VirF |
| ipaB | IpaB | Type III secretion effector; pore-forming; triggers macrophage apoptosis |
| ipaC | IpaC | Type III secretion effector; membrane insertion |
| ipaD | IpaD | Needle tip protein; controls effector secretion |

All six genes should be present in a complete, functional pINV plasmid. Isolates missing one or more genes may have undergone partial plasmid deletions, which is associated with attenuated virulence.

#### 3. IS element presence

Binary presence of individual IS element families:

| Element | Significance |
|---|---|
| IS1 | Ubiquitous; present in multiple copies across Shigella species |
| IS1A | IS1 variant |
| IS30 | Common in Shigella; associated with gene disruption events |
| IS186 | Detected in S. sonnei and S. flexneri |
| IS600 | High copy number in S. sonnei pandemic clade; lineage marker |
| IS629 | Associated with O-antigen modification in S. flexneri |

A simple binary (present/absent) is shown here; for copy-number detail, refer to Fig 1 Panel B (IS element landscape heatmap) which shows copy numbers.

### Interpreting Fig 9

**All pINV genes present + ipaH present**: fully invasive strain with intact virulence plasmid.

**ipaH present but icsA/virG absent**: pINV is present but icsA has been disrupted — common in clinical passaged or attenuated derivatives.

**ipaH absent in an otherwise *Shigella*-typed isolate**: virulence plasmid may have been cured during laboratory passage. Treat with caution; may not be pathogenic.

**IS600 in high copy number (Fig 1B) + SS serotype**: hallmark of the global *S. sonnei* pandemic lineage (ST152). High IS600 copy number is not associated with specific virulence or resistance per se, but is a strong lineage marker.

---

## Output files

- `{prefix}_fig8_shigella_serotypes.pdf` and `.png`
- `{prefix}_fig9_shigella_features.pdf` and `.png`
