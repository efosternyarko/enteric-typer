# Shigella Typing Tools

*Shigella* assemblies are run through five dedicated typing and screening steps in addition to MLST, AMRFinder, and PlasmidFinder. All results are aggregated into `shigella_typer_results.tsv`.

---

## 1. MLST

**Scheme**: `ecoli_achtman_4`  
*Shigella* species are phylogenetically *E. coli*, so the *E. coli* Achtman MLST scheme is used. This correctly assigns STs to *Shigella* isolates.

Notable STs:
- **ST152**: global *S. sonnei* pandemic lineage
- **ST148**, **ST145**: common *S. flexneri* lineages
- *S. dysenteriae* and *S. boydii* have no dominant single ST

---

## 2. ShigEiFinder — Species and serotype

**Tool**: [ShigEiFinder](https://github.com/LanLab/ShigEiFinder) v1.3.5  
**Output columns**: `shigeifinder_ipaH`, `shigeifinder_virulence_plasmid`, `shigeifinder_cluster`, `shigeifinder_serotype`, `shigeifinder_o_antigen`, `shigeifinder_h_antigen`

ShigEiFinder identifies *Shigella* species and serotype using a curated marker gene database. It covers all four *Shigella* species and their major serotypes, including all *S. flexneri* subtypes (1a, 1b, 2a, 2b, X, Y, and variants).

### Key output fields

| Column | Description |
|---|---|
| `shigeifinder_ipaH` | ipaH status: `pos` (present) or `neg` (absent) — critical virulence marker |
| `shigeifinder_virulence_plasmid` | Virulence plasmid classification (pINV presence/type) |
| `shigeifinder_cluster` | ShigEiFinder cluster ID |
| `shigeifinder_serotype` | Serotype call (e.g. `SS` = *S. sonnei*, `SF2a` = *S. flexneri* 2a) |
| `shigeifinder_o_antigen` | O-antigen formula |
| `shigeifinder_h_antigen` | H-antigen (non-motile in most *Shigella*) |

### Serotype prefix conventions

| Prefix | Species |
|---|---|
| `SS` | *S. sonnei* |
| `SF1a`, `SF2a`, etc. | *S. flexneri* (subtype follows) |
| `SB1`, `SB2`, etc. | *S. boydii* |
| `SD1`, `SD7`, etc. | *S. dysenteriae* |

### ipaH as a diagnostic marker

**ipaH is the gold standard molecular marker for *Shigella* and EIEC.** It encodes a leucine-rich repeat effector protein secreted by the type III secretion system and is found on the pINV virulence plasmid. All pathogenic *Shigella* carry ipaH; its absence strongly suggests the virulence plasmid has been lost (plasmid curing occurs during serial laboratory passage and rarely in clinical isolates).

---

## 3. Mykrobe — *S. sonnei* genotyping

**Tool**: [Mykrobe](https://github.com/Mykrobe-tools/mykrobe) with sonnei panel `20210201` (Hawkey et al. 2021)  
**Output columns**: `mykrobe_genotype`, `mykrobe_lineage`, `mykrobe_clade`, `mykrobe_subclade`, `mykrobe_genotype_name`, `mykrobe_confidence`

Mykrobe assigns *S. sonnei* isolates to the Hawkey 2021 global phylogenetic framework (lineages 1–4, with sub-clades). This framework was developed specifically for the global *S. sonnei* pandemic and correlates genotype with:
- Geographic origin
- Azithromycin resistance (lineage 3 clade III carries ICEShSon1 encoding *mph(A)* and *ermB*)
- Fluoroquinolone resistance (multiple independent acquisitions mapped to lineages)
- Ceftriaxone resistance emergence (lineage 3.6.1, XDR clade)

### Mykrobe output for non-sonnei isolates

Mykrobe is run on all *Shigella* assemblies because species is not always known in advance. For non-*S. sonnei* assemblies, Mykrobe will return a low-confidence or failed result; these are stored as NA and do not affect downstream figures.

### Mykrobe columns

| Column | Description |
|---|---|
| `mykrobe_lineage` | Top-level lineage (e.g. lineage3) |
| `mykrobe_clade` | Sub-lineage clade (e.g. lineage3.6) |
| `mykrobe_subclade` | Further resolution (e.g. lineage3.6.1) |
| `mykrobe_genotype` | Full genotype string |
| `mykrobe_genotype_name` | Human-readable name for this genotype |
| `mykrobe_confidence` | Confidence of the genotype call |

---

## 4. pINV screen — Invasion plasmid marker genes

**Method**: BLASTn against curated pINV marker gene database  
**Thresholds**: ≥ 80% identity, ≥ 80% query coverage  
**Output columns**: `pinv_present`, `pinv_genes`

The pINV screen detects six key pINV invasion genes individually:

| Gene | Protein | Role |
|---|---|---|
| icsA / virG | IcsA | Actin nucleation; intracellular motility and cell-to-cell spread |
| virF | VirF | Master virulence transcription factor (AraC family) |
| virB | VirB | Activates ipa operon; activated by VirF |
| ipaB | IpaB | TTSS effector; macrophage apoptosis; pore formation |
| ipaC | IpaC | TTSS effector; membrane insertion |
| ipaD | IpaD | TTSS needle tip; controls effector secretion |

`pinv_present`: `TRUE` if any gene detected at threshold; `FALSE` otherwise.  
`pinv_genes`: semicolon-separated list of detected genes with percentage identity and coverage.

A complete, functional pINV plasmid should carry all six genes. Partial gene sets suggest plasmid deletions.

---

## 5. IS element screen — Insertion sequence copy numbers

**Method**: BLASTn against IS element reference sequences  
**Thresholds**: ≥ 85% identity  
**Output column**: `is_elements` (format: `IS1(3);IS600(12);IS629(1)`)

Screens for six IS element families:

| Element | Copies in a typical *S. sonnei* pandemic isolate | Significance |
|---|---|---|
| IS1 | 5–15 | Ubiquitous; involved in gene disruption |
| IS1A | 0–5 | IS1 variant |
| IS30 | 0–10 | Gene disruption; present in multiple species |
| IS186 | 0–5 | Common |
| IS600 | 8–25 | High copy number marks *S. sonnei* pandemic clade ST152 |
| IS629 | 0–10 | O-antigen modification in *S. flexneri* |

Copy numbers are reported as counts of distinct BLAST hits above threshold on distinct contig positions. The `is_elements` column encodes copy numbers as `ISname(n)` tuples separated by semicolons.

IS copy numbers are visualised in Fig 1 Panel B (IS element landscape heatmap) and Fig 9 (Shigella features panel).

---

## 6. AMRFinder Plus

**Organism flag**: `--organism Escherichia` (correct for *Shigella* — same species complex)  
Applies the same point mutation database as for *E. coli*.

Key resistance mechanisms in *Shigella*:
- **Azithromycin (macrolide)**: *mph(A)* (macrolide phosphotransferase) and *erm(B)* (ribosomal methylase) — increasingly common in *S. sonnei* lineage 3, often on ICE elements
- **Ciprofloxacin (quinolone)**: *gyrA* mutations (S83L, D87N) + *parC* mutations — near-universal in global pandemic *S. sonnei*
- **Ceftriaxone (beta-lactam)**: *blaCTX-M-27* and *blaCTX-M-15* — emerging in XDR clades
- **Trimethoprim + sulfonamide**: *dfrA* + *sul* genes — historically very common; often chromosomally integrated via class 1 integrons

---

## 7. PlasmidFinder

**Database**: `enterobacteriaceae`  
**Output column**: `plasmidfinder_replicons`

*Shigella* plasmid content beyond pINV:
- **Col156**: cryptic small plasmid ubiquitous in *S. sonnei* ST152 pandemic clade (present in ~100% of pandemic isolates; lineage marker)
- **IncFII**: carries resistance genes in some lineages
- **IncI1**: broad host range; occasional ESBL carrier

pINV itself is a large (~214 kb) plasmid but its replicon is not reliably detected by PlasmidFinder — use the pINV screen (column `pinv_present`) to confirm its presence.
