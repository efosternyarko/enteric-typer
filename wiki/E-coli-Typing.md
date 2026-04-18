# E. coli Typing Tools

*E. coli* assemblies are run through six typing tools that together characterise sequence type, serotype, capsule type, phylogroup, pathotype, and virulence potential. All results are aggregated into the `ecoli_typer_results.tsv` output file.

---

## 1. MLST — Multi-locus sequence typing

**Tool**: [mlst](https://github.com/tseemann/mlst) (Torsten Seemann)  
**Scheme**: `ecoli_achtman_4` — the Warwick Achtman 7-locus scheme (PubMLST)  
**Output columns**: `mlst_scheme`, `mlst_st`, `mlst_st_complex`

The Achtman scheme uses seven housekeeping gene loci (*adk*, *fumC*, *gyrB*, *icd*, *mdh*, *purA*, *recA*). Each unique combination of alleles defines a sequence type (ST). STs are assigned a sequence type complex (ST complex) where relevant, based on a curated lookup table.

### Interpreting MLST results

- **Numbered ST** (e.g. ST131, ST73): confident assignment; all 7 alleles found
- **ST with tilde** (e.g. ~ST131): closest match but not an exact allele match — treat as provisional
- **NA**: fewer than 7 loci found; assembly may be fragmented or contaminated
- **Novel**: allele combination not yet in the PubMLST database

Common clinically important *E. coli* STs:

| ST | Phylogroup | Notable features |
|---|---|---|
| ST131 | B2 | Global ESBL/fluoroquinolone-resistant ExPEC pandemic clone |
| ST73 | B2 | Common ExPEC; often carries *iutA*, *hlyA* |
| ST95 | B2 | ExPEC; neonatal meningitis association |
| ST1193 | B2 | Emerging fluoroquinolone-resistant clade |
| ST10 | A | Broad environmental niche; diverse plasmids |
| ST58 | B1 | Livestock-associated; common ESBL carrier in food animals |
| ST410 | C | Carbapenemase-associated (NDM, OXA-48) |
| ST167 | A | Associated with ESBL and MCR |

---

## 2. ECTyper — O:H serotyping

**Tool**: [ECTyper](https://github.com/phac-nml/ecoli_serotyping) v1.0  
**Output columns**: `ectyper_O`, `ectyper_H`, `ectyper_serotype`, `ectyper_qc`, `ectyper_evidence`

ECTyper predicts O (somatic lipopolysaccharide) and H (flagellar) antigen types by BLASTing the assembly against a curated reference database of O-antigen biosynthesis genes (*wzx/wzy*, *wbgA-Z*) and H-antigen genes (*fliC*). Thresholds used: 90% identity, 50% coverage (for both O and H).

### Output fields

| Column | Description |
|---|---|
| `ectyper_O` | O-antigen type (e.g. O25, O6) |
| `ectyper_H` | H-antigen type (e.g. H4, H31) |
| `ectyper_serotype` | Combined O:H call (e.g. O25:H4) |
| `ectyper_qc` | QC flag: PASS / WARN / FAIL |
| `ectyper_evidence` | Number of marker genes found supporting the call |

### Common serotypes and their significance

| Serotype | Notes |
|---|---|
| O25:H4 | ST131 — dominant ESBL/ExPEC clone globally |
| O6:H1 | ST69 — cystitis-associated |
| O1:H7 / O2:H6 | ST95, ST73 — ExPEC |
| O157:H7 | EHEC/STEC — Shiga toxin-producing; notifiable |
| O104:H4 | STEC — caused 2011 Germany outbreak |

---

## 3. EzClermont — Clermont phylogroup

**Tool**: [EzClermont](https://github.com/nickp60/EzClermont)  
**Output column**: `clermont_phylogroup`

EzClermont applies the Clermont 2013 PCR-based phylotyping scheme *in silico* to assembled genomes. It determines which Clermont phylogroup (A, B1, B2, C, D, E, F, or cryptic clades) an isolate belongs to by checking for the presence or absence of PCR amplicon sequences: *chuA*, *yjaA*, *TspE4.C2*, and *arpA*.

Minimum contig length is set to 200 bp (rather than the default 500 bp) to retain diagnostic loci in short-contig assemblies.

### Phylogroup colour coding

Phylogroup colours used throughout the pipeline figures:

| Colour | Phylogroup | Typical niche |
|---|---|---|
| Blue | A | Commensal; diverse environments |
| Green | B1 | Commensal; livestock; water |
| Orange | B2 | ExPEC; UTI, bacteraemia |
| Teal | C | Clinical; less common |
| Red | D | Clinical; fluoroquinolone-resistant clones |
| Purple | E | O157:H7 and relatives |
| Pink | F | Rare; includes some ExPEC |
| Grey | Unknown | Assignment failed |

---

## 4. Kleborate — Pathotyping + virulence markers

**Tool**: [Kleborate](https://github.com/klebgenomics/Kleborate) v3, preset `escherichia`  
**Output columns**: `kleborate_pathovar`, `kleborate_Stx1`, `kleborate_Stx2`, `kleborate_eae`, `kleborate_ipaH`, `kleborate_LT`, `kleborate_ST_toxin`

Kleborate v3 with the `escherichia` preset runs several sub-modules:

- **Pathotype detection**: classifies each isolate as STEC, EHEC, EPEC, ETEC, or EIEC based on virulence markers
- **Toxin typing**: detailed Shiga toxin subtype (Stx1a, Stx1c, Stx2a, Stx2b, etc. via StxTyper)
- **Virulence markers**: individual presence/absence calls for Stx1, Stx2, eae (intimin), ipaH, LT, ST-toxin

### Pathotype logic

| Pathotype | Required markers | Disease |
|---|---|---|
| EHEC | Stx + eae | Haemorrhagic colitis, HUS |
| STEC | Stx (any), no eae required | Bloody/watery diarrhoea |
| EPEC | eae, no Stx | Infantile diarrhoea |
| EIEC | ipaH | Dysentery-like illness |
| ETEC | LT and/or ST-toxin | Travellers' diarrhoea |
| − | None of the above | ExPEC or commensal |

Note: Kleborate can produce compound calls (e.g. `STEC/EPEC`) when multiple marker combinations are present. Each component is counted independently in Fig 5.

---

## 5. Kaptive — Capsular K-locus typing

**Tool**: [Kaptive](https://github.com/klebsimonics/Kaptive) v3  
**Databases**: Group 2/3 (Gladstone et al. 2024) + Group 1/4 (custom, efosternyarko/EC-K-typing-G1G4)  
**Output columns**: `k_group`, `k_locus`, `k_type`, `k_confidence`

Kaptive types the *E. coli* capsular K-locus — the genomic region encoding the polysaccharide capsule biosynthesis machinery. The pipeline implements a two-step workflow:

### Step 1 — Group 2/3 database

All assemblies are typed against the Kaptive Group 2/3 reference database. Samples that receive a confident assignment are classified as Group 2 (K-loci with *kpsC/kpsS* gene cluster) or Group 3.

### Step 2 — Group 1/4 database (untypeables only)

Assemblies that are "Untypeable" under the Group 2/3 database are re-typed using the Group 1/4 database in `--scores` mode. Raw alignment scores are normalised by expected locus CDS length to achieve full typeability on G1/G4 strains.

**Why two separate databases?** Group 1/4 and Group 2/3 loci share some synteny and database cross-hits can cause misclassification (*wzy* interference). Running both simultaneously gives incorrect calls on a subset of samples.

### K-locus groups and clinical significance

| Group | Biosynthesis pathway | Notable associations |
|---|---|---|
| G1 | Wzy-dependent lipopolysaccharide-linked | Invasive ExPEC; neonatal meningitis; K1, K5 |
| G2 | ABC transporter-dependent | Diverse; K2 associated with virulent ExPEC |
| G3 | ABC transporter-dependent | Diverse |
| G4 | Wzy-dependent | Uropathogenic *E. coli*; K30, K40 |

K1 (Group 1) is the most common capsule type in neonatal meningitis isolates. K2 (Group 2) is associated with virulent liver abscess *E. coli* (Klebsiella-like presentation). K5 (Group 1) is common in UTI isolates.

### Output fields

| Column | Description |
|---|---|
| `k_group` | Broad group: G1, G2, G3, G4, or Unknown |
| `k_locus` | Specific locus designation (e.g. KL1, KL2) |
| `k_type` | Antigenic K-type where known (e.g. K1, K2) |
| `k_confidence` | Kaptive match confidence: Perfect / Good / Low / Untypeable |

---

## 6. PlasmidFinder — Replicon typing

**Tool**: [PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/) (CGE)  
**Database**: `enterobacteriaceae`  
**Output column**: `plasmidfinder_replicons`

Screens the assembly for known plasmid replicon sequences (Inc types). Results are used directly in Fig 4 (plasmid overview). See [Fig 4 — Plasmid Overview](Fig4-Plasmid-Overview.md) for interpretation.

---

## Tool versions (default containers)

| Tool | Version |
|---|---|
| mlst | 2.23.0 |
| ECTyper | 1.0.0 |
| EzClermont | current conda |
| Kleborate | 3.1.3 |
| Kaptive | 3.0.0b5 |
| AMRFinder Plus | see `amrfinder --version` |
| PlasmidFinder | current conda |
