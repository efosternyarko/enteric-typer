# AMR Classification

enteric-typer uses [AMRFinder Plus](https://github.com/ncbi/amr) to detect antimicrobial resistance (AMR) genes and mutations, then applies [AMRrules](https://github.com/AMRverse/AMRrules) to classify every detected gene as either **acquired** or **intrinsic**. This distinction is fundamental to how AMR results are presented across all figures and output tables.

---

## AMRFinder Plus: what it detects

AMRFinder Plus screens each assembly against the NCBI AMR Reference Gene Database and reports:

- **AMR genes** — genes conferring resistance by any mechanism (enzyme, efflux, target modification, etc.)
- **Point mutations** — chromosomal mutations known to confer resistance (e.g. *gyrA* S83L for fluoroquinolones)
- **Virulence genes** — pathogenicity-associated genes (reported separately; used in Fig 5 and the phylogeny annotation)
- **Stress response genes** — genes conferring resistance to biocides, metals, heat (reported but not plotted by default)

The pipeline runs AMRFinder with `--plus` (enables virulence and stress gene detection) and passes the appropriate `--organism` flag so that species-specific point mutation databases are applied.

---

## AMRrules: intrinsic vs acquired

Not all AMR genes detected by AMRFinder are clinically or epidemiologically meaningful. Many are **intrinsic** — present in the wild-type genome of every isolate of that species and not acquired by horizontal gene transfer. Reporting these alongside genuinely acquired resistance creates noise that obscures real resistance signals.

AMRrules provides a curated, species-specific list of genes that are intrinsic to each enteric pathogen. The pipeline applies these rules to every sample and splits the AMRFinder output into two lists:

| Column | Definition |
|---|---|
| `amrfinder_acquired_genes` | Genes and mutations **not** on the intrinsic list — clinically relevant acquired resistance |
| `amrfinder_intrinsic_genes` | Genes on the intrinsic list — wild-type background, recorded but not plotted |
| `amrfinder_drug_classes` | Drug classes represented by acquired genes only |
| `amrfinder_genes` | All raw AMRFinder output (for reference) |

**All AMR figures use acquired genes only.** Intrinsic genes are retained in the TSV for completeness and transparency but are never included in prevalence bars, bubble matrices, heatmaps, or MDR counts.

### Examples of intrinsic genes excluded

| Species | Gene | Drug class | Why intrinsic |
|---|---|---|---|
| *E. coli* / *Shigella* | `blaEC` | Beta-lactam | Chromosomal AmpC; present in all wild-type *E. coli* |
| *E. coli* / *Shigella* | `acrB`, `acrF`, `mdtM`, `emrD` | Efflux | Constitutive efflux pumps; chromosomally encoded |
| *Salmonella* | `mdsA`, `mdsB` | Efflux | Salmonella-specific efflux complex |
| *Salmonella* | `acrB` | Efflux | As for *E. coli* |

---

## Drug class abbreviations

All figures use the following short codes on axes and in legends:

| Code | Full class | Clinical relevance |
|---|---|---|
| AMG | Aminoglycoside | Gentamicin, amikacin, tobramycin — bacteraemia, complicated UTI |
| BLA | Beta-lactam | Ampicillin, cephalosporins, carbapenems — first-line for many infections |
| COL | Colistin | Last-resort for carbapenem-resistant organisms |
| FOS | Fosfomycin | UTI treatment; increasing use for ESBL isolates |
| FOSM | Fosmidomycin | Investigational; rarely reported in enteric pathogens |
| LIN | Lincosamide | Clindamycin; uncommon in enteric pathogens |
| MAC | Macrolide | Azithromycin — key drug for *Shigella sonnei* treatment |
| NIT | Nitrofuran | Nitrofurantoin — uncomplicated UTI |
| PHE | Phenicol | Chloramphenicol — historically important for typhoid/NTS |
| QNL | Quinolone | Ciprofloxacin, nalidixic acid — critically important for enteric infections |
| SGM | Streptogramin | Quinupristine–dalfopristine; rarely relevant for enteric pathogens |
| STR | Streptothricin | Nourseothricin; rare but carried on plasmids |
| SUL | Sulfonamide | Sulfamethoxazole — component of co-trimoxazole |
| TET | Tetracycline | Broad resistance marker; common on IncFII plasmids |
| TMP | Trimethoprim | Component of co-trimoxazole — widely used for NTS in Africa |

---

## Multi-drug resistance (MDR) definition

An isolate is classified as **multidrug resistant (MDR)** if it carries acquired resistance genes in **≥ 3** of the above drug classes.

- The **EFFLUX** class is excluded from MDR counting. Efflux-mediated resistance to multiple classes via constitutive pumps is not considered acquired MDR because it is background for the species.
- Only `amrfinder_acquired_genes` (AMRrules-filtered) are counted — intrinsic genes do not contribute to MDR status.

This definition follows WHO/ECDC categorical definitions for MDR in Gram-negative bacteria.

---

## QUINOLONE/TRICLOSAN

AMRFinder may report some mutations (e.g. *gyrA* point mutations at certain positions) under a combined class label `QUINOLONE/TRICLOSAN`. These are counted as **QUINOLONE** resistance for MDR classification and figure display, since clinical relevance in enteric pathogens is driven by fluoroquinolone use, not triclosan.
