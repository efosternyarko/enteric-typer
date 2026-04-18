# Output File Reference

All typing results for each species are aggregated into a single tab-separated results file. This page documents every column in these files. Missing data is represented as `NA`.

---

## Files produced

| File | Species |
|---|---|
| `{outdir}/ecoli_typer_results.tsv` | *E. coli* |
| `{outdir}/salmonella_typer_results.tsv` | *Salmonella enterica* |
| `{outdir}/shigella_typer_results.tsv` | *Shigella* spp. |

---

## Columns common to all species

| Column | Source | Description |
|---|---|---|
| `sample` | Input | Sample identifier |
| `mlst_scheme` | mlst | MLST scheme used |
| `mlst_st` | mlst | Predicted sequence type (e.g. ST131) |
| `mlst_st_complex` | mlst + lookup | ST complex assignment (e.g. ST131 complex) |
| `amrfinder_acquired_genes` | AMRFinder + AMRrules | Semicolon-separated acquired AMR genes and mutations only (intrinsic excluded) |
| `amrfinder_intrinsic_genes` | AMRFinder + AMRrules | Semicolon-separated intrinsic/wild-type genes (reference only; not in figures) |
| `amrfinder_virulence_genes` | AMRFinder | Semicolon-separated virulence-associated genes |
| `amrfinder_genes` | AMRFinder | Raw unfiltered gene list (all types, for reference) |
| `amrfinder_drug_classes` | AMRFinder + AMRrules | Semicolon-separated drug class list (acquired genes only) |
| `amrfinder_gene_classes` | AMRFinder + AMRrules | `gene=CLASS;gene2=CLASS2` mapping for acquired genes |
| `plasmidfinder_replicons` | PlasmidFinder | Semicolon-separated replicon types detected |

---

## *E. coli*-specific columns

| Column | Source | Description |
|---|---|---|
| `clermont_phylogroup` | EzClermont | Clermont phylogroup (A, B1, B2, C, D, E, F, or Unknown) |
| `kleborate_pathovar` | Kleborate | Pathotype assignment: STEC, EPEC, ETEC, EIEC, EHEC, or `-` (ExPEC/commensal) |
| `kleborate_Stx1` | Kleborate | Shiga toxin 1 call (subtype or `-`) |
| `kleborate_Stx2` | Kleborate | Shiga toxin 2 call (subtype or `-`) |
| `kleborate_eae` | Kleborate | Intimin (*eae*) call (`+` or `-`) |
| `kleborate_ipaH` | Kleborate | ipaH call (`+` or `-`) |
| `kleborate_LT` | Kleborate | Heat-labile toxin call (`+` or `-`) |
| `kleborate_ST_toxin` | Kleborate | Heat-stable toxin call (`+` or `-`) |
| `ectyper_O` | ECTyper | O antigen (e.g. O25) |
| `ectyper_H` | ECTyper | H antigen (e.g. H4) |
| `ectyper_serotype` | ECTyper | Combined O:H serotype (e.g. O25:H4) |
| `ectyper_qc` | ECTyper | QC status: PASS / WARN / FAIL |
| `ectyper_evidence` | ECTyper | Number of marker genes supporting the call |
| `k_group` | Kaptive | K-locus group: G1, G2, G3, G4, or Unknown |
| `k_locus` | Kaptive | K-locus designation (e.g. KL1) |
| `k_type` | Kaptive | Antigenic K-type (e.g. K1) |
| `k_confidence` | Kaptive | Match confidence: Perfect / Good / Low / Untypeable |

### PathogenWatch columns (*E. coli* and *Salmonella*)

These columns are populated when a PathogenWatch API query is made (optional feature).

| Column | Description |
|---|---|
| `pw_status` | PathogenWatch query status (matched / not_found) |
| `pw_species` | Species confirmed by PathogenWatch |
| `pw_genome_uuid` | PathogenWatch genome UUID |
| `pw_collection_url` | URL to PathogenWatch collection |
| `pw_cgmlst_st` | cgMLST ST assigned by PathogenWatch |
| `pw_cluster5_count` | Number of genomes within SNP distance 5 in PathogenWatch |
| `pw_cluster5_labels` | Sample labels within SNP distance 5 |
| `pw_cluster10_count` | Number within SNP distance 10 |
| `pw_cluster10_labels` | Labels within SNP distance 10 |
| `pw_cluster20_count` | Number within SNP distance 20 |
| `pw_cluster20_labels` | Labels within SNP distance 20 |
| `pw_cluster50_count` | Number within SNP distance 50 |
| `pw_cluster50_labels` | Labels within SNP distance 50 |
| `pw_tree_available` | Whether a PathogenWatch tree is available for this genome |

---

## *Salmonella*-specific columns

| Column | Source | Description |
|---|---|---|
| `sistr_serovar` | SISTR | Consensus serovar prediction (Kauffmann-White name) |
| `sistr_serovar_antigen` | SISTR | Serovar from antigen typing method alone |
| `sistr_serovar_cgmlst` | SISTR | Serovar from cgMLST method alone |
| `sistr_O` | SISTR | O-antigen formula (e.g. `4,5,12`) |
| `sistr_H1` | SISTR | Phase 1 H antigen |
| `sistr_H2` | SISTR | Phase 2 H antigen (`-` = monophasic) |
| `sistr_cgmlst_ST` | SISTR | cgMLST sequence type |
| `sistr_qc` | SISTR | QC status: PASS / WARNING / FAIL |

---

## *Shigella*-specific columns

| Column | Source | Description |
|---|---|---|
| `shigeifinder_ipaH` | ShigEiFinder | ipaH marker: `pos` or `neg` |
| `shigeifinder_virulence_plasmid` | ShigEiFinder | Virulence plasmid classification |
| `shigeifinder_cluster` | ShigEiFinder | Cluster ID |
| `shigeifinder_serotype` | ShigEiFinder | Serotype (e.g. `SS`, `SF2a`) |
| `shigeifinder_o_antigen` | ShigEiFinder | O-antigen |
| `shigeifinder_h_antigen` | ShigEiFinder | H-antigen |
| `mykrobe_genotype` | Mykrobe | Full genotype string |
| `mykrobe_lineage` | Mykrobe | Top-level lineage (e.g. lineage3) — *S. sonnei* only |
| `mykrobe_clade` | Mykrobe | Clade within lineage |
| `mykrobe_subclade` | Mykrobe | Sub-clade |
| `mykrobe_genotype_name` | Mykrobe | Human-readable genotype name |
| `mykrobe_confidence` | Mykrobe | Confidence of genotype call |
| `pinv_present` | pINV screen | `TRUE` if any pINV marker gene detected |
| `pinv_genes` | pINV screen | Semicolon-separated detected pINV genes with identity/coverage |
| `is_elements` | IS screen | IS element copy numbers: `IS1(3);IS600(12);...` |

---

## Parsing the semicolon-separated list columns

Most multi-value columns use `;` as separator. To parse in R:

```r
library(dplyr)
library(tidyr)
df |>
  separate_rows(amrfinder_acquired_genes, sep = ";") |>
  filter(!is.na(amrfinder_acquired_genes), amrfinder_acquired_genes != "NA")
```

To parse in Python:

```python
df["amrfinder_acquired_genes"].str.split(";").explode().dropna()
```

The `amrfinder_gene_classes` column uses `=` within each semicolon-separated token (e.g. `blaOXA-1=BETA-LACTAM;aph(3'')-Ib=AMINOGLYCOSIDE`). Parse with:

```python
{g.split("=")[0]: g.split("=")[1]
 for gene_cls in df["amrfinder_gene_classes"].dropna()
 for g in gene_cls.split(";") if "=" in g}
```
