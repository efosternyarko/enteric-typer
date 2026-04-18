# Salmonella Typing Tools

*Salmonella enterica* assemblies are run through MLST and SISTR for lineage and serovar assignment, plus AMRFinder and PlasmidFinder for resistance and plasmid characterisation. All results are aggregated into `salmonella_typer_results.tsv`.

---

## 1. MLST — Multi-locus sequence typing

**Tool**: [mlst](https://github.com/tseemann/mlst) (Torsten Seemann)  
**Scheme**: `senterica_achtman_2` — the Warwick Achtman 7-locus scheme for *S. enterica*  
**Output columns**: `mlst_scheme`, `mlst_st`, `mlst_st_complex`

The *Salmonella* Achtman scheme uses seven housekeeping loci (*aroC*, *dnaN*, *hemD*, *hisD*, *purE*, *sucA*, *thrA*). STs are mapped to ST complexes (e.g. ST19 complex, ST313 complex) using a curated lookup table.

### Common clinically important *Salmonella* STs

| ST | Serovar | Notable features |
|---|---|---|
| ST19 | Typhimurium | Global dominance; diverse plasmids; often MDR |
| ST34 | I 4,[5],12:i:- | Monophasic Typhimurium; pandemic MDR lineage |
| ST11 | Enteritidis | Dominant in poultry; often relatively susceptible |
| ST1 | Typhi | Typhoid fever; *Salmonella* Typhi exclusively |
| ST2 | Typhi | Typhoid fever |
| ST313 | Typhimurium ST313 | African invasive NTS (iNTS); MDR |
| ST198 | Kentucky | Globally disseminated fluoroquinolone-resistant clone |

---

## 2. SISTR — Salmonella in silico serotyping

**Tool**: [SISTR](https://github.com/peterk87/sistr_cmd) v1.1.2  
**Output columns**: `sistr_serovar`, `sistr_serovar_antigen`, `sistr_serovar_cgmlst`, `sistr_O`, `sistr_H1`, `sistr_H2`, `sistr_cgmlst_ST`, `sistr_qc`

SISTR predicts *Salmonella* serovar from assembled genomes by combining three independent approaches and returning a consensus call:

### Three prediction methods

**1. Antigen sequence typing** (primary)  
BLASTs the assembly against reference sequences for the flagellar antigens (*fliC*, *fljB*) and O-antigen biosynthesis genes (*wzx/wzy*, *wbgA*). Assigns H1, H2, and O antigens according to the Kauffmann-White scheme.

**2. cgMLST** (corroborative)  
Uses a 330-locus *Salmonella*-specific core genome MLST scheme to predict serovar. Provides an independent check and produces a cgMLST sequence type.

**3. mash sketch comparison** (fallback)  
For assemblies where antigen typing fails, compares assembly MinHash sketches to a sketch database of reference serovars.

### Output fields

| Column | Description |
|---|---|
| `sistr_serovar` | Consensus predicted serovar (Kauffmann-White name, e.g. Typhimurium) |
| `sistr_serovar_antigen` | Serovar predicted from antigen typing alone |
| `sistr_serovar_cgmlst` | Serovar predicted from cgMLST alone |
| `sistr_O` | O-antigen formula (e.g. `4,5,12`) |
| `sistr_H1` | Phase 1 H antigen (e.g. `i`) |
| `sistr_H2` | Phase 2 H antigen (e.g. `1,2` or `-` for monophasic) |
| `sistr_cgmlst_ST` | cgMLST sequence type |
| `sistr_qc` | QC status: PASS / WARNING / FAIL |

### The monophasic Typhimurium variant

*S.* Typhimurium variant **I 4,[5],12:i:-** is a monophasic form (H2 antigen = `-`) that has spread globally as a pandemic MDR lineage (predominantly ST34). SISTR correctly assigns this as a distinct serovar from classical Typhimurium (I 4,[5],12:i:1,2). Both appear prominently in European and African datasets.

### When SISTR disagrees with MLST

SISTR serovar and MLST ST are correlated but not equivalent. Multiple STs can share a serovar (e.g. Typhimurium includes ST19, ST34, ST313), and a single ST can produce different serovars if the O/H antigen genes are on mobile elements. If SISTR and MLST give apparently contradictory results, check:

1. The QC status (`sistr_qc`) — a FAIL or WARNING may indicate assembly fragmentation affecting antigen gene detection
2. The `sistr_serovar_antigen` vs `sistr_serovar_cgmlst` columns — agreement between the two methods increases confidence
3. Whether the discrepant serovar-ST combination is documented (some are; e.g. ST19 is associated with both Typhimurium and Agona)

---

## 3. AMRFinder Plus

**Organism flag**: `--organism Salmonella`  
Activates species-specific point mutation databases for AMRFinder, including chromosomal resistance mutations in *gyrA*, *parC* (quinolones), *pmrA/B*, *mgrB* (colistin — limited), and others.

See [AMR Classification](AMR-Classification.md) for full documentation.

---

## 4. PlasmidFinder

**Database**: `enterobacteriaceae`  
**Output column**: `plasmidfinder_replicons`

See [Fig 4 — Plasmid Overview](Fig4-Plasmid-Overview.md) for interpretation of replicon results.

Common *Salmonella* replicon types and their resistance associations:

| Replicon | Common cargo | Notes |
|---|---|---|
| IncFII(S) | MDR cassettes (blaTEM-1, tetAB, sul1, cat) | Classic DT104-related MDR plasmids |
| IncFIB(S) | Often co-resident with IncFII(S) | |
| IncI1 | blaCTX-M ESBLs, tetA | Broad host range; common in NTS |
| IncA/C2 | blaCMY, aac, sul | Cephalosporin resistance |
| IncX1 | blaNDM carbapenemase | Carbapenem resistance — alert |
| IncHI2 | mcr-1 colistin resistance | Last-resort resistance — alert |
| Col156 | Cryptic; no resistance | Background plasmid |
