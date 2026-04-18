# Fig 7 — AMR by Serovar / Phylogroup

Fig 7 is a second AMRnet-style tile heatmap, using the same colour scale and tile format as [Fig 6](Fig6-AMR-by-ST.md), but stratified by a biologically meaningful grouping that differs by species.

---

## Rows by species

| Species | Row grouping | Source column |
|---|---|---|
| *E. coli* | Clermont phylogroup | `clermont_phylogroup` |
| *Salmonella enterica* | Serovar | `sistr_serovar` (SISTR) |
| *Shigella* spp. | Serotype, or species if serotypes unavailable | `shigeifinder_serotype` |

### *E. coli*: Clermont phylogroups

Rows are the Clermont phylogroups (A, B1, B2, C, D, E, F) present in the dataset, plus Unknown. Up to 10 rows shown.

Phylogroup is a coarser grouping than ST but captures the major ecological and clinical niches:

- **B2**: ExPEC — urinary tract infections, bacteraemia, neonatal meningitis; dominates clinical settings in high-income countries
- **A**: ubiquitous commensal; many broad-host-range plasmids circulate
- **B1**: livestock and environmental; many produce ESBL in zoonotic contexts
- **D**: clinical; ST131 and other fluoroquinolone-resistant clones
- **C/E**: rarer; C includes ST410 (carbapenem-resistant lineage)

### *Salmonella*: Serovars

Top 15 serovars by isolate count. Salmonella resistance is strongly serovar-associated:

- **Typhimurium** / **I 4,[5],12:i:-** (monophasic Typhimurium): MDR common via DT104 chromosomal cassette or IncFII plasmids
- **Enteritidis**: often more susceptible, but increasing ESBLs in some regions
- **Kentucky**: globally disseminated fluoroquinolone-resistant lineage
- **Typhi**: azithromycin (macrolide) resistance — monitor MAC column

### *Shigella*: Serotypes or species

If ShigEiFinder serotype calls are available, the top 15 serotypes are used as rows. If serotype calls are absent (all NA), rows fall back to species-level grouping (*S. sonnei*, *S. flexneri*, *S. boydii*, *S. dysenteriae*).

---

## Columns

Drug classes sorted by overall prevalence (same logic as Fig 6). EFFLUX excluded.

---

## Interpreting Fig 7 vs Fig 6

Figs 6 and 7 answer complementary questions:

| Question | Fig 6 (by ST) | Fig 7 (by group) |
|---|---|---|
| Which lineage drives resistance? | Yes — ST-level resolution | No — coarser grouping |
| Is resistance serovar/phylogroup-wide or ST-specific? | Partially (multiple STs per serovar visible only in raw data) | Yes — shows group-level picture |
| Are all isolates of serovar X resistant to Y? | Not directly | Yes |
| How does my commensal collection compare to clinical? | No | Yes (by phylogroup) |

Use both figures together: if a drug class tile is dark in Fig 7 but the same class shows dark tiles only in one or two STs in Fig 6, resistance is concentrated in specific clones within that group. If multiple STs show dark tiles in Fig 6, resistance is lineage-independent (likely a mobile element).

---

## Output file

`{prefix}_fig7_amr_by_group.pdf` and `.png`
