# Fig 5 — Virulence Profiling

Fig 5 summarises virulence gene carriage across isolates. The figure layout differs by species because the biology of virulence differs substantially between *E. coli*, *Salmonella*, and *Shigella*.

---

## *E. coli* — three-panel layout

Fig 5 for *E. coli* has up to three panels depending on which virulence data are available in the results TSV.

### Panel A — Pathotype distribution (Kleborate)

A horizontal bar chart showing the number of isolates assigned to each diarrhoeagenic *E. coli* pathotype, as classified by Kleborate from the `kleborate_pathovar` column.

| Colour | Pathotype | Defining marker(s) |
|---|---|---|
| Dark red | EHEC | Stx + eae (intimin) |
| Red | STEC | Stx without eae |
| Orange | EPEC | eae without Stx |
| Purple | EIEC | ipaH (invasion plasmid antigen H) |
| Blue | ETEC | LT and/or ST enterotoxin |
| Grey | ExPEC | None of the above; presumed extra-intestinal pathogen |

Note: some isolates may carry multiple markers (e.g. Stx + LT) producing compound calls separated by `/`. Each component pathotype is counted independently.

**Isolates labelled "ExPEC"** carry none of the defined diarrhoeagenic markers but may still carry ExPEC virulence determinants (iron acquisition, toxins, adhesins) — these are shown in the AMRFinder panel (Panel C).

### Panel B — InPEC virulence marker prevalence (Kleborate)

A dot-plot showing the percentage of isolates with each intestinal pathogenic *E. coli* (InPEC) marker. Only markers present in at least one isolate are shown.

| Marker | Pathotype association | Notes |
|---|---|---|
| Stx1 | STEC / EHEC | Heat-stable Shiga toxins; renal damage |
| Stx2 | STEC / EHEC | More strongly associated with HUS than Stx1 |
| eae | EPEC / EHEC | Intimin; required for attaching-and-effacing lesion |
| ipaH | EIEC | Invasion plasmid antigen H; crosses epithelium |
| LT | ETEC | Heat-labile cholera-like toxin |
| ST-toxin | ETEC | Heat-stable toxin |

### Panel C — AMRFinder virulence genes

A horizontal bar chart of the top 25 most prevalent virulence-associated genes detected by AMRFinder Plus (`--plus` flag), from the `amrfinder_virulence_genes` column.

These include ExPEC virulence factors that Kleborate does not specifically classify:
- Iron acquisition: *iutA* (aerobactin receptor), *fyuA* (yersiniabactin receptor), *iucABCD* (aerobactin synthesis)
- Toxins: *hlyA* (alpha-haemolysin), *cnf1* (cytotoxic necrotising factor), *sat* (secreted autotransporter toxin)
- Adhesins: *papA/papC* (P fimbriae), *afa*, *sfa*
- Protectins: *kpsMII*, *traT* (serum resistance)

---

## *Salmonella* — single-panel layout

A horizontal bar chart of the top 25 most prevalent virulence genes, sourced from:
1. Abricate VFDB results (`abricate_vfdb_genes`) when available
2. AMRFinder virulence genes as fallback

For *Salmonella*, the most epidemiologically relevant virulence features are typically the presence of pathogenicity islands (SPI-1 to SPI-6) and plasmid-encoded virulence genes. These are reported by VFDB but are often present in virtually all isolates (near-100% prevalence). The figure therefore functions primarily as a confirmation that expected virulence genes are present and to flag any unusual patterns.

**Note**: VFDB-based virulence gene detection for *Salmonella* is less standardised than for *E. coli*, and the informative content of Fig 5 for *Salmonella* is more limited. The species-level comparison (which serovar/ST carries which additional factors) is better addressed by custom analysis of the raw AMRFinder TSV.

---

## *Shigella* — Fig 5 is redirected

For *Shigella*, the virulence profiling is handled by the dedicated *Shigella*-specific figures:
- Invasion gene and IS element details are in [Fig 9 — Shigella Virulence Features](Shigella-Specific-Figures.md)
- Serotype composition is in [Fig 8 — Shigella Serotypes](Shigella-Specific-Figures.md)

Fig 5 is therefore not produced for *Shigella* runs.

---

## Output file

`{prefix}_fig5_virulence.pdf` and `.png`
