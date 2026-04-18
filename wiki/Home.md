# enteric-typer Wiki

**enteric-typer** is a species-gated genotyping pipeline for enteric pathogens. It accepts genome assemblies, automatically identifies the species of each sample, runs the appropriate typing tools in parallel, builds a whole-genome SNP phylogeny, and produces publication-ready summary figures.

---

## Navigation

### Pipeline logic
- [Species Detection](Species-Detection) — how Mash-based species assignment works and why
- [Assembly QC](Assembly-QC) — genome size filters, contamination screening, what fails and why
- [AMR Classification](AMR-Classification) — AMRrules, intrinsic vs acquired resistance, drug class definitions
- [SNP Phylogenetics](SNP-Phylogenetics) — SKA2 alignment, IQ-TREE inference, what distances mean

### Per-species typing tools
- [E. coli Typing](E-coli-Typing) — MLST, ECTyper, EzClermont, Kleborate, Kaptive
- [Salmonella Typing](Salmonella-Typing) — MLST, SISTR, serovar assignment
- [Shigella Typing](Shigella-Typing) — ShigEiFinder, Mykrobe, pINV screen, IS element screen

### Output figures — interpretation guide
- [Fig 1 — Population Summary](Fig1-Population-Summary) — ST distribution, serovar/serotype bars, AMR prevalence, MDR
- [Fig 2 — Phylogeny and AMR](Fig2-Phylogeny-and-AMR) — whole-genome SNP tree annotated with typing and resistance
- [Fig 3 — AMR Genes](Fig3-AMR-Genes) — acquired resistance gene prevalence
- [Fig 4 — Plasmid Overview](Fig4-Plasmid-Overview) — replicon prevalence, AMR co-occurrence, auto-layout logic
- [Fig 5 — Virulence](Fig5-Virulence) — virulence factor prevalence
- [Fig 6 — AMR by ST](Fig6-AMR-by-ST) — drug class prevalence stratified by MLST sequence type
- [Fig 7 — AMR by Group](Fig7-AMR-by-Group) — drug class prevalence by serovar / phylogroup / serotype
- [Shigella-Specific Figures](Shigella-Specific-Figures) — serotype composition (Fig 8), virulence/invasion panel (Fig 9)
- [SNP Distance Heatmap](SNP-Distance-Heatmap) — pairwise whole-genome SNP distances

### Reference
- [Output File Reference](Output-File-Reference) — all TSV columns explained
- [Parameters](Parameters) — every `--flag` and its default
- [FAQ and Troubleshooting](FAQ-and-Troubleshooting)

---

## Quick orientation

```
Input assemblies
       │
       ▼
  Species detection (Mash)
       │
  ┌────┴────────────────┐
  ▼         ▼           ▼
E. coli  Salmonella  Shigella  → other: logged, skipped
  │         │           │
  └────┬────┘           │
       ▼                ▼
  Typing tools     Typing tools
  (in parallel)    (in parallel)
       │
       ▼
  SNP phylogeny (SKA2 + IQ-TREE)
       │
       ▼
  Aggregate → TSV results
       │
       ▼
  Summary figures (PDF + PNG)
```

Each stage is described in detail in the linked pages above.
