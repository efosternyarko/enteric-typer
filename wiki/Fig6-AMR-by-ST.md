# Fig 6 — AMR by Sequence Type

Fig 6 is an AMRnet-style tile heatmap showing the percentage of isolates carrying acquired resistance in each drug class, broken down by MLST sequence type. It is the primary figure for understanding which lineages drive resistance in your dataset.

---

## What is shown

A grid with:
- **Rows** — the top 12 most common MLST sequence types (STs) by isolate count, plus "Unknown" if present. Sample counts are shown in the row label (e.g. "ST131  (n=47)").
- **Columns** — all drug classes with any detected acquired resistance, sorted by overall prevalence (highest prevalence = leftmost column). Drug class abbreviations follow the [AMR Classification](AMR-Classification.md) table (BLA, AMG, SUL, etc.).
- **Tiles** — the percentage of isolates in that ST that carry acquired resistance to that drug class.

---

## Colour scale

The colour scale mirrors the AMRnet tool used by WHO/GLASS:

| Tile colour | Prevalence |
|---|---|
| Mid-grey | 0% — no resistance detected in this ST |
| Cream / pale yellow | 1–25% |
| Orange | ~50% |
| Dark red | ~75–80% |
| Dark red/purple | 100% |

The percentage is printed inside tiles where non-zero. Values ≥ 1% are rounded to the nearest integer; values between 0 and 1% are shown to one decimal place.

---

## Reading the figure

**A uniformly grey row**: the ST has no detected acquired resistance in any class. This is biologically possible (e.g. a truly susceptible lineage) but may also reflect a small sample size for that ST — check the n count.

**A uniformly dark row**: the ST is broadly MDR, with high resistance prevalence across many classes. ST131 in *E. coli* (ESBL + aminoglycoside + fluoroquinolone) and ST313 in *S.* Typhimurium (multidrug resistant African isolates) are canonical examples.

**A column with one very dark tile and all others grey or light**: a specific ST drives resistance to that drug class, with other STs showing low prevalence. This pattern suggests clonal dissemination of a resistance determinant in that lineage.

**All tiles in a column are light to moderate**: the resistance is distributed across multiple STs, suggesting the resistance gene is on mobile elements circulating across lineages rather than being clonally disseminated.

---

## Important caveats

### Sample size per ST

Each row tile represents 100% × (resistant isolates in that ST / total isolates in that ST). For an ST with n=2, one resistant isolate gives 50% and two gives 100%. Always read the n= label in the row name before interpreting a tile value.

### Missing STs

Only the top 12 STs are shown. STs not in the top 12 are excluded from this figure regardless of their resistance profile. To investigate resistance in rare STs, use the raw TSV.

### Drug classes shown

EFFLUX is excluded from all AMRnet figures. Classes with zero prevalence across all isolates are also omitted (empty columns add no information).

---

## Output file

`{prefix}_fig6_amr_by_st.pdf` and `.png`
