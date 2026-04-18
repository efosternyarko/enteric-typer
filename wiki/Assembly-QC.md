# Assembly QC

Before any typing tools run, each assembly is evaluated against quality thresholds. Assemblies that fail are logged and excluded from all downstream analysis — they do not cause the pipeline to abort.

---

## QC checks applied

### 1. Genome size

The total assembly length (sum of all contigs) must fall within species-specific bounds. Assemblies outside this range are flagged as likely contaminated, misassembled, or mislabelled.

| Species | Minimum (bp) | Maximum (bp) | Source |
|---|---|---|---|
| *E. coli* / *Shigella* | 4,300,000 | 5,900,000 | BIGSdb *E. coli* genome quality criteria |
| *Salmonella enterica* | 4,100,000 | 6,600,000 | PATH-SAFE consortium recommendations (Table 3) |

These thresholds can be overridden per-run:
```
--ecoli_min_length 4000000
--ecoli_max_length 6200000
--shigella_min_length 4000000
--salmonella_min_length 4000000
```

**Why size matters:**
- A genome much smaller than expected typically indicates incomplete assembly, high contamination, or an incorrect species assignment
- A genome much larger than expected suggests contamination with a second organism or a misassembly creating duplicated regions

### 2. Contamination screening (optional)

If `--kraken2_db` is provided, each assembly is screened with [Kraken2](https://github.com/DerrickWood/kraken2) against a standard bacterial database. The proportion of contigs classified as a secondary species is computed.

An assembly fails if the secondary species makes up > 3% of total contigs (`--max_contamination 3.0`).

This check is **skipped** when `--kraken2_db` is not provided and does not block the run.

---

## Decision tree

```
Assembly received
       │
       ▼
Is total length within species range?
       │
       ├── NO  → mark FAIL (size), log, skip all typing
       │
       └── YES
               │
               ▼
       Is --kraken2_db provided?
               │
               ├── NO  → proceed to typing
               │
               └── YES
                       │
                       ▼
               Kraken2: secondary species < 3%?
                       │
                       ├── NO  → mark FAIL (contamination), log, skip
                       │
                       └── YES → proceed to typing
```

---

## What counts as "failed"

Failed samples are:
- Written to the pipeline log with the reason for failure
- Excluded from all typing outputs (no row in the results TSV)
- Excluded from all figures (not plotted)
- Counted in the pipeline summary so you know how many samples did not make it through

A failed assembly does **not** stop the pipeline from running for other samples.

---

## Assembly metrics figure

Even for samples that pass QC, the **assembly metrics figure** (`{species}_assembly_metrics.{pdf,png}`) summarises per-assembly statistics across all samples of that species:

| Panel | Metric | Expected range |
|---|---|---|
| A / B | Genome length (bp) | *E. coli*: 4.5–5.5 Mb; *Salmonella*: 4.5–5.2 Mb |
| C / D | N50 (bp) | Higher is better; depends heavily on sequencing technology |
| E / F | Number of contigs | Short-read: 50–300; Long-read/hybrid: 1–20 |
| G / H | GC content (%) | *E. coli* / *Shigella*: ~51%; *Salmonella*: ~52% |

Each metric is shown as both a histogram and a box plot. Dashed reference lines indicate the expected values for the species. Assemblies that are outliers in multiple metrics simultaneously warrant manual inspection.

---

## Notes on N50 and contig count

**N50** is the contig length such that 50% of the total assembled sequence is in contigs of that length or longer. It is a summary of assembly contiguity:
- Short-read (Illumina) assemblies typically have N50 of 50–500 kb and 50–400 contigs
- Long-read (Nanopore/PacBio) assemblies typically have N50 of 2–6 Mb and 1–5 contigs
- Hybrid assemblies typically fall between these ranges

A very low N50 (< 20 kb) or very high contig count (> 1000) may indicate a poor-quality sequencing run or failed assembly step. These samples may still pass the size QC check but should be interpreted with caution.
