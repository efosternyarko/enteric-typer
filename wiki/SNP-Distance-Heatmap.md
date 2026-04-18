# SNP Distance Heatmap

The SNP distance heatmap visualises pairwise whole-genome SNP distances between all isolates of a species as a clustered square matrix. It is produced alongside the phylogenetic tree and complements it by giving exact numerical distances where the tree gives topological relationships.

---

## What is shown

A symmetric n × n matrix where:
- Each row and column is one isolate
- Cell colour encodes the pairwise SNP distance: **white/pale blue = close; dark blue = distant**
- The colour scale runs from 0 to the 95th percentile of non-zero distances (to prevent outliers from compressing the useful colour range)
- A colour bar on the right side shows the mapping from colour to SNP count

Rows and columns are reordered by **average-linkage hierarchical clustering** of the distance matrix, so isolates that are closest to each other are adjacent. This brings outbreak clusters and clonal groups into visually obvious blocks.

---

## Scale and labels

| n isolates | Labels shown |
|---|---|
| ≤ 200 | Sample names on both axes (font size scales down with n) |
| > 200 | No labels; x-axis message states count and threshold |

When labels are shown, x-axis labels are rotated 90° for readability.

---

## Interpreting the heatmap

### Dark diagonal

The main diagonal is always 0 (each sample compared to itself). If the diagonal does not appear as the darkest region in each row/column, check that sample names were not duplicated.

### Off-diagonal dark blocks

A dark blue square block off the diagonal indicates two groups of isolates that are all distant from each other. This is expected when the dataset contains multiple STs.

### Off-diagonal pale blocks

A pale/white square block means a group of isolates that are all very close to each other — a cluster. Cross-reference these clusters with the ST strip in Fig 2 to see if they correspond to a single ST, and with the AMR heatmap to see if they share resistance genes.

### Outbreak investigation thresholds

These are rough guides; appropriate thresholds depend on the organism and context:

| SNP distance | Typical interpretation |
|---|---|
| 0–5 | Essentially identical; same strain, potentially outbreak-linked |
| 5–20 | Probable outbreak cluster (apply species-specific guidance) |
| 20–50 | Recent common ancestor; short-term transmission possible |
| 50–500 | Same ST, unrelated; no transmission link implied |
| > 500 | Different lineages or STs |

**For *S.* Typhi**: thresholds are tighter; ≤ 3 SNPs is the standard cut-off for outbreak clusters in some WHO guidance.

**For *Shigella sonnei***: within-clade isolates from different outbreaks can share 5–15 SNPs; use the Mykrobe sub-clade assignment alongside the SNP distance.

---

## Source data

The heatmap is computed from the SKA2 pairwise SNP distance matrix (`ska2_snp_matrix.tsv`), which contains raw SNP counts across the whole-genome split k-mer alignment. See [SNP Phylogenetics](SNP-Phylogenetics.md) for how this matrix is produced.

---

## Output file

`{species}_snp_heatmap.pdf` and `.png`

For very large datasets (n > 200), the PNG may be large. The PDF provides vector graphics for zoom without pixelation.
