# SNP Phylogenetics

enteric-typer builds a whole-genome SNP phylogeny for each species when ≥ 3 samples are present. The phylogeny is used to annotate Fig 2 and to order samples in Fig 4 (plasmid heatmap) and the SNP distance heatmap.

Two tools are run sequentially: **SKA2** to compute a SNP alignment and pairwise distance matrix, then **IQ-TREE 2** to infer a maximum-likelihood phylogenetic tree.

---

## Step 1 — SKA2: split k-mer alignment

[SKA2](https://github.com/bacpop/ska.rust) (split k-mer analysis) aligns genome assemblies **without a reference genome**. Instead of mapping reads or contigs to a fixed reference, it represents every assembly as a set of split k-mers (two flanking sequences of length *k* with a variable central base) and then identifies positions where the central base differs between samples.

### Why no reference genome?

- Reference-free alignment avoids reference bias: positions absent from the chosen reference are not excluded.
- For a species-gated pipeline covering multiple serovars or phylogroups, no single reference adequately represents all input diversity.
- SKA2 is fast enough to run genome-wide on a typical workstation even for datasets of hundreds of samples.

### Parameters

| Parameter | Value | Rationale |
|---|---|---|
| k-mer length (`-k`) | 31 | Balances sensitivity and specificity; 31 is the SKA2 default and performs well for bacterial genomes |
| Minimum frequency | default | Retains all variant positions seen in ≥ 2 samples |

### Outputs

- `ska2_alignment.fasta` — multi-sample SNP alignment (one sequence per sample, gap = not called at that position)
- `snp_matrix.tsv` — pairwise SNP distance matrix (raw count of differing positions between each pair)

---

## Step 2 — IQ-TREE 2: maximum-likelihood tree

[IQ-TREE 2](http://www.iqtree.org/) infers a maximum-likelihood phylogeny from the SKA2 SNP alignment.

### Default model: GTR+G

The generalised time-reversible model with gamma-distributed rate variation (GTR+G) is used by default. GTR is appropriate because:

- It is the most parameter-rich standard substitution model, making no assumption that substitution rates are equal across base pairs.
- Gamma rate variation accounts for the fact that some SNP positions evolve faster than others (e.g. synonymous sites vs. regulatory regions).
- For whole-genome SNP data from bacterial genomes, GTR+G gives robust topology estimates in most settings.

To use automatic model selection (ModelFinder Plus), pass `--iqtree_model MFP`. This tests a large set of substitution models and selects the one with the best BIC score. It is slower but more rigorous for publication-quality analyses.

### Bootstrap support

By default, **100 ultrafast bootstrap replicates** (UFBoot) are computed. UFBoot values ≥ 95 are considered well-supported in IQ-TREE's framework (UFBoot values are not directly comparable to standard bootstrap — 95 UFBoot ≈ 70–80 standard bootstrap in terms of topological confidence).

### Midpoint rooting

IQ-TREE produces an unrooted tree. enteric-typer roots it at the midpoint (the point on the longest branch that makes both resulting subtrees as equal in total branch length as possible) before display. Midpoint rooting is the standard approach when no outgroup is available, and is appropriate for within-species datasets where all samples have similar evolutionary distance from a common ancestor.

---

## Interpreting branch lengths and SNP distances

Branch lengths in the tree represent **substitutions per site** in the SNP alignment. Because the alignment contains only variable positions (not the full genome), branch lengths are expressed relative to the number of variant sites detected, not the full chromosome length.

The **SNP distance matrix** gives raw pairwise distances in number of differing positions across the whole-genome SKA2 alignment. These values can be used for:

- **Cluster definition**: isolates within ≤ 10–20 SNPs are considered potentially linked in outbreak investigations (the threshold depends on the organism and study context)
- **Quality checking**: unexpectedly large distances between phenotypically similar isolates may indicate mislabelling or contamination

### Typical SNP distance ranges (within-species)

| Comparison | Expected SNP range |
|---|---|
| Same patient, same admission | 0–5 SNPs |
| Putative outbreak cluster | 5–20 SNPs |
| Same ST, unrelated | 50–500 SNPs |
| Different STs, same species | 500–50,000 SNPs |
| Different *Shigella* species | 20,000–100,000 SNPs |

These are rough guides. Ranges vary substantially by species, ST, and the evolutionary rate of the lineage in question.

---

## When phylogenetics is skipped

Phylogenetics runs per species when ≥ 3 samples of that species are present. It is automatically skipped when:

- Fewer than 3 samples are assigned to a species (`--ska2_min_samples`, default 3)
- `--skip_local_phylo` is passed

When skipped, Figs 2, 4 (tree panel), and the SNP heatmap are not produced. All other figures are unaffected.
