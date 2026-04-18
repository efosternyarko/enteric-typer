# FAQ and Troubleshooting

---

## Running the pipeline

### My run completed but I see no results for some samples. What happened?

Check the pipeline log and the `results/species_check/` directory. Common causes:

1. **Assembly failed QC** — size outside the expected range for the species. Check `results/assembly_qc/` for the fail log and the genome size of the offending sample.
2. **Species not recognised** — Mash distance > 0.05 to all references. The sample may be a different organism entirely, a very poor-quality assembly, or a contaminated sample.
3. **Kraken2 contamination** — if `--kraken2_db` was provided, the sample failed the secondary-species threshold. Reduce `--max_contamination` threshold or inspect the assembly.

### The pipeline runs to completion but SKA2/IQ-TREE outputs are missing

Phylogenetics requires ≥ 3 samples per species (`--ska2_min_samples`). If only 1 or 2 assemblies of a given species are present, phylogenetics is silently skipped. To force a run with 2 samples, set `--ska2_min_samples 2` (not recommended — trees with 2 tips are uninformative).

### I get an error about missing conda environments

Run the pipeline once with `-profile conda` to allow Nextflow to build the environments. First-time conda setup can take 15–30 minutes. Alternatively, use `-profile docker` if Docker is available.

### I see "No ST predicted" for many samples in MLST

This usually means the assembly is too fragmented for the MLST allele database search to find all 7 loci. Assemblies with N50 < 10 kb or > 500 contigs commonly fail MLST. If this affects many samples, check the assembly metrics figure for the batch.

---

## Species identification

### An assembly I know is *E. coli* is being classified as *Shigella*

*Shigella* and *E. coli* are genomically very close. If the assembly is borderline, the nearest reference may flip. Check:
- The Mash distance output in `results/species_check/`
- Whether the assembly carries ipaH or pINV-like sequences (true *Shigella* would)

If it is genuinely *E. coli* with atypically high *Shigella* similarity (e.g. EIEC), the typing will still proceed through the *Shigella* pipeline, which will produce `shigeifinder_ipaH = pos` for EIEC. This is biologically correct — EIEC and *Shigella* share the pINV plasmid.

### *Klebsiella* samples appear in my batch but are excluded

Correct behaviour. *K. pneumoniae* is detected by the Mash sketch and logged to `results/species_check/` but not typed. To type *Klebsiella*, use [Kleborate](https://github.com/klebgenomics/Kleborate) directly or a dedicated Klebsiella pipeline.

---

## AMR results

### Why does AMRFinder report different genes to what I expect from the phenotype?

AMRFinder screens for *known* resistance mechanisms that are in its database. Novel alleles, unusual mechanisms, or genes not yet curated may be missed. Phenotypic resistance without a genotypic explanation warrants:
1. Manual inspection of the AMRFinder output TSV in `results/amrfinder/`
2. Checking whether the gene may be present but below the default identity/coverage thresholds
3. Considering regulatory mutations not covered by the current AMRFinder database

### Why do some genes appear as "intrinsic" even though they are resistance genes?

AMRrules classifies certain chromosomally-encoded genes as intrinsic because they are present in virtually all wild-type isolates of that species and do not represent *acquired* resistance. Examples include *blaEC* (chromosomal AmpC in *E. coli*) and constitutive efflux pumps (*acrB*, *emrD*). These are present regardless of selective pressure and are not useful resistance markers. They are retained in the `amrfinder_intrinsic_genes` column for reference.

If you believe a gene classified as intrinsic in your dataset is actually meaningful (e.g. because AMRrules rules don't apply to your species context), you can use the `amrfinder_genes` (unfiltered) column for custom analysis.

### EFFLUX class is not shown in any figures

Efflux-mediated resistance is excluded from all AMR figures and MDR counts. Constitutive efflux pump overexpression contributes to resistance in wild-type backgrounds and its presence is not a marker of horizontally acquired resistance. To see efflux genes, use the raw `amrfinder_genes` column.

---

## Plasmid results

### There are no or very few AMR genes shown in the Fig 4 bubble matrix for my E. coli / Salmonella dataset. Is something wrong?

Not necessarily. If your assemblies are short-read (Illumina), plasmid fragmentation means that replicon sequences and AMR cassettes typically land on separate contigs even when they are on the same plasmid. The same-contig co-location criterion used by the pipeline correctly reflects this limitation — it does not fabricate links that don't exist in the assembly.

If the same-contig fraction is below 10%, the pipeline automatically switches to the simplified 2-panel layout (replicon prevalence by Inc family + tree heatmap), which is appropriate for short-read data. See [Fig 4 — Plasmid Overview](Fig4-Plasmid-Overview.md) for full explanation.

For reliable plasmid-AMR attribution, use long-read or hybrid assemblies.

### A replicon is at 100% prevalence in my dataset

This is real. Common causes:
- **Col-type small cryptic plasmids**: Col156 is found in nearly 100% of *S. sonnei* ST152 pandemic isolates.
- **Large virulence plasmid**: pINV in *Shigella* would appear at 100% if the dataset is all pathogenic isolates (but pINV replicon is not reliably detected by PlasmidFinder; use the pINV screen column instead).
- **Chromosome-integrated plasmid sequence**: rare but possible; the replicon sequence has integrated into the chromosome. Confirm by checking whether it co-localises with chromosomal genes in the assembly.

---

## Figures

### Fig 2 phylogenetic tree is not produced

Requires ≥ 3 samples with a passed phylogenetic analysis. Check whether `--skip_local_phylo` was set, whether too few samples passed QC, and whether the SKA2/IQ-TREE process completed in `results/`.

### Fig 4 shows the simplified layout when I expected the full layout

The simplified layout is triggered automatically when < 10% of plasmid-carrying isolates have same-contig AMR co-occurrence. For short-read assemblies this is expected. For long-read assemblies it may indicate that AMRFinder found no acquired genes at all (check `amrfinder_acquired_genes`), or that plasmid fragmentation is worse than expected (check N50 in the assembly metrics figure).

### Some samples are missing from the figures

Samples that fail assembly QC (size or contamination) are excluded from all figures. Samples that pass QC but have no data for a specific tool (e.g. MLST failed) will appear in the tree but show as "Unknown" in ST strips.

---

## Performance

### The pipeline is very slow on my laptop

- Set `--max_cpus 4` to reduce per-process CPU requests
- Set `--skip_local_phylo` to skip the computationally expensive SKA2 + IQ-TREE steps
- Use `-profile docker` if Docker is faster on your system than conda

### IQ-TREE is slow with many samples

IQ-TREE scales roughly O(n²) with sample count for bootstrap calculations. For > 200 samples, consider:
- Reducing `--iqtree_bootstraps` (minimum 1000 for UFBoot — cannot go lower)
- Using `-profile docker` for better CPU utilisation
- Running on a cluster with `--max_cpus 16` or more
