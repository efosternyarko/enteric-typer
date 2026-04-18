# Parameters

All parameters can be passed on the command line as `--param_name value`. Boolean flags (like `--skip_local_phylo`) are set by simply including the flag (no value needed).

---

## Input / output

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | `null` | Path to a two-column TSV (no header): sample name, assembly path. Use this when assembly files are in different locations. |
| `--input_dir` | `null` | Path to a directory containing FASTA assemblies (`.fasta`, `.fa`, `.fna`, `.fas`, `.fsa`). All FASTA files in the directory are processed. |
| `--outdir` | `results` | Output directory. Created if it does not exist. |

**One of `--samplesheet` or `--input_dir` is required.**

---

## Assembly QC — genome size thresholds

Assemblies outside the size range for their species are excluded from all typing. Size is the sum of all contig lengths.

| Parameter | Default | Species |
|---|---|---|
| `--ecoli_min_length` | `4300000` | *E. coli* |
| `--ecoli_max_length` | `5900000` | *E. coli* |
| `--shigella_min_length` | `4300000` | *Shigella* |
| `--shigella_max_length` | `5900000` | *Shigella* |
| `--salmonella_min_length` | `4100000` | *Salmonella enterica* |
| `--salmonella_max_length` | `6600000` | *Salmonella enterica* |

See [Assembly QC](Assembly-QC.md) for why these thresholds are set where they are.

---

## Contamination screening (Kraken2)

| Parameter | Default | Description |
|---|---|---|
| `--kraken2_db` | `null` | Path to a Kraken2 database directory. When provided, each assembly is screened for secondary species contamination. Skipped when null. |
| `--max_contamination` | `3.0` | Maximum percentage of contigs classified as a secondary species before the assembly is flagged as contaminated and excluded. |

**Recommended database**: `k2_standard_08gb` (~8 GB, covers all major bacterial genera). Download from [https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2).

---

## Phylogenetics

| Parameter | Default | Description |
|---|---|---|
| `--ska2_min_samples` | `3` | Minimum number of samples per species required to run SKA2 and IQ-TREE. If fewer samples, phylogenetics is skipped for that species. |
| `--ska2_prop_filter` | `0.95` | Core-genome filter: only SNP positions present in ≥ this fraction of samples are included in the alignment. Positions absent from many samples (e.g. due to assembly gaps) are excluded to prevent inflated distances. |
| `--iqtree_model` | `GTR+G` | Substitution model for IQ-TREE 2. `GTR+G` is appropriate for most bacterial WGS SNP alignments. Use `MFP` to enable ModelFinder Plus (automatic model selection; slower but more rigorous). |
| `--iqtree_bootstraps` | `1000` | Number of ultrafast bootstrap replicates (`-B`). Must be ≥ 1000 for UFBoot; lower values are clamped to 1000. |
| `--skip_local_phylo` | `false` | Set to skip SKA2 and IQ-TREE entirely. When set, Figs 2, 4 (tree panel), and the SNP distance heatmap are not produced. |

---

## Resource limits

| Parameter | Default | Description |
|---|---|---|
| `--max_cpus` | `16` | Maximum CPU cores a single process can request. Set lower on shared machines (e.g. `--max_cpus 4`). |
| `--max_memory` | `128.GB` | Maximum memory a single process can request. |

---

## Full usage example

```bash
nextflow run main.nf \
  --input_dir /path/to/assemblies \
  --outdir    my_results \
  --max_cpus  8 \
  --kraken2_db /data/k2_standard_08gb \
  --skip_local_phylo
```

For a samplesheet:

```bash
nextflow run main.nf \
  --samplesheet samples.tsv \
  --outdir      my_results \
  --iqtree_model MFP
```

Where `samples.tsv` has no header and two columns: sample name and absolute path to FASTA.

---

## Nextflow execution options

These are Nextflow options, not pipeline parameters:

| Flag | Description |
|---|---|
| `-profile docker` | Run tools in Docker containers |
| `-profile conda` | Run tools in conda environments (default) |
| `-resume` | Resume a failed run from cached results |
| `-work-dir /path` | Change the Nextflow work directory (default: `./work`) |
| `-N email@domain` | Send an email on run completion |
