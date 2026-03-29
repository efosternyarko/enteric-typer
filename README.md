# enteric-typer

A species-gated genotyping workflow for enteric pathogens. Given a folder of
genome assemblies, the pipeline speciates each sample and deploys the
appropriate species-specific typing tools, generates a core-SNP phylogeny, and
(optionally) uploads results to Pathogenwatch and Microreact.

## Supported species

| Species | Typing tools | Pathogenwatch |
|---|---|---|
| *Escherichia coli* | MLST (Achtman), AMRFinder, ECTyper (serotype), PlasmidFinder | Yes |
| *Salmonella enterica* | MLST, AMRFinder, SISTR (serovar), PlasmidFinder | Yes |
| Other / unclassified | Species ID only (logged and skipped) | — |

## Workflow overview

```
Input assemblies (folder or samplesheet)
        │
        ▼
┌───────────────────┐
│  1. SPECIES CHECK  │  Mash distance against 7-species reference sketch
└───────────────────┘
        │
   ┌────┴────┐
   ▼         ▼
E. coli   Salmonella   (other species logged and skipped)
   │         │
   ▼         ▼
┌────────────────────────────────────────────────────────────┐
│  2. SPECIES-SPECIFIC TYPING  (all tools run in parallel)   │
│                                                            │
│  E. coli              Salmonella                           │
│  ─────────────────    ──────────────────────               │
│  MLST (achtman_4)     MLST (salmonella)                    │
│  AMRFinder            AMRFinder                            │
│  ECTyper (O:H)        SISTR (serovar)                      │
│  PlasmidFinder        PlasmidFinder                        │
└────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────┐
│  3. PHYLOGENETICS  (per species, ≥ 3 samples needed)  │
│  SKA2 build → core-SNP alignment → IQ-TREE ML tree   │
└──────────────────────────────────────────────────────┘
        │
        ▼  (if --run_pathogenwatch)
┌───────────────────────────────────────────────────────────────┐
│  4. PATHOGENWATCH                                             │
│  Upload → wait → create collection → download tree           │
│  → cgMLST cluster search at thresholds 5, 10, 20, 50 alleles │
└───────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────┐
│  5. AGGREGATE  (one TSV per species)                  │
│  ecoli_typer_results.tsv                              │
│  salmonella_typer_results.tsv                         │
└──────────────────────────────────────────────────────┘
        │
        ▼  (if --upload_microreact)
┌──────────────────────────────────────────────────────┐
│  6. MICROREACT                                        │
│  Results TSV + IQ-TREE Newick → Microreact project   │
└──────────────────────────────────────────────────────┘
```

## Installation

### Prerequisites

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/)
  (or Docker / Singularity)

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/   # or wherever you keep executables

# Install mamba (optional but faster for conda env solving)
conda install -n base -c conda-forge mamba
```

### Build the Mash reference sketch (required before first run)

```bash
# Install ncbi-datasets-cli if needed
conda install -c conda-forge ncbi-datasets-cli

# Build the 7-species reference sketch (~1–5 min depending on bandwidth)
bash assets/build_references.sh
```

This downloads 7 reference genomes from NCBI and produces
`assets/enteric_species_refs.msh`.

## Quick start

### 1. From a folder of assemblies

```bash
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda
```

### 2. From a samplesheet (CSV: id,fasta)

```bash
# Auto-generate samplesheet from a folder
python bin/make_samplesheet.py \
    --input /path/to/assemblies/ \
    --output samples.csv

# Run with samplesheet
nextflow run main.nf \
    --samplesheet samples.csv \
    --outdir      results/ \
    -profile conda
```

### 3. Full run with Pathogenwatch + Microreact

```bash
# Set API keys once
nextflow secrets set PW_API_KEY       <your_pathogenwatch_key>
nextflow secrets set MICROREACT_TOKEN <your_microreact_token>

nextflow run main.nf \
    --input_dir            /path/to/assemblies/ \
    --outdir               results/ \
    --run_pathogenwatch \
    --upload_microreact \
    --microreact_project   "My enteric outbreak" \
    -profile conda
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_dir` | `null` | Folder of FASTA assemblies (`.fasta/.fa/.fna/.fas`) |
| `--samplesheet` | `null` | CSV with `id,fasta` columns |
| `--outdir` | `results` | Output directory |
| `--run_pathogenwatch` | `false` | Upload to Pathogenwatch (requires `PW_API_KEY`) |
| `--upload_microreact` | `false` | Create Microreact project (requires `MICROREACT_TOKEN`) |
| `--microreact_project` | `enteric-typer run` | Microreact project name prefix |
| `--pathogenwatch_cluster_thresholds` | `5,10,20,50` | cgMLST allele-difference thresholds |
| `--ska2_min_samples` | `3` | Minimum samples to attempt SKA2/IQ-TREE |
| `--ska2_prop_filter` | `0.95` | Core-genome proportion filter for SKA2 |
| `--iqtree_model` | `GTR+G` | IQ-TREE substitution model |
| `--iqtree_bootstraps` | `1000` | IQ-TREE ultrafast bootstrap replicates |

## Output files

```
results/
├── pipeline_info/
│   ├── timeline.html
│   ├── report.html
│   └── dag.svg
│
├── species_check/
│   └── *_species.txt          ← best species + Mash distance per sample
│
├── mlst_ecoli/
│   └── *_ecoli_achtman_4_mlst.tsv
├── amrfinder_ecoli/
│   └── *_amrfinder.tsv
├── ectyper/
│   └── *_ectyper.tsv
├── plasmidfinder_ecoli/
│   └── *_plasmidfinder.tsv
│
├── mlst_salmonella/
│   └── *_salmonella_mlst.tsv
├── amrfinder_salmonella/
│   └── *_amrfinder.tsv
├── sistr/
│   └── *_sistr.tsv
├── plasmidfinder_salmonella/
│   └── *_plasmidfinder.tsv
│
├── ska2_ecoli/
│   └── ska2_alignment.fasta
├── iqtree_ecoli/
│   ├── iqtree.treefile        ← Newick ML tree
│   └── iqtree.iqtree          ← IQ-TREE log
├── ska2_salmonella/
│   └── ska2_alignment.fasta
├── iqtree_salmonella/
│   ├── iqtree.treefile
│   └── iqtree.iqtree
│
├── pathogenwatch/             ← (if --run_pathogenwatch)
│   ├── ecoli_pathogenwatch_samples.tsv
│   ├── ecoli_pathogenwatch_collection.json
│   ├── ecoli_pathogenwatch_tree.nwk
│   ├── salmonella_pathogenwatch_samples.tsv
│   └── salmonella_pathogenwatch_collection.json
│
├── ecoli_typer_results.tsv    ← Master results table (E. coli)
├── salmonella_typer_results.tsv ← Master results table (Salmonella)
│
└── microreact_url_*.txt       ← (if --upload_microreact) Microreact project URLs
```

## Execution profiles

| Profile | Use case |
|---|---|
| `conda` | Local workstation with conda |
| `mamba` | Same as conda but faster env solving |
| `docker` | Local with Docker Desktop |
| `singularity` | HPC cluster with Singularity/Apptainer |
| `slurm` | SLURM HPC executor (combine with singularity: `-profile singularity,slurm`) |
| `pbs` | PBS/Torque HPC executor |
| `test` | Quick test run (no uploads) |

## HPC example (Monash M3 / MASSIVE)

```bash
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile singularity,slurm \
    -c custom.config   # optional: override queue, project account, etc.
```

## Updating the reference sketch

To add additional species or update reference genomes, edit
`assets/build_references.sh` and add the NCBI accession + label. Reference
filenames must start with a recognisable prefix:

| Prefix | Mapped species label |
|---|---|
| `Ecoli_` or `Escherichia_` | `E_coli` |
| `Salmonella_` | `Salmonella_enterica` |
| `Shigella_` | `Shigella` |
| `Klebsiella_` | `Klebsiella` |
| `Enterobacter_` | `Enterobacter` |

## Citation

If you use enteric-typer, please cite the tools it wraps:

- **Mash**: Ondov et al. (2016) Genome Biology 17:132
- **MLST**: Seemann (2016) github.com/tseemann/mlst
- **AMRFinder Plus**: Feldgarden et al. (2021) Scientific Reports 11:12728
- **ECTyper**: Laing et al. (2019) Microbial Genomics 5(12)
- **SISTR**: Yoshida et al. (2016) PLOS ONE 11(1):e0147101
- **PlasmidFinder**: Carattoli et al. (2014) Antimicrobial Agents and Chemotherapy
- **SKA2**: github.com/bacpop/ska.rust
- **IQ-TREE 2**: Minh et al. (2020) Molecular Biology and Evolution 37(5)
- **Pathogenwatch**: Argimón et al. (2021) Nature Communications 12:2732
- **Microreact**: Argimón et al. (2016) Microbial Genomics 2(11)
