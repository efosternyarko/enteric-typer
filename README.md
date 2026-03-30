# enteric-typer

A species-gated genotyping workflow for enteric pathogens. Given a folder of
genome assemblies, the pipeline speciates each sample and deploys the
appropriate species-specific typing tools, generates a core-SNP phylogeny, and
(optionally) uploads results to Pathogenwatch and Microreact.

## Supported species

| Species | Typing tools | Pathogenwatch |
|---|---|---|
| *Escherichia coli* | MLST (Achtman), AMRFinder, ECTyper (serotype), PlasmidFinder, Kaptive (K-locus) | Yes |
| *Salmonella enterica* | MLST, AMRFinder, SISTR (serovar), PlasmidFinder | Yes |
| Other / unclassified | Species ID only (logged and skipped) | — |

## Workflow overview

```
Input assemblies (folder or samplesheet)
        │
        ▼
┌────────────────────┐
│  1. SPECIES CHECK  │  Mash distance against 7-species reference sketch
└────────────────────┘
        │
   ┌────┴────┐
   ▼         ▼
E. coli   Salmonella   (other species logged and skipped)
   │         │
   ▼         ▼
┌────────────────────────────────────────────────────────────┐
│  2. SPECIES-SPECIFIC TYPING  (all tools run in parallel)   │
│                                                            │
│  E. coli                  Salmonella                       │
│  ─────────────────────    ───────────────────────          │
│  MLST (achtman_4)         MLST (salmonella)                │
│  AMRFinder                AMRFinder                        │
│  ECTyper (O:H serotype)   SISTR (serovar)                  │
│  PlasmidFinder            PlasmidFinder                    │
│  Kaptive K-locus (G2/G3                                    │
│    → G1/G4 on untypeables)                                 │
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
        │
        ▼
┌───────────────────────────────────────────────────────────────────┐
│  7. SUMMARY PLOTS  (always runs — no API keys required)           │
│  Fig 1: ST distribution · serotypes · AMR prevalence · MDR       │
│  Fig 2: Core-SNP tree + phylogroup strip + AMR heatmap           │
│  Fig 3: Top AMR genes (intrinsic genes excluded via AMRrules)     │
│  Fig 4: Plasmid replicon types                                    │
└───────────────────────────────────────────────────────────────────┘
```

---

## Installation

### Step 1 — Clone the repository

```bash
git clone https://github.com/efosternyarko/enteric-typer
cd enteric-typer
```

---

### Step 2 — Install Java

Nextflow requires Java 17 or later (Java 21 recommended).

**macOS**
```bash
brew install --cask temurin
```
> If Homebrew is not yet installed, run:
> ```bash
> /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
> ```

**Linux (Debian/Ubuntu)**
```bash
sudo apt update && sudo apt install -y default-jdk
```

**Windows (PowerShell)**
```powershell
winget install EclipseAdoptium.Temurin.21.JDK   # Java 21 (recommended)
# or
winget install EclipseAdoptium.Temurin.25.JDK   # Java 25 (latest)
```

Verify: `java -version`

---

### Step 3 — Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
```

The installer creates a `nextflow` executable in the current directory. Move it
somewhere on your `$PATH` so it can be run from anywhere:

**macOS**
```bash
# If Homebrew is installed, move to the Homebrew bin (already in $PATH):
mv nextflow /opt/homebrew/bin/        # Apple Silicon (M1/M2/M3)
# or
mv nextflow /usr/local/bin/           # Intel Mac

# If Homebrew is not installed, create ~/bin and add it to your PATH:
mkdir -p ~/bin
mv nextflow ~/bin/
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

**Linux**
```bash
mv nextflow ~/bin/
# If ~/bin/ does not exist yet:
mkdir -p ~/bin
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

**Windows (PowerShell)**
```powershell
mkdir "$HOME\bin"
Move-Item nextflow.exe "$HOME\bin\"
[Environment]::SetEnvironmentVariable("Path", $env:Path + ";$HOME\bin", "User")
```
Then restart PowerShell.

Verify: `nextflow -version`

---

### Step 4 — Install conda + mamba

If conda is not already installed, the quickest route is **Miniforge**, which
bundles both conda and mamba together:

**macOS (Apple Silicon / arm64)**
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
bash Miniforge3-MacOSX-arm64.sh
```

**macOS (Intel / x86_64)**
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh"
bash Miniforge3-MacOSX-x86_64.sh
```

**Linux (x86_64)**
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh
```

When the installer asks:
```
Do you wish the installer to initialize Miniforge3
by running conda init? [yes|no]
```
Answer **yes**. Then reload your shell:

```bash
source ~/.zshrc    # macOS (zsh)
# or
source ~/.bashrc   # Linux (bash)
```

**Windows**
Download the [Miniforge3 Windows installer](https://github.com/conda-forge/miniforge/releases/latest)
and run it. When prompted, select *Add Miniforge to my PATH*. Then reinitialise:
```powershell
conda init powershell
```
Restart PowerShell.

> **Already have conda but not mamba?**
> ```bash
> conda install -n base -c conda-forge mamba
> ```

---

### Step 5 — Install ncbi-datasets-cli and mash

Both are required by `build_references.sh` to download reference genomes and build the species sketch:

```bash
mamba install -c conda-forge ncbi-datasets-cli
mamba install -c bioconda mash
```

---

### Step 6 — Build the Mash reference sketch

This downloads 7 reference genomes from NCBI and builds the Mash sketch used
for species identification. Run once before the first pipeline execution:

```bash
bash assets/build_references.sh
```

Output: `assets/enteric_species_refs.msh`
Expected runtime: 1–5 minutes depending on bandwidth.

---

## Quick start

### For a given folder of assemblies

```bash
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda
```

### Or to use a samplesheet

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

### Full run with Pathogenwatch + Microreact

```bash
# Set API keys once
nextflow secrets set PW_API_KEY       <your_pathogenwatch_key>
nextflow secrets set MICROREACT_TOKEN <your_microreact_token>
```

```bash
# Run the workflow
nextflow run main.nf \
    --input_dir            /path/to/assemblies/ \
    --outdir               results/ \
    --run_pathogenwatch \
    --upload_microreact \
    --microreact_project   "My enteric outbreak" \
    -profile conda
```

---

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

---

## Output files

```
results/
├── pipeline_info/
│   ├── timeline.html
│   ├── report.html
│   └── dag.svg
│
├── species_check/
│   └── *_species.txt              ← best species + Mash distance per sample
│
├── mlst_ecoli/
│   └── *_ecoli_achtman_4_mlst.tsv
├── amrfinder_ecoli/
│   └── *_amrfinder.tsv
├── ectyper/
│   └── *_ectyper.tsv
├── kaptive/
│   └── *_ktype.tsv                ← K-locus group, locus, type, confidence
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
│   ├── iqtree.treefile             ← Newick ML tree
│   └── iqtree.iqtree               ← IQ-TREE log
├── ska2_salmonella/
│   └── ska2_alignment.fasta
├── iqtree_salmonella/
│   ├── iqtree.treefile
│   └── iqtree.iqtree
│
├── pathogenwatch/                  ← (if --run_pathogenwatch)
│   ├── ecoli_pathogenwatch_samples.tsv
│   ├── ecoli_pathogenwatch_collection.json
│   ├── ecoli_pathogenwatch_tree.nwk
│   ├── salmonella_pathogenwatch_samples.tsv
│   └── salmonella_pathogenwatch_collection.json
│
├── ecoli_typer_results.tsv         ← Master results table (E. coli)
├── salmonella_typer_results.tsv    ← Master results table (Salmonella)
│
├── microreact_url_*.txt            ← (if --upload_microreact) Microreact URLs
│
├── ecoli_fig1_population_summary.{pdf,png}
├── ecoli_fig3_amr_genes.{pdf,png}
├── ecoli_fig4_plasmid_replicons.{pdf,png}
├── ecoli_tree_amr.{pdf,png}        ← Phylogeny + AMR heatmap
├── salmonella_fig1_population_summary.{pdf,png}
├── salmonella_fig3_amr_genes.{pdf,png}
├── salmonella_fig4_plasmid_replicons.{pdf,png}
└── salmonella_tree_amr.{pdf,png}
```

---

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

---

## HPC example

```bash
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile singularity,slurm \
    -c custom.config   # optional: override queue, project account, etc.
```

---

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

---

## Citation

If you use enteric-typer, please cite the tools it wraps:

- **Mash**: Ondov et al. (2016) Genome Biology 17:132
- **MLST**: Seemann (2016) github.com/tseemann/mlst
- **AMRFinder Plus**: Feldgarden et al. (2021) Scientific Reports 11:12728
- **ECTyper**: Laing et al. (2019) Microbial Genomics 5(12)
- **SISTR**: Yoshida et al. (2016) PLOS ONE 11(1):e0147101
- **PlasmidFinder**: Carattoli et al. (2014) Antimicrobial Agents and Chemotherapy
- **Kaptive**: Wyres et al. (2016) Microbial Genomics 2(10); Lam et al. (2022) Nature Protocols
- **SKA2**: github.com/bacpop/ska.rust
- **IQ-TREE 2**: Minh et al. (2020) Molecular Biology and Evolution 37(5)
- **Pathogenwatch**: Argimón et al. (2021) Nature Communications 12:2732
- **Microreact**: Argimón et al. (2016) Microbial Genomics 2(11)
- **AMRrules**: github.com/AMRverse/AMRrules
