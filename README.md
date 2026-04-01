# enteric-typer

A species-gated genotyping workflow for enteric pathogens. Given a folder of
genome assemblies, the pipeline speciates each sample and deploys the
appropriate species-specific typing tools, generates a whole-genome SNP phylogeny, and
produces publication-ready summary figures.

## Supported species

| Species | Typing tools |
|---|---|
| *Escherichia coli* | MLST (Achtman), AMRFinder, ECTyper (serotype), EzClermont (Clermont phylogroup), Kleborate (pathotype), PlasmidFinder, Kaptive (K-locus) |
| *Salmonella enterica* | MLST, AMRFinder, SISTR (serovar), PlasmidFinder, Abricate VFDB |
| Other / unclassified | Species ID only (logged and skipped) |

> **Kleborate:** runs full pathotype detection on all platforms. On macOS
> Apple Silicon (ARM64) without Rosetta 2, it automatically falls back to
> MLST-only mode — all other tools run at full capability regardless of
> platform (Linux, macOS Intel/ARM64, HPC).

> **AMR gene classification:** all AMRFinder hits are classified by
> [AMRrules](https://github.com/AMRverse/AMRrules) into *acquired* resistance
> genes and *intrinsic* (species-wildtype) genes. Intrinsic genes are retained
> in the results TSV for reference but are **excluded from all AMR plots** so
> that figures reflect only clinically relevant acquired resistance.

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
│  E. coli                       Salmonella                  │
│  ──────────────────────────    ───────────────────────     │
│  MLST (achtman_4)              MLST (senterica_achtman_2)  │
│  AMRFinder                     AMRFinder                   │
│  ECTyper (O:H serotype)        SISTR (serovar)             │
│  EzClermont (phylogroup)       PlasmidFinder               │
│  Kleborate (pathotype)         Abricate VFDB               │
│  PlasmidFinder                                             │
│  Kaptive K-locus (G2/G3                                    │
│    → G1/G4 on untypeables)                                 │
└────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────┐
│  3. PHYLOGENETICS  (per species, ≥ 3 samples; skip with         │
│                    --skip_local_phylo)                          │
│  SKA2 build (k=31) → whole-genome SNP alignment + SNP distance matrix  │
│  → IQ-TREE ML tree (ModelFinder Plus automatic model selection) │
└─────────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────┐
│  4. AGGREGATE  (one TSV per species)                         │
│  ecoli_typer_results.tsv                                     │
│  salmonella_typer_results.tsv                                │
│  AMRFinder hits classified by AMRrules into:                 │
│    amrfinder_acquired_genes  — clinically relevant acquired  │
│                                resistance genes              │
│    amrfinder_intrinsic_genes — species-intrinsic genes       │
│                                (flagged, retained in TSV)    │
└──────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────────────┐
│  5. SUMMARY PLOTS                                                    │
│  Intrinsic resistance genes (classified by AMRrules) are excluded   │
│  from all AMR figures — only acquired genes are plotted.            │
│                                                                      │
│  Fig 1: ST distribution · serotypes · AMR drug classes · MDR        │
│  Fig 2: Core-SNP tree + ST/phylogroup strips + virulence + AMR panel │
│  Fig 3: Top acquired AMR genes (intrinsic genes excluded)            │
│  Fig 4: Plasmid replicon types                                       │
│  Fig 5: Virulence genes / pathotype (E. coli) or VFDB (Salmonella)  │
│  Fig 6: Pairwise whole-genome SNP distance heatmap                   │
└──────────────────────────────────────────────────────────────────────┘
```

> **E. coli K-locus typing (Kaptive):** K-locus typing uses two curated databases run sequentially.
> Assemblies are first typed against the **G2/G3 database** (`EC-K-typing_group2and3_v3.0.0.gbk`;
> [Gladstone et al. 2026](https://www.nature.com/articles/s41564-026-02283-w)).
> Samples that remain untypeable are then re-typed against the **G1/G4 database**
> (`EC-K-typing_group1and4_v1.2.gbk`;
> [Foster-Nyarko et al.](https://github.com/efosternyarko/EC-K-typing-G1G4)).
> G2/G3 and G1/G4 loci are mutually exclusive in *E. coli*, so sequential typing
> ensures each sample is assigned to the correct group without double-counting.

> **Phylogenetics (SKA2 + IQ-TREE):** whole-genome SNP phylogenies are built using
> [SKA2](https://github.com/bacpop/ska.rust) (split k-mer alignment, k=31),
> which aligns assemblies without a reference genome and produces a genome-wide SNP
> alignment and pairwise SNP distance matrix. The alignment is then passed to
> IQ-TREE 2 for maximum-likelihood tree inference with automatic model
> selection (ModelFinder Plus). Phylogenetics runs per species when ≥ 3
> samples are present; skip with `--skip_local_phylo` for faster runs.

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

> **Linux / Ubuntu / Debian:** use the Linux commands below.
> **Windows:** we recommend running enteric-typer inside
> [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) (Windows
> Subsystem for Linux), which gives a full Ubuntu environment. Install WSL2,
> then follow the Linux instructions. Alternatively use the `docker` profile
> with Docker Desktop for Windows.

### For a given folder of assemblies

```bash
# Linux / Ubuntu / Debian / Intel Mac
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda

# Apple Silicon Mac (M1/M2/M3/M4) — prefix CONDA_SUBDIR and add arm64 profile
CONDA_SUBDIR=osx-64 nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda,arm64

# Windows (Docker Desktop) — no conda required
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile docker
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

### Skip local phylogenetics (faster, no SKA2/IQ-TREE)

```bash
nextflow run main.nf \
    --input_dir        /path/to/assemblies/ \
    --outdir           results/ \
    --skip_local_phylo \
    -profile conda
```

> When `--skip_local_phylo` is set, the SNP distance matrix, SNP heatmap, and
> tree annotation figures are not produced. Microreact upload (if enabled) will
> not include a tree.

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_dir` | `null` | Folder of FASTA assemblies (`.fasta/.fa/.fna/.fas`) |
| `--samplesheet` | `null` | CSV with `id,fasta` columns |
| `--outdir` | `results` | Output directory |
| `--skip_local_phylo` | `false` | Skip SKA2 + IQ-TREE (no tree, no SNP matrix/heatmap) |
| `--ska2_min_samples` | `3` | Minimum samples to attempt SKA2/IQ-TREE |
| `--iqtree_model` | `MFP` | IQ-TREE substitution model (`MFP` = ModelFinder Plus automatic selection) |
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
│   └── *_ecoli_achtman_4_mlst.tsv  ← scheme ID in mlst tool database (7-gene Achtman scheme)
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
│   ├── ska2_alignment.fasta       ← whole-genome SNP alignment (input to IQ-TREE)
│   └── snp_matrix.tsv             ← pairwise whole-genome SNP distance matrix
├── iqtree_ecoli/
│   ├── iqtree.treefile             ← Newick ML tree
│   └── iqtree.iqtree               ← IQ-TREE log + best-fit model
├── ska2_salmonella/
│   ├── ska2_alignment.fasta
│   └── snp_matrix.tsv
├── iqtree_salmonella/
│   ├── iqtree.treefile
│   └── iqtree.iqtree
│
├── ecoli_typer_results.tsv         ← Master results table (E. coli)
├── salmonella_typer_results.tsv    ← Master results table (Salmonella)
│
│   Both tables include AMRFinder columns:
│     amrfinder_acquired_genes   — resistance genes (intrinsic excluded)
│     amrfinder_intrinsic_genes  — wildtype/intrinsic genes (flagged, not excluded)
│     amrfinder_genes            — all raw AMRFinder hits
│
├── {species}_fig1_population_summary.{pdf,png}   ← ST · serotype · AMR drug classes · MDR
├── {species}_fig3_amr_genes.{pdf,png}            ← Acquired AMR gene frequencies
├── {species}_fig4_plasmid_replicons.{pdf,png}    ← Plasmid replicon types
├── {species}_fig5_virulence.{pdf,png}            ← Pathotype / virulence genes
├── {species}_tree_amr.{pdf,png}                  ← Phylogeny + phylogroup + AMR heatmaps
└── {species}_snp_heatmap.{pdf,png}               ← Pairwise whole-genome SNP distance heatmap
```

---

## Cleaning up after a run

All final outputs are written to `--outdir` (default: `results/`). Once you are
happy with the results you can free up disk space by removing Nextflow's
temporary files:

```bash
rm -rf work/ .nextflow/ .nextflow.log*
```

> **Tip:** keep `work/` if you want to use `-resume` to rerun with different
> parameters without repeating completed steps. Delete it only when the run is
> finalised.

---

## Optional: install Graphviz

Nextflow generates an execution DAG diagram (`pipeline_info/dag.svg`) which
requires Graphviz to render. Without it you will see a harmless warning — the
pipeline still runs fine. To suppress the warning:

```bash
# macOS
brew install graphviz

# Linux
sudo apt install graphviz

# conda
conda install -c conda-forge graphviz
```

## Execution profiles

| Profile | Use case |
|---|---|
| `conda` | Local workstation with conda |
| `mamba` | Same as conda but faster env solving |
| `arm64` | **Add on Apple Silicon (M1/M2/M3/M4)** — forces osx-64 conda envs via Rosetta 2 |
| `docker` | Local with Docker Desktop |
| `singularity` | HPC cluster with Singularity/Apptainer |
| `slurm` | SLURM HPC executor (combine with singularity: `-profile singularity,slurm`) |
| `pbs` | PBS/Torque HPC executor |
| `test` | Quick test run (no uploads) |

### macOS Apple Silicon (M1/M2/M3/M4)

Some Bioconda packages (e.g. `abricate`) have no native arm64 build. Add the `arm64`
profile to force Rosetta 2 emulation — tools run at near-native speed and the pipeline
produces identical results to Linux. Kleborate will automatically use MLST-only mode
on ARM64 if its full pathotype dependencies are unavailable, but all other tools run
at full capability:

```bash
CONDA_SUBDIR=osx-64 nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda,arm64
```

> **`CONDA_SUBDIR=osx-64` is required** on Apple Silicon. Modern versions of
> libmamba ignore the `--platform` flag inside Nextflow and fall back to native
> `osx-arm64`, where some Bioconda packages are unavailable. Prefixing
> `CONDA_SUBDIR=osx-64` forces Rosetta 2 emulation reliably. You only pay the
> environment-build cost once; environments are cached and reused on subsequent
> runs.

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

If you use enteric-typer or figures generated by it, please cite:

- **enteric-typer**: Foster-Nyarko E et al. github.com/efosternyarko/enteric-typer

And the tools it wraps:

- **Mash**: Ondov et al. (2016) Genome Biology 17:132
- **MLST**: Seemann (2016) github.com/tseemann/mlst
- **AMRFinder Plus**: Feldgarden et al. (2021) Scientific Reports 11:12728
- **ECTyper**: Laing et al. (2019) Microbial Genomics 5(12)
- **EzClermont**: Waters et al. (2020) Microbial Genomics 6(9)
- **SISTR**: Yoshida et al. (2016) PLOS ONE 11(1):e0147101
- **PlasmidFinder**: Carattoli et al. (2014) Antimicrobial Agents and Chemotherapy
- **Kaptive**: Wyres et al. (2016) Microbial Genomics 2(10); Lam et al. (2022) Nature Protocols
- **EC K-typing G2/G3**: Gladstone et al. (2026) Nature Microbiology; github.com/rgladstone/EC-K-typing
- **EC K-typing G1/G4**: Foster-Nyarko E et al. github.com/efosternyarko/EC-K-typing-G1G4
- **Kleborate**: Lam et al. (2021) Nature Communications; github.com/klebgenomics/Kleborate
- **Abricate**: github.com/tseemann/abricate
- **SKA2**: github.com/bacpop/ska.rust
- **IQ-TREE 2**: Minh et al. (2020) Molecular Biology and Evolution 37(5)
- **AMRrules**: github.com/AMRverse/AMRrules
