🌐 **English** | [Français](README.fr.md) | [Português](README.pt.md)

# enteric-typer

A species-gated genotyping workflow for enteric pathogens. Given a folder of
genome assemblies, the pipeline speciates each sample and deploys the
appropriate species-specific typing tools, generates a whole-genome SNP phylogeny, and
produces publication-ready summary figures.

## Supported species

| Species | Typing tools |
|---|---|
| *Escherichia coli* | MLST (Achtman), AMRFinder, ECTyper (serotype), EzClermont (Clermont phylogroup), Kleborate (pathotype), PlasmidFinder, Kaptive (K-locus) |
| *Salmonella enterica* | MLST, AMRFinder, SISTR (serovar), PlasmidFinder |
| *Shigella* spp. | MLST (Achtman), AMRFinder, ShigEiFinder (serotype/species), Mykrobe (S. sonnei genotyping), PlasmidFinder, pINV screen, IS element screen |
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
┌───────────────────────────────────────────────────────────────────────────┐
│  1. SPECIES & QC CHECK                                                    │
│                                                                           │
│  Species identification — Mash distance against 13-genome reference       │
│  sketch; closest-reference-wins (same approach as Kleborate). Shigella   │
│  isolates are correctly identified because genuine Shigella are closer    │
│  to Shigella references than to any E. coli reference.                    │
│                                                                           │
│  Assembly QC filters  (failures logged and excluded from typing)          │
│  ──────────────────────────────────────────────────────────────────────   │
│  Genome size    E. coli / Shigella  4.3 – 5.9 Mb  [1]                    │
│                 Salmonella          4.1 – 6.6 Mb  [2]                    │
│  Contamination  Kraken2 secondary species < 3 % of total contigs           │
│                 (optional — provide --kraken2_db; skipped if not set)     │
└───────────────────────────────────────────────────────────────────────────┘
        │
   ┌────┼────────┐
   ▼    ▼        ▼
E. coli  Salmonella  Shigella   (other species logged and skipped)
   │         │           │
   ▼         ▼           ▼
┌───────────────────────────────────────────────────────────────────────────┐
│  2. SPECIES-SPECIFIC TYPING  (all tools run in parallel)                  │
│                                                                           │
│  E. coli                       Salmonella          Shigella               │
│  ──────────────────────────    ──────────────────  ──────────────────     │
│  MLST (achtman_4)              MLST (senterica_    MLST (achtman_4)       │
│  AMRFinder                       achtman_2)        AMRFinder              │
│  ECTyper (O:H serotype)        AMRFinder           ShigEiFinder           │
│  EzClermont (phylogroup)       SISTR (serovar)       (serotype/species)   │
│  Kleborate (pathotype)         PlasmidFinder       Mykrobe                │
│  PlasmidFinder                                       (S. sonnei genotype) │
│  Kaptive K-locus (G2/G3                            PlasmidFinder          │
│    → G1/G4 on untypeables)                         pINV screen            │
│                                                    IS element screen      │
└───────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  3. PHYLOGENETICS  (per species, ≥ 3 samples; skip with                  │
│                    --skip_local_phylo)                                   │
│  SKA2 build (k=31) → whole-genome SNP alignment + SNP distance matrix   │
│  → IQ-TREE ML tree (GTR+G model; override with --iqtree_model MFP)       │
└──────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────┐
│  4. AGGREGATE  (one TSV per species)                         │
│  ecoli_typer_results.tsv                                     │
│  salmonella_typer_results.tsv                                │
│  shigella_typer_results.tsv                                  │
│  AMRFinder hits classified by AMRrules into:                 │
│    amrfinder_acquired_genes  — clinically relevant acquired  │
│                                resistance genes              │
│    amrfinder_intrinsic_genes — species-intrinsic genes       │
│                                (flagged, retained in TSV)    │
└──────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  5. SUMMARY PLOTS                                                        │
│  Intrinsic resistance genes (classified by AMRrules) are excluded       │
│  from all AMR figures — only acquired genes are plotted.                │
│                                                                          │
│  All species                                                             │
│  ─────────────────────────────────────────────────────────────────────  │
│  Fig 1  Population summary (4 panels):                                  │
│           A — MLST sequence type distribution                            │
│           B — Serotype bars / IS element landscape (Shigella)           │
│           C — AMR drug class prevalence                                  │
│           D — MDR burden per isolate                                     │
│  Fig 2  Phylogenetic tree + ST/phylogroup strips + AMR heatmap          │
│  Fig 3  Top acquired AMR genes (intrinsic genes excluded)               │
│  Fig 4  Plasmid overview (3 panels):                                    │
│           A — Replicon prevalence stacked by AMR drug class             │
│           B — Plasmid–AMR co-occurrence bubble matrix                   │
│           C — SNP phylogeny + per-sample replicon presence heatmap      │
│         Individual panels also saved to individual_plasmid_plots/       │
│  Fig 5  Virulence genes / pathotype                                     │
│  Fig 6  AMR drug class prevalence by MLST sequence type                 │
│  Fig 7  AMR drug class prevalence by serovar / Clermont phylogroup      │
│  (SNP)  Pairwise whole-genome SNP distance heatmap                      │
│                                                                          │
│  Shigella only                                                           │
│  ─────────────────────────────────────────────────────────────────────  │
│  Fig 8  Species + serotype composition (stacked bars)                   │
│  Fig 9  Virulence & invasion feature panel (ipaH, pINV genes, IS elems) │
│                                                                          │
│  Assembly QC (all species)                                               │
│  ─────────────────────────────────────────────────────────────────────  │
│  Assembly metrics  Per-species 8-panel figure (PDF + PNG):              │
│           A — Genome length histogram   B — Genome length boxplot        │
│           C — N50 histogram             D — N50 boxplot                  │
│           E — Contig count histogram    F — Contig count boxplot         │
│           G — GC% histogram             H — GC% boxplot                  │
└──────────────────────────────────────────────────────────────────────────┘
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
> IQ-TREE 2 for maximum-likelihood tree inference (default: GTR+G; set
> `--iqtree_model MFP` to enable ModelFinder Plus). Phylogenetics runs per species
> when ≥ 3 samples are present; skip with `--skip_local_phylo` for faster runs.

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

This downloads 13 reference genomes from NCBI and builds the Mash sketch used
for species identification. Run once before the first pipeline execution:

```bash
bash assets/build_references.sh
```

Output: `assets/enteric_species_refs.msh`
Expected runtime: 1–5 minutes depending on bandwidth.

Reference genomes included:

| Species / phylogroup | Strain | Accession |
|---|---|---|
| *E. coli* phylogroup A | K-12 MG1655 | GCF_000005845.2 |
| *E. coli* phylogroup B1 | SE11 | GCF_000010485.1 |
| *E. coli* phylogroup B2 | CFT073 | GCF_000007445.1 |
| *E. coli* phylogroup D | UMN026 | GCF_000026265.1 |
| *E. coli* phylogroup E | O157:H7 EDL933 | GCF_000006665.1 |
| *S. enterica* Typhimurium | LT2 | GCF_000006945.2 |
| *S. enterica* Typhi | CT18 | GCF_000195995.1 |
| *S. enterica* Enteritidis | P125109 | GCF_000009505.1 |
| *S. sonnei* | ATCC 29930 | GCF_002950395.1 |
| *S. flexneri* | 2a 2457T | GCF_000007405.1 |
| *S. boydii* | Sb227 | GCF_000012025.1 |
| *S. dysenteriae* | Sd197 | GCF_000012005.1 |
| *K. pneumoniae* | HS11286 | GCF_000240185.1 |

> Species identification uses **closest-reference-wins** (same approach as
> Kleborate): the species group with the lowest Mash distance to any reference
> in the sketch is assigned. Genuine *Shigella* isolates are correctly
> identified because their closest reference is always a *Shigella* genome
> (Mash dist < 0.012), while *E. coli* phylogroup B2 isolates — which are
> phylogenetically close to *Shigella* but are not *Shigella* — score lower
> distances to *E. coli* references (~0.007) than to any *Shigella* reference
> (≥ 0.014) and are correctly called *E. coli*.

---

### Step 7 — Download Kraken2 database (optional)

The contamination screen is **opt-in** and is skipped entirely if `--kraken2_db`
is not provided. To enable it, download the `k2_standard_08gb` database (~8 GB)
once and point the pipeline to it at run time:

```bash
# Create a directory for the database
mkdir -p ~/kraken2_db

# Download and extract (requires ~8 GB disk space)
cd ~/kraken2_db
curl -L -O https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz
tar -xzf k2_standard_08gb_20240904.tar.gz
```

Then pass the database path when running the pipeline:

```bash
nextflow run main.nf \
    -profile conda \
    --input_dir assemblies/ \
    --outdir    results/ \
    --kraken2_db ~/kraken2_db/
```

The database can be reused across projects. Check the
[Kraken2 index page](https://benlangmead.github.io/aws-indexes/k2) for the
latest release if a newer build is preferred.

> **Genome size filters** require no setup — thresholds are applied automatically
> based on species classification. Override defaults with e.g.
> `--ecoli_min_length 4000000` if needed.
>
> Thresholds are based on published QC criteria:
> \[1\] E. coli / Shigella — [BIGSdb *E. coli* genome quality criteria](https://bigsdb.pasteur.fr/ecoli/genomes-quality-criteria/)
> \[2\] Salmonella — [PATH-SAFE consortium recommendations for genomic surveillance of foodborne diseases](https://science.food.gov.uk/article/143833-path-safe-consortium-recommendations-for-genomic-surveillance-of-foodborne-diseases-using-salmonella-as-an-exemplar?attachment_id=300637) (Table 3)

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

# Apple Silicon Mac (M1 and above) — prefix CONDA_SUBDIR and add arm64 profile
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
| `--iqtree_model` | `GTR+G` | IQ-TREE substitution model (use `MFP` for automatic model selection) |
| `--iqtree_bootstraps` | `100` | IQ-TREE ultrafast bootstrap replicates |
| **Assembly QC** | | |
| `--ecoli_min_length` | `4300000` | E. coli minimum assembly length (bp) |
| `--ecoli_max_length` | `5900000` | E. coli maximum assembly length (bp) |
| `--shigella_min_length` | `4300000` | Shigella minimum assembly length (bp) |
| `--shigella_max_length` | `5900000` | Shigella maximum assembly length (bp) |
| `--salmonella_min_length` | `4100000` | Salmonella minimum assembly length (bp) |
| `--salmonella_max_length` | `6600000` | Salmonella maximum assembly length (bp) |
| `--kraken2_db` | `null` | Path to Kraken2 database directory; contamination screen skipped if not set |
| `--max_contamination` | `3.0` | Maximum % secondary species allowed (Kraken2) |

---

## Output files

```
results/
├── pipeline_info/
│   ├── timeline.html
│   ├── report.html
│   └── dag.svg
│
├── qc/
│   ├── assembly_size_qc_summary.tsv  ← Pass/fail for every input sample (size + contamination)
│   └── *_kraken2_summary.tsv         ← Per-sample Kraken2 contamination results (if enabled)
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
├── kaptive_g2g3/ kaptive_g1g4/
│   └── *_ktype.tsv                ← K-locus group, locus, type, confidence
├── kleborate/
│   └── *_kleborate.tsv            ← Pathotype, Clermont phylogroup, virulence markers
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
├── mlst_shigella/
│   └── *_shigella_mlst.tsv
├── amrfinder_shigella/
│   └── *_amrfinder.tsv
├── shigeifinder/
│   └── *_shigeifinder.tsv        ← ipaH, virulence plasmid, cluster, serotype, O/H antigens
├── mykrobe/
│   └── *_mykrobe_parsed.tsv      ← S. sonnei lineage/genotype (NA for other Shigella spp.)
├── pinv_screen/
│   └── *_pinv.tsv                ← pINV marker gene presence (icsA/virG, virF, ipaB/C/D)
├── is_screen/
│   └── *_is.tsv                  ← IS element copy counts (IS1, IS30, IS600, IS629…)
├── plasmidfinder_shigella/
│   └── *_plasmidfinder.tsv
│
├── ska2_{species}/
│   ├── ska2_alignment.fasta       ← whole-genome SNP alignment (input to IQ-TREE)
│   └── snp_matrix.tsv             ← pairwise whole-genome SNP distance matrix
│
├── iqtree_{species}/
│   ├── iqtree.treefile            ← Maximum-likelihood tree (Newick)
│   ├── iqtree.contree             ← Consensus tree with bootstrap support
│   ├── iqtree.mldist              ← ML pairwise distance matrix
│   └── iqtree.iqtree              ← Full IQ-TREE log and model summary
│
├── ecoli_typer_results.tsv         ← Master results table (E. coli)
├── salmonella_typer_results.tsv    ← Master results table (Salmonella)
├── shigella_typer_results.tsv      ← Master results table (Shigella)
│
├── plasmid_amr_map.tsv             ← Per-contig plasmid–AMR linkage table
│     One row per (sample, replicon). Columns include:
│       sample_id, replicon, contig_id, drug_classes, likely_location
│       (chromosome / plasmid inferred from contig length)
│
├── plasmid_amr_map/
│   └── *_plasmid_amr_map.tsv      ← Per-sample plasmid–AMR contig map
│
│   All master tables include:
│     amrfinder_acquired_genes   — resistance genes (intrinsic excluded)
│     amrfinder_intrinsic_genes  — wildtype/intrinsic genes (flagged, not plotted)
│     amrfinder_genes            — all raw AMRFinder hits
│
│ ── Summary figures (PDF + PNG) ──────────────────────────────────────────
│
├── {species}_fig1_population_summary.{pdf,png}
│     Panel A: MLST sequence type distribution
│               E. coli — bars stacked by Clermont phylogroup
│               Salmonella — single-hue gradient
│               Shigella — bars stacked by species (S. sonnei / flexneri / boydii / dysenteriae)
│     Panel B: Serotype/serovar prevalence (E. coli / Salmonella)
│               IS element copy-number heatmap (Shigella)
│     Panel C: AMR drug class prevalence
│     Panel D: MDR burden per isolate (≥ 3 acquired drug classes = MDR)
│
├── {species}_fig2_tree_amr.{pdf,png}
│     Phylogenetic tree (IQ-TREE ML) with:
│       — ST colour strip
│       — Clermont phylogroup strip (E. coli only)
│       — Virulence gene binary heatmap
│       — AMR gene heatmap grouped by drug class
│
├── {species}_fig3_amr_genes.{pdf,png}      ← Acquired AMR gene frequencies
│
├── {species}_fig4_plasmid_overview.{pdf,png}
│     3-panel plasmid figure:
│       Panel A — Horizontal stacked-bar chart: replicon prevalence broken
│                 down by dominant AMR drug class (clinical priority order)
│       Panel B — Bubble matrix: % of isolates carrying each replicon × drug
│                 class co-occurrence
│       Panel C — SNP phylogeny with ST and phylogroup strips + per-sample
│                 replicon presence/absence heatmap (blue = present)
│
├── individual_plasmid_plots/
│   ├── {species}_fig4_panel_a_replicon_bars.{pdf,png}
│   ├── {species}_fig4_panel_b_bubble_matrix.{pdf,png}
│   └── {species}_fig4_panel_c_tree_heatmap.{pdf,png}
│
├── {species}_fig4_plasmid_replicons.{pdf,png}
│     Standalone replicon prevalence bar chart (also produced by fig4 Panel A)
│
├── {species}_fig5_virulence.{pdf,png}       ← Virulence genes / pathotype
│
├── {species}_fig6_amr_by_st.{pdf,png}
│     AMRnet-style tile heatmap: % isolates carrying each drug class, by MLST ST
│
├── {species}_fig7_amr_by_group.{pdf,png}
│     AMRnet-style tile heatmap: % isolates carrying each drug class, by
│     serovar (Salmonella) or Clermont phylogroup (E. coli) or serotype (Shigella)
│
├── {species}_snp_heatmap.{pdf,png}
│     Pairwise whole-genome SNP distance heatmap (hierarchically clustered).
│     Sample labels shown at 90° on both axes for datasets ≤ 200 isolates;
│     hidden (with count) for larger datasets.
│
│ ── Shigella-specific figures ────────────────────────────────────────────
│
├── shigella_fig8_shigella_serotypes.{pdf,png}
│     Species composition + serotype breakdown (stacked bar chart)
│
├── shigella_fig9_shigella_features.{pdf,png}
│     Binary heatmap: ipaH · virulence plasmid · pINV invasion genes (icsA/virG,
│     virF, virB, ipaB, ipaC, ipaD) · IS elements per sample
│
│ ── Assembly QC figures (PDF + PNG) ──────────────────────────────────────
│
├── {species}_assembly_metrics.{pdf,png}
│     8-panel assembly QC figure (one per species detected):
│       A  Genome length histogram    B  Genome length boxplot
│       C  Assembly N50 histogram     D  Assembly N50 boxplot
│       E  Contig count histogram     F  Contig count boxplot
│       G  GC% histogram              H  GC% boxplot
│
└── {species}_assembly_metrics_summary.tsv
      Merged per-sample assembly stats (genome_length, num_contigs,
      assembly_N50, gc_pct) for all samples of that species.
```

---

## AMR drug class abbreviations

Figures 6 and 7 use short abbreviations for antimicrobial drug classes.
Full class names and clinical notes are given below.

| Abbreviation | Drug class | Representative agents |
|---|---|---|
| **AMG** | Aminoglycoside | Gentamicin, amikacin, tobramycin, streptomycin |
| **BLA** | Beta-lactam | Ampicillin, cephalosporins, carbapenems, penicillins |
| **COL** | Colistin | Colistin (polymyxin E), polymyxin B |
| **FOS** | Fosfomycin | Fosfomycin |
| **FOSM** | Fosmidomycin | Fosmidomycin |
| **LIN** | Lincosamide | Lincomycin, clindamycin |
| **MAC** | Macrolide | Azithromycin, erythromycin |
| **NIT** | Nitrofuran | Nitrofurantoin |
| **PHE** | Phenicol | Chloramphenicol, florfenicol |
| **QNL** | Quinolone | Ciprofloxacin, nalidixic acid, levofloxacin (all fluoroquinolones) |
| **SGM** | Streptogramin | Quinupristin–dalfopristin (streptogramin A/B combinations) |
| **STR** | Streptothricin | Nourseothricin |
| **SUL** | Sulfonamide | Sulfamethoxazole, sulfisoxazole |
| **TET** | Tetracycline | Tetracycline, doxycycline, tigecycline |
| **TMP** | Trimethoprim | Trimethoprim (often combined with sulfonamide as co-trimoxazole) |

> **MDR definition:** an isolate is classified as multidrug-resistant (MDR) when
> it carries acquired resistance genes in **≥ 3** of the above drug classes.
> Class EFFLUX is excluded from the MDR count as near-universal efflux pumps
> are intrinsic to the species and are removed by AMRrules before plotting.

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
>
> **Note on `-resume` and tree annotation:** if a previous run failed to produce
> a tree annotation figure (e.g. due to a software update), Nextflow may cache
> the failed state. Run once **without `-resume`** to force a clean re-run of
> the annotation step.

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
| `arm64` | **Add on Apple Silicon (M1 and above)** — forces osx-64 conda envs via Rosetta 2 |
| `docker` | Local with Docker Desktop |
| `singularity` | HPC cluster with Singularity/Apptainer |
| `slurm` | SLURM HPC executor (combine with singularity: `-profile singularity,slurm`) |
| `pbs` | PBS/Torque HPC executor |
| `test` | Quick test run (no uploads) |

### macOS Apple Silicon (M1 and above)

Some Bioconda packages have no native arm64 build. Add the `arm64` profile to force
Rosetta 2 emulation — tools run at near-native speed and the pipeline produces
identical results to Linux. Kleborate will automatically use MLST-only mode on ARM64
if its full pathotype dependencies are unavailable, but all other tools run at full
capability:

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

After editing, rebuild the sketch:
```bash
bash assets/build_references.sh
```

---

## Vignettes

Example outputs from real datasets:

| Dataset | Species | Samples | Vignette |
|---|---|---|---|
| NHP gut isolates, The Gambia (Foster-Nyarko et al. 2020) | *Escherichia coli* | 98 | [View vignette](vignettes/ecoli_vignette.md) |
| Clinical NTS isolates, The Gambia (Darboe et al. 2022) | *Salmonella enterica* | 99 | [View vignette](vignettes/salmonella_vignette.md) |
| Diverse *Shigella* — all four species (reference genomes) | *Shigella* spp. | 15 | [View vignette](vignettes/shigella_vignette.md) |

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
- **ShigEiFinder**: LanLab (2022) Microbial Genomics; github.com/LanLab/ShigEiFinder
- **Mykrobe**: Hunt et al. (2015) Genome Biology 16:239; Hawkey et al. (2022) Microbial Genomics 8(6)
- **SKA2**: github.com/bacpop/ska.rust
- **IQ-TREE 2**: Minh et al. (2020) Molecular Biology and Evolution 37(5)
- **AMRrules**: github.com/AMRverse/AMRrules

