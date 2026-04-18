🌐 [English](README.md) | **Français** | [Português](README.pt.md)

# enteric-typer

Un processus de génotypage basé sur l'identification des espèces pour les agents pathogènes entériques. À partir d'un dossier contenant des assemblages de génomes, ce pipeline identifie l'espèce de chaque échantillon et utilise les outils de typage spécifiques à cette espèce, génère une phylogénie SNP à l'échelle du génome entier et produit des figures récapitulatives prêtes à être publiées.

## Espèces prises en charge

| Espèce | Outils de typage |
|---|---|
| *Escherichia coli* | MLST (Achtman), AMRFinder, ECTyper (sérotype), EzClermont (phylogroupe Clermont), Kleborate (pathotype), PlasmidFinder, Kaptive (locus K) |
| *Salmonella enterica* | MLST, AMRFinder, SISTR (sérovar), PlasmidFinder |
| *Shigella* spp. | MLST (Achtman), AMRFinder, ShigEiFinder (sérotype/espèce), Mykrobe (génotypage *S. sonnei*), PlasmidFinder, criblage pINV, criblage éléments IS |
| Autre / non classé | Identification de l'espèce uniquement (enregistré et ignoré) |

> **Kleborate :** effectue une détection complète du pathotype sur toutes les plateformes. Sur macOS Apple Silicon (ARM64) sans Rosetta 2, il passe automatiquement en mode MLST uniquement — tous les autres outils fonctionnent à pleine capacité quelle que soit la plateforme (Linux, macOS Intel/ARM64, HPC).

> **Classification des gènes de résistance :** tous les résultats AMRFinder sont classés par [AMRrules](https://github.com/AMRverse/AMRrules) en gènes de résistance *acquis* et gènes *intrinsèques* (propres à l'espèce sauvage). Les gènes intrinsèques sont conservés dans les fichiers TSV de résultats à titre de référence mais sont **exclus de tous les graphiques AMR** afin que les figures ne reflètent que la résistance acquise cliniquement pertinente.

## Vue d'ensemble du pipeline

```
Assemblages d'entrée (dossier ou feuille d'échantillons)
        │
        ▼
┌───────────────────────────────────────────────────────────────────────────┐
│  1. VÉRIFICATION D'ESPÈCE ET CONTRÔLE QUALITÉ                             │
│                                                                           │
│  Identification de l'espèce — distance Mash par rapport à une esquisse   │
│  de référence à 13 génomes ; la référence la plus proche l'emporte        │
│  (même approche que Kleborate).                                           │
│                                                                           │
│  Filtres de contrôle qualité des assemblages (échecs enregistrés et       │
│  exclus du typage)                                                        │
│  ──────────────────────────────────────────────────────────────────────   │
│  Taille du génome  E. coli / Shigella   4,3 – 5,9 Mb  [1]                │
│                    Salmonella           4,1 – 6,6 Mb  [2]                │
│  Contamination     Espèce secondaire Kraken2 < 3 % des contigs totaux     │
│                    (optionnel — fournir --kraken2_db ; ignoré si absent)  │
└───────────────────────────────────────────────────────────────────────────┘
        │
   ┌────┼────────┐
   ▼    ▼        ▼
E. coli  Salmonella  Shigella   (autres espèces enregistrées et ignorées)
   │         │           │
   ▼         ▼           ▼
┌───────────────────────────────────────────────────────────────────────────┐
│  2. TYPAGE SPÉCIFIQUE À L'ESPÈCE  (tous les outils en parallèle)          │
│                                                                           │
│  E. coli                       Salmonella          Shigella               │
│  ──────────────────────────    ──────────────────  ──────────────────     │
│  MLST (achtman_4)              MLST (senterica_    MLST (achtman_4)       │
│  AMRFinder                       achtman_2)        AMRFinder              │
│  ECTyper (sérotype O:H)        AMRFinder           ShigEiFinder           │
│  EzClermont (phylogroupe)      SISTR (sérovar)       (sérotype/espèce)   │
│  Kleborate (pathotype)         PlasmidFinder       Mykrobe                │
│  PlasmidFinder                                       (génotype S. sonnei) │
│  Kaptive locus K (G2/G3                            PlasmidFinder          │
│    → G1/G4 si non typable)                         criblage pINV          │
│                                                    criblage éléments IS   │
└───────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  3. PHYLOGÉNÉTIQUE  (par espèce, ≥ 3 échantillons ; à désactiver via     │
│                    --skip_local_phylo)                                   │
│  SKA2 build (k=31) → alignement SNP génome entier + matrice de           │
│  distances SNP → arbre ML IQ-TREE (modèle GTR+G ;                        │
│  remplacer par --iqtree_model MFP)                                        │
└──────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────┐
│  4. AGRÉGATION  (un TSV par espèce)                          │
│  ecoli_typer_results.tsv                                     │
│  salmonella_typer_results.tsv                                │
│  shigella_typer_results.tsv                                  │
│  Résultats AMRFinder classés par AMRrules en :               │
│    amrfinder_acquired_genes  — gènes de résistance acquise   │
│                                cliniquement pertinents        │
│    amrfinder_intrinsic_genes — gènes intrinsèques à l'espèce │
│                                (signalés, conservés dans TSV) │
└──────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  5. FIGURES RÉCAPITULATIVES                                               │
│  Les gènes de résistance intrinsèques (classés par AMRrules) sont        │
│  exclus de toutes les figures AMR — seuls les gènes acquis sont tracés.  │
│                                                                           │
│  Toutes espèces                                                           │
│  ─────────────────────────────────────────────────────────────────────   │
│  Fig 1  Résumé de la population (4 panneaux) :                           │
│           A — Distribution des types de séquence MLST                    │
│           B — Barres de sérotype / paysage éléments IS (Shigella)        │
│           C — Prévalence par classe d'antibiotiques                       │
│           D — Charge de multirésistance par isolat                        │
│  Fig 2  Arbre phylogénétique + bandes ST/phylogroupe + carte thermique de  │
│         la résistance aux antimicrobiens (AMR)                            │
│  Fig 3  Principaux gènes AMR acquis (gènes intrinsèques exclus)           │
│  Fig 4  Types de réplicons plasmidiques                                   │
│  Fig 5  Gènes de virulence / pathotype                                    │
│  Fig 6  Prévalence des classes AMR par type de séquence MLST              │
│  Fig 7  Prévalence des classes AMR par sérovar / phylogroupe Clermont     │
│  (SNP)  Carte thermique des distances SNP génome entier par paire         │
│                                                                           │
│  Shigella uniquement                                                      │
│  ─────────────────────────────────────────────────────────────────────   │
│  Fig 8  Composition en espèces + répartition des sérotypes (barres empil.)│
│  Fig 9  Panneau de virulence et invasion (ipaH, gènes pINV, éléments IS) │
│                                                                           │
│  Contrôle qualité des assemblages (toutes espèces)                        │
│  ─────────────────────────────────────────────────────────────────────   │
│  Métriques d'assemblage  Figure à 8 panneaux par espèce (PDF + PNG) :    │
│           A — Histogramme longueur de génome   B — Diagramme en boîte longueur génome  │
│           C — Histogramme N50                  D — Diagramme en boîte N50              │
│           E — Histogramme nombre de contigs    F — Diagramme en boîte nombre de contigs│
│           G — Histogramme GC%                  H — Diagramme en boîte GC%              │
└──────────────────────────────────────────────────────────────────────────┘
```

> **Typage du locus K d'*E. coli* (Kaptive) :** le typage du locus K utilise deux bases de données curées exécutées séquentiellement. Les assemblages sont d'abord typés sur la **base de données G2/G3** (`EC-K-typing_group2and3_v3.0.0.gbk` ; [Gladstone et al. 2026](https://www.nature.com/articles/s41564-026-02283-w)). Les échantillons restant non typables sont ensuite retypés sur la **base de données G1/G4** (`EC-K-typing_group1and4_v1.2.gbk` ; [Foster-Nyarko et al.](https://github.com/efosternyarko/EC-K-typing-G1G4)). Les loci G2/G3 et G1/G4 sont mutuellement exclusifs chez *E. coli*, donc le typage séquentiel garantit que chaque échantillon est assigné au bon groupe sans double comptage.

> **Phylogénétique (SKA2 + IQ-TREE) :** les phylogénies SNP génome entier sont construites avec [SKA2](https://github.com/bacpop/ska.rust) (alignement par k-mers divisés, k=31), qui aligne les assemblages sans génome de référence et produit un alignement SNP à l'échelle du génome ainsi qu'une matrice de distances SNP par paires. L'alignement est ensuite transmis à IQ-TREE 2 pour l'inférence d'arbre par maximum de vraisemblance (par défaut : GTR+G ; utiliser `--iqtree_model MFP` pour activer ModelFinder Plus). La phylogénétique s'exécute par espèce lorsque ≥ 3 échantillons sont présents ; désactivable via `--skip_local_phylo` pour des exécutions plus rapides.

---

## Installation

### Étape 1 — Cloner le dépôt

```bash
git clone https://github.com/efosternyarko/enteric-typer
cd enteric-typer
```

---

### Étape 2 — Installer Java

Nextflow nécessite Java 17 ou une version plus récente (Java 21 recommandé).

**macOS**
```bash
brew install --cask temurin
```
> Si Homebrew n'est pas encore installé, exécuter :
> ```bash
> /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
> ```

**Linux (Debian/Ubuntu)**
```bash
sudo apt update && sudo apt install -y default-jdk
```

**Windows (PowerShell)**
```powershell
winget install EclipseAdoptium.Temurin.21.JDK   # Java 21 (recommandé)
# ou
winget install EclipseAdoptium.Temurin.25.JDK   # Java 25 (dernière version)
```

Vérifier : `java -version`

---

### Étape 3 — Installer Nextflow

```bash
curl -s https://get.nextflow.io | bash
```

L'installeur crée un exécutable `nextflow` dans le répertoire courant. Le déplacer dans un dossier présent dans le `$PATH` pour pouvoir l'exécuter depuis n'importe où :

**macOS**
```bash
# Si Homebrew est installé, déplacer dans le dossier bin Homebrew (déjà dans $PATH) :
mv nextflow /opt/homebrew/bin/        # Apple Silicon (M1/M2/M3)
# ou
mv nextflow /usr/local/bin/           # Mac Intel

# Si Homebrew n'est pas installé, créer ~/bin et l'ajouter au PATH :
mkdir -p ~/bin
mv nextflow ~/bin/
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

**Linux**
```bash
mv nextflow ~/bin/
# Si ~/bin/ n'existe pas encore :
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
Puis redémarrer PowerShell.

Vérifier : `nextflow -version`

---

### Étape 4 — Installer conda + mamba

Si conda n'est pas encore installé, la voie la plus rapide est **Miniforge**, qui inclut à la fois conda et mamba :

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

Lorsque l'installeur demande :
```
Do you wish the installer to initialize Miniforge3
by running conda init? [yes|no]
```
Répondre **yes**. Puis recharger le shell :

```bash
source ~/.zshrc    # macOS (zsh)
# ou
source ~/.bashrc   # Linux (bash)
```

**Windows**
Télécharger l'[installeur Miniforge3 pour Windows](https://github.com/conda-forge/miniforge/releases/latest) et l'exécuter. Lorsque demandé, sélectionner *Add Miniforge to my PATH*. Puis réinitialiser :
```powershell
conda init powershell
```
Redémarrer PowerShell.

> **Vous avez déjà conda mais pas mamba ?**
> ```bash
> conda install -n base -c conda-forge mamba
> ```

---

### Étape 5 — Installer ncbi-datasets-cli et mash

Ces deux outils sont requis par `build_references.sh` pour télécharger les génomes de référence et construire l'esquisse d'espèces.

> **Bonne pratique :** installer ces outils dans un environnement dédié plutôt que dans `base`.
> Cela évite les conflits de versions et de dépendances avec d'autres logiciels déjà installés.

```bash
# Créer un environnement dédié aux utilitaires de téléchargement (à exécuter une seule fois)
mamba create -n enteric-tools -c conda-forge -c bioconda ncbi-datasets-cli mash -y
```

L'environnement `enteric-tools` n'est nécessaire que le temps d'exécuter `build_references.sh` (étape 6).
Après cela, le désactiver — le pipeline gère ensuite ses propres environnements automatiquement :

```bash
conda activate enteric-tools
bash assets/build_references.sh
conda deactivate
```

---

### Étape 6 — Construire l'esquisse de référence Mash

Cette étape télécharge 13 génomes de référence depuis NCBI et construit l'esquisse Mash utilisée pour l'identification des espèces. À exécuter une seule fois avant la première exécution du pipeline :

```bash
bash assets/build_references.sh
```

Sortie : `assets/enteric_species_refs.msh`
Durée estimée : 1–5 minutes selon la bande passante.

Génomes de référence inclus :

| Espèce / phylogroupe | Souche | Accession |
|---|---|---|
| *E. coli* phylogroupe A | K-12 MG1655 | GCF_000005845.2 |
| *E. coli* phylogroupe B1 | SE11 | GCF_000010485.1 |
| *E. coli* phylogroupe B2 | CFT073 | GCF_000007445.1 |
| *E. coli* phylogroupe D | UMN026 | GCF_000026265.1 |
| *E. coli* phylogroupe E | O157:H7 EDL933 | GCF_000006665.1 |
| *S. enterica* Typhimurium | LT2 | GCF_000006945.2 |
| *S. enterica* Typhi | CT18 | GCF_000195995.1 |
| *S. enterica* Enteritidis | P125109 | GCF_000009505.1 |
| *S. sonnei* | ATCC 29930 | GCF_002950395.1 |
| *S. flexneri* | 2a 2457T | GCF_000007405.1 |
| *S. boydii* | Sb227 | GCF_000012025.1 |
| *S. dysenteriae* | Sd197 | GCF_000012005.1 |
| *K. pneumoniae* | HS11286 | GCF_000240185.1 |

> L'identification des espèces repose sur le principe du « référence la plus proche » (même approche que Kleborate) : le groupe d'espèces présentant la distance Mash la plus faible par rapport à n'importe quelle référence figurant dans le schéma est attribué. Les isolats authentiques de *Shigella* sont correctement identifiés car leur référence la plus proche est toujours un génome de *Shigella* (distance Mash < 0,012), tandis que les isolats *E. coli* du phylogroupe B2 — phylogénétiquement proches de *Shigella* mais n'en étant pas — obtiennent des distances plus faibles aux références *E. coli* (~0,007) qu'à toute référence *Shigella* (≥ 0,014) et sont correctement classés comme *E. coli*.

---

### Étape 7 — Télécharger la base de données Kraken2 (optionnel)

Le criblage de contamination est **optionnel** et est ignoré si `--kraken2_db` n'est pas fourni. Pour l'activer, télécharger la base de données `k2_standard_08gb` (~8 Go) une fois et indiquer son chemin au pipeline lors de l'exécution :

```bash
# Créer un répertoire pour la base de données
mkdir -p ~/kraken2_db

# Télécharger et extraire (nécessite ~8 Go d'espace disque)
cd ~/kraken2_db
curl -L -O https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz
tar -xzf k2_standard_08gb_20240904.tar.gz
```

Puis indiquer le chemin de la base de données lors de l'exécution du pipeline :

```bash
nextflow run main.nf \
    -profile conda \
    --input_dir assemblies/ \
    --outdir    results/ \
    --kraken2_db ~/kraken2_db/
```

La base de données peut être réutilisée pour plusieurs projets. Consulter la [page d'index Kraken2](https://benlangmead.github.io/aws-indexes/k2) pour la version la plus récente.

> **Les filtres de taille de génome** ne nécessitent aucune configuration — les seuils sont appliqués automatiquement selon la classification de l'espèce. Les valeurs par défaut peuvent être remplacées avec par exemple `--ecoli_min_length 4000000`.
>
> Les seuils sont basés sur des critères de contrôle qualité publiés :
> \[1\] *E. coli* / *Shigella* — [Critères de qualité des génomes *E. coli* de BIGSdb](https://bigsdb.pasteur.fr/ecoli/genomes-quality-criteria/)
> \[2\] *Salmonella* — [Recommandations du consortium PATH-SAFE pour la surveillance génomique des maladies d'origine alimentaire](https://science.food.gov.uk/article/143833-path-safe-consortium-recommendations-for-genomic-surveillance-of-foodborne-diseases-using-salmonella-as-an-exemplar?attachment_id=300637) (Tableau 3)

---

## Démarrage rapide

> **Linux / Ubuntu / Debian :** utiliser les commandes Linux ci-dessous.
> **Windows :** nous recommandons d'exécuter enteric-typer dans [WSL2](https://learn.microsoft.com/fr-fr/windows/wsl/install) (Sous-système Windows pour Linux), qui fournit un environnement Ubuntu complet. Installer WSL2, puis suivre les instructions Linux. Vous pouvez également utiliser le profil `docker` avec Docker Desktop pour Windows.

### Pour un dossier d'assemblages

```bash
# Linux / Ubuntu / Debian / Mac Intel
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda

# Mac Apple Silicon (M1 et supérieur) — préfixer CONDA_SUBDIR et ajouter le profil arm64
CONDA_SUBDIR=osx-64 nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda,arm64

# Windows (Docker Desktop) — conda non requis
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile docker
```

### Ou avec une feuille d'échantillons

```bash
# Générer automatiquement une feuille d'échantillons depuis un dossier
python bin/make_samplesheet.py \
    --input /path/to/assemblies/ \
    --output samples.csv

# Exécuter avec la feuille d'échantillons
nextflow run main.nf \
    --samplesheet samples.csv \
    --outdir      results/ \
    -profile conda
```

### Ignorer la phylogénétique locale (plus rapide, sans SKA2/IQ-TREE)

```bash
nextflow run main.nf \
    --input_dir        /path/to/assemblies/ \
    --outdir           results/ \
    --skip_local_phylo \
    -profile conda
```

> Lorsque `--skip_local_phylo` est activé, la matrice de distances SNP, la carte thermique SNP et les figures d'annotation de l'arbre ne sont pas produites.

---

## Paramètres

| Paramètre | Défaut | Description |
|---|---|---|
| `--input_dir` | `null` | Dossier d'assemblages FASTA (`.fasta/.fa/.fna/.fas`) |
| `--samplesheet` | `null` | CSV avec les colonnes `id,fasta` |
| `--outdir` | `results` | Répertoire de sortie |
| `--skip_local_phylo` | `false` | Ignorer SKA2 + IQ-TREE (sans arbre, sans matrice SNP/carte thermique) |
| `--ska2_min_samples` | `3` | Nombre minimum d'échantillons pour tenter SKA2/IQ-TREE |
| `--iqtree_model` | `GTR+G` | Modèle de substitution IQ-TREE (utiliser `MFP` pour sélection automatique) |
| `--iqtree_bootstraps` | `100` | Réplicats bootstrap ultra-rapide IQ-TREE |
| **Contrôle qualité des assemblages** | | |
| `--ecoli_min_length` | `4300000` | Longueur minimale d'assemblage *E. coli* (pb) |
| `--ecoli_max_length` | `5900000` | Longueur maximale d'assemblage *E. coli* (pb) |
| `--shigella_min_length` | `4300000` | Longueur minimale d'assemblage *Shigella* (pb) |
| `--shigella_max_length` | `5900000` | Longueur maximale d'assemblage *Shigella* (pb) |
| `--salmonella_min_length` | `4100000` | Longueur minimale d'assemblage *Salmonella* (pb) |
| `--salmonella_max_length` | `6600000` | Longueur maximale d'assemblage *Salmonella* (pb) |
| `--kraken2_db` | `null` | Chemin vers le répertoire de la base de données Kraken2 ; criblage ignoré si absent |
| `--max_contamination` | `3.0` | Pourcentage maximal d'espèce secondaire autorisé (Kraken2) |

---

## Fichiers de sortie

```
results/
├── pipeline_info/
│   ├── timeline.html
│   ├── report.html
│   └── dag.svg
│
├── species_check/
│   └── *_species.txt              ← meilleure espèce + distance Mash par échantillon
│
├── mlst_ecoli/
│   └── *_ecoli_achtman_4_mlst.tsv
├── amrfinder_ecoli/
│   └── *_amrfinder.tsv
├── ectyper/
│   └── *_ectyper.tsv
├── kaptive_g2g3/ kaptive_g1g4/
│   └── *_ktype.tsv                ← groupe locus K, locus, type, niveau de confiance
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
│   └── *_shigeifinder.tsv        ← ipaH, plasmide de virulence, cluster, sérotype, antigènes O/H
├── mykrobe/
│   └── *_mykrobe_parsed.tsv      ← lignée/génotype S. sonnei (NA pour les autres Shigella spp.)
├── pinv_screen/
│   └── *_pinv.tsv                ← présence des gènes marqueurs pINV (icsA/virG, virF, ipaB/C/D)
├── is_screen/
│   └── *_is.tsv                  ← comptage des copies d'éléments IS (IS1, IS30, IS600, IS629…)
├── plasmidfinder_shigella/
│   └── *_plasmidfinder.tsv
│
├── ska2_{species}/
│   ├── ska2_alignment.fasta       ← alignement SNP génome entier (entrée pour IQ-TREE)
│   └── snp_matrix.tsv             ← matrice de distances SNP génome entier par paires
│
├── ecoli_typer_results.tsv         ← Tableau de résultats principal (E. coli)
├── salmonella_typer_results.tsv    ← Tableau de résultats principal (Salmonella)
├── shigella_typer_results.tsv      ← Tableau de résultats principal (Shigella)
│
│   Tous les tableaux principaux incluent :
│     amrfinder_acquired_genes   — gènes de résistance (intrinsèques exclus)
│     amrfinder_intrinsic_genes  — gènes intrinsèques sauvages (signalés, non tracés)
│     amrfinder_genes            — tous les résultats AMRFinder bruts
│
│ ── Figures récapitulatives (PDF + PNG) ──────────────────────────────────
│
├── {species}_fig1_population_summary.{pdf,png}
│     Panneau A : Distribution des types de séquence MLST
│               E. coli — barres empilées par phylogroupe Clermont
│               Salmonella — gradient monochrome
│               Shigella — barres empilées par espèce (S. sonnei / flexneri / boydii / dysenteriae)
│     Panneau B : Prévalence des sérotypes/sérovars (E. coli / Salmonella)
│               Carte thermique du nombre de copies d'éléments IS (Shigella)
│     Panneau C : Prévalence par classe d'antibiotiques
│     Panneau D : Charge de multirésistance par isolat (≥ 3 classes acquises = MDR)
│
├── {species}_fig2_tree_amr.{pdf,png}
│     Arbre phylogénétique (ML IQ-TREE) avec :
│       — bande de couleur ST
│       — bande de couleur phylogroupe Clermont (E. coli uniquement)
│       — carte thermique binaire gènes de virulence
│       — carte thermique gènes de résistance (AMR) regroupés par classe d'antibiotiques
│
├── {species}_fig3_amr_genes.{pdf,png}      ← Fréquences des gènes AMR acquis
├── {species}_fig4_plasmid_replicons.{pdf,png}
├── {species}_fig5_virulence.{pdf,png}       ← Gènes de virulence / pathotype
│
├── {species}_fig6_amr_by_st.{pdf,png}
│     Carte thermique AMRnet : % d'isolats portant chaque classe, par ST MLST
│
├── {species}_fig7_amr_by_group.{pdf,png}
│     Carte thermique AMRnet : % d'isolats portant chaque classe, par
│     sérovar (Salmonella) ou phylogroupe Clermont (E. coli) ou sérotype (Shigella)
│
├── {species}_snp_heatmap.{pdf,png}          ← Carte thermique des distances SNP par paires
│
│ ── Figures spécifiques à Shigella ──────────────────────────────────────
│
├── shigella_fig8_shigella_serotypes.{pdf,png}
│     Composition en espèces + répartition des sérotypes (barres empilées)
│
├── shigella_fig9_shigella_features.{pdf,png}
│     Carte thermique binaire : ipaH · plasmide de virulence · gènes d'invasion pINV
│     (icsA/virG, virF, virB, ipaB, ipaC, ipaD) · éléments IS par échantillon
│
│ ── Figures de contrôle qualité des assemblages (PDF + PNG) ─────────────
│
├── {species}_assembly_metrics.{pdf,png}
│     Figure CQ assemblage à 8 panneaux (une par espèce détectée) :
│       A  Histogramme longueur de génome    B  Diagramme en boîte longueur de génome
│       C  Histogramme N50                   D  Diagramme en boîte N50
│       E  Histogramme nombre de contigs     F  Diagramme en boîte nombre de contigs
│       G  Histogramme GC%                   H  Diagramme en boîte GC%
│
└── {species}_assembly_metrics_summary.tsv
      Statistiques d'assemblage par échantillon agrégées (genome_length, num_contigs,
      assembly_N50, gc_pct) pour tous les échantillons de cette espèce.
```

---

## Abréviations des classes d'antibiotiques

Les figures 6 et 7 utilisent des abréviations courtes pour les classes d'antibiotiques.
Les noms complets et les notes cliniques sont donnés ci-dessous.

| Abréviation | Classe d'antibiotiques | Agents représentatifs |
|---|---|---|
| **AMG** | Aminoside | Gentamicine, amikacine, tobramycine, streptomycine |
| **BLA** | Bêta-lactamine | Ampicilline, céphalosporines, carbapénèmes, pénicillines |
| **COL** | Colistine | Colistine (polymyxine E), polymyxine B |
| **FOS** | Fosfomycine | Fosfomycine |
| **FOSM** | Fosmidomycine | Fosmidomycine |
| **LIN** | Lincosamide | Lincomycine, clindamycine |
| **MAC** | Macrolide | Azithromycine, érythromycine |
| **NIT** | Nitrofurane | Nitrofurantoïne |
| **PHE** | Phénicol | Chloramphénicol, florfénicol |
| **QNL** | Quinolone | Ciprofloxacine, acide nalidixique, lévofloxacine (toutes fluoroquinolones) |
| **SGM** | Streptogramine | Quinupristine–dalfopristine (combinaisons streptogramine A/B) |
| **STR** | Streptothricine | Nourséothricine |
| **SUL** | Sulfonamide | Sulfaméthoxazole, sulfisoxazole |
| **TET** | Tétracycline | Tétracycline, doxycycline, tigécycline |
| **TMP** | Triméthoprime | Triméthoprime (souvent combiné avec sulfonamide en co-trimoxazole) |

> **Définition de la multirésistance (MDR) :** un isolat est classé multirésistant lorsqu'il porte des gènes de résistance acquise dans **≥ 3** des classes ci-dessus. La classe EFFLUX est exclue du comptage MDR car les pompes à efflux quasi-universelles sont intrinsèques à l'espèce et supprimées par AMRrules avant la production des figures.

---

## Nettoyage après une exécution

Tous les résultats définitifs sont écrits dans `--outdir` (par défaut : `results/`). Une fois les résultats vérifiés, les fichiers temporaires de Nextflow peuvent être supprimés pour libérer de l'espace disque :

```bash
rm -rf work/ .nextflow/ .nextflow.log*
```

> **Conseil :** conserver `work/` si vous souhaitez utiliser `-resume` pour relancer avec des paramètres différents sans répéter les étapes déjà terminées. Ne le supprimer que lorsque l'exécution est finalisée.
>
> **Note sur `-resume` et l'annotation de l'arbre :** si une exécution précédente n'a pas produit de figure d'annotation (par exemple suite à une mise à jour logicielle), Nextflow peut mettre en cache l'état d'échec. Exécuter une fois **sans `-resume`** pour forcer une réexécution propre de l'étape d'annotation.

---

## Optionnel : installer Graphviz

Nextflow génère un diagramme DAG d'exécution (`pipeline_info/dag.svg`) qui nécessite Graphviz pour le rendu. Sans lui, un avertissement sans conséquence s'affiche — le pipeline s'exécute normalement. Pour supprimer l'avertissement :

```bash
# macOS
brew install graphviz

# Linux
sudo apt install graphviz

# conda
conda install -c conda-forge graphviz
```

## Profils d'exécution

| Profil | Cas d'utilisation |
|---|---|
| `conda` | Station de travail locale avec conda |
| `mamba` | Identique à conda mais résolution d'environnement plus rapide |
| `arm64` | **À ajouter sur Apple Silicon (M1 et supérieur)** — force les environnements conda osx-64 via Rosetta 2 |
| `docker` | Local avec Docker Desktop |
| `singularity` | Cluster HPC avec Singularity/Apptainer |
| `slurm` | Exécuteur HPC SLURM (combiner avec singularity : `-profile singularity,slurm`) |
| `pbs` | Exécuteur HPC PBS/Torque |
| `test` | Exécution de test rapide |

### macOS Apple Silicon (M1 et supérieur)

Certains paquets Bioconda n'ont pas de version native arm64. Ajouter le profil `arm64` pour forcer l'émulation Rosetta 2 — les outils s'exécutent à une vitesse quasi-native et le pipeline produit des résultats identiques à Linux. Kleborate utilisera automatiquement le mode MLST uniquement sur ARM64 si ses dépendances complètes de pathotype sont indisponibles, mais tous les autres outils fonctionnent à pleine capacité :

```bash
CONDA_SUBDIR=osx-64 nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda,arm64
```

> **`CONDA_SUBDIR=osx-64` est requis** sur Apple Silicon. Les versions modernes de libmamba ignorent le paramètre `--platform` dans Nextflow et reviennent à `osx-arm64` natif, où certains paquets Bioconda sont indisponibles. Préfixer `CONDA_SUBDIR=osx-64` force l'émulation Rosetta 2 de manière fiable. Le coût de construction de l'environnement n'est payé qu'une seule fois ; les environnements sont mis en cache et réutilisés lors des exécutions suivantes.

---

## Exemple HPC

```bash
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile singularity,slurm \
    -c custom.config   # optionnel : remplacer la file d'attente, le compte de projet, etc.
```

---

## Mise à jour de l'esquisse de référence

Pour ajouter des espèces supplémentaires ou mettre à jour les génomes de référence, éditer `assets/build_references.sh` et ajouter l'accession NCBI et l'étiquette. Les noms de fichiers de référence doivent commencer par un préfixe reconnaissable :

| Préfixe | Étiquette d'espèce associée |
|---|---|
| `Ecoli_` ou `Escherichia_` | `E_coli` |
| `Salmonella_` | `Salmonella_enterica` |
| `Shigella_` | `Shigella` |
| `Klebsiella_` | `Klebsiella` |
| `Enterobacter_` | `Enterobacter` |

Après modification, reconstruire l'esquisse :
```bash
bash assets/build_references.sh
```

---

## Vignettes

Exemples de résultats à partir de vrais jeux de données :

| Jeu de données | Espèce | Échantillons | Vignette |
|---|---|---|---|
| Isolats intestinaux de PNH, Gambie (Foster-Nyarko et al. 2020) | *Escherichia coli* | 98 | [Voir la vignette](vignettes/ecoli_vignette.md) |
| Isolats NTS cliniques, Gambie (Darboe et al. 2022) | *Salmonella enterica* | 99 | [Voir la vignette](vignettes/salmonella_vignette.md) |
| *Shigella* diversifiée — les quatre espèces (génomes de référence) | *Shigella* spp. | 15 | [Voir la vignette](vignettes/shigella_vignette.md) |

---

## Citation

Si vous utilisez enteric-typer ou les figures qu'il génère, veuillez citer :

- **enteric-typer** : Foster-Nyarko E et al. github.com/efosternyarko/enteric-typer

Ainsi que les outils qu'il intègre :

- **Mash** : Ondov et al. (2016) Genome Biology 17:132
- **MLST** : Seemann (2016) github.com/tseemann/mlst
- **AMRFinder Plus** : Feldgarden et al. (2021) Scientific Reports 11:12728
- **ECTyper** : Laing et al. (2019) Microbial Genomics 5(12)
- **EzClermont** : Waters et al. (2020) Microbial Genomics 6(9)
- **SISTR** : Yoshida et al. (2016) PLOS ONE 11(1):e0147101
- **PlasmidFinder** : Carattoli et al. (2014) Antimicrobial Agents and Chemotherapy
- **Kaptive** : Wyres et al. (2016) Microbial Genomics 2(10) ; Lam et al. (2022) Nature Protocols
- **EC K-typing G2/G3** : Gladstone et al. (2026) Nature Microbiology ; github.com/rgladstone/EC-K-typing
- **EC K-typing G1/G4** : Foster-Nyarko E et al. github.com/efosternyarko/EC-K-typing-G1G4
- **Kleborate** : Lam et al. (2021) Nature Communications ; github.com/klebgenomics/Kleborate
- **ShigEiFinder** : LanLab (2022) Microbial Genomics ; github.com/LanLab/ShigEiFinder
- **Mykrobe** : Hunt et al. (2015) Genome Biology 16:239 ; Hawkey et al. (2022) Microbial Genomics 8(6)
- **SKA2** : github.com/bacpop/ska.rust
- **IQ-TREE 2** : Minh et al. (2020) Molecular Biology and Evolution 37(5)
- **AMRrules** : github.com/AMRverse/AMRrules
