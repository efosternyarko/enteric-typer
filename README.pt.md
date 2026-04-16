🌐 [English](README.md) | [Français](README.fr.md) | **Português**

# enteric-typer

Um pipeline de genotipagem com portão de espécie para agentes patogénicos entéricos. A partir de uma pasta de montagens genómicas, o pipeline identifica a espécie de cada amostra e executa as ferramentas de tipagem adequadas, gera uma filogenia SNP de genoma completo e produz figuras de resumo prontas para publicação.

## Espécies suportadas

| Espécie | Ferramentas de tipagem |
|---|---|
| *Escherichia coli* | MLST (Achtman), AMRFinder, ECTyper (serotipo), EzClermont (filogrupo Clermont), Kleborate (patotipo), PlasmidFinder, Kaptive (locus K) |
| *Salmonella enterica* | MLST, AMRFinder, SISTR (serovar), PlasmidFinder |
| *Shigella* spp. | MLST (Achtman), AMRFinder, ShigEiFinder (serotipo/espécie), Mykrobe (genotipagem *S. sonnei*), PlasmidFinder, rastreio pINV, rastreio de elementos IS |
| Outra / não classificada | Identificação de espécie apenas (registado e ignorado) |

> **Kleborate:** executa a detecção completa do patotipo em todas as plataformas. No macOS Apple Silicon (ARM64) sem Rosetta 2, muda automaticamente para o modo apenas MLST — todas as outras ferramentas funcionam com plena capacidade independentemente da plataforma (Linux, macOS Intel/ARM64, HPC).

> **Classificação de genes de resistência:** todos os resultados AMRFinder são classificados pelo [AMRrules](https://github.com/AMRverse/AMRrules) em genes de resistência *adquirida* e genes *intrínsecos* (característicos da espécie selvagem). Os genes intrínsecos são mantidos nos ficheiros TSV de resultados como referência, mas são **excluídos de todos os gráficos AMR** para que as figuras reflictam apenas a resistência adquirida clinicamente relevante.

## Visão geral do pipeline

```
Montagens de entrada (pasta ou folha de amostras)
        │
        ▼
┌───────────────────────────────────────────────────────────────────────────┐
│  1. VERIFICAÇÃO DE ESPÉCIE E CONTROLO DE QUALIDADE                        │
│                                                                           │
│  Identificação de espécie — distância Mash em relação a um esboço de     │
│  referência de 13 genomas; a referência mais próxima vence (mesma         │
│  abordagem que o Kleborate). Os isolados de Shigella são correctamente    │
│  identificados porque estão mais próximos das referências Shigella do que │
│  de qualquer referência E. coli.                                          │
│                                                                           │
│  Filtros de controlo de qualidade das montagens (falhas registadas e      │
│  excluídas da tipagem)                                                    │
│  ──────────────────────────────────────────────────────────────────────   │
│  Tamanho do genoma  E. coli / Shigella   4,3 – 5,9 Mb  [1]               │
│                     Salmonella           4,1 – 6,6 Mb  [2]               │
│  Contaminação       Espécie secundária Kraken2 < 3 % dos contigs totais   │
│                     (opcional — fornecer --kraken2_db; ignorado se ausente)│
└───────────────────────────────────────────────────────────────────────────┘
        │
   ┌────┼────────┐
   ▼    ▼        ▼
E. coli  Salmonella  Shigella   (outras espécies registadas e ignoradas)
   │         │           │
   ▼         ▼           ▼
┌───────────────────────────────────────────────────────────────────────────┐
│  2. TIPAGEM ESPECÍFICA DA ESPÉCIE  (todas as ferramentas em paralelo)     │
│                                                                           │
│  E. coli                       Salmonella          Shigella               │
│  ──────────────────────────    ──────────────────  ──────────────────     │
│  MLST (achtman_4)              MLST (senterica_    MLST (achtman_4)       │
│  AMRFinder                       achtman_2)        AMRFinder              │
│  ECTyper (serotipo O:H)        AMRFinder           ShigEiFinder           │
│  EzClermont (filogrupo)        SISTR (serovar)       (serotipo/espécie)   │
│  Kleborate (patotipo)          PlasmidFinder       Mykrobe                │
│  PlasmidFinder                                       (genótipo S. sonnei) │
│  Kaptive locus K (G2/G3                            PlasmidFinder          │
│    → G1/G4 se não tipável)                         rastreio pINV          │
│                                                    rastreio elementos IS   │
└───────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  3. FILOGENÉTICA  (por espécie, ≥ 3 amostras; ignorar com                │
│                  --skip_local_phylo)                                     │
│  SKA2 build (k=31) → alinhamento SNP de genoma completo + matriz de      │
│  distâncias SNP → árvore ML IQ-TREE (modelo GTR+G;                       │
│  substituir com --iqtree_model MFP)                                       │
└──────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────┐
│  4. AGREGAÇÃO  (um TSV por espécie)                          │
│  ecoli_typer_results.tsv                                     │
│  salmonella_typer_results.tsv                                │
│  shigella_typer_results.tsv                                  │
│  Resultados AMRFinder classificados pelo AMRrules em:        │
│    amrfinder_acquired_genes  — genes de resistência adquirida│
│                                clinicamente relevantes        │
│    amrfinder_intrinsic_genes — genes intrínsecos da espécie  │
│                                (sinalizados, mantidos no TSV) │
└──────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  5. FIGURAS DE RESUMO                                                     │
│  Os genes de resistência intrínsecos (classificados pelo AMRrules) são   │
│  excluídos de todas as figuras AMR — apenas os genes adquiridos são       │
│  representados.                                                           │
│                                                                           │
│  Todas as espécies                                                        │
│  ─────────────────────────────────────────────────────────────────────   │
│  Fig 1  Resumo da população (4 painéis):                                  │
│           A — Distribuição dos tipos de sequência MLST                    │
│           B — Barras de serotipo / paisagem de elementos IS (Shigella)    │
│           C — Prevalência por classe de antibióticos                      │
│           D — Carga de multirresistência por isolado                      │
│  Fig 2  Árvore filogenética + faixas ST/filogrupo + mapa de calor AMR    │
│  Fig 3  Principais genes AMR adquiridos (genes intrínsecos excluídos)     │
│  Fig 4  Tipos de replicões plasmídicos                                    │
│  Fig 5  Genes de virulência / patotipo                                    │
│  Fig 6  Prevalência de classes AMR por tipo de sequência MLST             │
│  Fig 7  Prevalência de classes AMR por serovar / filogrupo Clermont       │
│  (SNP)  Mapa de calor de distâncias SNP de genoma completo por pares      │
│                                                                           │
│  Apenas Shigella                                                          │
│  ─────────────────────────────────────────────────────────────────────   │
│  Fig 8  Composição de espécies + distribuição de serotipos (barras empi.) │
│  Fig 9  Painel de virulência e invasão (ipaH, genes pINV, elementos IS)  │
│                                                                           │
│  Controlo de qualidade das montagens (todas as espécies)                  │
│  ─────────────────────────────────────────────────────────────────────   │
│  Métricas de montagem  Figura de 8 painéis por espécie (PDF + PNG):      │
│           A — Histograma comprimento do genoma   B — Caixa comprimento    │
│           C — Histograma N50                     D — Caixa N50            │
│           E — Histograma número de contigs       F — Caixa contigs        │
│           G — Histograma GC%                     H — Caixa GC%            │
└──────────────────────────────────────────────────────────────────────────┘
```

> **Tipagem do locus K de *E. coli* (Kaptive):** a tipagem do locus K utiliza duas bases de dados curadas executadas sequencialmente. As montagens são primeiro tipadas na **base de dados G2/G3** (`EC-K-typing_group2and3_v3.0.0.gbk`; [Gladstone et al. 2026](https://www.nature.com/articles/s41564-026-02283-w)). As amostras que permanecem não tipáveis são depois retipadas na **base de dados G1/G4** (`EC-K-typing_group1and4_v1.2.gbk`; [Foster-Nyarko et al.](https://github.com/efosternyarko/EC-K-typing-G1G4)). Os loci G2/G3 e G1/G4 são mutuamente exclusivos em *E. coli*, pelo que a tipagem sequencial garante que cada amostra é atribuída ao grupo correcto sem dupla contagem.

> **Filogenética (SKA2 + IQ-TREE):** as filogenias SNP de genoma completo são construídas com [SKA2](https://github.com/bacpop/ska.rust) (alinhamento por k-meros divididos, k=31), que alinha montagens sem um genoma de referência e produz um alinhamento SNP à escala do genoma e uma matriz de distâncias SNP por pares. O alinhamento é depois passado ao IQ-TREE 2 para a inferência da árvore por máxima verosimilhança (por defeito: GTR+G; usar `--iqtree_model MFP` para activar o ModelFinder Plus). A filogenética é executada por espécie quando estão presentes ≥ 3 amostras; ignorar com `--skip_local_phylo` para execuções mais rápidas.

---

## Instalação

### Passo 1 — Clonar o repositório

```bash
git clone https://github.com/efosternyarko/enteric-typer
cd enteric-typer
```

---

### Passo 2 — Instalar o Java

O Nextflow requer o Java 17 ou posterior (Java 21 recomendado).

**macOS**
```bash
brew install --cask temurin
```
> Se o Homebrew ainda não estiver instalado, executar:
> ```bash
> /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
> ```

**Linux (Debian/Ubuntu)**
```bash
sudo apt update && sudo apt install -y default-jdk
```

**Windows (PowerShell)**
```powershell
winget install EclipseAdoptium.Temurin.21.JDK   # Java 21 (recomendado)
# ou
winget install EclipseAdoptium.Temurin.25.JDK   # Java 25 (mais recente)
```

Verificar: `java -version`

---

### Passo 3 — Instalar o Nextflow

```bash
curl -s https://get.nextflow.io | bash
```

O instalador cria um executável `nextflow` no directório actual. Movê-lo para uma localização no `$PATH` para que possa ser executado a partir de qualquer lado:

**macOS**
```bash
# Se o Homebrew estiver instalado, mover para o directório bin do Homebrew (já no $PATH):
mv nextflow /opt/homebrew/bin/        # Apple Silicon (M1/M2/M3)
# ou
mv nextflow /usr/local/bin/           # Mac Intel

# Se o Homebrew não estiver instalado, criar ~/bin e adicioná-lo ao PATH:
mkdir -p ~/bin
mv nextflow ~/bin/
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

**Linux**
```bash
mv nextflow ~/bin/
# Se ~/bin/ ainda não existir:
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
Depois reiniciar o PowerShell.

Verificar: `nextflow -version`

---

### Passo 4 — Instalar conda + mamba

Se o conda ainda não estiver instalado, a forma mais rápida é o **Miniforge**, que inclui tanto o conda como o mamba:

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

Quando o instalador perguntar:
```
Do you wish the installer to initialize Miniforge3
by running conda init? [yes|no]
```
Responder **yes**. Depois recarregar o shell:

```bash
source ~/.zshrc    # macOS (zsh)
# ou
source ~/.bashrc   # Linux (bash)
```

**Windows**
Descarregar o [instalador Miniforge3 para Windows](https://github.com/conda-forge/miniforge/releases/latest) e executá-lo. Quando solicitado, seleccionar *Add Miniforge to my PATH*. Depois reinicializar:
```powershell
conda init powershell
```
Reiniciar o PowerShell.

> **Já tem conda mas não mamba?**
> ```bash
> conda install -n base -c conda-forge mamba
> ```

---

### Passo 5 — Instalar ncbi-datasets-cli e mash

Ambos são necessários para o `build_references.sh` descarregar genomas de referência e construir o esboço de espécie:

```bash
mamba install -c conda-forge ncbi-datasets-cli
mamba install -c bioconda mash
```

---

### Passo 6 — Construir o esboço de referência Mash

Este passo descarrega 13 genomas de referência do NCBI e constrói o esboço Mash utilizado para a identificação de espécies. Executar uma vez antes da primeira execução do pipeline:

```bash
bash assets/build_references.sh
```

Saída: `assets/enteric_species_refs.msh`
Duração estimada: 1–5 minutos dependendo da largura de banda.

Genomas de referência incluídos:

| Espécie / filogrupo | Estirpe | Acessão |
|---|---|---|
| *E. coli* filogrupo A | K-12 MG1655 | GCF_000005845.2 |
| *E. coli* filogrupo B1 | SE11 | GCF_000010485.1 |
| *E. coli* filogrupo B2 | CFT073 | GCF_000007445.1 |
| *E. coli* filogrupo D | UMN026 | GCF_000026265.1 |
| *E. coli* filogrupo E | O157:H7 EDL933 | GCF_000006665.1 |
| *S. enterica* Typhimurium | LT2 | GCF_000006945.2 |
| *S. enterica* Typhi | CT18 | GCF_000195995.1 |
| *S. enterica* Enteritidis | P125109 | GCF_000009505.1 |
| *S. sonnei* | ATCC 29930 | GCF_002950395.1 |
| *S. flexneri* | 2a 2457T | GCF_000007405.1 |
| *S. boydii* | Sb227 | GCF_000012025.1 |
| *S. dysenteriae* | Sd197 | GCF_000012005.1 |
| *K. pneumoniae* | HS11286 | GCF_000240185.1 |

> A identificação de espécies utiliza a abordagem **a referência mais próxima vence** (mesma abordagem que o Kleborate): o grupo de espécies com a menor distância Mash em relação a qualquer referência no esboço é atribuído. Os isolados genuínos de *Shigella* são correctamente identificados porque a sua referência mais próxima é sempre um genoma de *Shigella* (distância Mash < 0,012), enquanto os isolados *E. coli* do filogrupo B2 — filogeneticamente próximos de *Shigella* mas não sendo *Shigella* — obtêm distâncias menores para as referências *E. coli* (~0,007) do que para qualquer referência *Shigella* (≥ 0,014) e são correctamente classificados como *E. coli*.

---

### Passo 7 — Descarregar a base de dados Kraken2 (opcional)

O rastreio de contaminação é **opcional** e é ignorado se `--kraken2_db` não for fornecido. Para o activar, descarregar a base de dados `k2_standard_08gb` (~8 GB) uma vez e indicar o seu caminho ao pipeline no momento da execução:

```bash
# Criar um directório para a base de dados
mkdir -p ~/kraken2_db

# Descarregar e extrair (requer ~8 GB de espaço em disco)
cd ~/kraken2_db
curl -L -O https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz
tar -xzf k2_standard_08gb_20240904.tar.gz
```

Depois indicar o caminho da base de dados ao executar o pipeline:

```bash
nextflow run main.nf \
    -profile conda \
    --input_dir assemblies/ \
    --outdir    results/ \
    --kraken2_db ~/kraken2_db/
```

A base de dados pode ser reutilizada entre projectos. Consultar a [página de índices Kraken2](https://benlangmead.github.io/aws-indexes/k2) para a versão mais recente.

> **Os filtros de tamanho de genoma** não requerem configuração — os limiares são aplicados automaticamente com base na classificação da espécie. Os valores por defeito podem ser substituídos com, por exemplo, `--ecoli_min_length 4000000`.
>
> Os limiares baseiam-se em critérios de controlo de qualidade publicados:
> \[1\] *E. coli* / *Shigella* — [Critérios de qualidade de genomas *E. coli* do BIGSdb](https://bigsdb.pasteur.fr/ecoli/genomes-quality-criteria/)
> \[2\] *Salmonella* — [Recomendações do consórcio PATH-SAFE para a vigilância genómica de doenças de origem alimentar](https://science.food.gov.uk/article/143833-path-safe-consortium-recommendations-for-genomic-surveillance-of-foodborne-diseases-using-salmonella-as-an-exemplar?attachment_id=300637) (Tabela 3)

---

## Início rápido

> **Linux / Ubuntu / Debian:** utilizar os comandos Linux abaixo.
> **Windows:** recomendamos executar o enteric-typer dentro do [WSL2](https://learn.microsoft.com/pt-pt/windows/wsl/install) (Subsistema Windows para Linux), que fornece um ambiente Ubuntu completo. Instalar o WSL2 e depois seguir as instruções Linux. Em alternativa, utilizar o perfil `docker` com o Docker Desktop para Windows.

### Para uma pasta de montagens

```bash
# Linux / Ubuntu / Debian / Mac Intel
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda

# Mac Apple Silicon (M1 e superior) — prefixar CONDA_SUBDIR e adicionar perfil arm64
CONDA_SUBDIR=osx-64 nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda,arm64

# Windows (Docker Desktop) — conda não necessário
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile docker
```

### Ou com uma folha de amostras

```bash
# Gerar automaticamente uma folha de amostras a partir de uma pasta
python bin/make_samplesheet.py \
    --input /path/to/assemblies/ \
    --output samples.csv

# Executar com a folha de amostras
nextflow run main.nf \
    --samplesheet samples.csv \
    --outdir      results/ \
    -profile conda
```

### Ignorar a filogenética local (mais rápido, sem SKA2/IQ-TREE)

```bash
nextflow run main.nf \
    --input_dir        /path/to/assemblies/ \
    --outdir           results/ \
    --skip_local_phylo \
    -profile conda
```

> Quando `--skip_local_phylo` está activo, a matriz de distâncias SNP, o mapa de calor SNP e as figuras de anotação da árvore não são produzidos.

---

## Parâmetros

| Parâmetro | Defeito | Descrição |
|---|---|---|
| `--input_dir` | `null` | Pasta de montagens FASTA (`.fasta/.fa/.fna/.fas`) |
| `--samplesheet` | `null` | CSV com as colunas `id,fasta` |
| `--outdir` | `results` | Directório de saída |
| `--skip_local_phylo` | `false` | Ignorar SKA2 + IQ-TREE (sem árvore, sem matriz SNP/mapa de calor) |
| `--ska2_min_samples` | `3` | Número mínimo de amostras para tentar SKA2/IQ-TREE |
| `--iqtree_model` | `GTR+G` | Modelo de substituição IQ-TREE (usar `MFP` para selecção automática de modelo) |
| `--iqtree_bootstraps` | `100` | Réplicas bootstrap ultra-rápidas IQ-TREE |
| **Controlo de qualidade das montagens** | | |
| `--ecoli_min_length` | `4300000` | Comprimento mínimo de montagem *E. coli* (pb) |
| `--ecoli_max_length` | `5900000` | Comprimento máximo de montagem *E. coli* (pb) |
| `--shigella_min_length` | `4300000` | Comprimento mínimo de montagem *Shigella* (pb) |
| `--shigella_max_length` | `5900000` | Comprimento máximo de montagem *Shigella* (pb) |
| `--salmonella_min_length` | `4100000` | Comprimento mínimo de montagem *Salmonella* (pb) |
| `--salmonella_max_length` | `6600000` | Comprimento máximo de montagem *Salmonella* (pb) |
| `--kraken2_db` | `null` | Caminho para o directório da base de dados Kraken2; rastreio ignorado se ausente |
| `--max_contamination` | `3.0` | Percentagem máxima de espécie secundária permitida (Kraken2) |

---

## Ficheiros de saída

```
results/
├── pipeline_info/
│   ├── timeline.html
│   ├── report.html
│   └── dag.svg
│
├── species_check/
│   └── *_species.txt              ← melhor espécie + distância Mash por amostra
│
├── mlst_ecoli/
│   └── *_ecoli_achtman_4_mlst.tsv
├── amrfinder_ecoli/
│   └── *_amrfinder.tsv
├── ectyper/
│   └── *_ectyper.tsv
├── kaptive_g2g3/ kaptive_g1g4/
│   └── *_ktype.tsv                ← grupo locus K, locus, tipo, nível de confiança
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
│   └── *_shigeifinder.tsv        ← ipaH, plasmídeo de virulência, cluster, serotipo, antigénios O/H
├── mykrobe/
│   └── *_mykrobe_parsed.tsv      ← linhagem/genótipo S. sonnei (NA para outras Shigella spp.)
├── pinv_screen/
│   └── *_pinv.tsv                ← presença de genes marcadores pINV (icsA/virG, virF, ipaB/C/D)
├── is_screen/
│   └── *_is.tsv                  ← contagens de cópias de elementos IS (IS1, IS30, IS600, IS629…)
├── plasmidfinder_shigella/
│   └── *_plasmidfinder.tsv
│
├── ska2_{species}/
│   ├── ska2_alignment.fasta       ← alinhamento SNP de genoma completo (entrada para IQ-TREE)
│   └── snp_matrix.tsv             ← matriz de distâncias SNP de genoma completo por pares
│
├── ecoli_typer_results.tsv         ← Tabela de resultados principal (E. coli)
├── salmonella_typer_results.tsv    ← Tabela de resultados principal (Salmonella)
├── shigella_typer_results.tsv      ← Tabela de resultados principal (Shigella)
│
│   Todas as tabelas principais incluem:
│     amrfinder_acquired_genes   — genes de resistência (intrínsecos excluídos)
│     amrfinder_intrinsic_genes  — genes intrínsecos selvagens (sinalizados, não representados)
│     amrfinder_genes            — todos os resultados AMRFinder em bruto
│
│ ── Figuras de resumo (PDF + PNG) ───────────────────────────────────────
│
├── {species}_fig1_population_summary.{pdf,png}
│     Painel A: Distribuição dos tipos de sequência MLST
│               E. coli — barras empilhadas por filogrupo Clermont
│               Salmonella — gradiente monocromático
│               Shigella — barras empilhadas por espécie (S. sonnei / flexneri / boydii / dysenteriae)
│     Painel B: Prevalência de serotipos/serovares (E. coli / Salmonella)
│               Mapa de calor do número de cópias de elementos IS (Shigella)
│     Painel C: Prevalência por classe de antibióticos
│     Painel D: Carga de multirresistência por isolado (≥ 3 classes adquiridas = MDR)
│
├── {species}_fig2_tree_amr.{pdf,png}
│     Árvore filogenética (ML IQ-TREE) com:
│       — faixa de cor ST
│       — faixa de cor do filogrupo Clermont (apenas E. coli)
│       — mapa de calor binário de genes de virulência
│       — mapa de calor de genes AMR agrupados por classe de antibióticos
│
├── {species}_fig3_amr_genes.{pdf,png}      ← Frequências de genes AMR adquiridos
├── {species}_fig4_plasmid_replicons.{pdf,png}
├── {species}_fig5_virulence.{pdf,png}       ← Genes de virulência / patotipo
│
├── {species}_fig6_amr_by_st.{pdf,png}
│     Mapa de calor AMRnet: % de isolados com cada classe, por ST MLST
│
├── {species}_fig7_amr_by_group.{pdf,png}
│     Mapa de calor AMRnet: % de isolados com cada classe, por
│     serovar (Salmonella) ou filogrupo Clermont (E. coli) ou serotipo (Shigella)
│
├── {species}_snp_heatmap.{pdf,png}          ← Mapa de calor de distâncias SNP por pares
│
│ ── Figuras específicas de Shigella ─────────────────────────────────────
│
├── shigella_fig8_shigella_serotypes.{pdf,png}
│     Composição de espécies + distribuição de serotipos (gráfico de barras empilhadas)
│
├── shigella_fig9_shigella_features.{pdf,png}
│     Mapa de calor binário: ipaH · plasmídeo de virulência · genes de invasão pINV
│     (icsA/virG, virF, virB, ipaB, ipaC, ipaD) · elementos IS por amostra
│
│ ── Figuras de controlo de qualidade das montagens (PDF + PNG) ──────────
│
├── {species}_assembly_metrics.{pdf,png}
│     Figura CQ de montagem com 8 painéis (uma por espécie detectada):
│       A  Histograma comprimento do genoma    B  Caixa comprimento do genoma
│       C  Histograma N50                      D  Caixa N50
│       E  Histograma número de contigs        F  Caixa número de contigs
│       G  Histograma GC%                      H  Caixa GC%
│
└── {species}_assembly_metrics_summary.tsv
      Estatísticas de montagem por amostra agregadas (genome_length, num_contigs,
      assembly_N50, gc_pct) para todas as amostras dessa espécie.
```

---

## Abreviaturas das classes de antibióticos

As figuras 6 e 7 utilizam abreviaturas curtas para as classes de antibióticos.
Os nomes completos e as notas clínicas são apresentados abaixo.

| Abreviatura | Classe de antibióticos | Agentes representativos |
|---|---|---|
| **AMG** | Aminoglicosídeo | Gentamicina, amicacina, tobramicina, estreptomicina |
| **BLA** | Beta-lactâmico | Ampicilina, cefalosporinas, carbapenemes, penicilinas |
| **COL** | Colistina | Colistina (polimixina E), polimixina B |
| **FOS** | Fosfomicina | Fosfomicina |
| **FOSM** | Fosmidomicina | Fosmidomicina |
| **LIN** | Lincosamida | Lincomicina, clindamicina |
| **MAC** | Macrólido | Azitromicina, eritromicina |
| **NIT** | Nitrofurano | Nitrofurantoína |
| **PHE** | Fenicol | Cloranfenicol, florfenicol |
| **QNL** | Quinolona | Ciprofloxacina, ácido nalidíxico, levofloxacina (todas as fluoroquinolonas) |
| **SGM** | Estreptogramina | Quinupristina–dalfopristina (combinações estreptogramina A/B) |
| **STR** | Estreptotricina | Nourseotricina |
| **SUL** | Sulfonamida | Sulfametoxazol, sulfisoxazol |
| **TET** | Tetraciclina | Tetraciclina, doxiciclina, tigeciclina |
| **TMP** | Trimetoprim | Trimetoprim (frequentemente combinado com sulfonamida como cotrimoxazol) |

> **Definição de MDR:** um isolado é classificado como multirresistente (MDR) quando possui genes de resistência adquirida em **≥ 3** das classes acima. A classe EFFLUX é excluída da contagem MDR, uma vez que as bombas de efluxo quase universais são intrínsecas à espécie e removidas pelo AMRrules antes da produção das figuras.

---

## Limpeza após uma execução

Todos os resultados finais são escritos em `--outdir` (por defeito: `results/`). Após verificar os resultados, pode libertar espaço em disco removendo os ficheiros temporários do Nextflow:

```bash
rm -rf work/ .nextflow/ .nextflow.log*
```

> **Sugestão:** manter `work/` se quiser utilizar `-resume` para executar novamente com parâmetros diferentes sem repetir os passos já concluídos. Apagar apenas quando a execução estiver finalizada.
>
> **Nota sobre `-resume` e anotação da árvore:** se uma execução anterior não produziu uma figura de anotação (por exemplo, devido a uma actualização de software), o Nextflow pode armazenar em cache o estado de falha. Executar uma vez **sem `-resume`** para forçar uma nova execução limpa do passo de anotação.

---

## Opcional: instalar o Graphviz

O Nextflow gera um diagrama DAG de execução (`pipeline_info/dag.svg`) que requer o Graphviz para ser renderizado. Sem ele, aparece um aviso inofensivo — o pipeline continua a funcionar normalmente. Para suprimir o aviso:

```bash
# macOS
brew install graphviz

# Linux
sudo apt install graphviz

# conda
conda install -c conda-forge graphviz
```

## Perfis de execução

| Perfil | Caso de utilização |
|---|---|
| `conda` | Estação de trabalho local com conda |
| `mamba` | Igual ao conda mas com resolução de ambiente mais rápida |
| `arm64` | **Adicionar no Apple Silicon (M1 e superior)** — força ambientes conda osx-64 via Rosetta 2 |
| `docker` | Local com Docker Desktop |
| `singularity` | Cluster HPC com Singularity/Apptainer |
| `slurm` | Executor HPC SLURM (combinar com singularity: `-profile singularity,slurm`) |
| `pbs` | Executor HPC PBS/Torque |
| `test` | Execução de teste rápida |

### macOS Apple Silicon (M1 e superior)

Alguns pacotes Bioconda não têm compilação nativa arm64. Adicionar o perfil `arm64` para forçar a emulação Rosetta 2 — as ferramentas funcionam a velocidade quase-nativa e o pipeline produz resultados idênticos ao Linux. O Kleborate utilizará automaticamente o modo apenas MLST no ARM64 se as suas dependências completas de patotipo estiverem indisponíveis, mas todas as outras ferramentas funcionam com plena capacidade:

```bash
CONDA_SUBDIR=osx-64 nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile conda,arm64
```

> **`CONDA_SUBDIR=osx-64` é necessário** no Apple Silicon. As versões modernas do libmamba ignoram o parâmetro `--platform` dentro do Nextflow e recorrem ao `osx-arm64` nativo, onde alguns pacotes Bioconda estão indisponíveis. Prefixar `CONDA_SUBDIR=osx-64` força a emulação Rosetta 2 de forma fiável. O custo de construção do ambiente é pago apenas uma vez; os ambientes são armazenados em cache e reutilizados nas execuções seguintes.

---

## Exemplo HPC

```bash
nextflow run main.nf \
    --input_dir /path/to/assemblies/ \
    --outdir    results/ \
    -profile singularity,slurm \
    -c custom.config   # opcional: substituir fila, conta de projecto, etc.
```

---

## Actualizar o esboço de referência

Para adicionar espécies adicionais ou actualizar os genomas de referência, editar `assets/build_references.sh` e adicionar a acessão NCBI e a etiqueta. Os nomes dos ficheiros de referência devem começar com um prefixo reconhecível:

| Prefixo | Etiqueta de espécie associada |
|---|---|
| `Ecoli_` ou `Escherichia_` | `E_coli` |
| `Salmonella_` | `Salmonella_enterica` |
| `Shigella_` | `Shigella` |
| `Klebsiella_` | `Klebsiella` |
| `Enterobacter_` | `Enterobacter` |

Após editar, reconstruir o esboço:
```bash
bash assets/build_references.sh
```

---

## Vinhetas

Exemplos de resultados a partir de conjuntos de dados reais:

| Conjunto de dados | Espécie | Amostras | Vinheta |
|---|---|---|---|
| Isolados intestinais de PNH, Gâmbia (Foster-Nyarko et al. 2020) | *Escherichia coli* | 98 | [Ver vinheta](vignettes/ecoli_vignette.md) |
| Isolados NTS clínicos, Gâmbia (Darboe et al. 2022) | *Salmonella enterica* | 99 | [Ver vinheta](vignettes/salmonella_vignette.md) |
| *Shigella* diversificada — as quatro espécies (genomas de referência) | *Shigella* spp. | 15 | [Ver vinheta](vignettes/shigella_vignette.md) |

---

## Citação

Se utilizar o enteric-typer ou as figuras por ele geradas, por favor cite:

- **enteric-typer**: Foster-Nyarko E et al. github.com/efosternyarko/enteric-typer

E as ferramentas que integra:

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
