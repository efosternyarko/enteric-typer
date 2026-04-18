# Species Detection

enteric-typer identifies the species of every input assembly before any typing tool runs. All downstream tool selection, QC thresholds, and figure labelling depend on this step.

---

## Method: Mash distance to a reference sketch

Species assignment uses [Mash](https://github.com/marbl/Mash) — a MinHash-based tool that estimates genomic distance extremely rapidly without alignment. Each input assembly is compared against a pre-built sketch containing 13 representative reference genomes (one per species or phylogroup).

The principle is **nearest-reference wins**: whichever reference genome has the lowest Mash distance to the query assembly determines the species call. This is the same approach used by [Kleborate](https://github.com/klebgenomics/Kleborate) for *K. pneumoniae* species complex assignment.

### Reference genomes in the sketch

| Species / phylogroup | Reference strain | NCBI accession |
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

Multiple *E. coli* phylogroup references are included because phylogroups B2 and D are genomically closer to *Shigella* than phylogroup A or B1 isolates are. Without phylogroup representation, B2 isolates might be miscalled as Shigella.

---

## Decision logic

```
For each input assembly:
  Compute Mash distance to every reference in the sketch
       │
       ▼
  Find minimum distance reference
       │
       ├── distance ≤ 0.05 and reference is E. coli phylogroup?
       │        └── call: E_coli
       │
       ├── distance ≤ 0.05 and reference is Shigella spp.?
       │        └── call: Shigella
       │
       ├── distance ≤ 0.05 and reference is Salmonella?
       │        └── call: Salmonella_enterica
       │
       ├── distance ≤ 0.05 and reference is K. pneumoniae?
       │        └── call: Klebsiella  (logged only, not typed)
       │
       └── distance > 0.05 or no confident match?
                └── call: Unknown  (logged, skipped)
```

The output for each sample is written to `results/species_check/<sample>_species.txt` and contains the called species label and the Mash distance to the winning reference.

---

## Why Shigella is correctly distinguished from E. coli

*Shigella* species are phylogenetically nested within the *E. coli* species complex — they are, genomically, *E. coli* that have acquired a virulence plasmid and undergone pathoadaptive evolution. This creates a species assignment challenge.

The pipeline handles this correctly because the sketch contains dedicated *Shigella* references from all four species (*S. sonnei*, *S. flexneri*, *S. boydii*, *S. dysenteriae*). Genuine *Shigella* isolates have their lowest Mash distance to one of these four references, not to any *E. coli* reference — the unique accessory genome content (pINV virulence plasmid, IS element repertoire, loss of genes like *cadA*) shifts the overall genomic distance enough that nearest-reference wins correctly.

*E. coli* phylogroup B2 isolates — which are the *E. coli* most genomically similar to *Shigella* — are still called *E. coli* because they retain sufficient *E. coli*-specific core genome content to be slightly closer to the CFT073 (B2) reference than to any *Shigella* reference.

---

## What happens to unrecognised species

Assemblies with no confident species match (Mash distance > 0.05 to all references, or nearest reference is *K. pneumoniae*) are:

1. Written to `results/species_check/` with their best-match label and distance
2. Excluded from all typing steps
3. Still counted in sample totals and noted in the pipeline log

They do not cause the pipeline to fail. This is intentional — mixed-species datasets (e.g. a batch that accidentally includes a *Klebsiella* sample) run to completion for the recognised species without manual intervention.

---

## Updating the reference sketch

To add a new species or update a reference genome, edit `assets/build_references.sh` and add the NCBI accession. Reference FASTA filenames must begin with a recognised prefix that maps to an internal species label:

| Filename prefix | Species label |
|---|---|
| `Ecoli_` or `Escherichia_` | `E_coli` |
| `Salmonella_` | `Salmonella_enterica` |
| `Shigella_` | `Shigella` |
| `Klebsiella_` | `Klebsiella` |

Then rebuild the sketch: `bash assets/build_references.sh`
