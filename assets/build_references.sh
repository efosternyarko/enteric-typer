#!/usr/bin/env bash
# build_references.sh
# ─────────────────────────────────────────────────────────────────────────────
# Downloads reference genome assemblies from NCBI and builds the Mash sketch
# used by enteric-typer for species identification.
#
# Reference genomes included:
#   E. coli K-12 MG1655                  GCF_000005845.2
#   E. coli O157:H7 EDL933               GCF_000006665.3
#   Salmonella Typhimurium LT2           GCF_000006945.2
#   Salmonella Typhi CT18                GCF_000195995.1
#   Salmonella Enteritidis P125109       GCF_000009505.1
#   Shigella sonnei Ss046                GCF_000006925.2  (to flag Shigella)
#   Klebsiella pneumoniae HS11286        GCF_000240185.1  (common enteric)
#
# Requirements:
#   ncbi-datasets-cli  (mamba install -c conda-forge ncbi-datasets-cli)
#   mash               (mamba install -c bioconda mash)
#
# Usage (run from the repo root):
#   bash assets/build_references.sh
#
# Output:
#   assets/enteric_species_refs.msh   (Mash sketch)
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFS_DIR="${SCRIPT_DIR}/reference_genomes"
SKETCH="${SCRIPT_DIR}/enteric_species_refs.msh"
KMER_SIZE=21
SKETCH_SIZE=1000

# ── Dependency checks ─────────────────────────────────────────────────────────
if ! command -v datasets &>/dev/null; then
    echo "ERROR: 'datasets' (ncbi-datasets-cli) is not installed."
    echo "Install with:  mamba install -c conda-forge ncbi-datasets-cli"
    exit 1
fi

if ! command -v mash &>/dev/null; then
    echo "ERROR: 'mash' is not installed."
    echo "Install with:  mamba install -c bioconda mash"
    exit 1
fi

mkdir -p "${REFS_DIR}"
cd "${REFS_DIR}"

FAILED=0

# ── Download function ─────────────────────────────────────────────────────────
download_genome() {
    local accession="$1"
    local label="$2"
    local outfile="${label}.fasta"

    if [ -f "${outfile}" ] && [ -s "${outfile}" ]; then
        echo "  [SKIP] ${label} already present"
        return 0
    fi

    echo "  Downloading ${label} (${accession})..."

    # Download ZIP from NCBI datasets
    if ! datasets download genome accession "${accession}" \
            --include genome \
            --filename "${accession}.zip" 2>/dev/null; then
        echo "  WARNING: datasets download failed for ${accession} — skipping"
        FAILED=$((FAILED + 1))
        return 1
    fi

    # Extract all .fna files from the ZIP into a temp directory
    local tmpdir="${accession}_tmp"
    mkdir -p "${tmpdir}"
    unzip -o "${accession}.zip" -d "${tmpdir}" >/dev/null 2>&1

    # Find the assembly FASTA (exclude *_from_genomic.fna CDS files)
    local fasta
    fasta=$(find "${tmpdir}" -name "*.fna" \
            ! -name "*_cds_from_*" \
            ! -name "*_rna_from_*" \
            | head -1)

    if [ -n "${fasta}" ] && [ -s "${fasta}" ]; then
        cp "${fasta}" "${outfile}"
        echo "    → ${outfile} ($(grep -c '^>' "${outfile}") sequences)"
    else
        echo "  WARNING: could not extract FASTA for ${accession}"
        FAILED=$((FAILED + 1))
    fi

    rm -rf "${tmpdir}" "${accession}.zip"
}

# ── Reference genome downloads ────────────────────────────────────────────────
echo "Downloading reference genomes from NCBI..."
echo ""

# E. coli references — prefix 'Ecoli_' → species_check maps to 'E_coli'
download_genome GCF_000005845.2  "Ecoli_K12_MG1655"
download_genome GCF_000006665.3  "Ecoli_O157H7_EDL933"

# Salmonella references — prefix 'Salmonella_' → 'Salmonella_enterica'
download_genome GCF_000006945.2  "Salmonella_Typhimurium_LT2"
download_genome GCF_000195995.1  "Salmonella_Typhi_CT18"
download_genome GCF_000009505.1  "Salmonella_Enteritidis_P125109"

# Shigella — rejected by species gate but identified in logs
download_genome GCF_000006925.2  "Shigella_sonnei_Ss046"

# Klebsiella — rejected by species gate but identified in logs
download_genome GCF_000240185.1  "Klebsiella_pneumoniae_HS11286"

# ── Build Mash sketch ─────────────────────────────────────────────────────────
echo ""
echo "Collecting FASTA files for sketch..."
FASTAS=()
for f in Ecoli_*.fasta Salmonella_*.fasta Shigella_*.fasta Klebsiella_*.fasta; do
    [ -s "${f}" ] && FASTAS+=("${f}")
done

N="${#FASTAS[@]}"
if [ "${N}" -eq 0 ]; then
    echo "ERROR: No reference genomes downloaded successfully. Cannot build sketch."
    exit 1
fi

echo "  ${N} reference genome(s) found — building Mash sketch..."
echo "  (k=${KMER_SIZE}, sketch_size=${SKETCH_SIZE})"

mash sketch \
    -k "${KMER_SIZE}" \
    -s "${SKETCH_SIZE}" \
    -o "${SKETCH%.msh}" \
    "${FASTAS[@]}"

echo ""
if [ -f "${SKETCH}" ]; then
    echo "Done! Reference sketch written to:"
    echo "  ${SKETCH}"
    echo ""
    echo "Verify with:"
    echo "  mash info ${SKETCH}"
    if [ "${FAILED}" -gt 0 ]; then
        echo ""
        echo "NOTE: ${FAILED} genome(s) failed to download. The sketch was built"
        echo "with the remaining ${N} genome(s). Re-run to retry failed downloads."
    fi
else
    echo "ERROR: Mash sketch was not created. Check the output above for errors."
    exit 1
fi
