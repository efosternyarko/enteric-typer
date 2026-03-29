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
# Requirements: ncbi-datasets-cli (conda install -c conda-forge ncbi-datasets-cli)
#               mash (conda install -c bioconda mash)
#               OR: wget + gunzip
#
# Usage:
#   bash assets/build_references.sh
#
# Output:
#   assets/enteric_species_refs.msh   (Mash sketch)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFS_DIR="${SCRIPT_DIR}/reference_genomes"
SKETCH="${SCRIPT_DIR}/enteric_species_refs.msh"
KMER_SIZE=21
SKETCH_SIZE=1000

mkdir -p "${REFS_DIR}"
cd "${REFS_DIR}"

echo "Downloading reference genomes from NCBI..."

# ── Download function using ncbi-datasets or fallback to direct FTP ───────────
download_genome() {
    local accession="$1"
    local label="$2"
    local outfile="${label}.fasta"

    if [ -f "${outfile}" ]; then
        echo "  [SKIP] ${label} already exists"
        return
    fi

    echo "  Downloading ${label} (${accession})..."

    # Try datasets CLI first (most reliable)
    if command -v datasets &>/dev/null; then
        datasets download genome accession "${accession}" \
            --include genome \
            --filename "${accession}.zip" 2>/dev/null
        unzip -p "${accession}.zip" "ncbi_dataset/data/${accession}/"*".fna" > "${outfile}" 2>/dev/null \
            || unzip -p "${accession}.zip" "*/*.fna" > "${outfile}" 2>/dev/null
        rm -f "${accession}.zip"
    else
        # Fallback: direct NCBI FTP via esearch/efetch (requires entrez-direct)
        if command -v efetch &>/dev/null; then
            efetch -db assembly -id "${accession}" -format docsum \
                | grep -oP 'FtpPath_RefSeq":\s*"\K[^"]+' \
                | xargs -I{} wget -q "{}/{}_genomic.fna.gz" -O "${accession}.fna.gz" 2>/dev/null
            gunzip -c "${accession}.fna.gz" > "${outfile}"
            rm -f "${accession}.fna.gz"
        else
            echo "ERROR: Neither 'datasets' nor 'efetch' found."
            echo "Install with: conda install -c conda-forge ncbi-datasets-cli"
            exit 1
        fi
    fi

    if [ -s "${outfile}" ]; then
        echo "    → ${outfile} ($(wc -l < "${outfile}") lines)"
    else
        echo "  WARNING: ${outfile} is empty. Check accession ${accession}."
    fi
}

# ── Reference genome downloads ────────────────────────────────────────────────
# E. coli references (both labelled Ecoli_ so species_check maps them to E_coli)
download_genome GCF_000005845.2  "Ecoli_K12_MG1655"
download_genome GCF_000006665.3  "Ecoli_O157H7_EDL933"

# Salmonella references (labelled Salmonella_ → Salmonella_enterica)
download_genome GCF_000006945.2  "Salmonella_Typhimurium_LT2"
download_genome GCF_000195995.1  "Salmonella_Typhi_CT18"
download_genome GCF_000009505.1  "Salmonella_Enteritidis_P125109"

# Shigella (labelled Shigella_ → rejected by species gate but logged)
download_genome GCF_000006925.2  "Shigella_sonnei_Ss046"

# Klebsiella (labelled Klebsiella_ → rejected but identified)
download_genome GCF_000240185.1  "Klebsiella_pneumoniae_HS11286"

# ── Build Mash sketch ─────────────────────────────────────────────────────────
echo ""
echo "Building Mash sketch (k=${KMER_SIZE}, s=${SKETCH_SIZE})..."
FASTAS=(Ecoli_*.fasta Salmonella_*.fasta Shigella_*.fasta Klebsiella_*.fasta)
N="${#FASTAS[@]}"
echo "  ${N} reference genome(s) found"

mash sketch \
    -k "${KMER_SIZE}" \
    -s "${SKETCH_SIZE}" \
    -o "${SKETCH%.msh}" \
    "${FASTAS[@]}"

echo ""
echo "Done! Reference sketch written to:"
echo "  ${SKETCH}"
echo ""
echo "Verify with: mash info ${SKETCH}"
