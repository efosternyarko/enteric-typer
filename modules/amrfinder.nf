// ── AMRFINDER PLUS: AMR, stress, and virulence gene detection ────────────────
// organism parameter:
//   Escherichia  →  E. coli (enables species-specific point mutations)
//   Salmonella   →  Salmonella enterica

process AMRFINDER {
    tag "${sample_id} [${organism}]"
    label 'medium'

    conda     "${projectDir}/envs/amrfinder.yml"
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0'

    input:
    tuple val(sample_id), path(fasta)
    val(organism)

    output:
    tuple val(sample_id), path("${sample_id}_amrfinder.tsv")

    script:
    """
    # AMRFinder Plus does not bundle its database with the conda package.
    # Download it on first use; subsequent runs reuse the cached copy.
    DB_DIR="\$(dirname \$(which amrfinder))/../share/amrfinderplus/data"
    if [ ! -d "\${DB_DIR}/latest" ]; then
        amrfinder -u
    fi

    amrfinder \\
        --nucleotide ${fasta} \\
        --organism   ${organism} \\
        --plus \\
        --threads    ${task.cpus} \\
        --name       ${sample_id} \\
        --output     ${sample_id}_amrfinder.tsv
    """
}
