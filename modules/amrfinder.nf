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
    val(db_ready)   // sentinel from AMRFINDER_UPDATE — ensures DB exists before any sample runs

    output:
    tuple val(sample_id), path("${sample_id}_amrfinder.tsv")

    script:
    """
    amrfinder \\
        --nucleotide ${fasta} \\
        --organism   ${organism} \\
        --plus \\
        --threads    ${task.cpus} \\
        --name       ${sample_id} \\
        --output     ${sample_id}_amrfinder.tsv
    """
}
