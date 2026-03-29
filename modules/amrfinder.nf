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
    # Update database if not already present in the container/conda env
    amrfinder -u 2>/dev/null || true

    amrfinder \\
        --nucleotide ${fasta} \\
        --organism   ${organism} \\
        --plus \\
        --threads    ${task.cpus} \\
        --name       ${sample_id} \\
        --output     ${sample_id}_amrfinder.tsv \\
    || echo -e 'Name\\tProtein identifier\\tContig id\\tStart\\tStop\\tStrand\\tGene symbol\\tSequence name\\tScope\\tElement type\\tElement subtype\\tClass\\tSubclass\\tMethod\\tTarget length\\tReference sequence length\\t% Coverage of reference sequence\\t% Identity to reference sequence\\tAlignment length\\tAccession of closest sequence\\tName of closest sequence\\tHMM id\\tHMM description' \\
         > ${sample_id}_amrfinder.tsv
    """
}
