// ── MLST: 7-locus multi-locus sequence typing ────────────────────────────────
// scheme parameter controls which PubMLST scheme is used:
//   ecoli_achtman_4  →  E. coli (Warwick Achtman 7-locus)
//   salmonella       →  Salmonella enterica (Warwick)

process MLST {
    tag "${sample_id} [${scheme}]"
    label 'low'

    conda     "${projectDir}/envs/mlst.yml"
    container 'staphb/mlst:2.23.0'

    input:
    tuple val(sample_id), path(fasta)
    val(scheme)

    output:
    tuple val(sample_id), path("${sample_id}_${scheme}_mlst.tsv")

    script:
    """
    mlst \\
        --scheme ${scheme} \\
        --threads ${task.cpus} \\
        --label   ${sample_id} \\
        ${fasta} > ${sample_id}_${scheme}_mlst.tsv \\
    || printf '%s\\t%s\\tNA\\n' ${sample_id} ${scheme} > ${sample_id}_${scheme}_mlst.tsv
    """
}
