// ── PLASMID_AMR_MAP: link plasmid replicons to co-located AMR genes ───────────
// Joins PlasmidFinder and AMRFinder outputs on contig ID to show which
// AMR genes share a contig with each identified replicon type.
//
// Output TSV columns:
//   sample_id  replicon  contig  amr_genes  drug_classes  identity  coverage

process PLASMID_AMR_MAP {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/utils.yml"
    container 'python:3.11'

    publishDir "${params.outdir}/plasmid_amr_map", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(plasmidfinder_tsv), path(amrfinder_tsv)

    output:
    tuple val(sample_id), path("${sample_id}_plasmid_amr_map.tsv"), emit: map

    script:
    """
    plasmid_amr_map.py \\
        --plasmidfinder ${plasmidfinder_tsv} \\
        --amrfinder     ${amrfinder_tsv} \\
        --sample        ${sample_id} \\
        --output        ${sample_id}_plasmid_amr_map.tsv
    """
}
