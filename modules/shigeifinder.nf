// ── SHIGEIFINDER: Shigella serotyping and species identification ──────────────
// Predicts Shigella species and serotype from assembled genomes using
// cluster-specific gene markers (LanLab/ShigEiFinder).
// Covers S. sonnei, S. flexneri (including subtypes 1a/1b/2a/2b/X/Y),
// S. boydii, and S. dysenteriae. Also reports ipaH (invasion plasmid marker).
//
// Output columns: sample  ipaH  Cluster  Serotype  Notes

process SHIGEIFINDER {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/shigeifinder.yml"
    container 'quay.io/biocontainers/shigeifinder:1.3.5--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_shigeifinder.tsv")

    script:
    """
    shigeifinder \\
        -i  ${fasta} \\
        -t  ${task.cpus} \\
        2>${sample_id}_shigeifinder.log \\
    | awk -v sid="${sample_id}" '
        NR == 1 { print "sample\\t" \$0 }
        NR  > 1 { print sid "\\t" \$0 }
    ' > ${sample_id}_shigeifinder.tsv \\
    || (
        printf 'sample\\tipaH\\tCluster\\tSerotype\\tNotes\\n' \\
            > ${sample_id}_shigeifinder.tsv
        printf '%s\\tNA\\tNA\\tNA\\tNA\\n' "${sample_id}" \\
            >> ${sample_id}_shigeifinder.tsv
    )
    """
}
