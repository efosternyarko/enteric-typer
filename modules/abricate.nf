// ── ABRICATE: VFDB virulence gene screening ───────────────────────────────────
// Screens genome assemblies against VFDB using abricate.
// Used for Salmonella virulence profiling (E. coli uses Kleborate instead).
// Reference: Seemann T, github.com/tseemann/abricate

process ABRICATE {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/abricate.yml"
    container 'staphb/abricate:1.0.1'

    input:
    tuple val(sample_id), path(fasta)
    val(db)      // database name: 'vfdb', 'card', etc.
    val(mincov)  // minimum coverage (default 80)
    val(minid)   // minimum identity (default 90)

    output:
    tuple val(sample_id), path("${sample_id}_abricate_${db}.tsv")

    script:
    """
    abricate \\
        --db      ${db} \\
        --mincov  ${mincov} \\
        --minid   ${minid} \\
        --threads ${task.cpus} \\
        ${fasta} \\
        > abricate_raw.tsv \\
        2>${sample_id}_abricate_${db}.log \\
    || true

    if [ -f abricate_raw.tsv ] && [ \$(wc -l < abricate_raw.tsv) -ge 2 ]; then
        # abricate output has #FILE header; prepend sample column
        awk -v s="${sample_id}" \\
            'NR==1 {print "sample\\t" \$0; next} {print s "\\t" \$0}' \\
            abricate_raw.tsv > ${sample_id}_abricate_${db}.tsv
    else
        printf 'sample\\t#FILE\\tSEQUENCE\\tSTART\\tEND\\tSTRAND\\tGENE\\tCOVERAGE\\tCOVERAGE_MAP\\tGAPS\\t%%COVERAGE\\t%%IDENTITY\\tDATABASE\\tACCESSION\\tPRODUCT\\tRESISTANCE\\n' \\
            > ${sample_id}_abricate_${db}.tsv
    fi
    """
}
