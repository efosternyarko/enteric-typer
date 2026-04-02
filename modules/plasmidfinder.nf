// ── PLASMIDFINDER: plasmid replicon typing ───────────────────────────────────
// Identifies plasmid replicons using the CGE PlasmidFinder database.
// database parameter: 'enterobacteriaceae' or 'gram_positive'

process PLASMIDFINDER {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/plasmidfinder.yml"
    container 'staphb/plasmidfinder:2.1.6'

    input:
    tuple val(sample_id), path(fasta)
    val(database)

    output:
    tuple val(sample_id), path("${sample_id}_plasmidfinder.tsv")

    script:
    """
    # Locate blastn and PlasmidFinder database bundled with the conda env
    BLASTN_PATH=\$(which blastn)
    ENV_PREFIX=\$(dirname \$(dirname "\$BLASTN_PATH"))
    DB_PATH=\$(find "\$ENV_PREFIX/share" -name "database" -type d 2>/dev/null | head -1)
    [ -z "\$DB_PATH" ] && DB_PATH="/database"

    plasmidfinder.py \\
        -i ${fasta} \\
        -o ./ \\
        -mp "\$BLASTN_PATH" \\
        -p  "\$DB_PATH" \\
        -d  ${database} \\
        -l 0.60 \\
        -t 0.90 \\
        -q \\
        2>${sample_id}_plasmidfinder.log \\
    || true

    if [ -f results_tab.tsv ] && [ \$(wc -l < results_tab.tsv) -gt 1 ]; then
        # Prepend sample ID column
        awk -v s="${sample_id}" 'NR==1 {print "sample\\t" \$0; next} {print s "\\t" \$0}' \\
            results_tab.tsv > ${sample_id}_plasmidfinder.tsv
    else
        printf 'sample\\tPlasmid\\tIdentity\\tQuery / Template length\\tContig\\tPosition in contig\\tNote\\tAccession number\\n' \\
            > ${sample_id}_plasmidfinder.tsv
        printf '%s\\tNo replicons found\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\n' ${sample_id} \\
            >> ${sample_id}_plasmidfinder.tsv
    fi
    """
}
