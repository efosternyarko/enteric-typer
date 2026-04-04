// ── IS_SCREEN: insertion sequence element detection via BLAST ─────────────────
// Screens assemblies for IS elements relevant to Shigella (IS1, IS1A, IS30,
// IS600, IS186, IS629) using a curated reference FASTA in assets/.
// Reports copy number and contig positions per IS element per sample.

process IS_SCREEN {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/mlst.yml"
    container 'staphb/mlst:2.23.0'

    input:
    tuple val(sample_id), path(fasta)
    path(is_db)

    output:
    tuple val(sample_id), path("${sample_id}_is_screen.tsv")

    script:
    """
    makeblastdb -in ${is_db} -dbtype nucl -out is_db -logfile /dev/null 2>/dev/null || \\
    makeblastdb -in ${is_db} -dbtype nucl -out is_db 2>/dev/null

    blastn \\
        -query        ${fasta} \\
        -db           is_db \\
        -out          blast_is.tsv \\
        -outfmt       "6 qseqid sseqid pident length qstart qend slen" \\
        -perc_identity 85 \\
        -num_threads  ${task.cpus} \\
        -dust         no \\
        2>/dev/null || true

    printf 'sample\\tIS_element\\tcopies\\tlocations\\n' > ${sample_id}_is_screen.tsv
    if [ -s blast_is.tsv ]; then
        awk -v s="${sample_id}" '
        {
            is[\$2]++
            loc[\$2] = loc[\$2] (loc[\$2] == "" ? "" : ";") \$1 ":" \$5 "-" \$6
        }
        END {
            for (e in is) print s, e, is[e], loc[e]
        }
        ' OFS="\\t" blast_is.tsv >> ${sample_id}_is_screen.tsv
    fi
    """
}
