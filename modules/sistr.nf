// ── SISTR: Salmonella in silico typing ───────────────────────────────────────
// SISTR predicts Salmonella serovar from genome assemblies using cgMLST,
// antigen gene sequences (fliC, fljB, wzx/wzy) and the Kauffmann-White scheme.
// Reference: Yoshida et al., 2016

process SISTR {
    tag "${sample_id}"
    label 'medium'

    conda     "${projectDir}/envs/sistr.yml"
    container 'staphb/sistr:1.1.2'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_sistr.tsv")

    script:
    """
    sistr \\
        --input-fasta   ${fasta} \\
        --output-format tab \\
        --output-prediction ${sample_id}_sistr_raw.tsv \\
        --threads       ${task.cpus} \\
        --qc \\
        2>${sample_id}_sistr.log \\
    || true

    # Add a sample column for aggregation
    if [ -f "${sample_id}_sistr_raw.tsv" ] && [ \$(wc -l < "${sample_id}_sistr_raw.tsv") -gt 1 ]; then
        awk -v s="${sample_id}" 'NR==1 {print "sample\\t" \$0; next} {print s "\\t" \$0}' \\
            "${sample_id}_sistr_raw.tsv" > "${sample_id}_sistr.tsv"
    else
        echo -e "sample\\tgenome\\tserovar\\tserovar_antigen\\tserovar_cgmlst\\tO_antigen\\tH1\\tH2\\tqc_status\\tcgmlst_ST" \\
            > "${sample_id}_sistr.tsv"
        echo -e "${sample_id}\\t${sample_id}\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA" \\
            >> "${sample_id}_sistr.tsv"
    fi
    """
}
