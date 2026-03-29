// ── ECTYPER: E. coli O:H serotyping ──────────────────────────────────────────
// ECTyper predicts O and H antigen serotype from genome assemblies.
// Reference: Laing et al., 2019 (github.com/phac-nml/ecoli_serotyping)

process ECTYPER {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/ectyper.yml"
    container 'kbessonov/ectyper:1.0.0'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_ectyper.tsv")

    script:
    """
    ectyper \\
        --input   ${fasta} \\
        --output  ectyper_out \\
        --cores   ${task.cpus} \\
        --opid    90 \\
        --hpid    90 \\
        --opcov   50 \\
        --hpcov   50 \\
        2>${sample_id}_ectyper.log \\
    || true

    # ECTyper writes output.tsv inside the output directory
    if [ -f ectyper_out/output.tsv ]; then
        # Add sample ID prefix column for aggregation
        awk -v s="${sample_id}" 'NR==1 {print "sample\\t" \$0; next} {print s "\\t" \$0}' \\
            ectyper_out/output.tsv > ${sample_id}_ectyper.tsv
    else
        echo -e "sample\\tName\\tO-type\\tH-type\\tSerotype\\tQC\\tEvidence\\tGenomeSize\\tAlleleDatabase\\tVersion" \\
            > ${sample_id}_ectyper.tsv
        echo -e "${sample_id}\\t${sample_id}\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA" \\
            >> ${sample_id}_ectyper.tsv
    fi
    """
}
