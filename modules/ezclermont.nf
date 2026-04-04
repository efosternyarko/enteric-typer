// ── EZCLERMONT: Clermont phylotyping for E. coli ──────────────────────────────
// EzClermont classifies E. coli assemblies into Clermont phylogroups
// (A, B1, B2, C, D, E, F, G or cryptic clades).
//
// Output: {sample_id}_clermont.tsv  with columns: sample, clermont_phylogroup
//
// Reference: Waters et al. (2020) Microbial Genomics 6(9)
//            github.com/nickp60/EzClermont

process EZCLERMONT {
    tag "${sample_id}"
    label 'low'

    conda "${projectDir}/envs/ezclermont.yml"

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_clermont.tsv")

    script:
    """
    # ezclermont outputs: "<experiment_name>\\t<phylogroup>" to stdout
    # -e sets the experiment name to sample_id so the first column matches
    # -m 200: lower minimum contig length from the 500 bp default so that
    #   assemblies with many short contigs still have enough sequence context
    #   for the control PCR amplicons (~300-400 bp) to be detected
    printf 'sample\\tclermont_phylogroup\\n' > ${sample_id}_clermont.tsv
    ezclermont ${fasta} \
        -e ${sample_id} \
        -m 200 \
        2>${sample_id}_clermont.log \
    | awk '{print \$1 "\\t" \$2}' \
    >> ${sample_id}_clermont.tsv \
    || true

    # Fallback if ezclermont failed or produced no output
    if [ \$(wc -l < "${sample_id}_clermont.tsv") -lt 2 ]; then
        printf '%s\\tUnknown\\n' "${sample_id}" >> ${sample_id}_clermont.tsv
    fi
    """
}
