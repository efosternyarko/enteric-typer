// ── KLEBORATE: E. coli virulence pathotyping + Clermont phylogroup ────────────
// Kleborate v3 --preset escherichia runs:
//   • escherichia__ezclermont  → Clermont phylogroup (clermont_type)
//   • escherichia__pathovar    → pathotype (STEC/EPEC/ETEC/EIEC/EHEC/-)
//                                + individual markers: Stx1, Stx2, eae, ipaH, LT, ST
//   • escherichia__stxtyper    → detailed Stx subtyping
//   • escherichia__ectyper     → O:H serotype
//   • escherichia__mlst_achtman → MLST ST
//   • general__contig_stats    → assembly QC
//
// Output written to escherichia_output.txt inside --outdir.
// --trim_headers strips the module__ prefix so column names are plain (Pathotype, clermont_type…).
//
// Reference: Lam et al., 2021 / github.com/klebgenomics/Kleborate

process KLEBORATE {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/kleborate.yml"
    container 'staphb/kleborate:3.1.3'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_kleborate.tsv")

    script:
    """
    kleborate \\
        --assemblies  ${fasta} \\
        --preset      escherichia \\
        --outdir      kleborate_out \\
        --trim_headers \\
        2>${sample_id}_kleborate.log \\
    || true

    OUTFILE="kleborate_out/escherichia_output.txt"

    if [ -f "\$OUTFILE" ] && [ \$(wc -l < "\$OUTFILE") -ge 2 ]; then
        # Prepend a sample column so downstream aggregation can join on it
        awk -v s="${sample_id}" \\
            'NR==1 {print "sample\\t" \$0; next} {print s "\\t" \$0}' \\
            "\$OUTFILE" > ${sample_id}_kleborate.tsv
    else
        # Fallback empty row — keeps AGGREGATE from failing if Kleborate had no result
        printf 'sample\\tstrain\\tclermont_type\\tPathotype\\tStx1\\tStx2\\tST\\tLT\\teae\\tipaH\\n' \\
            > ${sample_id}_kleborate.tsv
        printf '%s\\t%s\\t-\\t-\\t-\\t-\\t-\\t-\\t-\\t-\\n' \\
            "${sample_id}" "${sample_id}" >> ${sample_id}_kleborate.tsv
    fi
    """
}
