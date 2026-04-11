// ── ASSEMBLY_QC: per-sample assembly statistics ────────────────────────────
// Computes genome length, contig count, N50, and GC% from an assembled FASTA
// using seqkit stats. Output feeds PLOT_ASSEMBLY_METRICS for per-species violin plots.
//
// Output TSV columns:
//   sample_id  num_contigs  genome_length  assembly_N50  gc_pct

process ASSEMBLY_QC {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/seqkit.yml"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_assembly_stats.tsv"), emit: stats

    script:
    """
    seqkit stats -a -T -j ${task.cpus} ${fasta} > seqkit_raw.txt
    parse_assembly_stats.py ${sample_id} seqkit_raw.txt \\
        > ${sample_id}_assembly_stats.tsv
    """
}
