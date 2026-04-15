// ── KRAKEN2_SCREEN: contamination screening on assembled FASTAs ──────────────
// Classifies contigs with Kraken2 against a user-supplied database and flags
// samples where secondary species exceed --max_contamination % of total sequences.
//
// Inputs:
//   fasta  — assembled FASTA (any source)
//   db     — Kraken2 database directory (e.g. k2_standard_08gb)
//
// Outputs:
//   report — full Kraken2 report (published for review)
//   qc     — one-row TSV: sample_id, top_species, top_pct, secondary_pct, status, reason

process KRAKEN2_SCREEN {
    tag "${sample_id}"
    label 'high'

    conda     "${projectDir}/envs/kraken2.yml"
    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0'

    input:
    tuple val(sample_id), path(fasta)
    path(db)

    output:
    tuple val(sample_id), path("${sample_id}_kraken2_report.txt"), emit: report
    tuple val(sample_id), path("${sample_id}_kraken2_qc.tsv"),    emit: qc

    script:
    """
    kraken2 \
        --db ${db} \
        --threads ${task.cpus} \
        --report ${sample_id}_kraken2_report.txt \
        --output /dev/null \
        ${fasta}

    parse_kraken2.py \
        --report  ${sample_id}_kraken2_report.txt \
        --sample  ${sample_id} \
        --max_secondary ${params.max_contamination} \
        --output  ${sample_id}_kraken2_qc.tsv
    """
}
