// ── PLOT_ASSEMBLY_METRICS: violin plots of assembly QC per species ─────────
// Generates a 4-panel figure (genome length, contig count, N50, GC%) styled
// to match the upper panels of the reference assembly+read metrics figure.
// One figure is produced per species.
//
// Inputs:
//   stats_files  — collected list of *_assembly_stats.tsv files for one species
//   species_id   — short identifier: ecoli | salmonella | shigella

process PLOT_ASSEMBLY_METRICS {
    label 'low'

    conda     "${projectDir}/envs/plots.yml"
    container 'quay.io/biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0'

    input:
    path(stats_files)
    val(species_id)

    output:
    path("${species_id}_assembly_metrics.png"),          emit: plot
    path("${species_id}_assembly_metrics.pdf"),          emit: plot_pdf
    path("${species_id}_assembly_metrics_summary.tsv"),  emit: summary

    script:
    """
    # Guard: skip if channel produced no files (species absent from this run)
    if [ "${stats_files}" = "[]" ] || [ -z "${stats_files}" ]; then
        touch ${species_id}_assembly_metrics.png
        touch ${species_id}_assembly_metrics_summary.tsv
        exit 0
    fi

    plot_assembly_metrics.py \\
        --stats ${stats_files} \\
        --species ${species_id} \\
        --output ${species_id}_assembly_metrics.png \\
        --summary ${species_id}_assembly_metrics_summary.tsv
    """
}
