// ── PLOT_SUMMARY: publication-ready summary figures ──────────────────────────
// Produces Figs 1, 3, 4 (population summary, AMR genes, plasmid replicons).
// The tree-annotated resistome heatmap is produced separately by TREE_ANNOTATION.

process PLOT_SUMMARY {
    label 'low'

    conda     "${projectDir}/envs/plots.yml"
    container 'quay.io/biocontainers/python:3.11'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(results_tsv)
    val(species_label)

    output:
    path("${species_label}_fig*.{pdf,png}"), emit: figures

    script:
    """
    plot_summary.py \\
        --input   ${results_tsv} \\
        --format  enteric-typer \\
        --outdir  . \\
        --prefix  ${species_label}
    """
}
