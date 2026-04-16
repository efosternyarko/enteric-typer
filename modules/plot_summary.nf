// ── PLOT_SUMMARY: publication-ready summary figures ──────────────────────────
// Produces Figs 1, 3–7 (population summary, AMR genes, plasmid replicons,
// virulence, AMR-by-ST, AMR-by-group).
// Fig 2 (tree-annotated resistome heatmap) is produced by TREE_ANNOTATION.
// plasmid_amr_map is optional; pass 'NO_FILE' path to skip drug-class colouring.

process PLOT_SUMMARY {
    label 'low'

    conda     "${projectDir}/envs/plots.yml"
    container 'quay.io/biocontainers/python:3.11'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(results_tsv)
    val(species_label)
    path(plasmid_amr_map)

    output:
    path("${species_label}_fig*.{pdf,png}"), emit: figures

    script:
    def map_arg = plasmid_amr_map.name != 'NO_FILE' ? "--plasmid_map ${plasmid_amr_map}" : ""
    """
    plot_summary.py \\
        --input   ${results_tsv} \\
        --format  enteric-typer \\
        --outdir  . \\
        --prefix  ${species_label} \\
        ${map_arg}
    """
}
