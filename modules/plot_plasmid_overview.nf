// ── PLOT_PLASMID_OVERVIEW: multi-panel plasmid figure ────────────────────────
// Produces a three-panel figure:
//   Panel A (top, full width): SNP tree | ST strip | PG strip
//             | plasmid replicon heatmap (coloured by dominant AMR class)
//   Panel B (bottom-left):  horizontal stacked-bar chart — replicon prevalence
//             by drug class (dominant-class allocation, non-overlapping)
//   Panel C (bottom-right): plasmid–AMR co-occurrence bubble matrix
//
// Replaces: ecoli_fig4_plasmid_replicons, ecoli_plasmid_amr_map_bubble figures.
// Run for E. coli, Salmonella, and Shigella.

process PLOT_PLASMID_OVERVIEW {
    label 'low'

    conda     "${projectDir}/envs/plots.yml"
    container 'quay.io/biocontainers/python:3.11'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(treefile)
    path(results_tsv)
    path(plasmid_amr_map)
    val(species_label)

    output:
    path("${species_label}_fig4_plasmid_overview.{pdf,png}"), emit: figures
    path("individual_plasmid_plots/"),                         emit: panels, optional: true

    script:
    def tree_arg = (treefile && treefile.name != 'NO_FILE') ? "--tree ${treefile}" : ""
    """
    plot_plasmid_overview.py \\
        ${tree_arg} \\
        --metadata    ${results_tsv} \\
        --plasmid_map ${plasmid_amr_map} \\
        --outdir      . \\
        --prefix      ${species_label} \\
        --top_n       15
    """
}
