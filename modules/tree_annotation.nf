// ── TREE_ANNOTATION: phylogenetic tree annotated with AMR heatmap ─────────────
// Combines the IQ-TREE Newick output with the aggregated typing results to
// produce a single figure: tree | phylogroup strip | drug-class heatmap.
// Skipped gracefully if no tree is available (< 3 samples).

process TREE_ANNOTATION {
    label 'low'
    cache false   // fast plotting step — always re-run to avoid caching silent failures

    conda     "${projectDir}/envs/plots.yml"
    container 'quay.io/biocontainers/python:3.11'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(treefile)
    path(results_tsv)
    val(species_label)

    output:
    path("${species_label}_fig2_tree_amr.{pdf,png}"), emit: figures, optional: true

    script:
    """
    # Skip if tree was not produced (< 3 samples or SKA2 failed)
    if [ ! -f "${treefile}" ] || [ "${treefile}" = "NO_FILE" ]; then
        echo "INFO: No tree available for ${species_label} — skipping tree annotation." >&2
        exit 0
    fi

    plot_tree_annotation.py \\
        --tree     ${treefile} \\
        --metadata ${results_tsv} \\
        --outdir   . \\
        --prefix   ${species_label} \\
        --species  ${species_label}
    """
}
