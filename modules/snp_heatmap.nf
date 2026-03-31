// ── SNP_HEATMAP: pairwise SNP distance heatmap ───────────────────────────────
// Generates a standalone clustered heatmap from the SKA2 pairwise SNP matrix.

process SNP_HEATMAP {
    label 'low'

    conda     "${projectDir}/envs/plots.yml"
    container 'python:3.11'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(snp_matrix)         // snp_matrix.tsv from SKA2_BUILD
    val(species_label)       // 'ecoli' or 'salmonella'

    output:
    path("${species_label}_snp_heatmap.pdf"), optional: true
    path("${species_label}_snp_heatmap.png"), optional: true

    script:
    """
    # Skip if SNP matrix was not produced (< 3 samples or SKA2 failed)
    if [ ! -f "${snp_matrix}" ] || [ "${snp_matrix}" = "NO_FILE" ]; then
        echo "INFO: No SNP matrix available for ${species_label} — skipping heatmap." >&2
        exit 0
    fi

    plot_snp_heatmap.py \\
        --matrix  ${snp_matrix} \\
        --species ${species_label} \\
        --outdir  .
    """
}
