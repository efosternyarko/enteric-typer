// ── AGGREGATE: merge per-sample tool outputs into one results table ───────────
// Produces a single TSV per species with all typing columns merged on sample ID.
// Missing data for any tool is represented as 'NA'.

process AGGREGATE {
    label 'low'

    conda     "${projectDir}/envs/utils.yml"
    container 'python:3.11'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(mlst_files)
    path(amrfinder_files)
    path(serotyper_files)     // ectyper (E. coli) or sistr (Salmonella)
    path(plasmidfinder_files)
    val(pathogenwatch_tsv)    // path string or 'NO_FILE'
    path(st_complexes)
    val(species_label)        // 'ecoli' or 'salmonella'

    output:
    path("${species_label}_typer_results.tsv"), emit: results

    script:
    def pw_arg = (pathogenwatch_tsv && pathogenwatch_tsv != 'NO_FILE')
        ? "--pathogenwatch ${pathogenwatch_tsv}"
        : "--pathogenwatch NO_FILE"
    """
    aggregate_results.py \\
        --species       ${species_label} \\
        --mlst          ${mlst_files} \\
        --amrfinder     ${amrfinder_files} \\
        --serotyper     ${serotyper_files} \\
        --plasmidfinder ${plasmidfinder_files} \\
        ${pw_arg} \\
        --st-complexes  ${st_complexes} \\
        --output        ${species_label}_typer_results.tsv
    """
}
