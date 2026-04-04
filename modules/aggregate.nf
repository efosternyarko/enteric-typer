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
    path(ktype_files)         // parse_kaptive output (E. coli only; empty list for Salmonella)
    val(pathogenwatch_tsv)    // path string or 'NO_FILE'
    path(st_complexes)
    path(amrrules_tsv)        // AMRrules wildtype gene definitions for this species
    val(species_label)        // 'ecoli' or 'salmonella'
    path(kleborate_files)     // Kleborate output (E. coli only)
    path(clermont_files)      // EzClermont phylotyping (E. coli only)
    path(mykrobe_files)       // Mykrobe parsed TSVs (Shigella only)
    path(pinv_files)          // pINV screen results (Shigella only)
    path(is_files)            // IS element screen results (Shigella only)

    output:
    path("${species_label}_typer_results.tsv"), emit: results

    script:
    def pw_arg = (pathogenwatch_tsv && pathogenwatch_tsv != 'NO_FILE')
        ? "--pathogenwatch ${pathogenwatch_tsv}"
        : "--pathogenwatch NO_FILE"
    """
    KTYPE_ARG=""
    if [ -n "${ktype_files}" ] && [ "${ktype_files}" != "[]" ]; then
        KTYPE_ARG="--ktype ${ktype_files}"
    fi

    KLEB_ARG=""
    if [ -n "${kleborate_files}" ] && [ "${kleborate_files}" != "[]" ]; then
        KLEB_ARG="--kleborate ${kleborate_files}"
    fi

    CLERMONT_ARG=""
    if [ -n "${clermont_files}" ] && [ "${clermont_files}" != "[]" ]; then
        CLERMONT_ARG="--clermont ${clermont_files}"
    fi

    MYKROBE_ARG=""
    if [ -n "${mykrobe_files}" ] && [ "${mykrobe_files}" != "[]" ]; then
        MYKROBE_ARG="--mykrobe ${mykrobe_files}"
    fi

    PINV_ARG=""
    if [ -n "${pinv_files}" ] && [ "${pinv_files}" != "[]" ]; then
        PINV_ARG="--pinv ${pinv_files}"
    fi

    IS_ARG=""
    if [ -n "${is_files}" ] && [ "${is_files}" != "[]" ]; then
        IS_ARG="--is-screen ${is_files}"
    fi

    # If there are no MLST files, no samples of this species were detected — exit cleanly
    MLST_FILES="${mlst_files}"
    if [ -z "\$MLST_FILES" ] || [ "\$MLST_FILES" = "[]" ]; then
        echo "No samples for species ${species_label} — skipping aggregation" >&2
        touch ${species_label}_typer_results.tsv
        exit 0
    fi

    aggregate_results.py \\
        --species       ${species_label} \\
        --mlst          ${mlst_files} \\
        --amrfinder     ${amrfinder_files} \\
        --serotyper     ${serotyper_files} \\
        --plasmidfinder ${plasmidfinder_files} \\
        \$KTYPE_ARG \\
        \$KLEB_ARG \\
        \$CLERMONT_ARG \\
        \$MYKROBE_ARG \\
        \$PINV_ARG \\
        \$IS_ARG \\
        ${pw_arg} \\
        --st-complexes  ${st_complexes} \\
        --amrrules      ${amrrules_tsv} \\
        --output        ${species_label}_typer_results.tsv
    """
}
