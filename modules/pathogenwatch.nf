// ── PATHOGENWATCH: upload assemblies and retrieve cgMLST + cluster data ───────
// Uploads all confirmed species assemblies to Pathogenwatch, waits for
// processing, creates a collection, and runs cluster searches at multiple
// allele-difference thresholds (default: 5, 10, 20, 50 cgMLST alleles).
//
// Requires PW_API_KEY to be set as a Nextflow secret or environment variable.
//   Set with: nextflow secrets set PW_API_KEY <your_key>

process PATHOGENWATCH {
    label 'low'

    conda     "${projectDir}/envs/utils.yml"
    container 'python:3.11-slim'

    secret 'PW_API_KEY'

    publishDir "${params.outdir}/pathogenwatch", mode: 'copy', overwrite: true,
               saveAs: { fn -> fn.replaceAll('^pathogenwatch_', "${species_label}_pathogenwatch_") }

    input:
    val(samples)          // collected list of [id, fasta] tuples
    val(species_label)    // 'ecoli' or 'salmonella'
    val(collection_name)  // Pathogenwatch collection/folder name

    output:
    tuple val(species_label), path("pathogenwatch_samples.tsv"),    emit: results
    path "pathogenwatch_collection.json",                            emit: collection
    path "pathogenwatch_summary.json",                               emit: summary
    path "pathogenwatch_tree.nwk",                                   emit: tree, optional: true

    script:
    def rows = samples.collect { id, fasta -> "${id},${fasta}" }.join('\n')
    def thresholds = params.pathogenwatch_cluster_thresholds ?: '5,10,20,50'
    """
    # Write samplesheet for pathogenwatch_client.py
    cat > pw_samplesheet.csv << 'CSV'
id,fasta
${rows}
CSV

    pathogenwatch_client.py \\
        --samplesheet       pw_samplesheet.csv \\
        --collection-name   "${collection_name}" \\
        --thresholds        ${thresholds} \\
        --poll-seconds      ${params.pathogenwatch_poll_seconds} \\
        --max-wait-seconds  ${params.pathogenwatch_max_wait_seconds} \\
        --sample-output     pathogenwatch_samples.tsv \\
        --collection-output pathogenwatch_collection.json \\
        --summary-output    pathogenwatch_summary.json \\
        --tree-output       pathogenwatch_tree.nwk \\
        2>&1
    """
}
