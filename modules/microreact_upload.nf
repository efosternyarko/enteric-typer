// ── MICROREACT_UPLOAD: create a Microreact project with results + tree ────────
// POSTs the aggregated results table and (optionally) the IQ-TREE Newick tree
// to the Microreact v2 API.
//
// Requires MICROREACT_TOKEN to be set as a Nextflow secret or environment variable.
//   Set with: nextflow secrets set MICROREACT_TOKEN <your_token>

process MICROREACT_UPLOAD {
    label 'low'

    conda     "${projectDir}/envs/utils.yml"
    container 'python:3.11-slim'

    secret 'MICROREACT_TOKEN'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(results_tsv)
    path(tree_file)       // IQ-TREE .treefile (Newick) — may be 'NO_FILE' placeholder
    val(project_name)

    output:
    path("microreact_url_${results_tsv.baseName}.txt")

    script:
    def tree_arg = (tree_file.name != 'NO_FILE' && tree_file.size() > 0)
        ? "--tree ${tree_file}"
        : ""
    """
    upload_microreact.py \\
        --input   ${results_tsv} \\
        --project "${project_name}" \\
        --token   "\${MICROREACT_TOKEN}" \\
        ${tree_arg} \\
        --output  microreact_url_${results_tsv.baseName}.txt
    """
}
