// ── IQTREE: maximum-likelihood phylogenetic tree ─────────────────────────────
// Builds a ML tree from the SKA2 core-SNP alignment using IQ-TREE 2.
// Skips gracefully if the alignment is empty (< 2 sequences).

process IQTREE {
    label 'high'

    conda     "${projectDir}/envs/iqtree.yml"
    container 'quay.io/biocontainers/iqtree:2.3.6--h21ec9f0_0'

    input:
    path(alignment)

    output:
    path("iqtree.treefile"),   emit: treefile,  optional: true
    path("iqtree.iqtree"),     emit: log,       optional: true
    path("iqtree.*"),          emit: all,        optional: true

    script:
    def model       = params.iqtree_model      ?: 'GTR+G'
    def bootstraps  = params.iqtree_bootstraps ?: 1000
    """
    # Count sequences in alignment
    N_SEQ=\$(grep -c '^>' ${alignment} 2>/dev/null || echo 0)

    if [ "\$N_SEQ" -lt 2 ]; then
        echo "INFO: Alignment has \$N_SEQ sequence(s) — need ≥ 2 for IQ-TREE. Skipping." >&2
        exit 0
    fi

    iqtree2 \\
        -s       ${alignment} \\
        -m       ${model} \\
        -B       ${bootstraps} \\
        -T       ${task.cpus} \\
        --prefix iqtree \\
        --redo \\
        2>&1

    # Rename treefile for clarity if it exists
    [ -f iqtree.treefile ] || echo "WARNING: IQ-TREE did not produce a treefile." >&2
    """
}
