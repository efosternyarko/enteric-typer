// ── SKA2_BUILD: split k-mer core-SNP alignment ───────────────────────────────
// Takes all FASTA files for one species group, builds a split k-mer file (SKF),
// filters to core sites (≥ params.ska2_prop_filter fraction of samples),
// and produces a multi-FASTA alignment for IQ-TREE.
//
// Requires ≥ params.ska2_min_samples samples; produces an empty alignment otherwise.

process SKA2_BUILD {
    label 'medium'

    conda     "${projectDir}/envs/ska2.yml"
    container 'quay.io/biocontainers/ska2:0.5.1--h4ac6f70_0'

    input:
    path(fastas)    // collected list of FASTA files

    output:
    path("ska2_alignment.fasta"), emit: alignment
    path("*.skf"),                emit: skf, optional: true

    script:
    def min_samples   = params.ska2_min_samples  ?: 3
    def prop_filter   = params.ska2_prop_filter  ?: 0.95
    """
    # Count input FASTA files (filter out any sentinel 'NO_FASTAS' placeholder)
    FASTAS=( \$(ls *.fasta *.fa *.fna *.fas *.fsa 2>/dev/null | grep -v 'NO_FASTAS') )
    N=\${#FASTAS[@]}

    if [ "\$N" -lt ${min_samples} ]; then
        echo "INFO: Only \$N sample(s) — need ≥ ${min_samples} for phylogenetics. Producing empty alignment." >&2
        echo "" > ska2_alignment.fasta
        exit 0
    fi

    echo "INFO: Building SKA2 sketch for \$N samples..." >&2

    # Build merged split k-mer file
    ska2 build -o merged -k 31 "\${FASTAS[@]}" 2>&1

    # Filter to core sites (present in ≥ prop_filter fraction of samples)
    ska2 weed \\
        --proportion-filter ${prop_filter} \\
        --filter-ambiguous \\
        merged.skf \\
        2>&1

    # Generate FASTA alignment
    ska2 align merged.skf 2>&1 > ska2_alignment.fasta

    # Guard: if alignment is empty (no variant sites), write minimal placeholder
    if [ ! -s ska2_alignment.fasta ]; then
        echo "WARNING: SKA2 alignment is empty (no variant sites at proportion filter ${prop_filter})" >&2
    fi
    """
}
