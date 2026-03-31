// ── SKA2_BUILD: split k-mer core-SNP alignment ───────────────────────────────
// Takes all FASTA files for one species group, builds a split k-mer file (SKF),
// and produces a multi-FASTA alignment for IQ-TREE plus a pairwise SNP matrix.
//
// Parameters follow the lab's standard SKA protocol:
//   ska build  -k 31
//   ska distance --min-freq 1   (100% core — do not change)
//   ska align  --min-freq 1 --filter no-filter
//
// Note: the bioconda `ska2` package installs the binary as `ska` (not `ska2`).
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
    path("snp_matrix.tsv"),       emit: snp_matrix, optional: true
    path("*.skf"),                emit: skf, optional: true

    script:
    def min_samples = params.ska2_min_samples ?: 3
    """
    # Build tab-separated input list (sample_name<TAB>./file) required by ska build -f
    for f in *.fasta *.fa *.fna *.fas *.fsa; do
        [ -f "\$f" ] || continue
        name="\${f%.*}"
        printf "%s\\t./%s\\n" "\$name" "\$f"
    done > input_sequence.txt

    N=\$(wc -l < input_sequence.txt | tr -d ' ')

    if [ "\$N" -lt ${min_samples} ]; then
        echo "INFO: Only \$N sample(s) — need ≥ ${min_samples} for phylogenetics. Producing empty alignment." >&2
        echo "" > ska2_alignment.fasta
        exit 0
    fi

    echo "INFO: Building SKA2 sketch for \$N samples..." >&2

    # Build merged split k-mer file
    # (bioconda ska2 package installs the binary as 'ska')
    ska build \\
        -o seqs \\
        -k 31 \\
        --threads ${task.cpus} \\
        -f input_sequence.txt \\
        2>&1

    # Generate pairwise SNP distance matrix
    # --min-freq 1: require k-mer present in 100% of samples (core genome)
    ska distance \\
        --min-freq 1 \\
        -o snp_matrix.tsv \\
        seqs.skf \\
        2>&1

    # Generate FASTA alignment
    # --min-freq 1 --filter no-filter: core sites, retain all variant classes
    ska align \\
        --min-freq 1 \\
        --filter no-filter \\
        -o ska2_alignment.fasta \\
        seqs.skf \\
        2>&1

    # Guard: if alignment is empty (no variant sites), warn
    if [ ! -s ska2_alignment.fasta ]; then
        echo "WARNING: SKA2 alignment is empty (no variant sites in core genome)" >&2
    fi
    """
}
