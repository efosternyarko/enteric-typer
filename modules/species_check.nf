// ── SPECIES_CHECK: Mash species identification ───────────────────────────────
// Identifies the closest reference species for each assembly using Mash distance.
// Reference sketch must contain labelled sequences for:
//   E_coli, Salmonella_enterica, Shigella, Klebsiella (and others as desired).
//
// Output format (tab-separated, one line):
//   <species_label>\t<mash_distance>
// e.g.  E_coli\t0.0012

process SPECIES_CHECK {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/mash.yml"
    container 'quay.io/biocontainers/mash:2.3--hb105d93_10'

    input:
    tuple val(sample_id), path(fasta)
    path(reference_sketch)

    output:
    tuple val(sample_id), path("${sample_id}_species.txt"), emit: species

    script:
    """
    mash dist ${reference_sketch} ${fasta} 2>/dev/null \\
        | sort -k3,3n \\
        | head -1 \\
        | awk '{
            # Extract filename component from the reference path
            n = split(\$1, a, "/")
            ref = a[n]
            # Strip any FASTA extension
            gsub(/\\.(fasta|fa|fna|fas|fsa)\$/, "", ref)
            # Map reference name prefix to a canonical species label
            if      (ref ~ /^[Ee]_?[Cc]oli/ || ref ~ /^Escherichia/) species = "E_coli"
            else if (ref ~ /^[Ss]almonella/)                          species = "Salmonella_enterica"
            else if (ref ~ /^[Ss]higella/)                            species = "Shigella"
            else if (ref ~ /^[Kk]lebsiella/)                          species = "Klebsiella"
            else if (ref ~ /^[Ee]nterobacter/)                        species = "Enterobacter"
            else                                                       species = ref
            print species "\\t" \$3
        }' > ${sample_id}_species.txt

    # Guard: ensure file is non-empty
    [ -s "${sample_id}_species.txt" ] || echo -e "Unknown\\t1.0" > "${sample_id}_species.txt"
    """
}
