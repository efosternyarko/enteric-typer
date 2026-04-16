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
        | awk '
        {
            # Map each reference to a species group
            n = split(\$1, a, "/")
            ref = a[n]
            gsub(/\\.(fasta|fa|fna|fas|fsa)\$/, "", ref)
            if      (ref ~ /^[Ee]_?[Cc]oli/ || ref ~ /^Escherichia/) sp = "E_coli"
            else if (ref ~ /^[Ss]almonella/)                          sp = "Salmonella_enterica"
            else if (ref ~ /^[Ss]higella/)                            sp = "Shigella"
            else if (ref ~ /^[Kk]lebsiella/)                          sp = "Klebsiella"
            else if (ref ~ /^[Ee]nterobacter/)                        sp = "Enterobacter"
            else                                                       sp = ref

            dist = \$3 + 0
            # Track best (lowest) distance seen per species group
            if (!(sp in best) || dist < best[sp])
                best[sp] = dist
        }
        END {
            # Closest-reference-wins (same approach as Kleborate).
            # A previous Shigella-priority rule (< 0.025 → Shigella regardless
            # of E. coli distance) was removed because it misclassified E. coli
            # phylogroup B2 isolates whose Mash distance to Shigella references
            # (0.014–0.025) is below the threshold but whose distance to E. coli
            # references is clearly lower (~0.007). Genuine Shigella isolates
            # are closer to Shigella references than to any E. coli reference
            # and are correctly identified without a priority rule.
            best_sp   = ""
            best_dist = 1.0
            for (sp in best) {
                if (best[sp] < best_dist) {
                    best_dist = best[sp]
                    best_sp   = sp
                }
            }
            print best_sp "\\t" best_dist
        }' > ${sample_id}_species.txt

    # Guard: ensure file is non-empty
    [ -s "${sample_id}_species.txt" ] || echo -e "Unknown\\t1.0" > "${sample_id}_species.txt"
    """
}
