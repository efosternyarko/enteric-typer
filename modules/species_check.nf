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
            # Shigella-priority rule: Shigella and E. coli are phylogenetically
            # interleaved; if ANY Shigella reference is within 0.025 the sample
            # is Shigella, even when an E. coli reference is marginally closer.
            if ("Shigella" in best && best["Shigella"] < 0.025) {
                print "Shigella\\t" best["Shigella"]
            } else {
                best_sp   = ""
                best_dist = 1.0
                for (sp in best) {
                    if (best[sp] < best_dist) {
                        best_dist = best[sp]
                        best_sp   = sp
                    }
                }
                print best_sp "\\t" best_dist
            }
        }' > ${sample_id}_species.txt

    # Guard: ensure file is non-empty
    [ -s "${sample_id}_species.txt" ] || echo -e "Unknown\\t1.0" > "${sample_id}_species.txt"
    """
}
