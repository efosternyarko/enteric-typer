// ── PLASMIDFINDER: plasmid replicon typing ───────────────────────────────────
// Identifies plasmid replicons using the CGE PlasmidFinder database.
// database parameter: 'enterobacteriaceae' or 'gram_positive'

process PLASMIDFINDER {
    tag "${sample_id}"
    label 'low'

    conda     "${projectDir}/envs/plasmidfinder.yml"
    container 'staphb/plasmidfinder:2.1.6'

    input:
    tuple val(sample_id), path(fasta)
    val(database)

    output:
    tuple val(sample_id), path("${sample_id}_plasmidfinder.tsv")

    script:
    """
    # Locate blastn and PlasmidFinder database bundled with the conda env
    BLASTN_PATH=\$(which blastn)
    ENV_PREFIX=\$(dirname \$(dirname "\$BLASTN_PATH"))
    DB_PATH=\$(find "\$ENV_PREFIX/share" -name "database" -type d 2>/dev/null | head -1)
    [ -z "\$DB_PATH" ] && DB_PATH="/database"

    plasmidfinder.py \\
        -i ${fasta} \\
        -o ./ \\
        -mp "\$BLASTN_PATH" \\
        -p  "\$DB_PATH" \\
        -d  ${database} \\
        -l 0.60 \\
        -t 0.90 \\
        -q \\
        2>${sample_id}_plasmidfinder.log \\
    || true

    # v2.1.6 writes data.json instead of results_tab.tsv — parse JSON output
    python3 - <<'PYEOF'
import json, sys

sample = "${sample_id}"
header = "sample\tPlasmid\tIdentity\tQuery / Template length\tContig\tPosition in contig\tNote\tAccession number"

try:
    with open("data.json") as fh:
        data = json.load(fh)
    hits = []
    pf = data.get("plasmidfinder", {}).get("results", {})
    for db_group in pf.values():
        for db_hits in db_group.values():
            if not isinstance(db_hits, dict):
                continue
            for h in db_hits.values():
                hits.append("\t".join([
                    sample,
                    str(h.get("plasmid", "")),
                    str(h.get("identity", "")),
                    f"{h.get('HSP_length','')} / {h.get('template_length','')}",
                    str(h.get("contig_name", "")),
                    str(h.get("positions_in_contig", "")),
                    str(h.get("note", "")),
                    str(h.get("accession", "")),
                ]))
except Exception:
    hits = []

with open("${sample_id}_plasmidfinder.tsv", "w") as out:
    out.write(header + "\n")
    if hits:
        out.write("\n".join(hits) + "\n")
    else:
        out.write(f"{sample}\tNo replicons found\tNA\tNA\tNA\tNA\tNA\tNA\n")
PYEOF
    """
}
