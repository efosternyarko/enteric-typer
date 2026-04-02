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
import json

sample = "${sample_id}"
cols   = ["sample", "Plasmid", "Identity", "Query / Template length",
          "Contig", "Position in contig", "Note", "Accession number"]

try:
    with open("data.json") as fh:
        data = json.load(fh)
    hits = []
    for db_group in data.get("plasmidfinder", {}).get("results", {}).values():
        for db_hits in db_group.values():
            if not isinstance(db_hits, dict):
                continue
            for h in db_hits.values():
                row = [sample,
                       str(h.get("plasmid", "")),
                       str(h.get("identity", "")),
                       str(h.get("HSP_length", "")) + " / " + str(h.get("template_length", "")),
                       str(h.get("contig_name", "")),
                       str(h.get("positions_in_contig", "")),
                       str(h.get("note", "")),
                       str(h.get("accession", ""))]
                hits.append(chr(9).join(row))
except Exception:
    hits = []

with open("${sample_id}_plasmidfinder.tsv", "w") as out:
    print(chr(9).join(cols), file=out)
    if hits:
        for hit in hits:
            print(hit, file=out)
    else:
        print(chr(9).join([sample, "No replicons found", "NA", "NA", "NA", "NA", "NA", "NA"]), file=out)
PYEOF
    """
}
