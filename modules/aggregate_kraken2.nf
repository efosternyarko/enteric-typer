// ── AGGREGATE_KRAKEN2: merge per-sample Kraken2 QC results into one table ─────
// Produces a single TSV summarising contamination screening for all samples,
// sorted with FAIL samples first, then alphabetically by sample_id.

process AGGREGATE_KRAKEN2 {
    label 'low'

    conda     "${projectDir}/envs/utils.yml"
    container 'python:3.11'

    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    input:
    path(qc_files)

    output:
    path('kraken2_contamination_summary.tsv'), emit: summary

    script:
    """
    python3 - << 'EOF'
import glob, os

rows = []
for f in glob.glob('*_kraken2_qc.tsv'):
    with open(f) as fh:
        lines = fh.read().strip().split('\\n')
    if len(lines) < 2:
        continue
    rows.append(lines[1].split('\\t'))

# Sort: FAIL first, then by sample_id
rows.sort(key=lambda r: (0 if r[4] == 'FAIL' else 1, r[0]))

header = 'sample_id\\ttop_species\\ttop_pct\\tsecondary_pct\\tstatus\\treason'
with open('kraken2_contamination_summary.tsv', 'w') as out:
    out.write(header + '\\n')
    for r in rows:
        out.write('\\t'.join(r) + '\\n')
EOF
    """
}
