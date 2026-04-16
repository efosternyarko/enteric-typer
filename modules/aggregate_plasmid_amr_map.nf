// ── AGGREGATE_PLASMID_AMR_MAP: merge per-sample plasmid–AMR maps ──────────────
// Concatenates all per-sample plasmid_amr_map TSVs into one summary file,
// sorted by sample_id then replicon. Rows with no AMR genes are retained
// so the full plasmid inventory is visible.

process AGGREGATE_PLASMID_AMR_MAP {
    label 'low'

    conda     "${projectDir}/envs/utils.yml"
    container 'python:3.11'

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(map_files)

    output:
    path('plasmid_amr_map.tsv'), emit: summary

    script:
    """
    python3 - << 'EOF'
import glob

rows = []
for f in glob.glob('*_plasmid_amr_map.tsv'):
    with open(f) as fh:
        lines = fh.read().strip().split('\\n')
    for line in lines[1:]:       # skip per-file header
        if line.strip():
            rows.append(line.split('\\t'))

rows.sort(key=lambda r: (r[0], r[1], r[2]))  # sample, replicon, contig

with open('plasmid_amr_map.tsv', 'w') as out:
    out.write('sample_id\\treplicon\\tcontig\\tamr_genes\\tdrug_classes\\t'
              'identity\\tcoverage\\tlikely_location\\n')
    for r in rows:
        out.write('\\t'.join(r) + '\\n')
EOF
    """
}
