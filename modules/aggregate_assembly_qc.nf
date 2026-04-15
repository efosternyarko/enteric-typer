// ── AGGREGATE_ASSEMBLY_QC: combined assembly size QC summary ──────────────────
// Merges per-sample seqkit stats with species assignments and genome-size
// thresholds to produce a single summary TSV. FAIL samples are listed first.
//
// Output TSV columns:
//   sample_id  species  genome_length  status  reason

process AGGREGATE_ASSEMBLY_QC {
    label 'low'

    conda     "${projectDir}/envs/utils.yml"
    container 'python:3.11'

    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    input:
    path(stats_files)   // collected *_assembly_stats.tsv files
    path(species_map)   // two-column TSV: sample_id <TAB> species
    val(ecoli_min)
    val(ecoli_max)
    val(salmonella_min)
    val(salmonella_max)
    val(shigella_min)
    val(shigella_max)

    output:
    path('assembly_size_qc_summary.tsv'), emit: summary

    script:
    """
    python3 - << 'EOF'
import glob

# Load species map
sp_map = {}
with open('${species_map}') as f:
    for line in f:
        parts = line.rstrip('\\n').split('\\t')
        if len(parts) >= 2:
            sp_map[parts[0]] = parts[1]

thresholds = {
    'E. coli':             (${ecoli_min},      ${ecoli_max}),
    'Shigella':            (${shigella_min},    ${shigella_max}),
    'Salmonella enterica': (${salmonella_min},  ${salmonella_max}),
}

rows = []
for f in sorted(glob.glob('*_assembly_stats.tsv')):
    with open(f) as fh:
        lines = fh.read().strip().split('\\n')
    if len(lines) < 2:
        continue
    vals    = lines[1].split('\\t')
    sid     = vals[0]
    length  = int(vals[2])
    species = sp_map.get(sid, 'Unknown')
    lo, hi  = thresholds.get(species, (0, 999_999_999))
    if lo <= length <= hi:
        status = 'PASS'
        reason = 'ok'
    else:
        status = 'FAIL'
        reason = f'genome_length_{length}bp_outside_range_{lo}-{hi}bp'
    rows.append((sid, species, length, status, reason))

# FAILs first, then alphabetical by sample_id
rows.sort(key=lambda r: (0 if r[3] == 'FAIL' else 1, r[0]))

with open('assembly_size_qc_summary.tsv', 'w') as out:
    out.write('sample_id\\tspecies\\tgenome_length\\tstatus\\treason\\n')
    for r in rows:
        out.write('\\t'.join(str(x) for x in r) + '\\n')
EOF
    """
}
