#!/usr/bin/env python3
"""
parse_assembly_stats.py — extract assembly QC metrics from seqkit stats output.

Usage:
    parse_assembly_stats.py <sample_id> <seqkit_stats.txt>

Output (stdout): one-row TSV with columns:
    sample_id  num_contigs  genome_length  assembly_N50  gc_pct
"""

import sys


def safe_int(val: str, default: int = 0) -> int:
    try:
        return int(float(val.replace(',', '').strip()))
    except (ValueError, AttributeError):
        return default


def safe_float(val: str, default: float = 0.0) -> float:
    try:
        v = val.replace(',', '').strip()
        return float(v) if v not in ('', 'N/A', '-', 'NA') else default
    except (ValueError, AttributeError):
        return default


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <sample_id> <seqkit_stats.txt>")

    sample_id  = sys.argv[1]
    stats_file = sys.argv[2]

    with open(stats_file) as fh:
        lines = [l.rstrip('\n') for l in fh if l.strip()]

    if len(lines) < 2:
        sys.exit(f"ERROR: unexpected seqkit stats output in {stats_file} "
                 f"(expected ≥2 lines, got {len(lines)})")

    header = lines[0].split('\t')
    vals   = lines[1].split('\t')
    d      = dict(zip(header, vals))

    num_contigs    = safe_int(d.get('num_seqs', '0'))
    genome_length  = safe_int(d.get('sum_len',  '0'))
    assembly_n50   = safe_int(d.get('N50',      '0'))
    gc_pct         = safe_float(d.get('GC(%)', d.get('GC', '0')))

    print('\t'.join(['sample_id', 'num_contigs', 'genome_length',
                     'assembly_N50', 'gc_pct']))
    print('\t'.join([
        sample_id,
        str(num_contigs),
        str(genome_length),
        str(assembly_n50),
        f'{gc_pct:.2f}',
    ]))


if __name__ == '__main__':
    main()
