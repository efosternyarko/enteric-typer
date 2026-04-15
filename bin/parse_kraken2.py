#!/usr/bin/env python3
"""
parse_kraken2.py — parse a Kraken2 report and flag contaminated assemblies.

A sample fails if the combined percentage of all secondary species
(species-level hits other than the dominant species) exceeds --max_secondary
percent of total sequences (including unclassified).

Output TSV columns:
    sample_id  top_species  top_pct  secondary_pct  status  reason

Usage:
    parse_kraken2.py --report <kraken2_report.txt> \\
                     --sample <sample_id> \\
                     --max_secondary 3.0 \\
                     --output <sample_kraken2_qc.tsv>
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',        required=True,
                        help='Kraken2 report file')
    parser.add_argument('--sample',        required=True,
                        help='Sample ID')
    parser.add_argument('--max_secondary', type=float, default=3.0,
                        help='Maximum allowed secondary species %% of total sequences '
                             '(default: 3.0)')
    parser.add_argument('--output',        required=True,
                        help='Output TSV path')
    args = parser.parse_args()

    # ── Parse species-level (rank S) hits ────────────────────────────────────
    species = []  # list of (pct, name)

    with open(args.report) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            rank = parts[3].strip()
            if rank != 'S':
                continue
            try:
                pct = float(parts[0].strip())
            except ValueError:
                continue
            name = parts[5].strip()
            if pct > 0:
                species.append((pct, name))

    species.sort(reverse=True)

    # ── Compute metrics ───────────────────────────────────────────────────────
    if not species:
        top_species   = 'unclassified'
        top_pct       = 0.0
        secondary_pct = 0.0
        status        = 'PASS'
        reason        = 'no_species_classified'
    else:
        top_species   = species[0][1]
        top_pct       = species[0][0]
        secondary_pct = sum(p for p, _ in species[1:])

        if secondary_pct > args.max_secondary:
            status = 'FAIL'
            reason = (f'secondary_species_{secondary_pct:.2f}pct_'
                      f'exceeds_threshold_{args.max_secondary}pct')
        else:
            status = 'PASS'
            reason = 'ok'

    # ── Write output ──────────────────────────────────────────────────────────
    with open(args.output, 'w') as out:
        out.write('sample_id\ttop_species\ttop_pct\tsecondary_pct\tstatus\treason\n')
        out.write(
            f'{args.sample}\t{top_species}\t{top_pct:.2f}\t'
            f'{secondary_pct:.2f}\t{status}\t{reason}\n'
        )

    print(
        f'{args.sample}: top={top_species} ({top_pct:.1f}%), '
        f'secondary={secondary_pct:.1f}% → {status}',
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
