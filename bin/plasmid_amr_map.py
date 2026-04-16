#!/usr/bin/env python3
"""
plasmid_amr_map.py — link plasmid replicon types to co-located AMR genes.

Joins PlasmidFinder and AMRFinder outputs on contig ID to produce a table
showing which AMR genes share a contig with each identified replicon type.
AMR genes on contigs with no replicon are reported as 'no_replicon'.

Output TSV columns:
    sample_id  replicon  contig  amr_genes  drug_classes  identity  coverage

Usage:
    plasmid_amr_map.py --plasmidfinder <tsv> --amrfinder <tsv>
                       --sample <id> --output <tsv>
"""

import argparse
from collections import defaultdict


def normalise_contig(name):
    """Return the first whitespace-delimited token — handles 'ctg1 len=5000' etc."""
    return name.strip().split()[0] if name.strip() else name.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plasmidfinder', required=True)
    parser.add_argument('--amrfinder',     required=True)
    parser.add_argument('--sample',        required=True)
    parser.add_argument('--fasta',         required=True,
                        help='Assembly FASTA — used to determine contig lengths')
    parser.add_argument('--output',        required=True)
    # Size thresholds for likely_location inference
    parser.add_argument('--chrom_min_bp',  type=int, default=1_000_000,
                        help='Contigs >= this length are called chromosome (default 1 Mb)')
    parser.add_argument('--plasmid_max_bp', type=int, default=500_000,
                        help='No-replicon contigs < this length are called plasmid_fragment (default 500 kb)')
    args = parser.parse_args()

    # ── Contig lengths from FASTA ─────────────────────────────────────────────
    contig_lengths = {}
    current, length = None, 0
    with open(args.fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                if current is not None:
                    contig_lengths[normalise_contig(current)] = length
                current = line[1:].strip()
                length = 0
            else:
                length += len(line.strip())
    if current is not None:
        contig_lengths[normalise_contig(current)] = length

    def likely_location(contig, has_replicon):
        bp = contig_lengths.get(contig, 0)
        if has_replicon:
            return 'plasmid'
        if bp >= args.chrom_min_bp:
            return 'chromosome'
        if bp < args.plasmid_max_bp:
            return 'plasmid_fragment'
        return 'unknown'

    # ── PlasmidFinder: contig → list of (replicon, identity, coverage) ────────
    # Columns: sample Plasmid Identity "Query / Template length" Contig ...
    pf_contigs = defaultdict(list)   # normalised_contig -> [(replicon, pct_id, pct_cov)]

    with open(args.plasmidfinder) as fh:
        for line in fh:
            if line.startswith('sample') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            replicon  = parts[1].strip()
            pct_id    = parts[2].strip()
            qt_len    = parts[3].strip()          # "query / template length"
            raw_contig = parts[4].strip()
            contig    = normalise_contig(raw_contig)

            # Compute coverage from "query / template length"
            try:
                q, t = [int(x.strip()) for x in qt_len.split('/')]
                pct_cov = f'{100 * q / t:.1f}' if t > 0 else 'NA'
            except (ValueError, ZeroDivisionError):
                pct_cov = 'NA'

            pf_contigs[contig].append((replicon, pct_id, pct_cov))

    # ── AMRFinder: contig → list of (gene, drug_class) ────────────────────────
    # Columns: Name ProteinId ContigId Start Stop Strand GeneSymbol SeqName
    #          Scope ElementType ElementSubtype Class Subclass ...
    amr_contigs = defaultdict(list)   # normalised_contig -> [(gene, drug_class)]

    with open(args.amrfinder) as fh:
        header = None
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                # Build column index
                col = {h: i for i, h in enumerate(header)}
                continue
            if len(parts) < 12:
                continue

            element_type = parts[col.get('Element type', 9)].strip()
            if element_type != 'AMR':
                continue   # skip virulence / stress entries for this table

            gene        = parts[col.get('Gene symbol', 6)].strip()
            raw_contig  = parts[col.get('Contig id', 2)].strip()
            contig      = normalise_contig(raw_contig)
            drug_class  = parts[col.get('Class', 11)].strip()

            if gene:
                amr_contigs[contig].append((gene, drug_class))

    # ── Build output rows ─────────────────────────────────────────────────────
    rows = []

    # Contigs with a replicon hit
    for contig, pf_hits in sorted(pf_contigs.items()):
        amr_on_contig = amr_contigs.get(contig, [])
        genes  = ';'.join(sorted(set(g for g, _ in amr_on_contig))) or '-'
        classes = ';'.join(sorted(set(c for _, c in amr_on_contig))) or '-'

        for replicon, pct_id, pct_cov in pf_hits:
            rows.append((args.sample, replicon, contig, genes, classes,
                         pct_id, pct_cov, 'plasmid'))

    # AMR genes on contigs with no replicon identified
    replicon_contigs = set(pf_contigs.keys())
    for contig, amr_hits in sorted(amr_contigs.items()):
        if contig in replicon_contigs:
            continue
        genes   = ';'.join(sorted(set(g for g, _ in amr_hits)))
        classes = ';'.join(sorted(set(c for _, c in amr_hits)))
        loc     = likely_location(contig, has_replicon=False)
        rows.append((args.sample, 'no_replicon', contig, genes, classes,
                     '-', '-', loc))

    # ── Write output ──────────────────────────────────────────────────────────
    with open(args.output, 'w') as out:
        out.write('sample_id\treplicon\tcontig\tamr_genes\tdrug_classes\t'
                  'identity\tcoverage\tlikely_location\n')
        for r in rows:
            out.write('\t'.join(str(x) for x in r) + '\n')


if __name__ == '__main__':
    main()
