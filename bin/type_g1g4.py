#!/usr/bin/env python3
"""
type_g1g4.py — Normalised G1/G4 K-locus typing (independent of the Kaptive version).

Replaces the "kaptive --scores" + normalise_kaptive_scores.py two-step workflow.
Kaptive >= 3.3 removed --scores and ranks loci by gene coverage rather than
identity, which mis-assigns G1/G4 loci, so this script performs the scoring itself.

Every CDS of every reference locus is aligned to each assembly (minimap2 via
`mappy`, preset map-ont). For each locus the best alignment score of each gene is
summed (AS), then divided by the total expected CDS length of the locus
(Norm AS, ~2.0 = perfect). The locus with the highest Norm AS is the best match.
Output columns are identical to normalise_kaptive_scores.py.

Usage
-----
    python3 scripts/type_g1g4.py \
        --db  DB/EC-K-typing_group1and4_v1.3.gbk \
        --assemblies untypeable/*.fna \
        --out results_G14/kaptive_results_norm.tsv \
        [--scores results_G14/kaptive_scores.tsv] [-t 8]

Requirements: Python >= 3.8, biopython, mappy (pip install mappy biopython)
"""

import argparse
import csv
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import mappy
from Bio import SeqIO

MIN_GENE_COVERAGE = 0.50  # fraction of expected genes required for "Typeable"
MIN_GENE_QCOV = 0.50      # fraction of a gene that must align for the gene to count as found
MATCH, MISMATCH = 2, 4
GAP_PENALTIES = ((4, 2), (24, 1))  # minimap2 dual-affine gap costs (open, extend)

_GENES = None  # [(locus, gene_seq)] shared with worker processes
_LOCI = None   # {locus: total_expected_bp}


def get_locus_name(rec) -> str:
    for f in rec.features:
        if f.type == "source":
            for note in f.qualifiers.get("note", []):
                m = re.search(r"K locus:\s*(\S+)", note)
                if m:
                    return m.group(1)
    return rec.name.split("_")[0]


def load_db(db: Path):
    """Return ([(locus, gene_seq)], {locus: total_expected_CDS_bp})."""
    genes, loci = [], {}
    for rec in SeqIO.parse(str(db), "genbank"):
        name = get_locus_name(rec)
        total = 0
        for f in rec.features:
            if f.type != "CDS":
                continue
            seq = str(f.extract(rec.seq)).upper()
            genes.append((name, seq))
            total += len(seq)
        loci[name] = total
    return genes, loci


def alignment_score(hit) -> int:
    """Alignment score (minimap2 defaults: match +2, mismatch -4, dual-affine gaps) from a mappy hit."""
    score = 0
    for length, op in hit.cigar:
        if op == 0:  # M (match or mismatch)
            score += MATCH * length
        elif op in (1, 2):  # I / D
            score -= min(o + e * length for o, e in GAP_PENALTIES)
    # mappy's M includes mismatches; NM counts mismatches + gap bases
    gap_bases = sum(l for l, op in hit.cigar if op in (1, 2))
    mismatches = max(hit.NM - gap_bases, 0)
    return score - (MATCH + MISMATCH) * mismatches


def score_assembly(path: str):
    """Return (assembly_name, {locus: (AS, genes_found)})."""
    name = re.sub(r"\.(fa|fas|fasta|fna)(\.gz)?$", "", Path(path).name)
    aligner = mappy.Aligner(path, preset="map-ont", best_n=20)
    if not aligner:
        raise SystemExit(f"ERROR: could not index {path}")
    per_locus = {}
    cache = {}  # identical gene sequences are scored once
    for locus, seq in _GENES:
        if seq not in cache:
            best, best_as = None, 0
            for h in aligner.map(seq):
                a = alignment_score(h)
                if a > best_as:
                    best_as, best = a, h
            qcov = (best.q_en - best.q_st) / len(seq) if best else 0.0
            cache[seq] = (best_as, qcov)
        a, qcov = cache[seq]
        tot, found = per_locus.get(locus, (0, 0))
        per_locus[locus] = (tot + a, found + (qcov >= MIN_GENE_QCOV))
    return name, per_locus


def _init(genes):
    global _GENES
    _GENES = genes


def rank(name, per_locus, loci, genes_per_locus):
    best_locus, best_norm, best_as, best_found = "none", 0.0, 0, 0
    for locus, (a, found) in per_locus.items():
        norm = a / max(loci[locus], 1)
        if (norm, a) > (best_norm, best_as):
            best_locus, best_norm, best_as, best_found = locus, norm, a, found
    if best_as == 0:
        return {"Assembly": name, "Best match locus": "none", "Best match confidence": "Untypeable",
                "Genes found": 0, "Genes expected": 0, "Gene coverage": "0.0%", "Raw AS": 0, "Norm AS": 0.0}
    expected = genes_per_locus[best_locus]
    cov = best_found / expected if expected else 0.0
    return {"Assembly": name, "Best match locus": best_locus,
            "Best match confidence": "Typeable" if cov >= MIN_GENE_COVERAGE else "Untypeable",
            "Genes found": best_found, "Genes expected": expected, "Gene coverage": f"{100 * cov:.1f}%",
            "Raw AS": int(best_as), "Norm AS": round(best_norm, 4)}


def main():
    p = argparse.ArgumentParser(description="Normalised G1/G4 K-locus typing (no Kaptive --scores needed).")
    p.add_argument("--db", required=True, type=Path, help="G1/G4 GenBank database (.gbk)")
    p.add_argument("--assemblies", required=True, nargs="+", type=Path, help="Assembly FASTA files (optionally gzipped)")
    p.add_argument("--out", required=True, type=Path, help="Output TSV (best match per assembly)")
    p.add_argument("--scores", type=Path, help="Optional: also write the full locus x assembly score matrix")
    p.add_argument("-t", "--threads", type=int, default=1)
    args = p.parse_args()

    for f in [args.db] + args.assemblies:
        if not f.exists():
            sys.exit(f"ERROR: not found: {f}")

    genes, loci = load_db(args.db)
    genes_per_locus = {}
    for locus, _ in genes:
        genes_per_locus[locus] = genes_per_locus.get(locus, 0) + 1
    print(f"[1] {len(loci)} loci, {len(genes)} reference genes loaded from {args.db.name}")

    paths = [str(a) for a in args.assemblies]
    print(f"[2] Scoring {len(paths)} assemblies ...")
    with ProcessPoolExecutor(max_workers=args.threads, initializer=_init, initargs=(genes,)) as ex:
        results = list(ex.map(score_assembly, paths))

    rows = [rank(n, pl, loci, genes_per_locus) for n, pl in results]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    if args.scores:
        with open(args.scores, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["Assembly", "Locus", "AS", "genes_found", "genes_expected", "total_expected_bp"])
            for n, pl in results:
                for locus, (a, found) in pl.items():
                    w.writerow([n, locus, a, found, genes_per_locus[locus], loci[locus]])

    typeable = sum(r["Best match confidence"] == "Typeable" for r in rows)
    print(f"[3] Done — {typeable}/{len(rows)} assemblies typeable\n    Output: {args.out}")


if __name__ == "__main__":
    main()
