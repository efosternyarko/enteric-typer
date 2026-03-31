#!/usr/bin/env python3
"""
normalise_kaptive_scores.py — Re-rank Kaptive scores by AS / total_expected_gene_bp.

Kaptive's raw alignment score (AS) accumulates across all reference genes, so
large loci always score higher than small ones even at identical per-base identity.
This script corrects for that by dividing AS by the total expected CDS length of
each reference locus, producing a size-independent per-base identity metric
(Norm AS). A value of ~2.0 indicates a perfect 100% identity match.

Usage:
    python3 normalise_kaptive_scores.py \\
        --db  DB/EC-K-typing_group1and4_v0.9.gbk \\
        --in  results_G14/kaptive_scores.tsv \\
        --out results_G14/kaptive_results_norm.tsv

Arguments:
    --db            Path to the G1/G4 Kaptive GenBank database
    --in            Kaptive --scores output TSV (locus × assembly score matrix)
    --out           Output TSV path for normalised results
    --min-coverage  Minimum fraction of expected genes to call Typeable (default: 0.50)
"""

import argparse
import re

import pandas as pd
from Bio import SeqIO


def load_locus_bp(db_path: str) -> dict:
    """Return {locus_name: total_CDS_bp} from a Kaptive GenBank database."""
    bp = {}
    for rec in SeqIO.parse(db_path, "genbank"):
        lname = rec.name.split("_")[0]
        for f in rec.features:
            if f.type == "source":
                for note in f.qualifiers.get("note", []):
                    m = re.search(r"K locus:\s*(\S+)", note)
                    if m:
                        lname = m.group(1)
        cds = [f for f in rec.features if f.type == "CDS"]
        bp[lname] = sum(len(f) for f in cds)
    return bp


def normalise(scores_tsv: str, bp: dict, min_coverage: float) -> pd.DataFrame:
    df = pd.read_csv(scores_tsv, sep="\t")
    asm_col = df.columns[0]

    df["total_bp"] = df["Locus"].map(bp).fillna(df["q_len"])
    df["AS_norm"]  = df["AS"] / df["total_bp"].clip(lower=1)

    rows = []
    for asm, grp in df.groupby(asm_col):
        active = grp[grp["AS"] > 0]
        if active.empty:
            rows.append({
                "Assembly": asm, "Best match locus": "none",
                "Best match confidence": "Untypeable",
                "Genes found": 0, "Genes expected": 0,
                "Gene coverage": "0.0%", "Raw AS": 0, "Norm AS": 0.0,
            })
            continue
        # Rank by Norm AS descending; use raw AS as tiebreaker
        best = active.sort_values(["AS_norm", "AS"], ascending=[False, False]).iloc[0]
        gf   = int(best["genes_found"])
        ge   = int(best["genes_expected"])
        cov  = gf / ge if ge > 0 else 0
        conf = "Typeable" if cov >= min_coverage else "Untypeable"
        rows.append({
            "Assembly": asm,
            "Best match locus": best["Locus"],
            "Best match confidence": conf,
            "Genes found": gf,
            "Genes expected": ge,
            "Gene coverage": f"{100 * cov:.1f}%",
            "Raw AS": int(best["AS"]),
            "Norm AS": round(best["AS_norm"], 4),
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db",  required=True, metavar="PATH",
                    help="Kaptive GenBank database (.gbk)")
    ap.add_argument("--in",  required=True, dest="scores_tsv", metavar="PATH",
                    help="Kaptive --scores output TSV")
    ap.add_argument("--out", required=True, metavar="PATH",
                    help="Output TSV for normalised results")
    ap.add_argument("--min-coverage", type=float, default=0.50, metavar="FLOAT",
                    help="Min gene coverage to call Typeable (default: 0.50)")
    args = ap.parse_args()

    print(f"Loading locus stats from {args.db}...")
    bp = load_locus_bp(args.db)
    print(f"  {len(bp)} loci loaded")

    print(f"Normalising scores from {args.scores_tsv}...")
    results = normalise(args.scores_tsv, bp, args.min_coverage)

    n_typeable   = (results["Best match confidence"] == "Typeable").sum()
    n_untypeable = (results["Best match confidence"] == "Untypeable").sum()
    print(f"  {len(results)} assemblies: {n_typeable} Typeable, {n_untypeable} Untypeable")

    results.to_csv(args.out, sep="\t", index=False)
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
