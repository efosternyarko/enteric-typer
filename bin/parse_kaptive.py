#!/usr/bin/env python3
"""
parse_kaptive.py — Merge G2/G3 and G1/G4 Kaptive results into one K-type per sample.

Logic (mirrors the two-step workflow in efosternyarko/EC-K-typing-G1G4):
  - G2/G3 result is authoritative if Match confidence != "Untypeable"
  - Otherwise use G1/G4 normalised result
  - One sample can carry ONLY G2/G3 OR G1/G4 (they are mutually exclusive)
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


NA = "NA"

# Confidence rank: higher = more reliable
G2G3_CONFIDENCE_RANK = {
    "Typeable":   6,   # Kaptive v3 + Gladstone G2/G3 db vocabulary
    "Perfect":    5,
    "Very High":  4,
    "High":       3,
    "Good":       2,
    "Low":        1,
    "Untypeable": 0,
}


def read_tsv_first_data_row(path: Path) -> dict[str, str]:
    """Return first data row of a TSV as a dict (header → value)."""
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return dict(row)
    return {}


def parse_g2g3(tsv: Path) -> dict:
    row = read_tsv_first_data_row(tsv)
    locus      = (row.get("Best match locus") or "none").strip()
    ktype      = (row.get("Best match type")  or "none").strip()
    confidence = (row.get("Match confidence") or "Untypeable").strip()
    problems   = (row.get("Problems")         or "").strip()
    identity   = (row.get("Identity")         or NA).strip()
    coverage   = (row.get("Coverage")         or NA).strip()
    return {
        "k_group":          "G2/G3" if confidence != "Untypeable" else NA,
        "k_locus":          locus if confidence != "Untypeable" else NA,
        "k_type":           ktype if confidence != "Untypeable" else NA,
        "k_confidence":     confidence,
        "k_problems":       problems or NA,
        "k_identity":       identity,
        "k_coverage":       coverage,
        "kaptive_g2g3_locus":      locus,
        "kaptive_g2g3_type":       ktype,
        "kaptive_g2g3_confidence": confidence,
    }


def parse_g1g4(tsv: Path) -> dict:
    row = read_tsv_first_data_row(tsv)
    locus      = (row.get("Best match locus")      or "none").strip()
    confidence = (row.get("Best match confidence") or "Untypeable").strip()
    gene_cov   = (row.get("Gene coverage")         or NA).strip()
    norm_as    = (row.get("Norm AS")               or NA).strip()
    genes_f    = (row.get("Genes found")           or NA).strip()
    genes_e    = (row.get("Genes expected")        or NA).strip()
    return {
        "k_group":          "G1/G4" if confidence == "Typeable" else NA,
        "k_locus":          locus   if confidence == "Typeable" else NA,
        "k_type":           locus   if confidence == "Typeable" else NA,   # KL name IS the type for G1/G4
        "k_confidence":     confidence,
        "k_gene_coverage":  gene_cov,
        "k_norm_as":        norm_as,
        "k_genes_found":    genes_f,
        "k_genes_expected": genes_e,
        "kaptive_g1g4_locus":      locus,
        "kaptive_g1g4_confidence": confidence,
        "kaptive_g1g4_gene_cov":   gene_cov,
        "kaptive_g1g4_norm_as":    norm_as,
    }


COLUMNS = [
    "sample",
    # Final merged call
    "k_group",        # G2/G3 | G1/G4 | NA
    "k_locus",        # e.g. KL7, KL302
    "k_type",         # phenotypic K-type where known (e.g. K7); = KL for G1/G4
    "k_confidence",   # Kaptive confidence / Typeable
    "k_problems",
    "k_identity",     # G2/G3 only
    "k_coverage",     # G2/G3 only
    "k_gene_coverage",# G1/G4 only
    "k_norm_as",      # G1/G4 only
    "k_genes_found",
    "k_genes_expected",
    # Raw per-database results (kept for transparency)
    "kaptive_g2g3_locus", "kaptive_g2g3_type", "kaptive_g2g3_confidence",
    "kaptive_g1g4_locus", "kaptive_g1g4_confidence",
    "kaptive_g1g4_gene_cov", "kaptive_g1g4_norm_as",
]


def merge(sample_id: str, g2g3: dict, g1g4: dict) -> dict:
    """Choose the authoritative K-type call."""
    row: dict = {"sample": sample_id}

    g2g3_typed = G2G3_CONFIDENCE_RANK.get(g2g3.get("k_confidence", "Untypeable"), 0) > 0

    if g2g3_typed:
        # G2/G3 hit — use it; G1/G4 is not applicable
        row.update({
            "k_group":          g2g3["k_group"],
            "k_locus":          g2g3["k_locus"],
            "k_type":           g2g3["k_type"],
            "k_confidence":     g2g3["k_confidence"],
            "k_problems":       g2g3.get("k_problems", NA),
            "k_identity":       g2g3.get("k_identity", NA),
            "k_coverage":       g2g3.get("k_coverage", NA),
            "k_gene_coverage":  NA,
            "k_norm_as":        NA,
            "k_genes_found":    NA,
            "k_genes_expected": NA,
        })
    else:
        # G2/G3 untypeable → use G1/G4 result
        row.update({
            "k_group":          g1g4["k_group"],
            "k_locus":          g1g4["k_locus"],
            "k_type":           g1g4["k_type"],
            "k_confidence":     g1g4["k_confidence"],
            "k_problems":       NA,
            "k_identity":       NA,
            "k_coverage":       NA,
            "k_gene_coverage":  g1g4.get("k_gene_coverage", NA),
            "k_norm_as":        g1g4.get("k_norm_as", NA),
            "k_genes_found":    g1g4.get("k_genes_found", NA),
            "k_genes_expected": g1g4.get("k_genes_expected", NA),
        })

    # Always keep raw per-database columns
    for key in ("kaptive_g2g3_locus", "kaptive_g2g3_type", "kaptive_g2g3_confidence",
                "kaptive_g1g4_locus", "kaptive_g1g4_confidence",
                "kaptive_g1g4_gene_cov", "kaptive_g1g4_norm_as"):
        row[key] = g2g3.get(key) or g1g4.get(key) or NA

    return row


def main() -> None:
    p = argparse.ArgumentParser(description="Merge Kaptive G2/G3 + G1/G4 results")
    p.add_argument("--sample",  required=True)
    p.add_argument("--g2g3",    required=True, type=Path)
    p.add_argument("--g1g4",    required=True, type=Path)
    p.add_argument("--output",  required=True, type=Path)
    args = p.parse_args()

    g2g3 = parse_g2g3(args.g2g3) if args.g2g3.exists() else {}
    g1g4 = parse_g1g4(args.g1g4) if args.g1g4.exists() else {}

    row = merge(args.sample, g2g3, g1g4)

    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerow(row)

    print(f"{args.sample}: k_group={row['k_group']}  k_locus={row['k_locus']}  "
          f"k_type={row['k_type']}  confidence={row['k_confidence']}", file=sys.stderr)


if __name__ == "__main__":
    main()
