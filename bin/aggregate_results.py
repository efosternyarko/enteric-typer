#!/usr/bin/env python3
"""Aggregate per-sample tool outputs into a single results table.

Supports two species modes (--species ecoli | salmonella).
Missing data for any tool is represented as 'NA'.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys


# ── Loader helpers ─────────────────────────────────────────────────────────────

def load_mlst(files: list[str]) -> dict[str, dict]:
    """Parse mlst TSV (tab-sep: sample scheme ST allele...) → {sample: {mlst_st, mlst_scheme}}."""
    results: dict[str, dict] = {}
    for f in files:
        if not _valid(f):
            continue
        with open(f) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                sid    = parts[0].strip()
                scheme = parts[1].strip() if len(parts) > 1 else "NA"
                st     = parts[2].strip() if parts[2].strip() not in ("", "-") else "NA"
                if sid:
                    results[sid] = {"mlst_scheme": scheme, "mlst_st": st}
    return results


def load_st_complexes(path: str) -> dict[str, str]:
    """Load ST → ST-complex lookup table."""
    lookup: dict[str, str] = {}
    if not _valid(path):
        return lookup
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            st = (row.get("st") or "").strip()
            cx = (row.get("st_complex") or "").strip()
            if st:
                lookup[st] = cx or "NA"
    return lookup


def load_amrfinder(files: list[str]) -> dict[str, str]:
    """Collapse per-row AMRFinder output to {sample: semicolon-list-of-gene-symbols}."""
    results: dict[str, set] = {}
    for f in files:
        if not _valid(f):
            continue
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid  = (row.get("Name") or "").strip()
                gene = (row.get("Gene symbol") or row.get("Element symbol") or "").strip()
                cls  = (row.get("Class") or "").strip()
                if sid and gene:
                    results.setdefault(sid, set()).add(gene)
    return {sid: ";".join(sorted(genes)) if genes else "NA"
            for sid, genes in results.items()}


def load_amrfinder_classes(files: list[str]) -> dict[str, str]:
    """Collapse AMRFinder output to {sample: semicolon-list-of-unique-drug-classes}."""
    results: dict[str, set] = {}
    for f in files:
        if not _valid(f):
            continue
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid = (row.get("Name") or "").strip()
                cls = (row.get("Class") or "").strip()
                if sid and cls and cls not in ("NA", ""):
                    results.setdefault(sid, set()).add(cls)
    return {sid: ";".join(sorted(classes)) if classes else "NA"
            for sid, classes in results.items()}


def load_plasmidfinder(files: list[str]) -> dict[str, str]:
    """Return {sample: semicolon-list-of-replicons}."""
    results: dict[str, set] = {}
    for f in files:
        if not _valid(f):
            continue
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid     = (row.get("sample") or "").strip()
                plasmid = (row.get("Plasmid") or "").strip()
                if sid and plasmid and plasmid not in ("No replicons found", "NA", ""):
                    results.setdefault(sid, set()).add(plasmid)
    return {sid: ";".join(sorted(r)) if r else "NA" for sid, r in results.items()}


def load_ectyper(files: list[str]) -> dict[str, dict]:
    """Parse ECTyper output → {sample: {ectyper_serotype, ectyper_O, ectyper_H, ectyper_qc}}."""
    results: dict[str, dict] = {}
    for f in files:
        if not _valid(f):
            continue
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid = (row.get("sample") or "").strip()
                if not sid:
                    continue
                results[sid] = {
                    "ectyper_O":        row.get("O-type",   "NA").strip() or "NA",
                    "ectyper_H":        row.get("H-type",   "NA").strip() or "NA",
                    "ectyper_serotype": row.get("Serotype", "NA").strip() or "NA",
                    "ectyper_qc":       row.get("QC",       "NA").strip() or "NA",
                    "ectyper_evidence": row.get("Evidence", "NA").strip() or "NA",
                }
    return results


def load_sistr(files: list[str]) -> dict[str, dict]:
    """Parse SISTR output → {sample: {sistr_serovar, sistr_O, sistr_H1, sistr_H2, sistr_cgmlst_ST}}."""
    results: dict[str, dict] = {}
    for f in files:
        if not _valid(f):
            continue
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid = (row.get("sample") or "").strip()
                if not sid:
                    continue
                results[sid] = {
                    "sistr_serovar":       _na(row.get("serovar")),
                    "sistr_serovar_antigen": _na(row.get("serovar_antigen")),
                    "sistr_serovar_cgmlst": _na(row.get("serovar_cgmlst")),
                    "sistr_O":             _na(row.get("O_antigen")),
                    "sistr_H1":            _na(row.get("H1")),
                    "sistr_H2":            _na(row.get("H2")),
                    "sistr_cgmlst_ST":     _na(row.get("cgmlst_ST")),
                    "sistr_qc":            _na(row.get("qc_status")),
                }
    return results


def load_pathogenwatch(path: str) -> dict[str, dict]:
    """Load Pathogenwatch per-sample TSV."""
    results: dict[str, dict] = {}
    if not _valid(path):
        return results
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid = (row.get("sample") or "").strip()
            if not sid:
                continue
            # Collect all cluster threshold columns dynamically
            cluster_cols = {k: v for k, v in row.items() if k.startswith("pw_cluster")}
            results[sid] = {
                "pw_status":              _na(row.get("pw_status")),
                "pw_species":             _na(row.get("pw_species")),
                "pw_genome_uuid":         _na(row.get("pw_genome_uuid")),
                "pw_collection_url":      _na(row.get("pw_collection_url")),
                "pw_cgmlst_st":           _na(row.get("pw_cgmlst_st")),
                "pw_tree_available":      _na(row.get("pw_tree_available")),
                **cluster_cols,
            }
    return results


# ── Utilities ──────────────────────────────────────────────────────────────────

def _valid(path: str) -> bool:
    return bool(path) and path not in ("NO_FILE", "null", "") and os.path.isfile(path)


def _na(v) -> str:
    return str(v).strip() if v is not None and str(v).strip() not in ("", "-", "None") else "NA"


# ── Column schemas ─────────────────────────────────────────────────────────────

ECOLI_COLUMNS = [
    "sample",
    "mlst_scheme", "mlst_st", "mlst_st_complex",
    "ectyper_O", "ectyper_H", "ectyper_serotype", "ectyper_qc", "ectyper_evidence",
    "amrfinder_genes", "amrfinder_drug_classes",
    "plasmidfinder_replicons",
    "pw_status", "pw_species", "pw_genome_uuid", "pw_collection_url",
    "pw_cgmlst_st",
    "pw_cluster5_count",   "pw_cluster5_labels",
    "pw_cluster10_count",  "pw_cluster10_labels",
    "pw_cluster20_count",  "pw_cluster20_labels",
    "pw_cluster50_count",  "pw_cluster50_labels",
    "pw_tree_available",
]

SALMONELLA_COLUMNS = [
    "sample",
    "mlst_scheme", "mlst_st", "mlst_st_complex",
    "sistr_serovar", "sistr_serovar_antigen", "sistr_serovar_cgmlst",
    "sistr_O", "sistr_H1", "sistr_H2", "sistr_cgmlst_ST", "sistr_qc",
    "amrfinder_genes", "amrfinder_drug_classes",
    "plasmidfinder_replicons",
    "pw_status", "pw_species", "pw_genome_uuid", "pw_collection_url",
    "pw_cgmlst_st",
    "pw_cluster5_count",   "pw_cluster5_labels",
    "pw_cluster10_count",  "pw_cluster10_labels",
    "pw_cluster20_count",  "pw_cluster20_labels",
    "pw_cluster50_count",  "pw_cluster50_labels",
    "pw_tree_available",
]


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate enteric-typer results")
    parser.add_argument("--species",       required=True, choices=["ecoli", "salmonella"])
    parser.add_argument("--mlst",          nargs="+", default=[])
    parser.add_argument("--amrfinder",     nargs="+", default=[])
    parser.add_argument("--serotyper",     nargs="+", default=[],
                        help="ECTyper files (ecoli) or SISTR files (salmonella)")
    parser.add_argument("--plasmidfinder", nargs="+", default=[])
    parser.add_argument("--pathogenwatch", default="NO_FILE")
    parser.add_argument("--st-complexes",  default=None)
    parser.add_argument("--output",        required=True)
    args = parser.parse_args()

    mlst_data       = load_mlst(args.mlst)
    st_lookup       = load_st_complexes(args.st_complexes)
    amr_genes       = load_amrfinder(args.amrfinder)
    amr_classes     = load_amrfinder_classes(args.amrfinder)
    plasmid_data    = load_plasmidfinder(args.plasmidfinder)
    pw_data         = load_pathogenwatch(args.pathogenwatch)

    if args.species == "ecoli":
        sero_data   = load_ectyper(args.serotyper)
        columns     = ECOLI_COLUMNS
    else:
        sero_data   = load_sistr(args.serotyper)
        columns     = SALMONELLA_COLUMNS

    # Union of all sample IDs seen across any tool
    all_samples = sorted(set(
        list(mlst_data.keys()) +
        list(amr_genes.keys()) +
        list(sero_data.keys()) +
        list(plasmid_data.keys()) +
        list(pw_data.keys())
    ))

    if not all_samples:
        print(f"WARNING: No samples found for {args.species}. Writing empty output.", file=sys.stderr)

    # Gather all cluster columns seen in pathogenwatch data
    all_pw_cluster_cols: set[str] = set()
    for pw in pw_data.values():
        all_pw_cluster_cols.update(k for k in pw if k.startswith("pw_cluster"))

    # Extend column list with any dynamic cluster columns not already listed
    extra_cols = sorted(all_pw_cluster_cols - set(columns))
    final_columns = columns + extra_cols

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=final_columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()

        for sid in all_samples:
            ml = mlst_data.get(sid, {})
            st = ml.get("mlst_st", "NA")
            pw = pw_data.get(sid, {})
            sr = sero_data.get(sid, {})

            row: dict = {
                "sample":               sid,
                "mlst_scheme":          ml.get("mlst_scheme",       "NA"),
                "mlst_st":              st,
                "mlst_st_complex":      st_lookup.get(st,           "NA"),
                "amrfinder_genes":      amr_genes.get(sid,          "NA"),
                "amrfinder_drug_classes": amr_classes.get(sid,      "NA"),
                "plasmidfinder_replicons": plasmid_data.get(sid,    "NA"),
                "pw_status":            pw.get("pw_status",         "NA"),
                "pw_species":           pw.get("pw_species",        "NA"),
                "pw_genome_uuid":       pw.get("pw_genome_uuid",    "NA"),
                "pw_collection_url":    pw.get("pw_collection_url", "NA"),
                "pw_cgmlst_st":         pw.get("pw_cgmlst_st",      "NA"),
                "pw_tree_available":    pw.get("pw_tree_available",  "False"),
            }
            # Add cluster columns from pathogenwatch
            for col in all_pw_cluster_cols:
                row[col] = pw.get(col, "NA")

            # Species-specific serotyping columns
            if args.species == "ecoli":
                row.update({
                    "ectyper_O":        sr.get("ectyper_O",        "NA"),
                    "ectyper_H":        sr.get("ectyper_H",        "NA"),
                    "ectyper_serotype": sr.get("ectyper_serotype", "NA"),
                    "ectyper_qc":       sr.get("ectyper_qc",       "NA"),
                    "ectyper_evidence": sr.get("ectyper_evidence", "NA"),
                })
            else:
                row.update({
                    "sistr_serovar":          sr.get("sistr_serovar",          "NA"),
                    "sistr_serovar_antigen":  sr.get("sistr_serovar_antigen",  "NA"),
                    "sistr_serovar_cgmlst":   sr.get("sistr_serovar_cgmlst",   "NA"),
                    "sistr_O":                sr.get("sistr_O",                "NA"),
                    "sistr_H1":               sr.get("sistr_H1",               "NA"),
                    "sistr_H2":               sr.get("sistr_H2",               "NA"),
                    "sistr_cgmlst_ST":        sr.get("sistr_cgmlst_ST",        "NA"),
                    "sistr_qc":               sr.get("sistr_qc",               "NA"),
                })

            writer.writerow(row)

    print(f"Written: {args.output} ({len(all_samples)} {args.species} samples)", file=sys.stderr)


if __name__ == "__main__":
    main()
