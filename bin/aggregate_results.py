#!/usr/bin/env python3
"""Aggregate per-sample tool outputs into a single results table.

Supports two species modes (--species ecoli | salmonella).
Missing data for any tool is represented as 'NA'.

AMRrules integration
--------------------
When --amrrules is supplied, AMRFinder gene lists are split into:
  amrfinder_acquired_genes   — genes NOT flagged as wildtype by AMRrules
                               (these confer clinical resistance)
  amrfinder_intrinsic_genes  — genes flagged as wildtype/intrinsic by AMRrules
                               (present in most strains; do NOT confer resistance)
  amrfinder_genes            — raw unfiltered list (all genes, for reference)
  amrfinder_drug_classes     — drug classes from acquired genes only
"""

from __future__ import annotations

import argparse
import csv
import os
import sys


# ── AMRrules loader ────────────────────────────────────────────────────────────

def load_amrrules(path: str) -> frozenset[str]:
    """Return set of gene symbols classified as wildtype (phenotype=wildtype)."""
    if not _valid(path):
        return frozenset()
    genes: set[str] = set()
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("phenotype", "").strip().lower() == "wildtype":
                gene = row.get("gene", "").strip()
                if gene:
                    genes.add(gene)
    return frozenset(genes)


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


AMR_TYPES      = {"AMR", "POINT"}   # element types that count as resistance
VIRULENCE_TYPES = {"VIRULENCE"}      # element types that are virulence factors
STRESS_TYPES    = {"STRESS"}         # stress response — kept separate, excluded from AMR plot


def load_amrfinder(files: list[str],
                   wildtype_genes: frozenset[str] = frozenset()
                   ) -> dict[str, dict[str, str]]:
    """
    Parse AMRFinder output, splitting by element Type:

      amrfinder_acquired_genes   — AMR/POINT type, not in AMRrules wildtype set
      amrfinder_intrinsic_genes  — AMR/POINT type, in AMRrules wildtype set
      amrfinder_virulence_genes  — VIRULENCE type (separate column, not mixed into AMR)
      amrfinder_genes            — all genes raw (for reference)
    """
    raw: dict[str, dict[str, set]] = {}

    for f in files:
        if not _valid(f):
            continue
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid      = (row.get("Name") or "").strip()
                gene     = (row.get("Gene symbol") or row.get("Element symbol") or "").strip()
                etype    = (row.get("Type") or row.get("Element type") or "AMR").strip().upper()
                if not sid or not gene:
                    continue
                d = raw.setdefault(sid, {
                    "all": set(), "acquired": set(),
                    "intrinsic": set(), "virulence": set(),
                })
                d["all"].add(gene)
                if etype in VIRULENCE_TYPES:
                    d["virulence"].add(gene)
                elif etype in AMR_TYPES:
                    if gene in wildtype_genes:
                        d["intrinsic"].add(gene)
                    else:
                        d["acquired"].add(gene)
                # STRESS and other types are captured in 'all' only

    return {
        sid: {
            "all":       ";".join(sorted(d["all"]))       if d["all"]       else "NA",
            "acquired":  ";".join(sorted(d["acquired"]))  if d["acquired"]  else "NA",
            "intrinsic": ";".join(sorted(d["intrinsic"])) if d["intrinsic"] else "NA",
            "virulence": ";".join(sorted(d["virulence"])) if d["virulence"] else "NA",
        }
        for sid, d in raw.items()
    }


def load_amrfinder_classes(files: list[str],
                            wildtype_genes: frozenset[str] = frozenset()
                            ) -> dict[str, str]:
    """
    Collapse AMRFinder output to {sample: semicolon-list-of-unique-drug-classes}.
    Only acquired AMR/POINT-type (non-intrinsic, non-virulence) genes are counted.
    """
    results: dict[str, set] = {}
    for f in files:
        if not _valid(f):
            continue
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid   = (row.get("Name") or "").strip()
                gene  = (row.get("Gene symbol") or row.get("Element symbol") or "").strip()
                etype = (row.get("Type") or row.get("Element type") or "AMR").strip().upper()
                cls   = (row.get("Class") or "").strip()
                if sid and gene and cls and cls not in ("NA", ""):
                    if etype in AMR_TYPES and gene not in wildtype_genes:
                        results.setdefault(sid, set()).add(cls)
    return {sid: ";".join(sorted(classes)) if classes else "NA"
            for sid, classes in results.items()}


def load_amrfinder_gene_classes(files: list[str],
                                 wildtype_genes: frozenset[str] = frozenset()
                                 ) -> dict[str, str]:
    """
    Return {sample: 'gene=CLASS;gene2=CLASS2;...'} for acquired AMR/POINT genes.
    Used by plot_summary.py to draw gene-filled drug class bars.
    """
    raw: dict[str, dict[str, str]] = {}   # {sample: {gene: class}}
    for f in files:
        if not _valid(f):
            continue
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid   = (row.get("Name") or "").strip()
                gene  = (row.get("Gene symbol") or row.get("Element symbol") or "").strip()
                etype = (row.get("Type") or row.get("Element type") or "AMR").strip().upper()
                cls   = (row.get("Class") or "").strip().upper()
                if sid and gene and cls and cls not in ("NA", ""):
                    if etype in AMR_TYPES and gene not in wildtype_genes:
                        raw.setdefault(sid, {})[gene] = cls
    return {
        sid: ";".join(f"{g}={c}" for g, c in sorted(d.items())) if d else "NA"
        for sid, d in raw.items()
    }


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
    """Parse SISTR output → {sample: {sistr_serovar, ...}}."""
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
                    "sistr_serovar":         _na(row.get("serovar")),
                    "sistr_serovar_antigen": _na(row.get("serovar_antigen")),
                    "sistr_serovar_cgmlst":  _na(row.get("serovar_cgmlst")),
                    "sistr_O":               _na(row.get("o_antigen")),
                    "sistr_H1":              _na(row.get("h1")),
                    "sistr_H2":              _na(row.get("h2")),
                    "sistr_cgmlst_ST":       _na(row.get("cgmlst_ST")),
                    "sistr_qc":              _na(row.get("qc_status")),
                }
    return results


def load_clermont(files: list[str]) -> dict[str, str]:
    """
    Parse EzClermont per-sample TSV files.
    Expected columns: sample, clermont_phylogroup
    Returns {sample_id: phylogroup_string}.
    """
    result: dict[str, str] = {}
    for f in files:
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid = (row.get("sample") or "").strip()
                pg  = (row.get("clermont_phylogroup") or "Unknown").strip()
                if sid:
                    result[sid] = pg if pg not in ("", "-", "NA", "nan") else "Unknown"
    return result


def load_kleborate(files: list[str]) -> dict[str, dict]:
    """
    Parse Kleborate v3 --preset escherichia output (with --trim_headers).

    Extracts:
      clermont_phylogroup   — Clermont type (clermont_type column; also exposed as kleborate_phylogroup for legacy)
      kleborate_pathovar    — pathotype: STEC/EPEC/ETEC/EIEC/EHEC or '-'
      kleborate_Stx1        — Stx1 gene hits (';'-separated) or '-'
      kleborate_Stx2        — Stx2 gene hits
      kleborate_eae         — intimin '+'/'-'
      kleborate_ipaH        — ipaH invasion gene '+'/'-'
      kleborate_LT          — heat-labile enterotoxin '+'/'-'
      kleborate_ST_toxin    — heat-stable enterotoxin '+'/'-' (renamed to avoid ST ambiguity)
    """
    def _get(row: dict, *keys: str) -> str:
        """Try multiple possible column names (handles full and trimmed header formats)."""
        for k in keys:
            v = row.get(k)
            if v is not None:
                return v.strip()
        return ""

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
                    # Phylogroup (from Kleborate's internal EzClermont run; used only as fallback)
                    "clermont_phylogroup": _na(
                        _get(row, "clermont_type",
                             "escherichia__ezclermont__clermont_type") or None
                    ),
                    # Pathotype
                    "kleborate_pathovar": _na(
                        _get(row, "Pathotype",
                             "escherichia__pathovar__Pathotype") or None
                    ),
                    # Pathovar markers
                    "kleborate_Stx1": _na(
                        _get(row, "Stx1", "escherichia__pathovar__Stx1") or None
                    ),
                    "kleborate_Stx2": _na(
                        _get(row, "Stx2", "escherichia__pathovar__Stx2") or None
                    ),
                    "kleborate_eae": _na(
                        _get(row, "eae", "escherichia__pathovar__eae") or None
                    ),
                    "kleborate_ipaH": _na(
                        _get(row, "ipaH", "escherichia__pathovar__ipaH") or None
                    ),
                    "kleborate_LT": _na(
                        _get(row, "LT", "escherichia__pathovar__LT") or None
                    ),
                    # 'ST' in Kleborate pathovar module = heat-stable enterotoxin;
                    # do NOT fall back to the bare "ST" key, which is the MLST
                    # sequence-type column and would pollute this field with ST IDs.
                    "kleborate_ST_toxin": _na(
                        _get(row, "escherichia__pathovar__ST") or None
                    ),
                }
    return results


def load_abricate(files: list[str], mincov: float = 80.0, minid: float = 90.0) -> dict[str, str]:
    """
    Parse abricate output → {sample: semicolon-list-of-gene-names}.
    Applies optional mincov / minid thresholds (abricate's own --mincov/--minid
    already filter, but this provides a second check if defaults differ).
    """
    raw: dict[str, set] = {}
    for f in files:
        if not _valid(f):
            continue
        with open(f) as fh:
            # abricate writes a header starting with '#FILE'; strip leading #
            first = fh.readline().lstrip("#").strip()
            headers = first.split("\t")
            reader = csv.DictReader(fh, fieldnames=headers, delimiter="\t")
            for row in reader:
                sid  = (row.get("sample") or "").strip()
                gene = (row.get("GENE") or "").strip()
                try:
                    cov = float(row.get("%COVERAGE") or 0)
                    pid = float(row.get("%IDENTITY") or 0)
                except ValueError:
                    continue
                if sid and gene and cov >= mincov and pid >= minid:
                    raw.setdefault(sid, set()).add(gene)
    return {sid: ";".join(sorted(genes)) if genes else "NA"
            for sid, genes in raw.items()}


def load_ktype(files: list[str]) -> dict[str, dict]:
    """Load parse_kaptive per-sample TSV files → {sample: {k_group, k_locus, k_type, ...}}."""
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
                    "k_group":      _na(row.get("k_group")),
                    "k_locus":      _na(row.get("k_locus")),
                    "k_type":       _na(row.get("k_type")),
                    "k_confidence": _na(row.get("k_confidence")),
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
            cluster_cols = {k: v for k, v in row.items() if k.startswith("pw_cluster")}
            results[sid] = {
                "pw_status":         _na(row.get("pw_status")),
                "pw_species":        _na(row.get("pw_species")),
                "pw_genome_uuid":    _na(row.get("pw_genome_uuid")),
                "pw_collection_url": _na(row.get("pw_collection_url")),
                "pw_cgmlst_st":      _na(row.get("pw_cgmlst_st")),
                "pw_tree_available": _na(row.get("pw_tree_available")),
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
    "clermont_phylogroup",
    "kleborate_pathovar",
    "kleborate_Stx1", "kleborate_Stx2", "kleborate_eae",
    "kleborate_ipaH", "kleborate_LT", "kleborate_ST_toxin",
    "ectyper_O", "ectyper_H", "ectyper_serotype", "ectyper_qc", "ectyper_evidence",
    "k_group", "k_locus", "k_type", "k_confidence",
    "amrfinder_acquired_genes",
    "amrfinder_intrinsic_genes",
    "amrfinder_virulence_genes",
    "amrfinder_genes",
    "amrfinder_drug_classes",
    "amrfinder_gene_classes",
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
    "amrfinder_acquired_genes",
    "amrfinder_intrinsic_genes",
    "amrfinder_virulence_genes",
    "amrfinder_genes",
    "amrfinder_drug_classes",
    "amrfinder_gene_classes",
    "plasmidfinder_replicons",
    "pw_status", "pw_species", "pw_genome_uuid", "pw_collection_url",
    "pw_cgmlst_st",
    "pw_cluster5_count",   "pw_cluster5_labels",
    "pw_cluster10_count",  "pw_cluster10_labels",
    "pw_cluster20_count",  "pw_cluster20_labels",
    "pw_cluster50_count",  "pw_cluster50_labels",
    "pw_tree_available",
]

SHIGELLA_COLUMNS = [
    "sample",
    "mlst_scheme", "mlst_st", "mlst_st_complex",
    "shigeifinder_ipaH", "shigeifinder_virulence_plasmid",
    "shigeifinder_cluster", "shigeifinder_serotype",
    "shigeifinder_o_antigen", "shigeifinder_h_antigen",
    "mykrobe_genotype", "mykrobe_lineage", "mykrobe_clade",
    "mykrobe_subclade", "mykrobe_genotype_name", "mykrobe_confidence",
    "amrfinder_acquired_genes",
    "amrfinder_intrinsic_genes",
    "amrfinder_virulence_genes",
    "amrfinder_genes",
    "amrfinder_drug_classes",
    "amrfinder_gene_classes",
    "plasmidfinder_replicons",
    "pinv_present", "pinv_genes",
    "is_elements",
]


def load_shigeifinder(files: list) -> dict:
    """Parse ShigEiFinder output → per-sample dict of typed fields.

    ShigEiFinder v1.3.5 output columns (all uppercase):
      SAMPLE  ipaH  VIRULENCE_PLASMID  CLUSTER  SEROTYPE  O_ANTIGEN  H_ANTIGEN  NOTES
    """
    data: dict = {}
    for path in files:
        try:
            with open(path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    # ShigEiFinder v1.3.5 uses "#SAMPLE" as the header (comment-style)
                    sid = (str(row.get("#SAMPLE", "")
                               or row.get("SAMPLE", "")
                               or row.get("sample", "")).strip())
                    if not sid:
                        continue
                    def _get(col: str) -> str:
                        return str(row.get(col, "NA")).strip() or "NA"
                    data[sid] = {
                        "shigeifinder_ipaH":              _get("ipaH"),
                        "shigeifinder_virulence_plasmid": _get("VIRULENCE_PLASMID"),
                        "shigeifinder_cluster":           _get("CLUSTER"),
                        "shigeifinder_serotype":          _get("SEROTYPE"),
                        "shigeifinder_o_antigen":         _get("O_ANTIGEN"),
                        "shigeifinder_h_antigen":         _get("H_ANTIGEN"),
                    }
        except Exception:
            pass
    return data


def load_mykrobe_tsv(files: list) -> dict:
    """Parse parse_mykrobe.py TSV → {sample: mykrobe field dict}."""
    data: dict = {}
    for path in files:
        try:
            with open(path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    sid = str(row.get("sample", "")).strip()
                    if not sid:
                        continue
                    data[sid] = {k: (str(v).strip() or "NA") for k, v in row.items()}
        except Exception:
            pass
    return data


def load_pinv(files: list) -> dict:
    """Parse pINV screen output → {sample: {pinv_present, pinv_genes}}."""
    genes_by_sample: dict = {}
    for path in files:
        try:
            with open(path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    sid  = str(row.get("sample", "")).strip()
                    gene = str(row.get("gene",   "")).strip()
                    if sid and gene and gene != "NA":
                        genes_by_sample.setdefault(sid, set()).add(gene)
        except Exception:
            pass
    return {
        sid: {"pinv_present": "Y", "pinv_genes": ";".join(sorted(gs))}
        for sid, gs in genes_by_sample.items()
    }


def load_is_screen(files: list) -> dict:
    """Parse IS screen output → {sample: 'IS1(3);IS30(1)' string}."""
    counts: dict = {}
    for path in files:
        try:
            with open(path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    sid     = str(row.get("sample",     "")).strip()
                    elem    = str(row.get("IS_element", "")).strip()
                    copies  = str(row.get("copies",     "")).strip()
                    if sid and elem and elem != "NA":
                        counts.setdefault(sid, {})[elem] = copies
        except Exception:
            pass
    return {
        sid: ";".join(f"{e}({c})" for e, c in sorted(d.items()))
        for sid, d in counts.items()
    }


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate enteric-typer results")
    parser.add_argument("--species",       required=True,
                        choices=["ecoli", "salmonella", "shigella"])
    parser.add_argument("--mlst",          nargs="+", default=[])
    parser.add_argument("--amrfinder",     nargs="+", default=[])
    parser.add_argument("--serotyper",     nargs="+", default=[],
                        help="ECTyper (ecoli) / SISTR (salmonella) / ShigEiFinder (shigella)")
    parser.add_argument("--plasmidfinder", nargs="+", default=[])
    parser.add_argument("--ktype",         nargs="+", default=[],
                        help="parse_kaptive output TSV files (E. coli only)")
    parser.add_argument("--kleborate",     nargs="+", default=[],
                        help="Kleborate output TSV files (E. coli only)")
    parser.add_argument("--clermont",      nargs="+", default=[],
                        help="EzClermont phylotyping TSV files (E. coli only)")
    parser.add_argument("--mykrobe",       nargs="+", default=[],
                        help="parse_mykrobe.py TSV files (Shigella only)")
    parser.add_argument("--pinv",          nargs="+", default=[],
                        help="pINV screen TSV files (Shigella only)")
    parser.add_argument("--is-screen",     nargs="+", default=[], dest="is_screen",
                        help="IS element screen TSV files (Shigella only)")
    parser.add_argument("--pathogenwatch", default="NO_FILE")
    parser.add_argument("--st-complexes",  default=None)
    parser.add_argument("--amrrules",      default=None,
                        help="AMRrules TSV for this species (assets/amrrules/*.tsv)")
    parser.add_argument("--output",        required=True)
    args = parser.parse_args()

    # Load AMRrules wildtype gene set (empty frozenset if not provided)
    wildtype_genes = load_amrrules(args.amrrules) if args.amrrules else frozenset()
    if wildtype_genes:
        print(f"INFO: AMRrules loaded — {len(wildtype_genes)} intrinsic gene(s) will be flagged: "
              f"{', '.join(sorted(wildtype_genes))}", file=sys.stderr)
    else:
        print("INFO: No AMRrules file provided — amrfinder_intrinsic_genes will be NA for all samples.",
              file=sys.stderr)

    mlst_data      = load_mlst(args.mlst)
    st_lookup      = load_st_complexes(args.st_complexes)
    amr_data       = load_amrfinder(args.amrfinder, wildtype_genes)
    amr_classes    = load_amrfinder_classes(args.amrfinder, wildtype_genes)
    amr_gene_cls   = load_amrfinder_gene_classes(args.amrfinder, wildtype_genes)
    plasmid_data   = load_plasmidfinder(args.plasmidfinder)
    ktype_data     = load_ktype(args.ktype)         if args.species == "ecoli" else {}
    kleb_data      = load_kleborate(args.kleborate) if args.species == "ecoli" else {}
    clermont_data  = load_clermont(args.clermont)   if args.species == "ecoli" else {}
    mykrobe_data   = load_mykrobe_tsv(args.mykrobe) if args.species == "shigella" else {}
    pinv_data      = load_pinv(args.pinv)            if args.species == "shigella" else {}
    is_data        = load_is_screen(args.is_screen)  if args.species == "shigella" else {}
    pw_data        = load_pathogenwatch(args.pathogenwatch)

    if args.species == "ecoli":
        sero_data = load_ectyper(args.serotyper)
        columns   = ECOLI_COLUMNS
    elif args.species == "salmonella":
        sero_data = load_sistr(args.serotyper)
        columns   = SALMONELLA_COLUMNS
    else:
        sero_data = load_shigeifinder(args.serotyper)
        columns   = SHIGELLA_COLUMNS

    all_samples = sorted(set(
        list(mlst_data.keys()) +
        list(amr_data.keys()) +
        list(sero_data.keys()) +
        list(plasmid_data.keys()) +
        list(ktype_data.keys()) +
        list(kleb_data.keys()) +
        list(mykrobe_data.keys()) +
        list(pw_data.keys())
    ))

    if not all_samples:
        print(f"WARNING: No samples found for {args.species}. Writing empty output.", file=sys.stderr)

    all_pw_cluster_cols: set[str] = set()
    for pw in pw_data.values():
        all_pw_cluster_cols.update(k for k in pw if k.startswith("pw_cluster"))

    extra_cols   = sorted(all_pw_cluster_cols - set(columns))
    final_columns = columns + extra_cols

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=final_columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()

        for sid in all_samples:
            ml   = mlst_data.get(sid, {})
            st   = ml.get("mlst_st", "NA")
            pw   = pw_data.get(sid, {})
            sr   = sero_data.get(sid, {})
            kt   = ktype_data.get(sid, {})
            kleb = kleb_data.get(sid, {})
            # EzClermont phylogroup takes priority over Kleborate's (often absent) clermont_type
            clermont_pg = clermont_data.get(sid, None)
            amr  = amr_data.get(sid, {})

            row: dict = {
                "sample":                   sid,
                "mlst_scheme":              ml.get("mlst_scheme",       "NA"),
                "mlst_st":                  st,
                "mlst_st_complex":          st_lookup.get(st,           "NA"),
                "amrfinder_acquired_genes":   amr.get("acquired",   "NA"),
                "amrfinder_intrinsic_genes":  amr.get("intrinsic",  "NA"),
                "amrfinder_virulence_genes":  amr.get("virulence",  "NA"),
                "amrfinder_genes":            amr.get("all",        "NA"),
                "amrfinder_drug_classes":   amr_classes.get(sid,        "NA"),
                "amrfinder_gene_classes":   amr_gene_cls.get(sid,       "NA"),
                "plasmidfinder_replicons":  plasmid_data.get(sid,       "NA"),
                "pw_status":               pw.get("pw_status",          "NA"),
                "pw_species":              pw.get("pw_species",         "NA"),
                "pw_genome_uuid":          pw.get("pw_genome_uuid",     "NA"),
                "pw_collection_url":       pw.get("pw_collection_url",  "NA"),
                "pw_cgmlst_st":            pw.get("pw_cgmlst_st",       "NA"),
                "pw_tree_available":       pw.get("pw_tree_available",  "False"),
            }
            for col in all_pw_cluster_cols:
                row[col] = pw.get(col, "NA")

            if args.species == "ecoli":
                row.update({
                    "clermont_phylogroup":   (clermont_pg if clermont_pg
                                              else kleb.get("clermont_phylogroup", "NA")),
                    "kleborate_pathovar":    kleb.get("kleborate_pathovar",    "NA"),
                    "kleborate_Stx1":        kleb.get("kleborate_Stx1",        "NA"),
                    "kleborate_Stx2":        kleb.get("kleborate_Stx2",        "NA"),
                    "kleborate_eae":         kleb.get("kleborate_eae",         "NA"),
                    "kleborate_ipaH":        kleb.get("kleborate_ipaH",        "NA"),
                    "kleborate_LT":          kleb.get("kleborate_LT",          "NA"),
                    "kleborate_ST_toxin":    kleb.get("kleborate_ST_toxin",    "NA"),
                })
                row.update({
                    "k_group":      kt.get("k_group",      "NA"),
                    "k_locus":      kt.get("k_locus",      "NA"),
                    "k_type":       kt.get("k_type",       "NA"),
                    "k_confidence": kt.get("k_confidence", "NA"),
                })
                row.update({
                    "ectyper_O":        sr.get("ectyper_O",        "NA"),
                    "ectyper_H":        sr.get("ectyper_H",        "NA"),
                    "ectyper_serotype": sr.get("ectyper_serotype", "NA"),
                    "ectyper_qc":       sr.get("ectyper_qc",       "NA"),
                    "ectyper_evidence": sr.get("ectyper_evidence", "NA"),
                })
            elif args.species == "salmonella":
                row.update({
                    "sistr_serovar":         sr.get("sistr_serovar",         "NA"),
                    "sistr_serovar_antigen": sr.get("sistr_serovar_antigen", "NA"),
                    "sistr_serovar_cgmlst":  sr.get("sistr_serovar_cgmlst",  "NA"),
                    "sistr_O":               sr.get("sistr_O",               "NA"),
                    "sistr_H1":              sr.get("sistr_H1",              "NA"),
                    "sistr_H2":              sr.get("sistr_H2",              "NA"),
                    "sistr_cgmlst_ST":       sr.get("sistr_cgmlst_ST",       "NA"),
                    "sistr_qc":              sr.get("sistr_qc",              "NA"),
                })
            else:  # shigella
                mk = mykrobe_data.get(sid, {})
                pv = pinv_data.get(sid, {})
                row.update({
                    "shigeifinder_ipaH":              sr.get("shigeifinder_ipaH",              "NA"),
                    "shigeifinder_virulence_plasmid":  sr.get("shigeifinder_virulence_plasmid", "NA"),
                    "shigeifinder_cluster":            sr.get("shigeifinder_cluster",           "NA"),
                    "shigeifinder_serotype":           sr.get("shigeifinder_serotype",          "NA"),
                    "shigeifinder_o_antigen":          sr.get("shigeifinder_o_antigen",         "NA"),
                    "shigeifinder_h_antigen":          sr.get("shigeifinder_h_antigen",         "NA"),
                    "mykrobe_genotype":      mk.get("mykrobe_genotype",      "NA"),
                    "mykrobe_lineage":       mk.get("mykrobe_lineage",       "NA"),
                    "mykrobe_clade":         mk.get("mykrobe_clade",         "NA"),
                    "mykrobe_subclade":      mk.get("mykrobe_subclade",      "NA"),
                    "mykrobe_genotype_name": mk.get("mykrobe_genotype_name", "NA"),
                    "mykrobe_confidence":    mk.get("mykrobe_confidence",    "NA"),
                    "pinv_present":          pv.get("pinv_present",          "N"),
                    "pinv_genes":            pv.get("pinv_genes",            "NA"),
                    "is_elements":           is_data.get(sid,                "NA"),
                })

            writer.writerow(row)

    print(f"Written: {args.output} ({len(all_samples)} {args.species} samples)", file=sys.stderr)


if __name__ == "__main__":
    main()
