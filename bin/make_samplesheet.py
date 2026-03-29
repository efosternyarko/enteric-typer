#!/usr/bin/env python3
"""Generate an enteric-typer samplesheet CSV from a folder of FASTA assemblies.

Usage
-----
  make_samplesheet.py --input /path/to/assemblies/ --output samples.csv
  make_samplesheet.py --input /path/to/assemblies/ --output samples.csv --strip _contigs _assembly

The sample ID is derived from the filename by removing the extension and any
optional suffix strings passed via --strip.

Example
-------
  Filename:  ERR1234567_contigs.fasta
  --strip _contigs
  ID:        ERR1234567
"""

import argparse
import csv
import fnmatch
import os
import sys

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna", ".fas", ".fsa"}


def get_sample_id(filename: str, strip_suffixes: list[str]) -> str:
    name = filename
    for ext in FASTA_EXTENSIONS:
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break
    for suffix in strip_suffixes:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate an enteric-typer samplesheet from a FASTA folder",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--input",   "-i", required=True,  help="Folder containing FASTA assembly files")
    parser.add_argument("--output",  "-o", required=True,  help="Output samplesheet CSV path")
    parser.add_argument("--strip",   nargs="+", default=[], metavar="SUFFIX",
                        help="Filename suffixes to remove when deriving sample IDs")
    parser.add_argument("--pattern", default=None, metavar="GLOB",
                        help="Optional glob pattern to filter files (e.g. '*.fasta')")
    args = parser.parse_args()

    folder = os.path.abspath(args.input)
    if not os.path.isdir(folder):
        sys.exit(f"ERROR: Input folder not found: {folder}")

    if args.pattern:
        all_files = [f for f in os.listdir(folder) if fnmatch.fnmatch(f, args.pattern)]
    else:
        all_files = [f for f in os.listdir(folder)
                     if os.path.splitext(f)[1].lower() in FASTA_EXTENSIONS]

    if not all_files:
        sys.exit(f"ERROR: No FASTA files found in {folder}.\n"
                 f"Expected extensions: {', '.join(sorted(FASTA_EXTENSIONS))}")

    all_files.sort()

    rows: list[dict] = []
    seen_ids: dict[str, str] = {}
    for fname in all_files:
        fpath = os.path.join(folder, fname)
        sid   = get_sample_id(fname, args.strip)
        if sid in seen_ids:
            sys.exit(
                f"ERROR: Duplicate sample ID '{sid}' from:\n"
                f"  {seen_ids[sid]}\n  {fpath}\n"
                f"Use --strip to disambiguate."
            )
        seen_ids[sid] = fpath
        rows.append({"id": sid, "fasta": fpath})

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "fasta"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written: {args.output} ({len(rows)} samples)")
    for row in rows[:5]:
        print(f"  {row['id']}  →  {row['fasta']}")
    if len(rows) > 5:
        print(f"  ... and {len(rows) - 5} more")


if __name__ == "__main__":
    main()
