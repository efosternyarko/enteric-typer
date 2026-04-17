#!/usr/bin/env python3
"""
plot_snp_heatmap.py — standalone clustered heatmap from a SKA2 pairwise SNP matrix.

SKA2 `ska distance` produces a tab-separated file where each row is:
    Sample1  Sample2  SNP_distance  [other fields...]

This script pivots that into a symmetric matrix, clusters it hierarchically,
and draws a publication-quality heatmap.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_matrix(path: str) -> pd.DataFrame:
    """
    Load a SKA2 distance file (long format: seq1 seq2 dist …) into a
    symmetric square DataFrame.  Returns an empty DataFrame on failure.
    """
    try:
        # ska distance outputs a header: Sample1  Sample2  Distance  ...
        df = pd.read_csv(path, sep="\t", header=0, comment="#")
    except Exception as e:
        print(f"ERROR: cannot read matrix file '{path}': {e}", file=sys.stderr)
        return pd.DataFrame()

    if df.shape[1] < 3:
        print(f"WARNING: expected ≥3 columns, got {df.shape[1]}", file=sys.stderr)
        return pd.DataFrame()

    # Normalise column names regardless of ska version header wording
    cols = list(df.columns)
    df.columns = ["s1", "s2", "dist"] + [f"_c{i}" for i in range(len(cols) - 3)]

    # Collect all sample names
    samples = sorted(set(df["s1"].tolist() + df["s2"].tolist()))
    n = len(samples)
    idx = {s: i for i, s in enumerate(samples)}

    mat = np.zeros((n, n), dtype=float)
    for _, row in df.iterrows():
        i, j = idx[row["s1"]], idx[row["s2"]]
        mat[i, j] = row["dist"]
        mat[j, i] = row["dist"]

    return pd.DataFrame(mat, index=samples, columns=samples)


def cluster_matrix(mat: pd.DataFrame) -> pd.DataFrame:
    """Return a re-ordered DataFrame using average-linkage hierarchical clustering."""
    if mat.shape[0] < 2:
        return mat
    condensed = squareform(mat.values, checks=False)
    Z = linkage(condensed, method="average")
    order = leaves_list(Z)
    labels = mat.index.tolist()
    ordered = [labels[i] for i in order]
    return mat.loc[ordered, ordered]


# ── Main plot ─────────────────────────────────────────────────────────────────

LABEL_MAX = 200   # show individual sample labels up to this many samples


def plot_snp_heatmap(matrix_path: str, species_label: str, outdir: str) -> None:
    mat = load_matrix(matrix_path)
    if mat.empty:
        print("WARNING: Empty SNP matrix — no heatmap produced.", file=sys.stderr)
        return

    n = mat.shape[0]
    mat = cluster_matrix(mat)

    show_labels = n <= LABEL_MAX

    # Font size: shrink gracefully as n grows; floor at 4 pt
    label_fs = max(4, min(8, int(700 / n))) if show_labels else 0

    # Figure sizing
    # Cell size in inches: shrink with n but keep minimum legibility
    cell_in = max(0.055, min(0.20, 8.0 / n))
    heat_sz = n * cell_in                          # heatmap square side

    if show_labels:
        # Extra margin for rotated x-tick labels (below) and y-tick labels (left)
        # Approximate label length from a typical sample ID (≤ 20 chars at label_fs pt)
        label_margin = max(1.0, label_fs * 20 * 0.6 / 72)   # chars × pts-per-char / pts-per-inch
        fig_w = heat_sz + label_margin + 1.2        # +1.2 for colorbar
        fig_h = heat_sz + label_margin
    else:
        fig_w = heat_sz + 1.5
        fig_h = heat_sz

    # Cap figure dimensions to avoid gigantic files
    fig_w = min(fig_w, 30.0)
    fig_h = min(fig_h, 30.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Choose colourmap: white→dark-blue for SNP counts
    cmap = plt.get_cmap("Blues")
    vmax = np.percentile(mat.values[mat.values > 0], 95) if (mat.values > 0).any() else 1

    im = ax.imshow(mat.values, aspect="auto", cmap=cmap, vmin=0, vmax=vmax,
                   interpolation="nearest")

    # Axis tick labels
    if show_labels:
        ax.set_xticks(range(n))
        ax.set_xticklabels(mat.columns, rotation=90, fontsize=label_fs, ha="center")
        ax.set_yticks(range(n))
        ax.set_yticklabels(mat.index, fontsize=label_fs)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(f"{n} samples (labels hidden at this scale — reduce to ≤{LABEL_MAX} to show)",
                      fontsize=9)

    sp_display = "E. coli" if species_label == "ecoli" else "Salmonella enterica"
    ax.set_title(f"Pairwise whole-genome SNP distances — {sp_display}", fontsize=11, pad=8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Pairwise SNP distance", fontsize=9)

    fig.tight_layout()

    stem = f"{species_label}_snp_heatmap"
    for ext in ("pdf", "png"):
        out = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")

    plt.close(fig)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="SNP distance heatmap from SKA2 output")
    p.add_argument("--matrix",  required=True, help="SKA2 distance TSV (long format)")
    p.add_argument("--species", required=True, help="'ecoli' or 'salmonella'")
    p.add_argument("--outdir",  default=".",   help="Output directory")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    plot_snp_heatmap(args.matrix, args.species, args.outdir)


if __name__ == "__main__":
    main()
