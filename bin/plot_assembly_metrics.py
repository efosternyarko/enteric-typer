#!/usr/bin/env python3
"""
plot_assembly_metrics.py — assembly QC metrics figure per species.

Generates an 8-panel figure (2 rows × 4 columns):
  Row 1: Genome Length histogram | Genome Length boxplot | N50 histogram | N50 boxplot
  Row 2: Contig Count histogram  | Contig Count boxplot  | GC% histogram | GC% boxplot

Style matches panels E–H of the reference assembly+read metrics figure:
histogram with KDE overlay paired with a horizontal boxplot for each metric.

Usage:
    plot_assembly_metrics.py --stats sample1.tsv [sample2.tsv ...] \\
                             --species ecoli \\
                             --output ecoli_assembly_metrics.png \\
                             --summary ecoli_assembly_metrics_summary.tsv
"""

import argparse
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# ── Species display names ──────────────────────────────────────────────────
SPECIES_LABELS = {
    'ecoli':      'Escherichia coli',
    'salmonella': 'Salmonella enterica',
    'shigella':   'Shigella spp.',
}

# ── Colour palette (one colour per metric, matching reference palette) ─────
COLOURS = {
    'length':  '#6a9fca',  # steel blue
    'n50':     '#9b7cbc',  # purple
    'contigs': '#e08c6a',  # warm orange
    'gc':      '#7bbf7b',  # muted green
}

PANEL_LABELS = 'ABCDEFGH'


# ── Helpers ───────────────────────────────────────────────────────────────

def _scale_data(data: np.ndarray) -> tuple:
    """Return (scaled_data, suffix) so axis values are readable."""
    mx = np.max(data) if len(data) > 0 else 0
    if mx >= 1e6:
        return data / 1e6, 'Mb'
    if mx >= 1e3:
        return data / 1e3, 'kb'
    return data, ''


def _n_bins(data: np.ndarray) -> int:
    n = len(data)
    if n <= 5:  return n
    if n <= 20: return max(5, n // 2)
    return min(30, max(10, int(np.ceil(np.sqrt(n)))))


def plot_histogram(ax, data: np.ndarray, colour: str, xlabel: str,
                   title: str, panel_label: str, scale: bool = True) -> None:
    """Histogram + KDE overlay, matching the reference figure style (panels E, G)."""
    if len(data) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                transform=ax.transAxes, color='#999999')
        ax.set_title(f'{panel_label}   {title}', loc='left', fontweight='bold')
        return

    scaled, suffix = _scale_data(data) if scale else (data, '')
    bins = _n_bins(scaled)

    ax.hist(scaled, bins=bins, color=colour, alpha=0.70,
            edgecolor='white', linewidth=0.5)

    # KDE overlay (needs ≥ 2 distinct values)
    if len(scaled) >= 2 and np.std(scaled) > 0:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            kde = gaussian_kde(scaled)
        xs = np.linspace(scaled.min(), scaled.max(), 300)
        bin_width = (scaled.max() - scaled.min()) / bins
        ax.plot(xs, kde(xs) * len(scaled) * bin_width,
                color=colour, linewidth=2.0)

    x_label = f'{xlabel} ({suffix})' if suffix else xlabel
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel('Frequency', fontsize=9)
    ax.set_title(f'{panel_label}   {title}', loc='left', fontweight='bold', fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.7)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=8)


def plot_boxplot(ax, data: np.ndarray, colour: str, xlabel: str,
                 title: str, panel_label: str, scale: bool = True) -> None:
    """Horizontal boxplot matching the reference figure style (panels F, H)."""
    if len(data) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                transform=ax.transAxes, color='#999999')
        ax.set_title(f'{panel_label}   {title}', loc='left', fontweight='bold')
        return

    scaled, suffix = _scale_data(data) if scale else (data, '')

    ax.boxplot(
        scaled,
        vert=False,
        patch_artist=True,
        widths=0.5,
        flierprops=dict(marker='D', markerfacecolor='#888888',
                        markeredgecolor='none', markersize=5, alpha=0.7),
        medianprops=dict(color='#333333', linewidth=2),
        boxprops=dict(facecolor=colour, alpha=0.70),
        whiskerprops=dict(color='#555555', linewidth=1.2),
        capprops=dict(color='#555555', linewidth=1.2),
    )

    x_label = f'{xlabel} ({suffix})' if suffix else xlabel
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_yticks([])
    ax.set_title(f'{panel_label}   {title}', loc='left', fontweight='bold', fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.7, axis='x')
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.tick_params(labelsize=8)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stats',   nargs='+', required=True,
                        help='Per-sample *_assembly_stats.tsv files')
    parser.add_argument('--species', required=True,
                        help='Species ID: ecoli | salmonella | shigella')
    parser.add_argument('--output',  required=True,
                        help='Output PNG file path')
    parser.add_argument('--summary', required=True,
                        help='Output merged summary TSV path')
    args = parser.parse_args()

    # ── Load and merge per-sample stats ──────────────────────────────────────
    frames = []
    for f in args.stats:
        try:
            df = pd.read_csv(f, sep='\t')
            frames.append(df)
        except Exception as e:
            print(f"WARNING: could not read {f}: {e}", file=sys.stderr)

    if not frames:
        sys.exit("ERROR: no readable stats files provided")

    data = pd.concat(frames, ignore_index=True)
    data.to_csv(args.summary, sep='\t', index=False)

    n = len(data)
    species_label = SPECIES_LABELS.get(args.species, args.species)

    # ── Extract metric arrays ─────────────────────────────────────────────────
    genome_len  = data['genome_length'].to_numpy(dtype=float)
    num_contigs = data['num_contigs'].to_numpy(dtype=float)
    asm_n50     = data['assembly_N50'].to_numpy(dtype=float)
    gc_pct      = data['gc_pct'].to_numpy(dtype=float)

    # ── Build figure ──────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    fig.suptitle(f'{species_label} — Assembly Metrics  (n = {n} samples)',
                 fontsize=14, fontweight='bold', y=0.98)

    # Row 1 — Genome Length + N50
    plot_histogram(axes[0, 0], genome_len, COLOURS['length'],
                   'Genome Length', 'Genome Length Distribution', 'A')
    plot_boxplot  (axes[0, 1], genome_len, COLOURS['length'],
                   'Genome Length', 'Genome Length Boxplot',      'B')
    plot_histogram(axes[0, 2], asm_n50,   COLOURS['n50'],
                   'Assembly N50',  'N50 Distribution',           'C')
    plot_boxplot  (axes[0, 3], asm_n50,   COLOURS['n50'],
                   'Assembly N50',  'N50 Boxplot',                'D')

    # Row 2 — Contig Count + GC%
    plot_histogram(axes[1, 0], num_contigs, COLOURS['contigs'],
                   'Number of Contigs', 'Contig Count Distribution', 'E', scale=False)
    plot_boxplot  (axes[1, 1], num_contigs, COLOURS['contigs'],
                   'Number of Contigs', 'Contig Count Boxplot',      'F', scale=False)
    plot_histogram(axes[1, 2], gc_pct,      COLOURS['gc'],
                   'GC%',               'GC% Distribution',          'G', scale=False)
    plot_boxplot  (axes[1, 3], gc_pct,      COLOURS['gc'],
                   'GC%',               'GC% Boxplot',               'H', scale=False)

    plt.tight_layout(rect=[0, 0, 1, 0.97])

    # Save PNG and PDF
    pdf_path = args.output.replace('.png', '.pdf')
    fig.savefig(args.output,  dpi=150, bbox_inches='tight')
    fig.savefig(pdf_path,             bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {args.output} + {pdf_path}  ({species_label}, n={n})")


if __name__ == '__main__':
    main()
