#!/usr/bin/env python3
"""
plot_plasmid_amr_map.py — Visualise the plasmid–AMR gene co-occurrence map.

Reads the aggregate plasmid_amr_map.tsv produced by AGGREGATE_PLASMID_AMR_MAP
and generates a bubble-matrix plot:
  - Rows   = top-N most prevalent replicon types
  - Columns = AMR drug classes present in the dataset
  - Bubble  = % of unique isolates where that replicon co-occurs with that
              drug class on the same contig
  - Colour  = drug class (consistent with plot_summary.py fig4 palette)

Also generates a companion stacked-bar figure showing the fraction of isolates
carrying each replicon type broken down by dominant AMR burden
(no AMR | single class | 2 classes | ≥3 classes).

Output (PDF + PNG):
  {prefix}_plasmid_amr_map_bubble.pdf/png
  {prefix}_plasmid_amr_map_bars.pdf/png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Global aesthetics ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica Neue", "DejaVu Sans"],
    "font.size":         8,
    "axes.labelsize":    9,
    "axes.titlesize":    10,
    "axes.titleweight":  "bold",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "xtick.labelsize":   7.5,
    "ytick.labelsize":   8,
    "legend.fontsize":   7.5,
    "legend.frameon":    False,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "savefig.pad_inches": 0.08,
})

# Drug class display labels (same order as plot_summary.py)
CLASS_PRIORITY = [
    "BETA-LACTAM", "QUINOLONE", "COLISTIN", "AMINOGLYCOSIDE",
    "TETRACYCLINE", "MACROLIDE", "PHENICOL", "SULFONAMIDE",
    "TRIMETHOPRIM", "FOSFOMYCIN", "STREPTOTHRICIN", "NITROFURAN",
]
CLASS_LABEL = {
    "BETA-LACTAM":     "Beta-lactam",
    "QUINOLONE":       "Quinolone",
    "COLISTIN":        "Colistin",
    "AMINOGLYCOSIDE":  "Aminoglycoside",
    "TETRACYCLINE":    "Tetracycline",
    "MACROLIDE":       "Macrolide",
    "PHENICOL":        "Phenicol",
    "SULFONAMIDE":     "Sulfonamide",
    "TRIMETHOPRIM":    "Trimethoprim",
    "FOSFOMYCIN":      "Fosfomycin",
    "STREPTOTHRICIN":  "Streptothricin",
    "NITROFURAN":      "Nitrofuran",
}
CLASS_COLOR = {
    "BETA-LACTAM":     "#d62728",
    "QUINOLONE":       "#ff7f0e",
    "COLISTIN":        "#7f7f7f",
    "AMINOGLYCOSIDE":  "#1f77b4",
    "TETRACYCLINE":    "#9467bd",
    "MACROLIDE":       "#e377c2",
    "PHENICOL":        "#8c564b",
    "SULFONAMIDE":     "#2ca02c",
    "TRIMETHOPRIM":    "#bcbd22",
    "FOSFOMYCIN":      "#17becf",
    "STREPTOTHRICIN":  "#aec7e8",
    "NITROFURAN":      "#ffbb78",
}

_SENTINEL = {"NA", "nan", "None", "-", ""}

# Replicon values that are not real replicon types — filter from all visualisations
_REPLICON_SENTINELS = {
    "no_replicon", "No replicons found", "No replicon",
    "no_replicons", "No replicons", "No plasmids found",
}


def _save(fig: plt.Figure, outdir: Path, stem: str) -> None:
    for ext in ("pdf", "png"):
        dpi = 300 if ext == "pdf" else 150
        p   = outdir / f"{stem}.{ext}"
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        print(f"  Saved: {p}", file=sys.stderr)
    plt.close(fig)


def load_map(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    required = {"sample_id", "replicon", "drug_classes"}
    missing  = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: plasmid_amr_map.tsv missing columns: {missing}")
    return df


def _expand_classes(val) -> list[str]:
    """Split a ';'-delimited drug-class string into a list."""
    if pd.isna(val):
        return []
    s = str(val).strip()
    if not s or s in _SENTINEL:
        return []
    return [c.strip() for c in s.split(";") if c.strip() and c.strip() not in _SENTINEL]


# ── Bubble-matrix plot ────────────────────────────────────────────────────────

def plot_bubble_matrix(df: pd.DataFrame, outdir: Path, prefix: str, top_n: int = 20) -> None:
    """
    Bubble matrix: rows = replicon types, columns = drug classes.
    Bubble area ∝ % of isolates with that replicon carrying that drug class.
    """
    n_isolates = df["sample_id"].nunique()

    # Count how many isolates have each replicon (to rank top-N).
    # Exclude sentinel values that are not real replicon types.
    rep_isolates: dict[str, set] = {}
    for _, row in df.iterrows():
        rep = str(row["replicon"]).strip()
        sid = str(row["sample_id"]).strip()
        if rep and rep not in _SENTINEL and rep not in _REPLICON_SENTINELS:
            rep_isolates.setdefault(rep, set()).add(sid)

    if not rep_isolates:
        print("WARNING: no replicon data — skipping bubble matrix", file=sys.stderr)
        return

    # Top-N replicons by number of isolates carrying them
    top_reps = sorted(rep_isolates, key=lambda r: len(rep_isolates[r]), reverse=True)[:top_n]

    # For each replicon × drug class, count isolates
    # matrix[rep][cls] = set of sample_ids
    matrix: dict[str, dict[str, set]] = {r: {} for r in top_reps}
    for _, row in df.iterrows():
        rep = str(row["replicon"]).strip()
        sid = str(row["sample_id"]).strip()
        if rep not in matrix:
            continue
        for cls in _expand_classes(row.get("drug_classes")):
            if cls in CLASS_COLOR:
                matrix[rep].setdefault(cls, set()).add(sid)

    # Classes present in any top-N replicon, ordered by CLASS_PRIORITY
    all_classes = set()
    for rep in top_reps:
        all_classes |= set(matrix[rep].keys())
    cols = [c for c in CLASS_PRIORITY if c in all_classes]

    if not cols:
        print("WARNING: no drug class data in bubble matrix — skipping", file=sys.stderr)
        return

    n_rows = len(top_reps)
    n_cols = len(cols)
    fig_w  = max(6.5, n_cols * 0.75 + 2.5)
    fig_h  = max(4.0, n_rows * 0.45 + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Draw bubbles
    max_r = 0.42  # max bubble radius in data units
    for i, rep in enumerate(top_reps):
        for j, cls in enumerate(cols):
            count = len(matrix[rep].get(cls, set()))
            pct   = 100 * count / n_isolates
            if pct > 0:
                # Scale bubble area proportional to pct; cap at max_r
                r = max_r * np.sqrt(pct / 100)
                ax.scatter(j, i,
                           s=(r * 72) ** 2,   # convert axes-units to points²
                           color=CLASS_COLOR[cls],
                           alpha=0.80, linewidths=0.4,
                           edgecolors="white", zorder=3)
                if pct >= 5:
                    ax.text(j, i, f"{pct:.0f}%",
                            ha="center", va="center",
                            fontsize=5.5, color="white", fontweight="bold", zorder=4)

    # Axes formatting
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([CLASS_LABEL.get(c, c) for c in cols],
                       rotation=40, ha="right", fontsize=7.5)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(
        [f"{r}  ({len(rep_isolates[r])/n_isolates*100:.0f}%)" for r in top_reps],
        fontsize=8
    )
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(-0.6, n_rows - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("AMR drug class")
    ax.set_ylabel("Replicon type")
    ax.set_title(
        f"Plasmid-AMR co-occurrence matrix  (n = {n_isolates} isolates)\n"
        "Bubble area proportional to % of isolates; % shown when >= 5%",
        fontweight="bold"
    )

    # Add light grid
    for j in range(n_cols):
        ax.axvline(j, color="#eeeeee", lw=0.5, zorder=1)
    for i in range(n_rows):
        ax.axhline(i, color="#eeeeee", lw=0.5, zorder=1)

    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)

    fig.tight_layout()
    _save(fig, outdir, f"{prefix}_plasmid_amr_map_bubble")


# ── Stacked-bar: replicon × AMR burden ───────────────────────────────────────

def plot_replicon_amr_bars(df: pd.DataFrame, outdir: Path, prefix: str, top_n: int = 20) -> None:
    """
    Horizontal stacked bar chart:
    Rows = top-N replicon types.
    Segments = % of isolates carrying that replicon, split by number of
               distinct drug classes on the same contig (0 / 1 / 2 / ≥3).
    """
    n_isolates = df["sample_id"].nunique()

    # For each (sample, replicon), count distinct drug classes
    # Use the maximum class count seen for a given sample–replicon pair
    # (a replicon may map to multiple contigs)
    pair_classes: dict[tuple, set] = {}
    for _, row in df.iterrows():
        rep = str(row["replicon"]).strip()
        sid = str(row["sample_id"]).strip()
        if not rep or rep in _SENTINEL:
            continue
        key = (sid, rep)
        for cls in _expand_classes(row.get("drug_classes")):
            if cls in CLASS_COLOR:
                pair_classes.setdefault(key, set()).add(cls)

    if not pair_classes:
        print("WARNING: no data for replicon AMR bars — skipping", file=sys.stderr)
        return

    # Rank replicons
    from collections import Counter
    rep_count: Counter = Counter()
    for (sid, rep) in pair_classes:
        rep_count[rep] += 1
    # Also count reps with 0 AMR (present in df but no classes)
    for _, row in df.iterrows():
        rep = str(row["replicon"]).strip()
        sid = str(row["sample_id"]).strip()
        if rep and rep not in _SENTINEL:
            key = (sid, rep)
            if key not in pair_classes:
                rep_count[rep] += 1

    top_reps = [r for r, _ in rep_count.most_common(top_n)]

    # Build counts per burden category
    burden_cats = ["No AMR", "1 class", "2 classes", "≥3 classes"]
    burden_colors = ["#d3d3d3", "#aec7e8", "#ff9896", "#d62728"]
    counts = {cat: [] for cat in burden_cats}

    for rep in top_reps:
        n0 = n1 = n2 = n3 = 0
        # Collect unique isolates for this rep
        isolates_seen: set = set()
        for (sid, r), cls_set in pair_classes.items():
            if r == rep:
                isolates_seen.add(sid)
                n = len(cls_set)
                if n == 1:   n1 += 1
                elif n == 2: n2 += 1
                else:        n3 += 1
        # Count no-AMR isolates (in df but not in pair_classes for this rep)
        for _, row in df.iterrows():
            r2  = str(row["replicon"]).strip()
            sid = str(row["sample_id"]).strip()
            if r2 == rep and sid not in isolates_seen:
                n0 += 1
                isolates_seen.add(sid)
        total = n0 + n1 + n2 + n3
        if total == 0:
            total = 1  # avoid div/0
        counts["No AMR"].append(100 * n0 / n_isolates)
        counts["1 class"].append(100 * n1 / n_isolates)
        counts["2 classes"].append(100 * n2 / n_isolates)
        counts["≥3 classes"].append(100 * n3 / n_isolates)

    fig_h = max(4.0, len(top_reps) * 0.42 + 1.6)
    fig, ax = plt.subplots(figsize=(7.0, fig_h))
    y = np.arange(len(top_reps))
    left = np.zeros(len(top_reps))
    for cat, col in zip(burden_cats, burden_colors):
        vals = np.array(counts[cat])
        ax.barh(y, vals, left=left, color=col, label=cat,
                edgecolor="white", linewidth=0.3, height=0.72)
        left += vals

    ax.set_yticks(y)
    ax.set_yticklabels(top_reps, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (% of isolates)")
    ax.set_title("Plasmid replicons — AMR burden per replicon type\n"
                 "(PlasmidFinder + AMRFinder)", fontweight="bold")
    ax.set_xlim(0, 108)
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    ax.legend(title="Drug classes on same contig",
              title_fontsize=7.5, fontsize=7.5, loc="lower right")
    fig.tight_layout()
    _save(fig, outdir, f"{prefix}_plasmid_amr_map_bars")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(description="Visualise plasmid–AMR co-occurrence map")
    p.add_argument("--input",  "-i", required=True,
                   help="Aggregate plasmid_amr_map.tsv")
    p.add_argument("--outdir", "-o", default=".",
                   help="Output directory (default: current dir)")
    p.add_argument("--prefix", "-p", default="enteric_typer",
                   help="Output file prefix")
    p.add_argument("--top_n",  "-n", type=int, default=20,
                   help="Top N replicon types to display (default: 20)")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.input}…", file=sys.stderr)
    df = load_map(args.input)
    print(f"  {df['sample_id'].nunique()} isolates  ·  "
          f"{df['replicon'].nunique()} distinct replicons", file=sys.stderr)

    print("Generating plasmid AMR map figures…", file=sys.stderr)
    plot_bubble_matrix(df, outdir, args.prefix, top_n=args.top_n)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
