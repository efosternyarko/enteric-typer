#!/usr/bin/env python3
"""
plot_tree_annotation.py — Phylogenetic tree annotated with AMR resistance profile.

Reads:
  --tree      Newick tree file (IQ-TREE .treefile)
  --metadata  enteric-typer results TSV (ecoli_typer_results.tsv)
  --outdir    output directory
  --prefix    file name prefix
  --species   ecoli | salmonella

Produces (PDF + PNG at 300 dpi):
  {prefix}_{species}_tree_amr.pdf / .png

Layout (left → right):
  Phylogenetic tree  |  Phylogroup strip  |  AMR drug class heatmap
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from Bio import Phylo
    HAS_PHYLO = True
except ImportError:
    HAS_PHYLO = False

# ── Global aesthetics (keep in sync with plot_summary.py) ────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica Neue", "DejaVu Sans"],
    "font.size":         8,
    "axes.labelsize":    9,
    "axes.titlesize":    10,
    "axes.titleweight":  "bold",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "xtick.labelsize":   8,
    "ytick.labelsize":   8,
    "legend.fontsize":   7.5,
    "legend.frameon":    False,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "savefig.pad_inches": 0.08,
})

# Phylogroup colours (Clermont scheme)
PHYLOGROUP_COLORS = {
    "A":  "#4e79a7", "B1": "#59a14f", "B2": "#f28e2b",
    "C":  "#76b7b2", "D":  "#e15759", "E":  "#b07aa1",
    "F":  "#ff9da7", "Unknown": "#bab0ac",
}

ST_PHYLOGROUP: dict[str, str] = {
    "ST10": "A",  "ST12": "B2", "ST14": "D",  "ST23": "A",
    "ST34": "A",  "ST38": "D",  "ST52": "A",  "ST58": "B1",
    "ST69": "D",  "ST73": "B2", "ST88": "D",  "ST95": "B2",
    "ST117":"D",  "ST127":"B2", "ST131":"B2", "ST141":"B1",
    "ST144":"A",  "ST155":"B1", "ST167":"A",  "ST218":"A",
    "ST354":"B1", "ST393":"A",  "ST405":"D",  "ST410":"C",
    "ST443":"A",  "ST448":"A",  "ST453":"A",  "ST617":"B1",
    "ST636":"A",  "ST648":"F",  "ST1193":"B2",
}

CLINICAL_CLASSES = [
    "BETA-LACTAM", "AMINOGLYCOSIDE", "SULFONAMIDE", "TRIMETHOPRIM",
    "TETRACYCLINE", "FOSFOMYCIN", "COLISTIN", "QUINOLONE", "PHENICOL",
    "MACROLIDE", "STREPTOTHRICIN", "FOSMIDOMYCIN", "NITROFURAN",
]

CLASS_LABEL = {
    "AMINOGLYCOSIDE": "Aminoglycoside", "BETA-LACTAM":    "Beta-lactam",
    "COLISTIN":       "Colistin",       "FOSFOMYCIN":     "Fosfomycin",
    "FOSMIDOMYCIN":   "Fosmidomycin",   "MACROLIDE":      "Macrolide",
    "NITROFURAN":     "Nitrofuran",     "PHENICOL":       "Phenicol",
    "QUINOLONE":      "Quinolone",      "STREPTOTHRICIN": "Streptothricin",
    "SULFONAMIDE":    "Sulfonamide",    "TETRACYCLINE":   "Tetracycline",
    "TRIMETHOPRIM":   "Trimethoprim",
}


# ── Tree layout ───────────────────────────────────────────────────────────────

def _layout(root) -> tuple[dict, list[str]]:
    """
    Compute {id(clade): (x, y)} and ordered tip name list (top → bottom).
    x = cumulative branch length from root.
    y = tip index (0, 1, 2, ...) for terminals; mean of children for internals.
    """
    tips: list[str] = []
    pos:  dict[int, tuple[float, float]] = {}

    def _recurse(clade, x: float) -> float:
        x_here = x + (clade.branch_length or 0.0)
        if clade.is_terminal():
            y = float(len(tips))
            tips.append(clade.name or "")
        else:
            ys = [_recurse(c, x_here) for c in clade.clades]
            y  = (min(ys) + max(ys)) / 2.0
        pos[id(clade)] = (x_here, y)
        return y

    _recurse(root, 0.0)
    return pos, tips


def _draw_tree(ax: plt.Axes, root, pos: dict) -> None:
    """Recursively draw horizontal + vertical tree branches."""
    def _draw(clade):
        x, y = pos[id(clade)]
        if not clade.is_terminal():
            cys = [pos[id(c)][1] for c in clade.clades]
            ax.plot([x, x], [min(cys), max(cys)],
                    color="#333333", lw=0.75, solid_capstyle="butt")
            for c in clade.clades:
                cx, cy = pos[id(c)]
                ax.plot([x, cx], [cy, cy], color="#333333", lw=0.75)
                _draw(c)

    _draw(root)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _parse_classes(val) -> set[str]:
    if pd.isna(val) or not str(val).strip():
        return set()
    out: set[str] = set()
    for tok in str(val).split(","):
        for part in tok.strip().split("/"):
            p = part.strip()
            if p:
                out.add(p)
    return out


def _clean_st(raw) -> str:
    s = str(raw).strip() if pd.notna(raw) else ""
    if s in ("", "-", "NA", "nan", "None", "No ST predicted"):
        return "Unknown"
    return s if s.startswith("ST") else f"ST{s}"


def _phylogroup_color(st: str) -> str:
    pg = ST_PHYLOGROUP.get(st, "Unknown")
    return PHYLOGROUP_COLORS.get(pg, PHYLOGROUP_COLORS["Unknown"])


def _scale_bar_value(max_x: float) -> float:
    """Round down to a clean value for scale bar."""
    target = max_x / 5.0
    if target <= 0:
        return 0.001
    mag = 10 ** math.floor(math.log10(target))
    return round(target / mag) * mag


# ── Main figure ───────────────────────────────────────────────────────────────

def plot_tree_amr(
    treefile:     Path,
    metadata_tsv: Path,
    outdir:       Path,
    prefix:       str,
    species:      str,
) -> None:
    if not HAS_PHYLO:
        print("ERROR: biopython is not installed. Cannot draw annotated tree.", file=sys.stderr)
        sys.exit(1)

    # ── Load tree ─────────────────────────────────────────────────────────────
    tree = Phylo.read(str(treefile), "newick")
    tree.root_at_midpoint()
    pos, tips = _layout(tree.root)
    n = len(tips)

    if n < 2:
        print("INFO: Fewer than 2 tips in tree — skipping annotation figure.", file=sys.stderr)
        return

    # ── Load metadata ─────────────────────────────────────────────────────────
    df   = pd.read_csv(str(metadata_tsv), sep="\t", low_memory=False)
    df["sample"] = df["sample"].astype(str).str.strip()
    meta = df.set_index("sample")

    matched = set(tips) & set(meta.index)
    if len(matched) == 0:
        print("WARNING: No tree tip labels match metadata sample IDs. "
              "Check that sample names are consistent.", file=sys.stderr)
        return
    if len(matched) < n:
        print(f"INFO: {n - len(matched)} tree tip(s) have no metadata — shown as Unknown/absent.",
              file=sys.stderr)

    # ── Determine AMR column ──────────────────────────────────────────────────
    amr_col = next(
        (c for c in ("amrfinder_drug_classes", "amr_classes", "amrfinderplus_amr_classes")
         if c in meta.columns),
        None
    )

    # Order classes by descending prevalence; skip absent classes
    if amr_col:
        class_counts = {
            cls: sum(1 for v in meta[amr_col] if cls in _parse_classes(v))
            for cls in CLINICAL_CLASSES
        }
        classes_ordered = [c for c in CLINICAL_CLASSES
                           if class_counts.get(c, 0) > 0]
        classes_ordered.sort(key=lambda c: -class_counts[c])
    else:
        class_counts   = {}
        classes_ordered = []

    n_cls = len(classes_ordered)

    # ── Build binary AMR matrix (rows = tips in tree order) ───────────────────
    mat = np.zeros((n, max(n_cls, 1)), dtype=np.uint8)
    if amr_col and n_cls:
        for i, tip in enumerate(tips):
            if tip in meta.index:
                row_cls = _parse_classes(meta.loc[tip, amr_col])
                for j, cls in enumerate(classes_ordered):
                    mat[i, j] = 1 if cls in row_cls else 0

    # ── Figure geometry ───────────────────────────────────────────────────────
    row_h   = max(0.05, min(0.18, 10.0 / n))
    fig_h   = max(5.0, n * row_h + 2.8)
    tree_w  = 4.5
    strip_w = 0.35
    heat_w  = max(1.5, n_cls * 0.50) if n_cls else 0

    col_ratios = [tree_w, strip_w]
    n_axes = 2
    if n_cls:
        col_ratios.append(heat_w)
        n_axes = 3

    fig, axes = plt.subplots(
        1, n_axes,
        figsize=(sum(col_ratios) + 1.2, fig_h),
        gridspec_kw={"width_ratios": col_ratios, "wspace": 0.025},
    )
    ax_tree = axes[0]
    ax_pg   = axes[1]
    ax_heat = axes[2] if n_cls else None

    # ── Draw tree ─────────────────────────────────────────────────────────────
    _draw_tree(ax_tree, tree.root, pos)

    max_x = max(x for x, _ in pos.values())
    ax_tree.set_ylim(n - 0.5, -0.5)   # inverted: tip 0 at top
    ax_tree.set_yticks([])
    ax_tree.set_xlabel("Substitutions / site", fontsize=8)
    ax_tree.spines["left"].set_visible(False)
    ax_tree.spines["top"].set_visible(False)
    ax_tree.spines["right"].set_visible(False)

    # Tip labels (show when ≤ 60 isolates)
    if n <= 60:
        x_label = max_x * 1.04
        for i, tip in enumerate(tips):
            ax_tree.text(x_label, i, tip, fontsize=5.5,
                         va="center", ha="left", fontfamily="monospace")
        ax_tree.set_xlim(-max_x * 0.02, max_x * 1.5)
    else:
        ax_tree.set_xlim(-max_x * 0.02, max_x * 1.05)

    # Scale bar
    sb = _scale_bar_value(max_x)
    y_sb = n - 0.3
    ax_tree.plot([0, sb], [y_sb, y_sb], color="#333333", lw=1.2, clip_on=False)
    ax_tree.text(sb / 2, y_sb + 0.55, f"{sb:.3g}",
                 ha="center", va="top", fontsize=6.5)

    # ── Phylogroup strip ──────────────────────────────────────────────────────
    for i, tip in enumerate(tips):
        st  = _clean_st(meta.loc[tip, "mlst_st"]) if tip in meta.index and "mlst_st" in meta.columns else "Unknown"
        col = _phylogroup_color(st)
        ax_pg.barh(i, 1, color=col, height=1.0, linewidth=0, align="center")

    ax_pg.set_xlim(0, 1)
    ax_pg.set_ylim(n - 0.5, -0.5)
    ax_pg.set_yticks([])
    ax_pg.set_xticks([])
    ax_pg.set_title("PG", fontsize=7.5, fontweight="bold", pad=3)
    for sp in ax_pg.spines.values():
        sp.set_visible(False)

    # ── AMR heatmap ───────────────────────────────────────────────────────────
    if ax_heat is not None and n_cls:
        cmap = mcolors.ListedColormap(["#f5f5f5", "#c0392b"])
        ax_heat.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                       interpolation="nearest", origin="upper")

        col_labels = [
            f"{CLASS_LABEL.get(c, c)}\n({100 * class_counts[c] / n:.0f}%)"
            for c in classes_ordered
        ]
        ax_heat.set_xticks(np.arange(n_cls))
        ax_heat.set_xticklabels(col_labels, rotation=40, ha="right", fontsize=7.5)
        ax_heat.set_yticks([])
        ax_heat.set_title("AMR drug class", fontsize=8, fontweight="bold", pad=3)
        for sp in ax_heat.spines.values():
            sp.set_visible(False)

        # Colour legend for absent / present
        ax_heat.legend(
            handles=[
                mpatches.Patch(facecolor="#c0392b", label="Resistant"),
                mpatches.Patch(facecolor="#f5f5f5", edgecolor="#cccccc",
                               linewidth=0.5, label="Susceptible / absent"),
            ],
            fontsize=7, loc="upper right",
            bbox_to_anchor=(1.0, -0.18),
            frameon=True, framealpha=0.9, edgecolor="none",
        )

    # ── Phylogroup legend on tree panel ───────────────────────────────────────
    seen_pgs: set[str] = set()
    for tip in tips:
        if tip in meta.index and "mlst_st" in meta.columns:
            st = _clean_st(meta.loc[tip, "mlst_st"])
        else:
            st = "Unknown"
        seen_pgs.add(ST_PHYLOGROUP.get(st, "Unknown"))

    pg_patches = [
        mpatches.Patch(facecolor=PHYLOGROUP_COLORS.get(pg, "#bab0ac"), label=f"Phylogroup {pg}")
        for pg in sorted(seen_pgs)
    ]
    ax_tree.legend(
        handles=pg_patches, fontsize=6.5, title_fontsize=7,
        loc="lower right", handlelength=1,
        frameon=True, framealpha=0.85, edgecolor="none",
    )

    # ── Title ─────────────────────────────────────────────────────────────────
    sp_label = "E. coli" if species == "ecoli" else "Salmonella enterica"
    fig.suptitle(
        f"{sp_label} core-SNP phylogeny with AMR profile  (n = {n} isolates)",
        fontsize=11, fontweight="bold", y=1.01,
    )

    # ── Save ──────────────────────────────────────────────────────────────────
    for ext in ("pdf", "png"):
        out = outdir / f"{prefix}_{species}_tree_amr.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"  Saved: {out}", file=sys.stderr)
    plt.close(fig)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(
        description="Annotate IQ-TREE phylogeny with AMR drug class heatmap"
    )
    p.add_argument("--tree",     "-t", required=True,  help="Newick tree file")
    p.add_argument("--metadata", "-m", required=True,  help="enteric-typer results TSV")
    p.add_argument("--outdir",   "-o", default=".",    help="Output directory")
    p.add_argument("--prefix",   "-p", default="enteric_typer")
    p.add_argument("--species",  "-s", required=True,  choices=["ecoli", "salmonella"])
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    plot_tree_amr(
        treefile     = Path(args.tree),
        metadata_tsv = Path(args.metadata),
        outdir       = outdir,
        prefix       = args.prefix,
        species      = args.species,
    )


if __name__ == "__main__":
    main()
