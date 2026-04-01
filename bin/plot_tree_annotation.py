#!/usr/bin/env python3
"""
plot_tree_annotation.py — Phylogenetic tree annotated with virulence & AMR profiles.

Reads:
  --tree      Newick tree file (IQ-TREE .treefile)
  --metadata  enteric-typer results TSV (ecoli_typer_results.tsv)
  --outdir    output directory
  --prefix    file name prefix
  --species   ecoli | salmonella

Produces (PDF + PNG at 300 dpi):
  {prefix}_tree_amr.pdf / .png

Layout (left → right):
  Phylogenetic tree  |  ST strip  |  PG strip (E. coli only)  |  Virulence heatmap  |  AMR genes heatmap (grouped by class)
"""

from __future__ import annotations
from typing import Optional

import argparse
import math
import sys
from pathlib import Path

from collections import Counter, defaultdict

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

_SENTINEL = {"", "-", "NA", "nan", "None", "none"}

def _parse_classes(val) -> set[str]:
    """Parse a ';' or ','-separated drug-class string into a set of class names."""
    if pd.isna(val) or str(val).strip() in _SENTINEL:
        return set()
    out: set[str] = set()
    # aggregate_results.py writes ';'-separated lists; tolerate ',' too
    for tok in str(val).replace(";", ",").split(","):
        for part in tok.strip().split("/"):
            p = part.strip()
            if p and p not in _SENTINEL:
                out.add(p)
    return out


def _parse_genes(val) -> list[str]:
    """Parse a ';'-separated gene list, filtering sentinel values."""
    if pd.isna(val) or str(val).strip() in _SENTINEL:
        return []
    return [g.strip() for g in str(val).split(";")
            if g.strip() and g.strip() not in _SENTINEL]


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

    # ── Phylogroup per tip ────────────────────────────────────────────────────
    # Prefer kleborate_phylogroup directly; fall back to ST-lookup
    def _tip_phylogroup(tip: str) -> str:
        if tip not in meta.index:
            return "Unknown"
        if "kleborate_phylogroup" in meta.columns:
            pg = str(meta.loc[tip, "kleborate_phylogroup"]).strip()
            if pg not in _SENTINEL and pg not in {"NA", "nan", "None"}:
                return pg
        if "mlst_st" in meta.columns:
            st = _clean_st(meta.loc[tip, "mlst_st"])
            return ST_PHYLOGROUP.get(st, "Unknown")
        return "Unknown"

    # ── AMR genes: parse class info and build class-grouped column order ──────
    gene_col = next((c for c in ("amrfinder_acquired_genes", "amrfinder_genes")
                     if c in meta.columns), None)
    gc_col   = "amrfinder_gene_classes"   # gene=CLASS;... format

    gene_cls_map: dict[str, str] = {}
    if gc_col in meta.columns:
        for val in meta[gc_col]:
            s = str(val).strip()
            if s in _SENTINEL:
                continue
            for pair in s.split(";"):
                pair = pair.strip()
                if "=" in pair:
                    g, c = pair.split("=", 1)
                    gene_cls_map[g.strip()] = c.strip().upper()

    gene_freq: Counter = Counter()
    if gene_col:
        for v in meta[gene_col]:
            gene_freq.update(_parse_genes(v))

    MAX_GENES = 20
    raw_top_genes = [g for g, _ in gene_freq.most_common(MAX_GENES) if gene_freq[g] > 0]

    # Group genes by drug class; order classes by total prevalence
    class_genes_map: dict[str, list[str]] = defaultdict(list)
    for g in raw_top_genes:
        class_genes_map[gene_cls_map.get(g, "OTHER")].append(g)
    cls_prev = {cls: sum(gene_freq[g] for g in gs) for cls, gs in class_genes_map.items()}
    sorted_cls = sorted(class_genes_map, key=lambda c: -cls_prev[c])
    top_genes   = [g for cls in sorted_cls for g in class_genes_map[cls]]
    gene_to_cls = [gene_cls_map.get(g, "OTHER") for g in top_genes]
    n_genes = len(top_genes)

    # ── Virulence panel ───────────────────────────────────────────────────────
    KLEBORATE_VIR = [
        ("kleborate_Stx1",     "Stx1"),
        ("kleborate_Stx2",     "Stx2"),
        ("kleborate_eae",      "eae"),
        ("kleborate_ipaH",     "ipaH"),
        ("kleborate_LT",       "LT"),
        ("kleborate_ST_toxin", "ST-toxin"),
    ]
    vir_labels:  list[str]   = []
    vir_sources: list[tuple] = []   # ("kleb", col) or ("gene", gene_name, vir_gene_col)

    vir_gene_col = None
    if species == "ecoli":
        tip_set = set(tips)
        for col, label in KLEBORATE_VIR:
            if col in meta.columns:
                present = any(
                    str(meta.loc[t, col]).strip() not in {"-", "NA", "nan", "None", ""}
                    for t in tip_set if t in meta.index
                )
                if present:
                    vir_labels.append(label)
                    vir_sources.append(("kleb", col))
        if "amrfinder_virulence_genes" in meta.columns:
            vir_gene_col = "amrfinder_virulence_genes"
    else:
        vir_gene_col = next((c for c in ("amrfinder_virulence_genes", "abricate_vfdb_genes")
                             if c in meta.columns), None)

    MAX_VIR_GENES = 12
    vir_gene_freq: Counter = Counter()
    top_vir_genes: list[str] = []
    if vir_gene_col:
        for v in meta[vir_gene_col]:
            vir_gene_freq.update(_parse_genes(v))
        top_vir_genes = [g for g, _ in vir_gene_freq.most_common(MAX_VIR_GENES)
                         if vir_gene_freq[g] > 0]
        for g in top_vir_genes:
            vir_labels.append(g)
            vir_sources.append(("gene", g, vir_gene_col))

    n_vir = len(vir_labels)

    # ── Build binary matrices ─────────────────────────────────────────────────
    mat_gene = np.zeros((n, max(n_genes, 1)), dtype=np.uint8)
    if gene_col and n_genes:
        for i, tip in enumerate(tips):
            if tip in meta.index:
                row_genes = set(_parse_genes(meta.loc[tip, gene_col]))
                for j, g in enumerate(top_genes):
                    mat_gene[i, j] = 1 if g in row_genes else 0

    mat_vir = np.zeros((n, max(n_vir, 1)), dtype=np.uint8)
    if n_vir:
        for i, tip in enumerate(tips):
            if tip not in meta.index:
                continue
            for j, src in enumerate(vir_sources):
                if src[0] == "kleb":
                    val = str(meta.loc[tip, src[1]]).strip()
                    mat_vir[i, j] = 0 if val in {"-", "NA", "nan", "None", ""} else 1
                else:  # "gene"
                    vcol = src[2]
                    row_vg = set(_parse_genes(meta.loc[tip, vcol]))
                    mat_vir[i, j] = 1 if src[1] in row_vg else 0

    # ── ST colour map (qualitative palette, consistent across strip + legend) ────
    _ST_PALETTE = [
        "#e63946", "#457b9d", "#2a9d8f", "#e9c46a", "#f4a261",
        "#a8dadc", "#6a4c93", "#1982c4", "#8ac926", "#ff595e",
        "#ffca3a", "#6a994e",
    ]
    st_vals = [
        _clean_st(meta.loc[tip, "mlst_st"])
        if (tip in meta.index and "mlst_st" in meta.columns) else "Unknown"
        for tip in tips
    ]
    unique_sts = sorted(set(st_vals) - {"Unknown"})
    st_color_map: dict[str, str] = {
        st: _ST_PALETTE[i % len(_ST_PALETTE)] for i, st in enumerate(unique_sts)
    }
    st_color_map["Unknown"] = "#bab0ac"

    # ── Figure geometry ───────────────────────────────────────────────────────
    row_h   = max(0.05, min(0.18, 10.0 / n))
    fig_h   = max(5.0, n * row_h + 3.5)   # extra bottom margin for class sub-labels
    tree_w  = 4.5
    st_w    = 0.35
    strip_w = 0.35   # PG strip (E. coli only)
    vir_w   = max(1.5, n_vir   * 0.45) if n_vir   else 0
    gene_w  = max(2.0, n_genes * 0.35) if n_genes else 0

    # E. coli: tree | ST | PG | vir | AMR
    # Salmonella: tree | ST | vir | AMR  (PG strip hidden)
    col_ratios = [tree_w, st_w]
    n_axes = 2
    if species == "ecoli":
        col_ratios.append(strip_w); n_axes += 1
    if n_vir:
        col_ratios.append(vir_w);   n_axes += 1
    if n_genes:
        col_ratios.append(gene_w);  n_axes += 1

    fig, axes = plt.subplots(
        1, n_axes,
        figsize=(sum(col_ratios) + 1.2, fig_h),
        gridspec_kw={"width_ratios": col_ratios, "wspace": 0.03},
    )
    axes    = list(axes) if n_axes > 1 else [axes]
    ax_tree = axes[0]
    ax_st   = axes[1]
    _aidx   = 2
    ax_pg   = axes[_aidx] if species == "ecoli" else None
    if species == "ecoli": _aidx += 1
    ax_vir  = axes[_aidx] if n_vir   else None; _aidx += (1 if n_vir   else 0)
    ax_gene = axes[_aidx] if n_genes else None

    # ── Draw tree ─────────────────────────────────────────────────────────────
    _draw_tree(ax_tree, tree.root, pos)

    max_x = max(x for x, _ in pos.values())
    ax_tree.set_ylim(n - 0.5, -0.5)
    ax_tree.set_yticks([])
    ax_tree.set_xlabel("")   # scale bar is drawn at top; no bottom x-label needed
    ax_tree.spines["left"].set_visible(False)
    ax_tree.spines["top"].set_visible(False)
    ax_tree.spines["right"].set_visible(False)

    if n <= 60:
        x_label = max_x * 1.04
        for i, tip in enumerate(tips):
            ax_tree.text(x_label, i, tip, fontsize=5.5,
                         va="center", ha="left", fontfamily="monospace")
        ax_tree.set_xlim(-max_x * 0.02, max_x * 1.5)
    else:
        ax_tree.set_xlim(-max_x * 0.02, max_x * 1.05)

    # Scale bar — top-left, aligned with adjacent panel titles
    # Blended transform: x in data units, y in axes fraction (1.0 = top of axes)
    import matplotlib.transforms as mtransforms
    sb        = _scale_bar_value(max_x)
    blended   = mtransforms.blended_transform_factory(ax_tree.transData, ax_tree.transAxes)
    y_bar     = 1.01   # just above axes top — same level as set_title(pad=3) labels
    ax_tree.plot([0, sb], [y_bar, y_bar], color="#333333", lw=1.2,
                 transform=blended, clip_on=False)
    ax_tree.text(0, y_bar + 0.015, f"{sb:.3g}",
                 transform=blended, ha="left", va="bottom", fontsize=6.5)

    # ── ST strip (all species) ────────────────────────────────────────────────
    for i, st in enumerate(st_vals):
        ax_st.barh(i, 1, color=st_color_map[st], height=1.0, linewidth=0, align="center")
    ax_st.set_xlim(0, 1)
    ax_st.set_ylim(n - 0.5, -0.5)
    ax_st.set_yticks([])
    ax_st.set_xticks([])
    ax_st.set_title("ST", fontsize=7.5, fontweight="bold", pad=3)
    for sp in ax_st.spines.values():
        sp.set_visible(False)

    # ── Phylogroup strip (E. coli only) ──────────────────────────────────────
    seen_pgs: set[str] = set()
    if species == "ecoli" and ax_pg is not None:
        for i, tip in enumerate(tips):
            pg  = _tip_phylogroup(tip)
            col = PHYLOGROUP_COLORS.get(pg, PHYLOGROUP_COLORS["Unknown"])
            ax_pg.barh(i, 1, color=col, height=1.0, linewidth=0, align="center")
            seen_pgs.add(pg)
        ax_pg.set_xlim(0, 1)
        ax_pg.set_ylim(n - 0.5, -0.5)
        ax_pg.set_yticks([])
        ax_pg.set_xticks([])
        ax_pg.set_title("PG", fontsize=7.5, fontweight="bold", pad=3)
        for sp in ax_pg.spines.values():
            sp.set_visible(False)

    # ── Virulence heatmap ─────────────────────────────────────────────────────
    # Teal/green palette — distinct from the AMR red
    cmap_vir = mcolors.ListedColormap(["#f5f5f5", "#2a9d8f"])
    cmap_amr = mcolors.ListedColormap(["#f5f5f5", "#c0392b"])

    if ax_vir is not None and n_vir:
        ax_vir.imshow(mat_vir, aspect="auto", cmap=cmap_vir, vmin=0, vmax=1,
                      interpolation="nearest", origin="upper")

        # Build prevalence-annotated labels consistent with AMR genes panel
        vir_pct_labels = []
        for label, src in zip(vir_labels, vir_sources):
            if src[0] == "kleb":
                col = src[1]
                pct = 100 * sum(
                    1 for t in tips if t in meta.index and
                    str(meta.loc[t, col]).strip() not in {"-", "NA", "nan", "None", ""}
                ) / n
            else:
                pct = 100 * int(vir_gene_freq.get(src[1], 0)) / n
            vir_pct_labels.append(f"{label}\n({pct:.0f}%)")

        ax_vir.set_xticks(np.arange(n_vir))
        ax_vir.set_xticklabels(vir_pct_labels, rotation=40, ha="right", fontsize=7,
                                fontstyle="italic")
        ax_vir.set_yticks([])
        vir_title = ("Virulence genes" if species == "ecoli"
                     else "Virulence genes (VFDB)")
        ax_vir.set_title(vir_title, fontsize=8, fontweight="bold", pad=3)
        for sp in ax_vir.spines.values():
            sp.set_visible(False)
        # (virulence legend consolidated into tree-panel legend below)

    # ── AMR genes heatmap (columns grouped by drug class) ─────────────────────
    if ax_gene is not None and n_genes:
        ax_gene.imshow(mat_gene, aspect="auto", cmap=cmap_amr, vmin=0, vmax=1,
                       interpolation="nearest", origin="upper")

        gene_pct_labels = [
            f"{g}\n({100 * int(gene_freq[g]) / n:.0f}%)"
            for g in top_genes
        ]
        ax_gene.set_xticks(np.arange(n_genes))
        ax_gene.set_xticklabels(gene_pct_labels, rotation=40, ha="right",
                                 fontsize=6.5, fontstyle="italic")
        ax_gene.set_yticks([])
        ax_gene.set_title("Acquired AMR genes", fontsize=8, fontweight="bold", pad=3)
        for sp in ax_gene.spines.values():
            sp.set_visible(False)

        # Drug-class sub-labels below the gene tick labels
        if gene_to_cls:
            trans = ax_gene.get_xaxis_transform()
            cls_spans: list[tuple[str, int, int]] = []
            prev_cls, span_start = None, 0
            for j, cls in enumerate(gene_to_cls):
                if cls != prev_cls:
                    if prev_cls is not None:
                        cls_spans.append((prev_cls, span_start, j - 1))
                    span_start, prev_cls = j, cls
            if prev_cls is not None:
                cls_spans.append((prev_cls, span_start, len(gene_to_cls) - 1))

            # Scale bracket/label gap with axes height: shorter axes (fewer
            # samples) need more negative offset so labels clear the gene ticks
            cls_bracket_y = -max(0.07, 1.5 / n)
            cls_text_y    = cls_bracket_y - 0.02

            for cls_name, j0, j1 in cls_spans:
                x_mid = (j0 + j1) / 2.0
                label = CLASS_LABEL.get(cls_name, cls_name.capitalize())
                # Bracket line spanning the class columns
                ax_gene.plot([j0 - 0.4, j1 + 0.4], [cls_bracket_y, cls_bracket_y],
                             transform=trans, color="#555555", lw=0.9, clip_on=False)
                # Label rotated 40° to match gene tick labels
                ax_gene.text(x_mid, cls_text_y, label,
                             transform=trans, ha="right", va="top",
                             fontsize=6.5, fontweight="bold", color="#333333",
                             rotation=40, rotation_mode="anchor")

        # (AMR genes legend consolidated into tree-panel legend below)

    # ── Left legend: Phylogroup (E. coli) + Virulence + AMR present/absent ──────
    left_handles: list[mpatches.Patch] = []

    if species == "ecoli" and seen_pgs:
        left_handles.append(mpatches.Patch(color="none", label="Phylogroup"))
        for pg in sorted(seen_pgs):
            left_handles.append(
                mpatches.Patch(facecolor=PHYLOGROUP_COLORS.get(pg, "#bab0ac"),
                               label=f"  {pg}")
            )

    if n_vir:
        left_handles.append(mpatches.Patch(color="none", label="Virulence genes"))
        left_handles.append(mpatches.Patch(facecolor="#2a9d8f", label="  Present"))
        left_handles.append(mpatches.Patch(facecolor="#f5f5f5", edgecolor="#cccccc",
                                           linewidth=0.5, label="  Absent"))

    if n_genes:
        left_handles.append(mpatches.Patch(color="none", label="AMR genes"))
        left_handles.append(mpatches.Patch(facecolor="#c0392b", label="  Present"))
        left_handles.append(mpatches.Patch(facecolor="#f5f5f5", edgecolor="#cccccc",
                                           linewidth=0.5, label="  Absent"))

    if left_handles:
        ax_tree.legend(
            handles=left_handles,
            fontsize=6.5,
            loc="lower left",
            bbox_to_anchor=(0.0, -0.22),
            handlelength=1,
            frameon=True, framealpha=0.88, edgecolor="none",
            labelspacing=0.3,
        )

    # ── Bottom figure legend: ST colour blocks (wide multi-column) ────────────
    if unique_sts:
        st_handles = [mpatches.Patch(facecolor=st_color_map[st], label=st)
                      for st in unique_sts]
        if "Unknown" in set(st_vals):
            st_handles.append(mpatches.Patch(facecolor="#bab0ac", label="Unknown"))
        # ~5 columns gives 4–8 rows for typical datasets; scale with count
        ncols = max(5, math.ceil(len(st_handles) / 4))
        fig.text(0.5, -0.04, "Sequence type", ha="center", va="top",
                 fontsize=7, fontweight="bold")
        fig.legend(
            handles=st_handles,
            ncol=ncols,
            fontsize=6.5,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.12),
            frameon=True, framealpha=0.88, edgecolor="none",
            handlelength=1, handleheight=1,
            labelspacing=0.3, columnspacing=1.0,
        )

    # ── Title ─────────────────────────────────────────────────────────────────
    sp_label = "E. coli" if species == "ecoli" else "Salmonella enterica"
    fig.suptitle(
        f"{sp_label} core-SNP phylogeny with virulence & AMR profiles  (n = {n} isolates)",
        fontsize=11, fontweight="bold", y=1.01,
    )

    # ── Save ──────────────────────────────────────────────────────────────────
    for ext in ("pdf", "png"):
        out = outdir / f"{prefix}_tree_amr.{ext}"
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
