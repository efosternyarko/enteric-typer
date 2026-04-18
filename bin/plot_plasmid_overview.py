#!/usr/bin/env python3
"""
plot_plasmid_overview.py — Multi-panel plasmid overview figure.

Layout
------
  Panel A (top, full width):
      Midpoint-rooted SNP tree  |  ST strip  |  PG strip
      | Plasmid replicon heatmap (one column per replicon,
        coloured by dominant AMR drug class)

  Panel B (bottom-left):
      Horizontal stacked-bar chart — replicon prevalence broken down
      by AMR drug class.  Each isolate carrying a replicon is allocated
      to exactly one class segment using clinical priority order, so
      segments are non-overlapping and sum to total replicon prevalence.
      Grey = replicon present but no AMR genes on same contig.

  Panel C (bottom-right):
      Plasmid–AMR co-occurrence bubble matrix (same as standalone
      ecoli_plasmid_amr_map_bubble figure).

All three panels show the same top-N replicons in the same order.

Output: {prefix}_plasmid_overview.{pdf,png}
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from Bio import Phylo
    HAS_PHYLO = True
except ImportError:
    HAS_PHYLO = False

# ── Global aesthetics ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "Helvetica Neue", "DejaVu Sans"],
    "font.size":          8,
    "axes.labelsize":     9,
    "axes.titlesize":     9,
    "axes.titleweight":   "bold",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "xtick.labelsize":    7.5,
    "ytick.labelsize":    8,
    "legend.fontsize":    7,
    "legend.frameon":     False,
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.1,
})

# ── Drug-class palette (consistent across all enteric-typer figures) ──────────
CLASS_PRIORITY = [
    "BETA-LACTAM", "QUINOLONE", "COLISTIN", "AMINOGLYCOSIDE",
    "TETRACYCLINE", "MACROLIDE", "PHENICOL", "SULFONAMIDE",
    "TRIMETHOPRIM", "FOSFOMYCIN", "STREPTOTHRICIN", "NITROFURAN",
]
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
    "No AMR":          "#bab0ac",
}
CLASS_LABEL = {
    "BETA-LACTAM":    "Beta-lactam",   "QUINOLONE":      "Quinolone",
    "COLISTIN":       "Colistin",      "AMINOGLYCOSIDE": "Aminoglycoside",
    "TETRACYCLINE":   "Tetracycline",  "MACROLIDE":      "Macrolide",
    "PHENICOL":       "Phenicol",      "SULFONAMIDE":    "Sulfonamide",
    "TRIMETHOPRIM":   "Trimethoprim",  "FOSFOMYCIN":     "Fosfomycin",
    "STREPTOTHRICIN": "Streptothricin","NITROFURAN":      "Nitrofuran",
    "No AMR":         "No AMR",
}
# Pastel version of CLASS_COLOR used for the plasmid heatmap in Panel A so it
# is visually distinct from the saturated AMR-gene heatmap in fig2/fig7.
# Each colour is blended 50% toward white.
PLASMID_CLASS_COLOR = {
    k: mcolors.to_hex(tuple(c * 0.55 + 0.45 for c in mcolors.to_rgb(v)))
    for k, v in CLASS_COLOR.items()
}
# Override two classes whose pastels land too close to PG Unknown (#bab0ac, grey):
#   COLISTIN  pastel #b9b9b9  →  d=0.062 from PG Unknown  (nearly identical)
#   PHENICOL  pastel #c0a29c  →  d=0.087 from PG Unknown  (very close warm-grey)
PLASMID_CLASS_COLOR["COLISTIN"] = "#c5e8a0"  # light sage-green
PLASMID_CLASS_COLOR["PHENICOL"] = "#d4b8f0"  # light lavender

PHYLOGROUP_COLORS = {
    "A":  "#4e79a7", "B1": "#59a14f", "B2": "#f28e2b",
    "C":  "#76b7b2", "D":  "#e15759", "E":  "#b07aa1",
    "F":  "#ff9da7", "Unknown": "#bab0ac",
}
ST_PHYLOGROUP: dict[str, str] = {
    "ST10":"A",  "ST12":"B2", "ST14":"D",  "ST23":"A",  "ST34":"A",
    "ST38":"D",  "ST52":"A",  "ST58":"B1", "ST69":"D",  "ST73":"B2",
    "ST88":"D",  "ST95":"B2", "ST117":"D", "ST127":"B2","ST131":"B2",
    "ST141":"B1","ST144":"A", "ST155":"B1","ST167":"A",  "ST218":"A",
    "ST354":"B1","ST393":"A", "ST405":"D", "ST410":"C", "ST443":"A",
    "ST448":"A", "ST453":"A", "ST617":"B1","ST636":"A",  "ST648":"F",
    "ST1193":"B2",
}
_SENTINEL = {"", "-", "NA", "nan", "None", "none"}
# Single colour used for plasmid-present blocks in Panel C (tree heatmap).
# Plain presence/absence — drug-class breakdown is in Panels A and B.
_PLASMID_PRESENT_COLOR = "#2471a3"   # medium navy blue
_REPLICON_SENTINELS = {
    "no_replicon", "No replicons found", "No replicon",
    "no_replicons", "No replicons", "No plasmids found",
}
_ST_PALETTE = [
    "#e63946","#457b9d","#2a9d8f","#e9c46a","#f4a261",
    "#a8dadc","#6a4c93","#1982c4","#8ac926","#ff595e",
    "#ffca3a","#6a994e",
]

# ── Replicon-family palette (simplified layout) ───────────────────────────────
# Used when same-contig AMR co-occurrence is too sparse for the full drug-class
# breakdown (typical of fragmented short-read assemblies).
FAMILY_COLOR = {
    "IncF":   "#4e79a7",
    "IncI":   "#f28e2b",
    "IncN":   "#e15759",
    "IncA/C": "#76b7b2",
    "IncH":   "#59a14f",
    "IncL/M": "#edc948",
    "IncX":   "#b07aa1",
    "IncP":   "#ff9da7",
    "IncQ":   "#9c755f",
    "Col":    "#bab0ac",
    "Other":  "#d3d3d3",
}

# Threshold: if < this fraction of plasmid-carrying isolates have at least one
# AMR gene on the same contig as a replicon, use the simplified layout.
_SIMPLE_LAYOUT_THRESHOLD = 0.10


def _replicon_family(rep: str) -> str:
    """Return broad Inc/Col family for colour grouping."""
    r = rep.split("(")[0].rstrip("0123456789_")
    if r.startswith("IncF"):  return "IncF"
    if r.startswith("IncI"):  return "IncI"
    if r.startswith("IncN"):  return "IncN"
    if r.startswith("IncA") or r.startswith("IncC"):  return "IncA/C"
    if r.startswith("IncH"):  return "IncH"
    if r.startswith("IncL") or r.startswith("IncM"):  return "IncL/M"
    if r.startswith("IncX"):  return "IncX"
    if r.startswith("IncP"):  return "IncP"
    if r.startswith("IncQ"):  return "IncQ"
    if r.startswith("Col"):   return "Col"
    return "Other"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _save(fig: plt.Figure, outdir: Path, stem: str) -> None:
    for ext in ("pdf", "png"):
        p = outdir / f"{stem}.{ext}"
        fig.savefig(p, dpi=300 if ext == "pdf" else 150, bbox_inches="tight")
        print(f"  Saved: {p}", file=sys.stderr)
    plt.close(fig)


def _norm(s: str) -> str:
    return s.replace("_complete", "").replace("_incomplete", "")


def _clean_st(raw) -> str:
    s = str(raw).strip() if pd.notna(raw) else ""
    if s in ("", "-", "NA", "nan", "None", "No ST predicted"):
        return "Unknown"
    if "." in s:
        try:
            s = str(int(float(s)))
        except ValueError:
            pass
    return s if s.startswith("ST") else f"ST{s}"


def _dominant_class(classes_set: set) -> str:
    for cls in CLASS_PRIORITY:
        if cls in classes_set:
            return cls
    return "No AMR"


def _expand_classes(val) -> list[str]:
    if pd.isna(val):
        return []
    s = str(val).strip()
    if not s or s in _SENTINEL:
        return []
    return [c.strip() for c in s.split(";") if c.strip() and c.strip() not in _SENTINEL]


# ── Tree helpers ──────────────────────────────────────────────────────────────

def _layout(root):
    tips: list[str] = []
    pos:  dict      = {}

    def _rec(clade, x):
        xh = x + (clade.branch_length or 0.0)
        if clade.is_terminal():
            y = float(len(tips))
            tips.append(clade.name or "")
        else:
            ys = [_rec(c, xh) for c in clade.clades]
            y  = (min(ys) + max(ys)) / 2.0
        pos[id(clade)] = (xh, y)
        return y

    _rec(root, 0.0)
    return pos, tips


def _draw_tree(ax, root, pos):
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


# ── Data loading ──────────────────────────────────────────────────────────────

def _build_plasmid_data(
    df_map: pd.DataFrame,
    tips: list[str],
    top_n: int,
    sample_acq_classes: Optional[dict] = None,
):
    """
    Returns
    -------
    top_reps     : list[str]  — top_n replicons by isolate prevalence
    rep_dom_color: dict       — replicon → dominant-class hex colour
    mat          : ndarray    — (n_tips, n_reps) binary presence
    pair_dom     : dict       — (norm_sample_id, replicon) → dominant class

    sample_acq_classes (optional): dict mapping normalised sample IDs to their
    acquired drug class sets (from the metadata TSV, AMRrules-filtered).  When
    provided, each (sample, replicon) pair is supplemented with the sample's
    acquired classes so that fragmented short-read assemblies — where the
    replicon and AMR genes land on different contigs — are handled correctly.
    """
    tip_norm = {t: _norm(t) for t in tips}
    norm_set = set(tip_norm.values())

    # per-replicon: set of normalised sample IDs + all drug classes seen
    rep_samples: dict[str, set] = {}
    rep_classes: dict[str, set] = {}
    # per (sample, replicon): dominant drug class
    pair_dom: dict[tuple, str]  = {}

    for _, row in df_map.iterrows():
        rep = str(row["replicon"]).strip()
        sid = str(row["sample_id"]).strip()
        if not rep or rep in _SENTINEL or rep in _REPLICON_SENTINELS:
            continue
        sid_n = _norm(sid)
        if sid_n not in norm_set and sid not in tips:
            continue
        key = sid_n if sid_n in norm_set else sid
        rep_samples.setdefault(rep, set()).add(key)

        classes = set(_expand_classes(row.get("drug_classes")))
        rep_classes.setdefault(rep, set()).update(classes)
        dom = _dominant_class(classes)
        pair_key = (key, rep)
        # keep highest-priority class seen for this (sample, replicon) pair
        if pair_key not in pair_dom or (
            CLASS_PRIORITY.index(dom) < CLASS_PRIORITY.index(pair_dom[pair_key])
            if dom != "No AMR" and pair_dom[pair_key] != "No AMR"
            else dom != "No AMR"
        ):
            pair_dom[pair_key] = dom

    if not rep_samples:
        return [], {}, np.zeros((len(tips), 0), dtype=np.uint8), {}

    top_reps = sorted(rep_samples, key=lambda r: len(rep_samples[r]), reverse=True)[:top_n]

    # Use the pastel palette so the heatmap is visually distinct from
    # the saturated AMR-gene heatmap in fig2/fig7.
    rep_dom_color = {
        r: PLASMID_CLASS_COLOR[_dominant_class(rep_classes.get(r, set()))]
        for r in top_reps
    }

    n = len(tips)
    mat = np.zeros((n, len(top_reps)), dtype=np.uint8)
    for j, rep in enumerate(top_reps):
        present = rep_samples[rep]
        for i, tip in enumerate(tips):
            if tip_norm[tip] in present or tip in present:
                mat[i, j] = 1

    return top_reps, rep_dom_color, mat, pair_dom


# ── Panel A: mini-tree + ST + PG + plasmid heatmap ───────────────────────────

def _draw_panel_a(
    fig: plt.Figure,
    gs_top,
    tree_root,
    pos: dict,
    tips: list[str],
    meta: pd.DataFrame,
    top_reps: list[str],
    mat: np.ndarray,
) -> plt.Axes:
    """Draws Panel C (tree + plasmid heatmap). Returns ax_tree for label placement.
    Heatmap uses a single colour (_PLASMID_PRESENT_COLOR) for presence — drug-class
    breakdown is shown in Panels A and B instead."""
    n = len(tips)

    # ── sub-axes widths ───────────────────────────────────────────────────────
    n_reps  = len(top_reps)
    sub_ratios = [4.5, 0.32, 0.32, max(2.0, n_reps * 0.38)]
    sub_gs = gridspec.GridSpecFromSubplotSpec(
        1, 4,
        subplot_spec=gs_top,
        width_ratios=sub_ratios,
        wspace=0.02,
    )
    ax_tree = fig.add_subplot(sub_gs[0])
    ax_st   = fig.add_subplot(sub_gs[1])
    ax_pg   = fig.add_subplot(sub_gs[2])
    ax_plas = fig.add_subplot(sub_gs[3])

    # ── tree ──────────────────────────────────────────────────────────────────
    _draw_tree(ax_tree, tree_root, pos)
    max_x = max(x for x, _ in pos.values())
    ax_tree.set_ylim(n - 0.5, -0.5)
    ax_tree.set_yticks([])
    ax_tree.set_xticks([])
    for sp in ("left", "top", "right", "bottom"):
        ax_tree.spines[sp].set_visible(False)
    ax_tree.set_xlim(-max_x * 0.02, max_x * 1.05)
    ax_tree.set_title("SNP phylogeny", fontsize=7.5, fontweight="bold", pad=3)

    # scale bar
    import matplotlib.transforms as mtrans
    target = max_x / 5.0
    mag    = 10 ** math.floor(math.log10(target)) if target > 0 else 0.001
    sb     = round(target / mag) * mag
    blended = mtrans.blended_transform_factory(ax_tree.transData, ax_tree.transAxes)
    ax_tree.plot([0, sb], [1.02, 1.02], color="#333333", lw=1.0,
                 transform=blended, clip_on=False)
    ax_tree.text(0, 1.035, f"{sb:.3g}", transform=blended,
                 ha="left", va="bottom", fontsize=6)

    # ── ST strip ──────────────────────────────────────────────────────────────
    st_vals = [
        _clean_st(meta.loc[t, "mlst_st"]) if t in meta.index and "mlst_st" in meta.columns
        else "Unknown"
        for t in tips
    ]
    unique_sts = sorted(set(st_vals) - {"Unknown"})
    st_color_map = {st: _ST_PALETTE[i % len(_ST_PALETTE)]
                    for i, st in enumerate(unique_sts)}
    st_color_map["Unknown"] = "#bab0ac"

    for i, st in enumerate(st_vals):
        ax_st.barh(i, 1, color=st_color_map[st], height=1.0, linewidth=0, align="center")
    ax_st.set_xlim(0, 1); ax_st.set_ylim(n - 0.5, -0.5)
    ax_st.set_yticks([]); ax_st.set_xticks([])
    ax_st.set_title("ST", fontsize=7, fontweight="bold", pad=3)
    for sp in ax_st.spines.values(): sp.set_visible(False)

    # ── PG strip ──────────────────────────────────────────────────────────────
    seen_pgs: set = set()
    for i, tip in enumerate(tips):
        pg = "Unknown"
        if tip in meta.index:
            for col in ("clermont_phylogroup", "kleborate_phylogroup"):
                if col in meta.columns:
                    v = str(meta.loc[tip, col]).strip()
                    if v not in _SENTINEL and v not in {"NA", "nan", "None",
                                                         "EC_control_fail"}:
                        pg = v; break
            if pg == "Unknown" and "mlst_st" in meta.columns:
                pg = ST_PHYLOGROUP.get(_clean_st(meta.loc[tip, "mlst_st"]), "Unknown")
        col = PHYLOGROUP_COLORS.get(pg, PHYLOGROUP_COLORS["Unknown"])
        ax_pg.barh(i, 1, color=col, height=1.0, linewidth=0, align="center")
        seen_pgs.add(pg)
    ax_pg.set_xlim(0, 1); ax_pg.set_ylim(n - 0.5, -0.5)
    ax_pg.set_yticks([]); ax_pg.set_xticks([])
    ax_pg.set_title("PG", fontsize=7, fontweight="bold", pad=3)
    for sp in ax_pg.spines.values(): sp.set_visible(False)

    # ── Plasmid heatmap (single presence colour) ─────────────────────────────
    if n_reps > 0:
        absent_rgb  = np.array(mcolors.to_rgb("#f5f5f5"))
        present_rgb = np.array(mcolors.to_rgb(_PLASMID_PRESENT_COLOR))
        rgb_mat     = np.tile(absent_rgb, (n, n_reps, 1)).astype(float)
        for i in range(n):
            for j in range(n_reps):
                if mat[i, j]:
                    rgb_mat[i, j] = present_rgb

        ax_plas.imshow(rgb_mat, aspect="auto", interpolation="nearest", origin="upper")
        rep_pct = [100 * int(mat[:, j].sum()) / n for j in range(n_reps)]
        plas_labels = [f"{r}\n({rep_pct[j]:.0f}%)" for j, r in enumerate(top_reps)]
        ax_plas.set_xticks(np.arange(n_reps))
        ax_plas.set_xticklabels(plas_labels, rotation=40, ha="right", fontsize=6.5)
        ax_plas.set_yticks([])
        ax_plas.set_title("Plasmid replicons (presence/absence)",
                           fontsize=7.5, fontweight="bold", pad=3)
        for sp in ax_plas.spines.values(): sp.set_visible(False)

    # ── Legends — PG and ST side by side at the same y ───────────────────────
    # PG legend (left)
    pg_handles = [mpatches.Patch(color="none", label="Phylogroup")]
    for pg in sorted(seen_pgs):
        pg_handles.append(mpatches.Patch(
            facecolor=PHYLOGROUP_COLORS.get(pg, "#bab0ac"), label=f"  {pg}"))
    leg_pg = ax_tree.legend(handles=pg_handles, fontsize=6, loc="lower left",
                             bbox_to_anchor=(0.0, -0.28), handlelength=1,
                             frameon=True, framealpha=0.88, edgecolor="none",
                             labelspacing=0.25)
    ax_tree.add_artist(leg_pg)  # keep it when we add the ST legend

    # ST legend (right of PG, same y)
    if unique_sts:
        st_handles = [mpatches.Patch(facecolor=st_color_map[st], label=st)
                      for st in unique_sts]
        if "Unknown" in set(st_vals):
            st_handles.append(mpatches.Patch(facecolor="#bab0ac", label="Unknown"))
        ncols = max(5, math.ceil(len(st_handles) / 3))
        ax_tree.legend(
            handles=st_handles, ncol=ncols, fontsize=6,
            title="Sequence type (ST)", title_fontsize=6.5,
            loc="lower left", bbox_to_anchor=(0.22, -0.28),
            frameon=True, framealpha=0.88, edgecolor="none",
            handlelength=1, labelspacing=0.25, columnspacing=0.8,
        )

    return ax_tree


# ── Panel A: stacked-bar chart (all drug classes by dominant-class allocation) ─

def _draw_panel_b(
    ax: plt.Axes,
    df_map: pd.DataFrame,
    top_reps: list[str],
    n_isolates: int,
    pair_dom: dict,
) -> None:
    """
    For each replicon, stack segments coloured by drug class.
    Each (sample, replicon) pair contributes to exactly one class segment
    (dominant class by priority), so segments are non-overlapping.
    """
    # Build counts: replicon → {class → count}
    class_order = CLASS_PRIORITY + ["No AMR"]
    rep_cls_count: dict[str, dict[str, int]] = {r: {c: 0 for c in class_order}
                                                  for r in top_reps}
    for (sid, rep), dom in pair_dom.items():
        if rep in rep_cls_count:
            rep_cls_count[rep][dom] = rep_cls_count[rep].get(dom, 0) + 1

    # Which classes are actually present in any bar?
    used_classes = [c for c in class_order
                    if any(rep_cls_count[r].get(c, 0) > 0 for r in top_reps)]

    y      = np.arange(len(top_reps))
    left   = np.zeros(len(top_reps))
    for cls in used_classes:
        vals = np.array([100 * rep_cls_count[r].get(cls, 0) / n_isolates
                         for r in top_reps])
        ax.barh(y, vals, left=left, color=CLASS_COLOR[cls],
                label=CLASS_LABEL.get(cls, cls), edgecolor="white",
                linewidth=0.25, height=0.76)
        left += vals

    ax.set_yticks(y)
    ax.set_yticklabels(top_reps, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (% of isolates)", fontsize=8)
    ax.set_title("Replicon prevalence by AMR drug class\n"
                 "(each isolate assigned to highest-priority class on that replicon)",
                 fontweight="bold")
    ax.set_xlim(0, 108)
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))

    # Legend — outside the plot, below x-axis
    leg_handles = [mpatches.Patch(facecolor=CLASS_COLOR[c],
                                   label=CLASS_LABEL.get(c, c))
                   for c in used_classes]
    ncol_leg = max(2, math.ceil(len(leg_handles) / 2))
    ax.legend(handles=leg_handles, fontsize=6.5, ncol=ncol_leg,
              title="AMR drug class", title_fontsize=7,
              loc="upper left", bbox_to_anchor=(0.0, -0.18),
              frameon=True, framealpha=0.88, edgecolor="none",
              handlelength=1, columnspacing=0.8)


# ── Panel C: bubble matrix ────────────────────────────────────────────────────

def _draw_panel_c(
    ax: plt.Axes,
    df_map: pd.DataFrame,
    top_reps: list[str],
    n_isolates: int,
) -> None:
    # replicon → {class → set of sample_ids}
    matrix: dict[str, dict[str, set]] = {r: {} for r in top_reps}
    rep_isolates: dict[str, set] = {r: set() for r in top_reps}

    for _, row in df_map.iterrows():
        rep = str(row["replicon"]).strip()
        sid = str(row["sample_id"]).strip()
        if rep not in matrix:
            continue
        sid_n = _norm(sid)
        rep_isolates[rep].add(sid_n)
        for cls in [c.strip() for c in str(row.get("drug_classes","")).split(";")
                    if c.strip() and c.strip() not in _SENTINEL and c.strip() in CLASS_COLOR]:
            matrix[rep].setdefault(cls, set()).add(sid_n)

    all_classes = set()
    for r in top_reps:
        all_classes |= set(matrix[r].keys())
    cols = [c for c in CLASS_PRIORITY if c in all_classes]

    if not cols:
        ax.text(0.5, 0.5, "No drug class data available",
                ha="center", va="center", transform=ax.transAxes, color="grey")
        return

    max_r = 0.42
    for i, rep in enumerate(top_reps):
        for j, cls in enumerate(cols):
            count = len(matrix[rep].get(cls, set()))
            pct   = 100 * count / n_isolates
            if pct > 0:
                r = max_r * np.sqrt(pct / 100)
                ax.scatter(j, i, s=(r * 72)**2, color=CLASS_COLOR[cls],
                           alpha=0.80, linewidths=0.4, edgecolors="white", zorder=3)
                if pct >= 5:
                    ax.text(j, i, f"{pct:.0f}%",
                            ha="center", va="center", fontsize=5.5,
                            color="white", fontweight="bold", zorder=4)

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([CLASS_LABEL.get(c, c) for c in cols],
                       rotation=40, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(top_reps)))
    ax.set_yticklabels(
        [f"{r}  ({len(rep_isolates[r])/n_isolates*100:.0f}%)" for r in top_reps],
        fontsize=7.5
    )
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylim(-0.6, len(top_reps) - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("AMR drug class", fontsize=8)
    ax.set_title("Plasmid–AMR co-occurrence\n"
                 "(bubble area proportional to % of isolates)",
                 fontweight="bold")
    for j in range(len(cols)):
        ax.axvline(j, color="#eeeeee", lw=0.5, zorder=1)
    for i in range(len(top_reps)):
        ax.axhline(i, color="#eeeeee", lw=0.5, zorder=1)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)


# ── Panel A (simplified): plain replicon prevalence coloured by Inc family ────

def _draw_panel_simple_bars(
    ax: plt.Axes,
    top_reps: list[str],
    mat: np.ndarray,
    n_isolates: int,
) -> None:
    """
    Simple horizontal bar chart of replicon prevalence, coloured by Inc/Col
    family.  Used for short-read assemblies where same-contig AMR attribution
    is not meaningful.
    """
    y    = np.arange(len(top_reps))
    pcts = [100 * int(mat[:, j].sum()) / n_isolates for j in range(len(top_reps))]
    fams = [_replicon_family(r) for r in top_reps]
    cols = [FAMILY_COLOR.get(f, FAMILY_COLOR["Other"]) for f in fams]

    ax.barh(y, pcts, color=cols, height=0.76, edgecolor="white", linewidth=0.25)
    ax.set_yticks(y)
    ax.set_yticklabels(top_reps, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (% of isolates)", fontsize=8)
    ax.set_title("Replicon prevalence by Inc/Col family", fontweight="bold")
    ax.set_xlim(0, 108)
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))

    # Legend — unique families in the order they appear
    seen = list(dict.fromkeys(fams))
    handles = [mpatches.Patch(facecolor=FAMILY_COLOR.get(f, FAMILY_COLOR["Other"]),
                               label=f) for f in seen]
    ncol = max(1, math.ceil(len(handles) / 4))
    ax.legend(handles=handles, fontsize=6.5, ncol=ncol,
              title="Replicon family", title_fontsize=7,
              loc="upper left", bbox_to_anchor=(0.0, -0.18),
              frameon=True, framealpha=0.88, edgecolor="none",
              handlelength=1, columnspacing=0.8)


# ── Main ──────────────────────────────────────────────────────────────────────

def plot_plasmid_overview(
    treefile:     Optional[Path],
    metadata_tsv: Path,
    plasmid_map:  Path,
    outdir:       Path,
    prefix:       str,
    top_n:        int = 15,
) -> None:

    # ── Load data ─────────────────────────────────────────────────────────────
    df_meta = pd.read_csv(str(metadata_tsv), sep="\t", low_memory=False)
    df_meta["sample"] = df_meta["sample"].astype(str).str.strip()
    meta = df_meta.set_index("sample")
    n_isolates = len(meta)

    df_map = pd.read_csv(str(plasmid_map), sep="\t", low_memory=False)

    # ── Tree ──────────────────────────────────────────────────────────────────
    tree_root = pos = tips = None
    if treefile and treefile.is_file() and HAS_PHYLO:
        try:
            tree = Phylo.read(str(treefile), "newick")
            tree.root_at_midpoint()
            pos, tips = _layout(tree.root)
            tree_root = tree.root
            print(f"  Tree: {len(tips)} tips", file=sys.stderr)
        except Exception as e:
            print(f"WARNING: could not read tree — {e}", file=sys.stderr)

    if tips is None:
        # Fall back to metadata order if no tree
        tips = list(meta.index)
        pos  = None

    # ── Build plasmid data (same top_reps for all panels) ─────────────────────
    top_reps, rep_dom_color, mat, pair_dom = _build_plasmid_data(df_map, tips, top_n)
    if not top_reps:
        print("WARNING: no plasmid replicon data — aborting", file=sys.stderr)
        return
    print(f"  Top {len(top_reps)} replicons selected", file=sys.stderr)

    # ── Auto-detect layout ────────────────────────────────────────────────────
    # In fragmented short-read assemblies the replicon sequence and AMR gene
    # cassettes typically land on *different* contigs, so same-contig matching
    # yields almost no drug-class signal.  If fewer than 10 % of plasmid-
    # carrying isolates have any AMR on the same contig as a replicon, use the
    # simplified 2-panel layout (replicon prevalence bars + tree heatmap)
    # instead of the full 3-panel layout (stacked-drug-class bars + bubble
    # matrix + tree heatmap).
    n_rep_samples  = len({sid for sid, _ in pair_dom})
    n_with_amr     = len({sid for (sid, _), dom in pair_dom.items()
                          if dom != "No AMR"})
    frac_with_amr  = n_with_amr / max(1, n_rep_samples)
    use_simple     = frac_with_amr < _SIMPLE_LAYOUT_THRESHOLD
    print(f"  Same-contig AMR fraction: {frac_with_amr:.2f} "
          f"({'simple' if use_simple else 'full'} layout)",
          file=sys.stderr)

    # ── Figure dimensions ─────────────────────────────────────────────────────
    n           = len(tips)
    row_h       = max(0.05, min(0.16, 9.0 / n))
    panel_tree_h = max(6.0,  n * row_h + 2.0)
    panel_bars_h = max(5.0,  top_n * 0.40 + 2.5)
    fig_w        = 22.0

    # ── Draw ──────────────────────────────────────────────────────────────────
    if use_simple:
        # ── Simplified layout: 2 rows (bars top, tree+heatmap bottom) ─────────
        fig_h = panel_tree_h + panel_bars_h + 1.5
        fig   = plt.figure(figsize=(fig_w, fig_h))
        outer = gridspec.GridSpec(
            2, 1, figure=fig,
            height_ratios=[panel_bars_h, panel_tree_h],
            hspace=0.35,
        )
        ax_a = fig.add_subplot(outer[0])
        _draw_panel_simple_bars(ax_a, top_reps, mat, n_isolates)

        if tree_root is not None and pos is not None:
            ax_tree = _draw_panel_a(fig, outer[1], tree_root, pos, tips, meta,
                                    top_reps, mat)
        else:
            ax_tree = fig.add_subplot(outer[1])
            ax_tree.text(0.5, 0.5, "No tree available", ha="center",
                         va="center", transform=ax_tree.transAxes,
                         color="grey", fontsize=12)
            ax_tree.axis("off")

        for ax, lbl in [(ax_a, "A"), (ax_tree, "B")]:
            ax.text(-0.04, 1.03, lbl, transform=ax.transAxes,
                    fontsize=14, fontweight="bold", va="bottom", ha="left")

        _save(fig, outdir, f"{prefix}_fig4_plasmid_overview")

        # Individual panels
        indiv_dir = outdir / "individual_plasmid_plots"
        indiv_dir.mkdir(parents=True, exist_ok=True)

        fig_a, ax_ia = plt.subplots(figsize=(8, max(4.0, top_n * 0.40 + 2.0)))
        _draw_panel_simple_bars(ax_ia, top_reps, mat, n_isolates)
        ax_ia.text(-0.04, 1.03, "A", transform=ax_ia.transAxes,
                   fontsize=14, fontweight="bold", va="bottom", ha="left")
        _save(fig_a, indiv_dir, f"{prefix}_fig4_panel_a_replicon_bars")

        if tree_root is not None and pos is not None:
            fig_b = plt.figure(figsize=(fig_w, panel_tree_h + 1.5))
            gs_b  = gridspec.GridSpec(1, 1, figure=fig_b)
            ax_ib = _draw_panel_a(fig_b, gs_b[0], tree_root, pos, tips, meta,
                                   top_reps, mat)
            ax_ib.text(-0.04, 1.03, "B", transform=ax_ib.transAxes,
                       fontsize=14, fontweight="bold", va="bottom", ha="left")
            _save(fig_b, indiv_dir, f"{prefix}_fig4_panel_b_tree_heatmap")

    else:
        # ── Full layout: 3 panels (stacked-drug bars | bubble matrix | tree) ──
        fig_h = panel_tree_h + panel_bars_h + 1.5
        fig   = plt.figure(figsize=(fig_w, fig_h))
        outer = gridspec.GridSpec(
            2, 1, figure=fig,
            height_ratios=[panel_bars_h, panel_tree_h],
            hspace=0.35,
        )
        top_gs = gridspec.GridSpecFromSubplotSpec(
            1, 2,
            subplot_spec=outer[0],
            wspace=0.35,
            width_ratios=[0.42, 0.58],
        )
        ax_a = fig.add_subplot(top_gs[0])
        ax_b = fig.add_subplot(top_gs[1])

        _draw_panel_b(ax_a, df_map, top_reps, n_isolates, pair_dom)
        _draw_panel_c(ax_b, df_map, top_reps, n_isolates)

        if tree_root is not None and pos is not None:
            ax_tree = _draw_panel_a(fig, outer[1], tree_root, pos, tips, meta,
                                    top_reps, mat)
        else:
            ax_tree = fig.add_subplot(outer[1])
            ax_tree.text(0.5, 0.5, "No tree available", ha="center",
                         va="center", transform=ax_tree.transAxes,
                         color="grey", fontsize=12)
            ax_tree.axis("off")

        for ax, lbl in [(ax_a, "A"), (ax_b, "B"), (ax_tree, "C")]:
            ax.text(-0.04, 1.03, lbl, transform=ax.transAxes,
                    fontsize=14, fontweight="bold", va="bottom", ha="left")

        _save(fig, outdir, f"{prefix}_fig4_plasmid_overview")

        # Individual panels
        indiv_dir = outdir / "individual_plasmid_plots"
        indiv_dir.mkdir(parents=True, exist_ok=True)

        fig_a, ax_ia = plt.subplots(figsize=(8, max(4.0, top_n * 0.40 + 2.0)))
        _draw_panel_b(ax_ia, df_map, top_reps, n_isolates, pair_dom)
        ax_ia.text(-0.04, 1.03, "A", transform=ax_ia.transAxes,
                   fontsize=14, fontweight="bold", va="bottom", ha="left")
        _save(fig_a, indiv_dir, f"{prefix}_fig4_panel_a_replicon_bars")

        fig_b, ax_ib = plt.subplots(figsize=(9, max(4.0, top_n * 0.40 + 2.0)))
        _draw_panel_c(ax_ib, df_map, top_reps, n_isolates)
        ax_ib.text(-0.04, 1.03, "B", transform=ax_ib.transAxes,
                   fontsize=14, fontweight="bold", va="bottom", ha="left")
        _save(fig_b, indiv_dir, f"{prefix}_fig4_panel_b_bubble_matrix")

        if tree_root is not None and pos is not None:
            fig_c = plt.figure(figsize=(fig_w, panel_tree_h + 1.5))
            gs_c  = gridspec.GridSpec(1, 1, figure=fig_c)
            ax_ic = _draw_panel_a(fig_c, gs_c[0], tree_root, pos, tips, meta,
                                   top_reps, mat)
            ax_ic.text(-0.04, 1.03, "C", transform=ax_ic.transAxes,
                       fontsize=14, fontweight="bold", va="bottom", ha="left")
            _save(fig_c, indiv_dir, f"{prefix}_fig4_panel_c_tree_heatmap")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(description="Plasmid overview multi-panel figure")
    p.add_argument("--tree",        "-t", default=None,
                   help="Newick tree file (optional; panels B/C still produced without it)")
    p.add_argument("--metadata",    "-m", required=True,
                   help="enteric-typer results TSV (ecoli_typer_results.tsv)")
    p.add_argument("--plasmid_map", "-q", required=True,
                   help="Aggregate plasmid_amr_map.tsv")
    p.add_argument("--outdir",      "-o", default=".")
    p.add_argument("--prefix",      "-p", default="ecoli")
    p.add_argument("--top_n",       "-n", type=int, default=15,
                   help="Top N replicons to display across all panels (default: 15)")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    plot_plasmid_overview(
        treefile     = Path(args.tree) if args.tree else None,
        metadata_tsv = Path(args.metadata),
        plasmid_map  = Path(args.plasmid_map),
        outdir       = outdir,
        prefix       = args.prefix,
        top_n        = args.top_n,
    )


if __name__ == "__main__":
    main()
