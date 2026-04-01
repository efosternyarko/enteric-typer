#!/usr/bin/env python3
"""
Generate 4 AMRnet-style tile heatmaps for review:
  1. Salmonella — AMR drug class % by MLST ST
  2. Salmonella — AMR drug class % by serovar
  3. E. coli    — AMR drug class % by MLST ST
  4. E. coli    — AMR drug class % by phylogroup
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from collections import Counter

# ---------------------------------------------------------------------------
# Colour map  –  0 % = mid-grey; 0+% → cream → orange → dark-red → purple
# ---------------------------------------------------------------------------
GREY_0   = np.array([0.502, 0.502, 0.502, 1.0])   # #808080
CREAM    = np.array([1.000, 0.973, 0.882, 1.0])   # #FFF8E1
ORANGE   = np.array([1.000, 0.600, 0.302, 1.0])   # #FF994D
DARK_RED = np.array([0.741, 0.031, 0.051, 1.0])   # #BD080D
PURPLE   = np.array([0.290, 0.004, 0.549, 1.0])   # #4A018C

# Custom cmap only used for pct > 0; we colour 0-cells separately
_cmap_data = {
    "red":   [(0.0, CREAM[0], CREAM[0]),
              (0.5, ORANGE[0], ORANGE[0]),
              (0.8, DARK_RED[0], DARK_RED[0]),
              (1.0, PURPLE[0], PURPLE[0])],
    "green": [(0.0, CREAM[1], CREAM[1]),
              (0.5, ORANGE[1], ORANGE[1]),
              (0.8, DARK_RED[1], DARK_RED[1]),
              (1.0, PURPLE[1], PURPLE[1])],
    "blue":  [(0.0, CREAM[2], CREAM[2]),
              (0.5, ORANGE[2], ORANGE[2]),
              (0.8, DARK_RED[2], DARK_RED[2]),
              (1.0, PURPLE[2], PURPLE[2])],
}
POS_CMAP = LinearSegmentedColormap("amrnet_pos", _cmap_data)

DRUG_ABBREV = {
    "TETRACYCLINE":         "TET",
    "AMINOGLYCOSIDE":       "AMG",
    "SULFONAMIDE":          "SUL",
    "BETA-LACTAM":          "BLA",
    "TRIMETHOPRIM":         "TMP",
    "FOSFOMYCIN":           "FOS",
    "MACROLIDE":            "MAC",
    "PHENICOL":             "PHE",
    "QUINOLONE":            "QNL",
    "QUINOLONE/TRICLOSAN":  "QNL",   # merge with QUINOLONE
    "MULTIDRUG":            "MDR",
    "FOSMIDOMYCIN":         "FOSM",
    "EFFLUX":               "EFX",   # will be excluded for E. coli
}


def clean_st(val):
    """Convert float-ish ST strings (e.g. '19.0') to clean integers."""
    try:
        return str(int(float(val)))
    except (ValueError, TypeError):
        return str(val)


def build_matrix(df, row_col, drug_col, top_n=15, exclude_classes=None):
    """
    Returns (matrix_df, row_labels) where matrix_df has rows = categories,
    cols = drug classes, values = % of that category carrying each class.
    """
    if exclude_classes is None:
        exclude_classes = set()

    # Parse per-isolate drug classes into a set
    def parse_classes(cell):
        if pd.isna(cell) or cell == "" or cell == "-":
            return set()
        return {DRUG_ABBREV.get(x.strip(), x.strip())
                for x in cell.split(";")
                if x.strip() and x.strip() not in exclude_classes
                and DRUG_ABBREV.get(x.strip(), x.strip()) != "EFX"}

    df = df.copy()
    df["_row"] = df[row_col].fillna("Unknown").apply(
        lambda v: clean_st(v) if "st" in row_col.lower() else str(v)
    )
    df["_classes"] = df[drug_col].apply(parse_classes)

    # Top N rows by frequency (Unknown always last)
    counts = Counter(df["_row"])
    ordered = [r for r, _ in counts.most_common() if r != "Unknown"][:top_n]
    if "Unknown" in counts:
        ordered.append("Unknown")

    # Collect all drug class abbreviations (excluding excluded)
    all_classes_set = set()
    for s in df["_classes"]:
        all_classes_set |= s
    # Order by overall prevalence
    class_counts = Counter()
    for s in df["_classes"]:
        for c in s:
            class_counts[c] += 1
    all_classes = [c for c, _ in class_counts.most_common()]

    rows = []
    labels = []
    for cat in ordered:
        sub = df[df["_row"] == cat]
        n = len(sub)
        row = []
        for dc in all_classes:
            pct = 100 * sum(1 for s in sub["_classes"] if dc in s) / n
            row.append(pct)
        rows.append(row)
        labels.append(f"{cat}  (n={n})")

    mat = pd.DataFrame(rows, index=labels, columns=all_classes)
    return mat


def tile_heatmap(mat, title, outpath, figsize=None):
    nrows, ncols = mat.shape
    if figsize is None:
        figsize = (max(6, ncols * 0.85 + 3), max(4, nrows * 0.55 + 2))

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(-0.5, ncols - 0.5)
    ax.set_ylim(-0.5, nrows - 0.5)
    ax.invert_yaxis()

    TILE = 0.88  # fraction of cell filled

    for ri, row_label in enumerate(mat.index):
        for ci, col_label in enumerate(mat.columns):
            pct = mat.iloc[ri, ci]
            if pct == 0:
                face = GREY_0
            else:
                face = np.array(POS_CMAP(pct / 100))

            rect = mpatches.FancyBboxPatch(
                (ci - TILE / 2, ri - TILE / 2),
                TILE, TILE,
                boxstyle="square,pad=0",
                facecolor=face, edgecolor="white", linewidth=0.5,
            )
            ax.add_patch(rect)

            # Text colour: white on dark tiles, dark on light
            lum = 0.299 * face[0] + 0.587 * face[1] + 0.114 * face[2]
            txt_color = "white" if lum < 0.45 else "#333333"
            label = f"{pct:.0f}" if pct >= 1 else (f"{pct:.1f}" if pct > 0 else "")
            if label:
                ax.text(ci, ri, label, ha="center", va="center",
                        fontsize=7.5, color=txt_color, fontweight="normal")

    # Axes labels
    ax.set_xticks(range(ncols))
    ax.set_xticklabels(mat.columns, rotation=45, ha="left", fontsize=9)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.tick_params(axis="x", length=0, pad=3)

    ax.set_yticks(range(nrows))
    ax.set_yticklabels(mat.index, fontsize=9)
    ax.tick_params(axis="y", length=0, pad=3)

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_title(title, fontsize=12, fontweight="bold", pad=18)

    # Colour legend
    legend_items = [
        mpatches.Patch(facecolor=GREY_0,
                       label="0 %", linewidth=0),
        mpatches.Patch(facecolor=CREAM[:3],
                       label="1–25 %", linewidth=0),
        mpatches.Patch(facecolor=ORANGE[:3],
                       label="50 %", linewidth=0),
        mpatches.Patch(facecolor=DARK_RED[:3],
                       label="75–80 %", linewidth=0),
        mpatches.Patch(facecolor=PURPLE[:3],
                       label="100 %", linewidth=0),
    ]
    ax.legend(handles=legend_items, loc="lower right",
              bbox_to_anchor=(1.0, -0.02),
              ncol=len(legend_items), frameon=False,
              fontsize=7.5, handlelength=1.2, handletextpad=0.4,
              columnspacing=0.8,
              bbox_transform=fig.transFigure)

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
SAL_TSV = "/Users/lshef4/Documents/dropbox/enteric_typer_project/results_gambia_salmonella/salmonella_typer_results.tsv"
ECO_TSV = "/Users/lshef4/Documents/dropbox/enteric_typer_project/results/ecoli_typer_results.tsv"
OUT     = "/Users/lshef4/Documents/dropbox/enteric_typer_project"

sal = pd.read_csv(SAL_TSV, sep="\t")
eco = pd.read_csv(ECO_TSV, sep="\t")

# ---------------------------------------------------------------------------
# Plot 1 — Salmonella: AMR drug class % by MLST ST
# ---------------------------------------------------------------------------
mat1 = build_matrix(sal, "mlst_st", "amrfinder_drug_classes", top_n=12)
tile_heatmap(
    mat1,
    "Salmonella enterica — AMR drug class prevalence by sequence type (ST)",
    f"{OUT}/salmonella_amr_by_st.png",
)

# ---------------------------------------------------------------------------
# Plot 2 — Salmonella: AMR drug class % by serovar
# ---------------------------------------------------------------------------
mat2 = build_matrix(sal, "sistr_serovar", "amrfinder_drug_classes", top_n=15)
tile_heatmap(
    mat2,
    "Salmonella enterica — AMR drug class prevalence by serovar",
    f"{OUT}/salmonella_amr_by_serovar.png",
)

# ---------------------------------------------------------------------------
# Plot 3 — E. coli: AMR drug class % by MLST ST  (exclude EFFLUX)
# ---------------------------------------------------------------------------
mat3 = build_matrix(eco, "mlst_st", "amrfinder_drug_classes", top_n=12,
                    exclude_classes={"EFFLUX"})
tile_heatmap(
    mat3,
    "Escherichia coli — AMR drug class prevalence by sequence type (ST)",
    f"{OUT}/ecoli_amr_by_st.png",
)

# ---------------------------------------------------------------------------
# Plot 4 — E. coli: AMR drug class % by phylogroup  (exclude EFFLUX)
# ---------------------------------------------------------------------------
mat4 = build_matrix(eco, "kleborate_phylogroup", "amrfinder_drug_classes", top_n=10,
                    exclude_classes={"EFFLUX"})
tile_heatmap(
    mat4,
    "Escherichia coli — AMR drug class prevalence by Clermont phylogroup",
    f"{OUT}/ecoli_amr_by_phylogroup.png",
)
