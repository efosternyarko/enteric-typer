#!/usr/bin/env python3
"""
plot_summary.py — Publication-ready summary figures for enteric-typer results.

Reads either:
  (a) enteric-typer output TSV  --format enteric-typer
  (b) TheiaProk / duplex CSV    --format theiaprok

Produces (PDF + PNG at 300 dpi):
  {prefix}_fig1_population_summary   — 4-panel: ST · serotype · AMR prevalence · MDR burden
  {prefix}_fig2_resistome_heatmap    — sample × drug-class binary matrix ordered by ST
  {prefix}_fig3_amr_genes            — top-25 AMR genes/determinants
  {prefix}_fig4_plasmid_replicons    — top-15 plasmid replicons
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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
    "A": "#4e79a7", "B1": "#59a14f", "B2": "#f28e2b",
    "C": "#76b7b2", "D": "#e15759", "E": "#b07aa1",
    "F": "#ff9da7", "Unknown": "#bab0ac", "Other": "#d3d3d3",
}

# ST → phylogroup (curated from published data)
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

# Clinically relevant AMR classes to show (exclude EFFLUX/MULTIDRUG from heatmap)
CLINICAL_CLASSES = [
    "BETA-LACTAM", "AMINOGLYCOSIDE", "SULFONAMIDE", "TRIMETHOPRIM",
    "TETRACYCLINE", "FOSFOMYCIN", "COLISTIN", "QUINOLONE", "PHENICOL",
    "MACROLIDE", "STREPTOTHRICIN", "FOSMIDOMYCIN", "NITROFURAN",
]

CLASS_LABEL = {
    "AMINOGLYCOSIDE":  "Aminoglycoside",
    "BETA-LACTAM":     "Beta-lactam",
    "COLISTIN":        "Colistin",
    "FOSFOMYCIN":      "Fosfomycin",
    "FOSMIDOMYCIN":    "Fosmidomycin",
    "MACROLIDE":       "Macrolide",
    "NITROFURAN":      "Nitrofuran",
    "PHENICOL":        "Phenicol",
    "QUINOLONE":       "Quinolone",
    "STREPTOTHRICIN":  "Streptothricin",
    "SULFONAMIDE":     "Sulfonamide",
    "TETRACYCLINE":    "Tetracycline",
    "TRIMETHOPRIM":    "Trimethoprim",
}


# ── Data loading ──────────────────────────────────────────────────────────────

def load_data(path: str, fmt: str) -> pd.DataFrame:
    if fmt == "auto":
        with open(path) as fh:
            hdr = fh.readline()
        fmt = "theiaprok" if ("entity:kleb_ecoli_duplex_id" in hdr
                              or "gambit_predicted_taxon" in hdr) else "enteric-typer"

    df = pd.read_csv(path, low_memory=False, sep="\t" if path.endswith(".tsv") else ",")

    if fmt == "theiaprok":
        if "gambit_predicted_taxon" in df.columns:
            df = df[df["gambit_predicted_taxon"].str.contains("Escherichia coli", na=False)].copy()
        df = df.rename(columns={
            "entity:kleb_ecoli_duplex_id":  "sample",
            "ts_mlst_predicted_st":         "mlst_st",
            "ectyper_predicted_serotype":   "ectyper_serotype",
            "amrfinderplus_amr_classes":    "amr_classes",
            "amrfinderplus_amr_core_genes": "amr_genes",
            "plasmidfinder_plasmids":       "replicons",
        })
    else:
        df = df.rename(columns={
            "amrfinder_drug_classes":    "amr_classes",
            "amrfinder_genes":           "amr_genes",
            "plasmidfinder_replicons":   "replicons",
            "ectyper_serotype":          "ectyper_serotype",
        })

    return df.reset_index(drop=True)


# ── Parsing helpers ───────────────────────────────────────────────────────────

def clean_st(st) -> str:
    s = str(st).strip() if pd.notna(st) else ""
    if s in ("", "-", "NA", "nan", "No ST predicted"):
        return "Unknown"
    return s if s.startswith("ST") else f"ST{s}"


def parse_classes(val) -> set[str]:
    if pd.isna(val) or not str(val).strip():
        return set()
    out: set[str] = set()
    for tok in str(val).split(","):
        for part in tok.strip().split("/"):
            p = part.strip()
            if p:
                out.add(p)
    return out


def parse_list(val, sep=",") -> list[str]:
    if pd.isna(val) or not str(val).strip():
        return []
    return [x.strip() for x in str(val).split(sep) if x.strip()]


def get_phylogroup(st: str) -> str:
    return ST_PHYLOGROUP.get(st, "Unknown")


# ── Figure 1: 4-panel population summary ─────────────────────────────────────

def fig_population_summary(df: pd.DataFrame, outdir: Path, prefix: str) -> None:
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.40)
    _panel_st(df,   fig.add_subplot(gs[0, 0]))
    _panel_sero(df, fig.add_subplot(gs[0, 1]))
    _panel_amr(df,  fig.add_subplot(gs[1, 0]))
    _panel_mdr(df,  fig.add_subplot(gs[1, 1]))
    fig.suptitle(f"E. coli genomic surveillance summary  (n = {len(df)} isolates)",
                 fontsize=12, fontweight="bold", y=1.01)
    _save(fig, outdir, f"{prefix}_fig1_population_summary")


def _panel_st(df: pd.DataFrame, ax: plt.Axes, top_n: int = 15) -> None:
    sts = df["mlst_st"].apply(clean_st) if "mlst_st" in df.columns else pd.Series(["Unknown"] * len(df))
    ctr = Counter(sts)
    unk_n = ctr.pop("Unknown", 0)

    top   = ctr.most_common(top_n)
    other_n = sum(v for k, v in ctr.items() if k not in dict(top))

    labels, values = zip(*top) if top else ([], [])
    labels, values = list(labels), list(values)
    if other_n:  labels.append("Other");    values.append(other_n)
    if unk_n:    labels.append("Unknown");  values.append(unk_n)

    def _col(lbl):
        if lbl in ("Other", "Unknown"):
            return PHYLOGROUP_COLORS[lbl]
        return PHYLOGROUP_COLORS.get(get_phylogroup(lbl), PHYLOGROUP_COLORS["Unknown"])

    colors = [_col(l) for l in labels]
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.5, height=0.72)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Number of isolates")
    ax.set_title("A   MLST sequence types")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=5))

    # Phylogroup legend
    seen: dict[str, str] = {}
    for lbl in labels:
        pg = get_phylogroup(lbl) if lbl not in ("Other", "Unknown") else lbl
        seen.setdefault(pg, _col(lbl))
    patches = [mpatches.Patch(facecolor=c, label=pg) for pg, c in seen.items()]
    ax.legend(handles=patches, title="Phylogroup", fontsize=7, title_fontsize=7,
              loc="upper right", handlelength=1)


def _panel_sero(df: pd.DataFrame, ax: plt.Axes, top_n: int = 15) -> None:
    col = next((c for c in ["ectyper_serotype"] if c in df.columns), None)
    if col is None:
        ax.text(0.5, 0.5, "No serotype data", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_title("B   Serotypes"); return

    sero = df[col].fillna("Unknown").replace({"NA": "Unknown", "-:-": "Unknown", "": "Unknown"})
    ctr  = Counter(sero)
    unk  = ctr.pop("Unknown", 0)
    top  = ctr.most_common(top_n)
    other_n = sum(v for k, v in ctr.items() if k not in dict(top))

    labels = [k for k, _ in top]
    values = [v for _, v in top]
    if other_n + unk:
        labels.append("Other / Unknown")
        values.append(other_n + unk)

    palette = sns.color_palette("husl", len(labels))
    y = np.arange(len(labels))
    ax.barh(y, values, color=palette, edgecolor="white", linewidth=0.5, height=0.72)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Number of isolates")
    ax.set_title("B   O:H serotypes (ECTyper)")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=5))


def _panel_amr(df: pd.DataFrame, ax: plt.Axes) -> None:
    col = next((c for c in ["amrfinder_drug_classes", "amr_classes",
                            "amrfinderplus_amr_classes"] if c in df.columns), None)
    if col is None:
        ax.text(0.5, 0.5, "No AMR data", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_title("C   AMR prevalence"); return

    n = len(df)
    counts: dict[str, int] = defaultdict(int)
    for val in df[col]:
        for cls in parse_classes(val):
            if cls in CLINICAL_CLASSES:
                counts[cls] += 1

    items  = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    labels = [CLASS_LABEL.get(k, k) for k, _ in items]
    pcts   = [100 * v / n for _, v in items]

    y = np.arange(len(labels))
    ax.hlines(y, 0, pcts, colors="#d0d0d0", linewidth=1.5, zorder=1)
    dot_colors = ["#e15759" if p >= 70 else "#f28e2b" if p >= 35 else "#4e79a7" for p in pcts]
    ax.scatter(pcts, y, color=dot_colors, s=45, zorder=2)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (% of isolates)")
    ax.set_title("C   AMR drug class prevalence")
    ax.set_xlim(0, 108)
    ax.axvline(80, color="#cccccc", linestyle="--", linewidth=0.8)
    ax.text(81, -0.6, "80%", fontsize=6.5, color="#888888")
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))


def _panel_mdr(df: pd.DataFrame, ax: plt.Axes) -> None:
    col = next((c for c in ["amrfinder_drug_classes", "amr_classes",
                            "amrfinderplus_amr_classes"] if c in df.columns), None)
    if col is None:
        ax.set_title("D   MDR burden"); return

    burdens = [len(parse_classes(v) & set(CLINICAL_CLASSES)) for v in df[col]]
    ctr = Counter(burdens)
    xs  = sorted(ctr)
    ys  = [ctr[x] for x in xs]

    colors = ["#e15759" if x >= 3 else "#4e79a7" for x in xs]
    ax.bar(xs, ys, color=colors, edgecolor="white", linewidth=0.5, width=0.72)
    ax.set_xlabel("No. resistant drug classes")
    ax.set_ylabel("No. isolates")
    ax.set_title("D   MDR burden per isolate")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    mean_v = float(np.mean(burdens))
    ax.axvline(mean_v, color="black", linestyle="--", linewidth=1)
    ax.text(mean_v + 0.08, max(ys) * 0.96, f"mean = {mean_v:.1f}", fontsize=7.5, va="top")

    ax.legend(handles=[
        mpatches.Patch(color="#e15759", label="MDR  (≥ 3 classes)"),
        mpatches.Patch(color="#4e79a7", label="Non-MDR"),
    ], fontsize=7.5)


# ── Figure 2: Resistome heatmap ───────────────────────────────────────────────

def fig_resistome_heatmap(df: pd.DataFrame, outdir: Path, prefix: str) -> None:
    col = next((c for c in ["amrfinder_drug_classes", "amr_classes",
                            "amrfinderplus_amr_classes"] if c in df.columns), None)
    if col is None:
        print("WARNING: no AMR class column — skipping heatmap", file=sys.stderr)
        return

    # Order classes by prevalence
    classes_ordered = sorted(
        CLINICAL_CLASSES,
        key=lambda c: -sum(1 for v in df[col] if c in parse_classes(v))
    )
    display_cls = [CLASS_LABEL.get(c, c) for c in classes_ordered]

    # Sort rows by ST (most-common ST first, Unknown last)
    df2 = df.copy()
    df2["_st"] = df2["mlst_st"].apply(clean_st) if "mlst_st" in df2.columns else "Unknown"
    st_rank = {st: i for i, (st, _) in
               enumerate(Counter(df2["_st"]).most_common())}
    st_rank["Unknown"] = 9999
    df2 = df2.sort_values("_st", key=lambda s: s.map(st_rank))

    # Binary matrix
    mat = np.array([
        [1 if cls in parse_classes(row[col]) else 0 for cls in classes_ordered]
        for _, row in df2.iterrows()
    ], dtype=np.uint8)

    n_rows  = len(df2)
    row_h   = max(0.055, min(0.13, 12.0 / n_rows))
    fig_h   = min(n_rows * row_h + 2.0, 22)
    fig_w   = len(classes_ordered) * 0.60 + 2.8

    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h),
                             gridspec_kw={"width_ratios": [len(classes_ordered), 1],
                                          "wspace": 0.02})
    ax_heat, ax_st = axes

    # Heatmap
    cmap = matplotlib.colors.ListedColormap(["#f5f5f5", "#c0392b"])
    ax_heat.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                   interpolation="nearest")
    ax_heat.set_xticks(np.arange(len(display_cls)))
    ax_heat.set_xticklabels(display_cls, rotation=45, ha="right", fontsize=8)
    ax_heat.set_yticks([])
    ax_heat.set_ylabel(f"Isolates  (n = {n_rows})", fontsize=9)
    ax_heat.set_title("Resistome — AMR drug class presence/absence  (ordered by ST)",
                      fontsize=10, fontweight="bold", pad=8)

    # ST colour bar — unique color per ST so no two STs share a color
    _unique_sts = [s for s in dict.fromkeys(df2["_st"]) if s != "Unknown"]
    _rgb = sns.color_palette("husl", max(len(_unique_sts), 1))
    st_palette: dict[str, tuple] = {st: _rgb[i] for i, st in enumerate(_unique_sts)}
    st_palette["Unknown"] = (0.85, 0.85, 0.85)

    for i, st in enumerate(df2["_st"]):
        ax_st.barh(i, 1, color=st_palette[st], height=1, linewidth=0, align="edge")
    ax_st.set_xlim(0, 1); ax_st.set_ylim(0, n_rows)
    ax_st.invert_yaxis()
    ax_st.set_yticks([]); ax_st.set_xticks([])
    ax_st.set_title("ST", fontsize=8, fontweight="bold")
    for spine in ax_st.spines.values():
        spine.set_visible(False)

    # ST legend (top 12 most common)
    top_sts = [st for st, _ in Counter(df2["_st"].tolist()).most_common(12)]
    ax_st.legend(
        handles=[mpatches.Patch(facecolor=st_palette[s], label=s) for s in top_sts],
        fontsize=6.5, loc="lower right", frameon=True, framealpha=0.9,
        edgecolor="none", bbox_to_anchor=(2.0, 0),
    )

    _save(fig, outdir, f"{prefix}_fig2_resistome_heatmap")


# ── Figure 3: AMR genes ───────────────────────────────────────────────────────

def fig_amr_genes(df: pd.DataFrame, outdir: Path, prefix: str, top_n: int = 25) -> None:
    # Prefer amrfinder_acquired_genes (intrinsic already excluded upstream by AMRrules);
    # fall back to amr_genes / theiaprok column if running on pre-AMRrules data
    col = next((c for c in ["amrfinder_acquired_genes", "amr_genes",
                             "amrfinderplus_amr_core_genes"] if c in df.columns), None)
    if col is None:
        print("WARNING: no AMR gene column — skipping gene plot", file=sys.stderr)
        return

    gene_ctr: Counter = Counter()
    for val in df[col]:
        for g in parse_list(val):
            gene_ctr[g] += 1

    top    = gene_ctr.most_common(top_n)
    labels = [g for g, _ in top]
    pcts   = [100 * v / len(df) for _, v in top]

    # Acquired gene vs chromosomal point mutation
    # Heuristic: contains underscore with capital letter suffix = point mutation
    import re
    def is_mutation(g: str) -> bool:
        return bool(re.search(r'_[A-Z]\d+[A-Z]$', g))

    colors = ["#aec6cf" if is_mutation(g) else "#4e79a7" for g in labels]

    fig, ax = plt.subplots(figsize=(6.5, top_n * 0.30 + 1.4))
    y = np.arange(len(labels))
    ax.barh(y, pcts, color=colors, edgecolor="white", linewidth=0.3, height=0.76)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5, fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (% of isolates)")
    ax.set_title("AMR genes and resistance determinants", fontweight="bold")
    ax.set_xlim(0, 108)
    ax.axvline(80, color="#cccccc", linestyle="--", linewidth=0.8)
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    ax.legend(handles=[
        mpatches.Patch(color="#4e79a7", label="Acquired gene"),
        mpatches.Patch(color="#aec6cf", label="Point mutation"),
    ], fontsize=7.5)

    _save(fig, outdir, f"{prefix}_fig3_amr_genes")


# ── Figure 4: Plasmid replicons ───────────────────────────────────────────────

def fig_plasmid_replicons(df: pd.DataFrame, outdir: Path, prefix: str, top_n: int = 15) -> None:
    col = next((c for c in ["replicons", "plasmidfinder_plasmids"] if c in df.columns), None)
    if col is None:
        print("WARNING: no plasmid column — skipping replicon plot", file=sys.stderr)
        return

    rep_ctr: Counter = Counter()
    for val in df[col]:
        if pd.notna(val) and "No plasmids" not in str(val):
            for r in parse_list(val):
                rep_ctr[r] += 1

    top    = rep_ctr.most_common(top_n)
    labels = [r for r, _ in top]
    pcts   = [100 * v / len(df) for _, v in top]

    def _rcol(r: str) -> str:
        if r.startswith("IncF"):  return "#f28e2b"
        if r.startswith("Col"):   return "#bab0ac"
        if r.startswith("Inc"):   return "#4e79a7"
        return "#76b7b2"

    colors = [_rcol(l) for l in labels]

    fig, ax = plt.subplots(figsize=(6.5, top_n * 0.34 + 1.4))
    y = np.arange(len(labels))
    ax.barh(y, pcts, color=colors, edgecolor="white", linewidth=0.3, height=0.76)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (% of isolates)")
    ax.set_title("Plasmid replicon types (PlasmidFinder)", fontweight="bold")
    ax.set_xlim(0, 108)
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    ax.legend(handles=[
        mpatches.Patch(color="#f28e2b", label="IncF family"),
        mpatches.Patch(color="#4e79a7", label="Other Inc"),
        mpatches.Patch(color="#bab0ac", label="Col plasmid"),
    ], fontsize=7.5)

    _save(fig, outdir, f"{prefix}_fig4_plasmid_replicons")


# ── Save helper ───────────────────────────────────────────────────────────────

def _save(fig: plt.Figure, outdir: Path, stem: str) -> None:
    for ext in ("pdf", "png"):
        p = outdir / f"{stem}.{ext}"
        fig.savefig(p)
        print(f"  Saved: {p}", file=sys.stderr)
    plt.close(fig)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(description="enteric-typer summary plots")
    p.add_argument("--input",   "-i", required=True)
    p.add_argument("--format",  "-f", default="auto",
                   choices=["auto", "enteric-typer", "theiaprok"])
    p.add_argument("--outdir",  "-o", default=".")
    p.add_argument("--prefix",  "-p", default="enteric_typer")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data…", file=sys.stderr)
    df = load_data(args.input, args.format)
    print(f"  {len(df)} isolates", file=sys.stderr)
    if df.empty:
        sys.exit("ERROR: no data after loading — check --format")

    print("Generating figures…", file=sys.stderr)
    fig_population_summary(df, outdir, args.prefix)
    # Fig 2 (tree-annotated resistome heatmap) is generated by plot_tree_annotation.py
    fig_amr_genes(df, outdir, args.prefix)
    fig_plasmid_replicons(df, outdir, args.prefix)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
