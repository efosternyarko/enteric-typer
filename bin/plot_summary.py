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
    # Accept both ';' (enteric-typer) and ',' (theiaprok / legacy) separators
    raw = str(val).replace(";", ",")
    for tok in raw.split(","):
        for part in tok.strip().split("/"):
            p = part.strip()
            if p:
                out.add(p)
    return out


_SENTINEL = {"NA", "nan", "None", "-", ""}

def parse_list(val, sep=None) -> list[str]:
    """Split a semicolon- or comma-separated gene/replicon list."""
    if pd.isna(val):
        return []
    s = str(val).strip()
    if not s or s in _SENTINEL:
        return []
    if sep is None:
        sep = ";" if ";" in s else ","
    return [x.strip() for x in s.split(sep) if x.strip() and x.strip() not in _SENTINEL]


def get_phylogroup(st: str) -> str:
    return ST_PHYLOGROUP.get(st, "Unknown")


# ── Figure 1: 4-panel population summary ─────────────────────────────────────

def fig_population_summary(df: pd.DataFrame, outdir: Path, prefix: str) -> None:
    # Extra width (16 → 18) gives each panel room for its outside legend
    fig = plt.figure(figsize=(18, 11))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.50, wspace=0.55)
    _panel_st(df,   fig.add_subplot(gs[0, 0]))
    _panel_sero(df, fig.add_subplot(gs[0, 1]))
    _panel_amr(df,  fig.add_subplot(gs[1, 0]))
    _panel_mdr(df,  fig.add_subplot(gs[1, 1]))
    sp_label = ("Salmonella enterica" if "salmonella" in prefix.lower()
                else "E. coli")
    fig.suptitle(f"{sp_label} genomic surveillance summary  (n = {len(df)} isolates)",
                 fontsize=12, fontweight="bold", y=1.01)
    _save(fig, outdir, f"{prefix}_fig1_population_summary")


def _panel_st(df: pd.DataFrame, ax: plt.Axes, top_n: int = 15) -> None:
    sts = df["mlst_st"].apply(clean_st) if "mlst_st" in df.columns else pd.Series(
        ["Unknown"] * len(df), index=df.index)

    # Use real Kleborate phylogroup if present; fall back to static ST lookup
    has_real_pg = (
        "kleborate_phylogroup" in df.columns and
        df["kleborate_phylogroup"].notna().any() and
        df["kleborate_phylogroup"].astype(str).str.strip().ne("NA").any()
    )

    # Build ST → dominant phylogroup dict
    if has_real_pg:
        pg_col = df["kleborate_phylogroup"].astype(str).str.strip()
        pg_col = pg_col.where(~pg_col.isin({"NA", "nan", "None", "-", ""}), "Unknown")
        st_pg_ctr: dict[str, Counter] = {}
        for st, pg in zip(sts, pg_col):
            st_pg_ctr.setdefault(st, Counter())[pg] += 1
        st_pg_map = {st: ctr.most_common(1)[0][0] for st, ctr in st_pg_ctr.items()}
        title_suffix = " (Clermont)"
    else:
        st_pg_map = {}
        title_suffix = " (ST-inferred)"

    def _pg_for(lbl: str) -> str:
        if lbl in ("Other", "Unknown"):
            return lbl
        if has_real_pg:
            return st_pg_map.get(lbl, "Unknown")
        return get_phylogroup(lbl)

    def _col(lbl: str) -> str:
        return PHYLOGROUP_COLORS.get(_pg_for(lbl), PHYLOGROUP_COLORS["Unknown"])

    ctr = Counter(sts)
    unk_n = ctr.pop("Unknown", 0)
    top   = ctr.most_common(top_n)
    other_n = sum(v for k, v in ctr.items() if k not in dict(top))

    labels, values = zip(*top) if top else ([], [])
    labels, values = list(labels), list(values)
    if other_n:  labels.append("Other");    values.append(other_n)
    if unk_n:    labels.append("Unknown");  values.append(unk_n)

    colors = [_col(l) for l in labels]
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.6, height=0.78)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Number of isolates", fontsize=8.5)
    ax.set_title(f"A. MLST sequence types{title_suffix}", fontweight="bold")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=5))
    ax.tick_params(axis="both", length=3)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.yaxis.set_ticks_position("none")

    # Phylogroup legend — outside, right
    seen: dict[str, str] = {}
    for lbl in labels:
        pg = _pg_for(lbl)
        seen.setdefault(pg, _col(lbl))
    patches = [mpatches.Patch(facecolor=c, label=pg, linewidth=0) for pg, c in seen.items()]
    ax.legend(handles=patches, title="Phylogroup", fontsize=7, title_fontsize=7.5,
              bbox_to_anchor=(1.02, 1), loc="upper left",
              frameon=False, handlelength=1.2, handleheight=1.2)


# K-locus group colours (G1–G4 + Unknown) — consistent across all figures
KGROUP_COLORS = {
    "G1": "#4e79a7", "G2": "#f28e2b",
    "G3": "#59a14f", "G4": "#e15759",
    "Unknown": "#bab0ac",
}
KGROUP_ORDER = ["G1", "G2", "G3", "G4", "Unknown"]


def _panel_sero(df: pd.DataFrame, ax: plt.Axes, top_n: int = 15) -> None:
    """
    Horizontal bar chart of top serotypes / serovars.

    E. coli: bars are stacked and filled by K-locus group (G1–G4) when
             the k_group column is present, giving a capsule-type context.
    Salmonella: bars filled by MLST ST complex when available; plain otherwise.
    """
    _sero_cols = [("ectyper_serotype", "B. O:H serotypes (ECTyper)",  "ecoli"),
                  ("sistr_serovar",    "B. Serovars (SISTR)",         "salmonella")]
    match = next(((c, t, sp) for c, t, sp in _sero_cols if c in df.columns), None)
    if match is None:
        ax.text(0.5, 0.5, "No serotype data", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_title("B   Serotypes"); return

    col, title, species = match

    sero = df[col].fillna("Unknown").replace(
        {"NA": "Unknown", "-:-": "Unknown", "": "Unknown"})

    ctr = Counter(sero)
    unk = ctr.pop("Unknown", 0)
    top = ctr.most_common(top_n)
    other_n = sum(v for k, v in ctr.items() if k not in dict(top))

    labels = [k for k, _ in top]
    if other_n + unk:
        labels.append("Other / Unknown")

    # ── Decide fill strategy ──────────────────────────────────────────────────
    use_klocus    = (species == "ecoli"      and "k_locus"        in df.columns)
    use_stcomplex = (species == "salmonella" and "mlst_st_complex" in df.columns)

    # ── Build per-serotype fill counts ────────────────────────────────────────
    if use_klocus:
        raw_klocus = df["k_locus"].fillna("Unknown").replace(
            {"NA": "Unknown", "": "Unknown", "-": "Unknown"})
        # Build locus → group mapping for colour assignment
        locus_grp: dict[str, str] = {}
        if "k_group" in df.columns:
            for locus, grp in zip(raw_klocus,
                                  df["k_group"].fillna("Unknown").replace(
                                      {"NA": "Unknown", "": "Unknown"})):
                locus_grp.setdefault(str(locus), str(grp).strip())
        locus_freq = Counter(raw_klocus)
        top_loci   = [l for l, _ in locus_freq.most_common(12) if l != "Unknown"]
        has_other  = any(l not in top_loci and l != "Unknown" for l in locus_freq)
        fill_vals  = top_loci
        if has_other:              fill_vals = fill_vals + ["Other"]
        if "Unknown" in locus_freq: fill_vals = fill_vals + ["Unknown"]
        fill_col   = raw_klocus.apply(
            lambda l: l if l in top_loci
                      else ("Unknown" if l == "Unknown" else "Other"))

        def _locus_color(locus: str) -> str:
            if locus in ("Other", "Unknown"):
                return KGROUP_COLORS.get(locus, "#bab0ac")
            return KGROUP_COLORS.get(locus_grp.get(locus, "Unknown"), "#bab0ac")

        fill_label = {l: l for l in fill_vals}
        fill_color = {l: _locus_color(l) for l in fill_vals}
        fill_title = f"{title}  ·  K-locus type (coloured by group)"

    elif use_stcomplex:
        raw_stc    = df["mlst_st_complex"].fillna("Unknown").replace(
            {"NA": "Unknown", "": "Unknown"})
        fill_col   = raw_stc
        top_stc    = [s for s, _ in Counter(raw_stc).most_common(8)]
        fill_vals  = top_stc + (["Other"] if len(Counter(raw_stc)) > 8 else [])
        stc_pal    = sns.color_palette("husl", len(top_stc))
        fill_color = {s: stc_pal[i] for i, s in enumerate(top_stc)}
        fill_color.update({"Other": "#bab0ac", "Unknown": "#d3d3d3"})
        fill_label = {s: s for s in fill_vals}
        fill_title = f"{title}  ·  filled by ST complex"

    else:
        # Plain bars: one colour per serotype
        palette = sns.color_palette("husl", len(labels))
        y = np.arange(len(labels))
        vals = [dict(ctr).get(l, other_n + unk) if l != "Other / Unknown"
                else other_n + unk for l in labels]
        ax.barh(y, vals, color=palette, edgecolor="white", linewidth=0.5, height=0.72)
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
        ax.invert_yaxis()
        ax.set_xlabel("Number of isolates")
        ax.set_title(title)
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=5))
        return

    # ── Pre-compute per-serotype fill breakdown ───────────────────────────────
    sero_fill: dict[str, Counter] = {}
    sero_series = sero.reset_index(drop=True)
    fill_series = fill_col.reset_index(drop=True)

    for idx in range(len(df)):
        sv = str(sero_series.iloc[idx]).strip()
        fv = str(fill_series.iloc[idx]).strip()

        # Normalise fill value to known categories
        valid = set(fill_vals)
        if fv not in valid:
            fv = "Other" if use_stcomplex else "Unknown"

        # Bin serotype into labelled buckets
        if sv in labels:
            key = sv
        else:
            key = "Other / Unknown"

        sero_fill.setdefault(key, Counter())[fv] += 1

    # ── Draw stacked bars ─────────────────────────────────────────────────────
    BAR_H = 0.78
    y_pos = np.arange(len(labels))
    for y_idx, sero_label in enumerate(labels):
        fill_ctr = sero_fill.get(sero_label, Counter())
        x = 0
        for fv in fill_vals:
            cnt = fill_ctr.get(fv, 0)
            if cnt:
                ax.barh(y_idx, cnt, left=x, height=BAR_H,
                        color=fill_color.get(fv, "#bab0ac"),
                        edgecolor="white", linewidth=0.6, zorder=2)
                x += cnt

    ax.set_yticks(y_pos); ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Number of isolates", fontsize=8.5)
    ax.set_title(fill_title, fontsize=8.5, fontweight="bold")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=5))
    ax.tick_params(axis="both", length=3)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.yaxis.set_ticks_position("none")

    patches = [mpatches.Patch(facecolor=fill_color.get(fv, "#bab0ac"),
                               label=fill_label.get(fv, fv), linewidth=0)
               for fv in fill_vals if fv in fill_color]
    ax.legend(handles=patches, fontsize=7, ncol=1,
              bbox_to_anchor=(1.02, 1), loc="upper left",
              frameon=False, handlelength=1.2, handleheight=1.2,
              title=("K-locus type" if use_klocus else "ST complex"),
              title_fontsize=7.5)


def _parse_gene_classes(val) -> dict[str, str]:
    """Parse 'gene=CLASS;gene2=CLASS2' → {gene: class}."""
    result: dict[str, str] = {}
    for pair in parse_list(val):
        if "=" in pair:
            gene, cls = pair.split("=", 1)
            result[gene.strip()] = cls.strip().upper()
    return result


# Bold, saturated qualitative palette — modelled on the figure style.
# 12 distinct hues before cycling; "Other" always gets the grey slot.
_BOLD_PALETTE = [
    "#e63946",  # vivid red
    "#457b9d",  # steel blue
    "#2a9d8f",  # teal
    "#e9c46a",  # amber
    "#264653",  # dark slate
    "#f4a261",  # sandy orange
    "#6a4c93",  # purple
    "#1982c4",  # sky blue
    "#8ac926",  # lime green
    "#ff595e",  # coral
    "#6a994e",  # olive green
    "#ffca3a",  # yellow
]
_OTHER_COLOR = "#d0d0d0"


def _make_gene_palette(global_gene_ctr: Counter) -> dict[str, str]:
    """Assign a consistent bold colour to each gene by global frequency rank."""
    top_genes = [g for g, _ in global_gene_ctr.most_common(len(_BOLD_PALETTE))]
    return {g: _BOLD_PALETTE[i] for i, g in enumerate(top_genes)}


def _panel_amr(df: pd.DataFrame, ax: plt.Axes) -> None:
    """
    Stacked horizontal bar chart styled after the reference figure:
      • One thick bar per clinical drug class, sorted by total prevalence
      • Each bar filled with clean colour-coded gene segments (no text inside)
      • Bold distinct palette; legend to the right of/below the chart
      • No gridlines, minimal spines
    Falls back to simple bars if the gene_classes column is absent.
    """
    gc_col  = "amrfinder_gene_classes"
    cls_col = next((c for c in ["amrfinder_drug_classes", "amr_classes",
                                "amrfinderplus_amr_classes"] if c in df.columns), None)

    if gc_col not in df.columns and cls_col is None:
        ax.text(0.5, 0.5, "No AMR data", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_title("C   AMR prevalence"); return

    n = len(df)

    # ── Build {drug_class: Counter({gene: n_samples})} ────────────────────────
    class_gene: dict[str, Counter] = defaultdict(Counter)

    if gc_col in df.columns:
        for val in df[gc_col]:
            for gene, cls in _parse_gene_classes(val).items():
                if cls in CLINICAL_CLASSES:
                    class_gene[cls][gene] += 1
    else:
        for val in df[cls_col]:
            for cls in parse_classes(val):
                if cls in CLINICAL_CLASSES:
                    class_gene[cls]["–"] += 1

    if not class_gene:
        ax.text(0.5, 0.5, "No clinical AMR detected", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_title("C   AMR drug class prevalence"); return

    # ── Colour palette: rank genes globally by total frequency ────────────────
    global_ctr: Counter = Counter()
    for ctr in class_gene.values():
        global_ctr.update(ctr)
    gene_color = _make_gene_palette(global_ctr)

    # ── Sort classes by total sample prevalence (descending) ──────────────────
    class_prev = {cls: sum(1 for val in df.get(gc_col, df.get(cls_col, []))
                           for _ in [1] if cls in str(val))
                  for cls in class_gene}
    # Simpler: use gene hit count as proxy for prevalence ordering
    class_prev = {cls: len({g for val in df.get(gc_col, pd.Series([])).dropna()
                             for g, c in _parse_gene_classes(val).items() if c == cls})
                  for cls in class_gene}
    # Most reliable: count samples that have ≥1 gene in this class
    class_prev2: dict[str, int] = defaultdict(int)
    if gc_col in df.columns:
        for val in df[gc_col]:
            seen_cls: set[str] = set()
            for gene, cls in _parse_gene_classes(val).items():
                if cls in CLINICAL_CLASSES and cls not in seen_cls:
                    class_prev2[cls] += 1
                    seen_cls.add(cls)
    else:
        for val in df[cls_col]:
            for cls in parse_classes(val):
                if cls in CLINICAL_CLASSES:
                    class_prev2[cls] += 1

    sorted_cls = sorted(class_gene, key=lambda c: class_prev2.get(c, 0), reverse=True)

    # ── Draw ──────────────────────────────────────────────────────────────────
    TOP_N = 6
    legend_entries: dict[str, str] = {}   # {gene: colour}, insertion-ordered

    y_pos = np.arange(len(sorted_cls))
    BAR_H = 0.78   # thick bars like the reference figure

    for y_idx, cls in enumerate(sorted_cls):
        ctr   = class_gene[cls]
        top   = ctr.most_common(TOP_N)
        other = sum(ctr.values()) - sum(c for _, c in top)

        x = 0.0
        for gene, count in top:
            w     = 100 * count / n
            color = gene_color.get(gene, _OTHER_COLOR)
            ax.barh(y_idx, w, left=x, height=BAR_H,
                    color=color, edgecolor="white", linewidth=0.6, zorder=2)
            legend_entries.setdefault(gene, color)
            x += w

        if other > 0:
            w = 100 * other / n
            ax.barh(y_idx, w, left=x, height=BAR_H,
                    color=_OTHER_COLOR, edgecolor="white", linewidth=0.6, zorder=2)
            legend_entries.setdefault("Other", _OTHER_COLOR)

    # ── Axes styling — match the clean figure aesthetic ───────────────────────
    ax.set_yticks(y_pos)
    ax.set_yticklabels([CLASS_LABEL.get(c, c) for c in sorted_cls], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (% of isolates)", fontsize=8.5)
    ax.set_title("C. AMR drug class prevalence", fontweight="bold")
    ax.set_xlim(0, 105)
    ax.xaxis.set_major_locator(plt.MultipleLocator(25))
    ax.tick_params(axis="both", length=3)
    # Remove all spines except bottom
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.set_axisbelow(False)
    ax.yaxis.set_ticks_position("none")

    # ── Legend (outside, below or right) ─────────────────────────────────────
    shown = list(legend_entries.items())
    patches = [mpatches.Patch(facecolor=c, label=g, linewidth=0) for g, c in shown]
    ax.legend(handles=patches, fontsize=7, ncol=1,
              bbox_to_anchor=(1.02, 1), loc="upper left",
              frameon=False, handlelength=1.2, handleheight=1.2,
              title="Gene", title_fontsize=7.5)


def _panel_mdr(df: pd.DataFrame, ax: plt.Axes) -> None:
    col = next((c for c in ["amrfinder_drug_classes", "amr_classes",
                            "amrfinderplus_amr_classes"] if c in df.columns), None)
    if col is None:
        ax.set_title("D   MDR burden"); return

    burdens = [len(parse_classes(v) & set(CLINICAL_CLASSES)) for v in df[col]]
    ctr = Counter(burdens)
    xs  = sorted(ctr)
    ys  = [ctr[x] for x in xs]

    colors = ["#e63946" if x >= 3 else "#457b9d" for x in xs]
    ax.bar(xs, ys, color=colors, edgecolor="white", linewidth=0.6, width=0.76)
    ax.set_xlabel("No. resistant drug classes", fontsize=8.5)
    ax.set_ylabel("No. isolates", fontsize=8.5)
    ax.set_title("D. MDR burden per isolate", fontweight="bold")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.tick_params(axis="both", length=3)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.8)

    mean_v = float(np.mean(burdens))
    ax.axvline(mean_v, color="#333333", linestyle="--", linewidth=1)
    ax.text(mean_v + 0.1, max(ys) * 0.96, f"mean = {mean_v:.1f}", fontsize=7.5, va="top")

    ax.legend(handles=[
        mpatches.Patch(color="#e63946", label="MDR  (≥ 3 classes)", linewidth=0),
        mpatches.Patch(color="#457b9d", label="Non-MDR",             linewidth=0),
    ], fontsize=7.5, frameon=False)


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

    if not gene_ctr:
        ax = plt.subplots(1, 1, figsize=(6.5, 2))[1]
        ax.text(0.5, 0.5, "No AMR genes detected in this dataset",
                ha="center", va="center", transform=ax.transAxes, color="grey", fontsize=10)
        ax.set_title("AMR genes and resistance determinants", fontweight="bold")
        ax.axis("off")
        _save(plt.gcf(), outdir, f"{prefix}_fig3_amr_genes")
        return

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
        s = str(val).strip() if pd.notna(val) else ""
        if s and s not in ("NA", "nan") and "No plasmids" not in s:
            for r in parse_list(s):
                if r not in ("NA", "nan", ""):
                    rep_ctr[r] += 1

    if not rep_ctr:
        ax = plt.subplots(1, 1, figsize=(6.5, 2))[1]
        ax.text(0.5, 0.5, "No plasmid replicons detected in this dataset",
                ha="center", va="center", transform=ax.transAxes, color="grey", fontsize=10)
        ax.set_title("Plasmid replicon types (PlasmidFinder)", fontweight="bold")
        ax.axis("off")
        _save(plt.gcf(), outdir, f"{prefix}_fig4_plasmid_replicons")
        return

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


# ── Figure 5: Virulence profiling ─────────────────────────────────────────────
# E. coli:  3-panel — pathotype distribution · pathovar marker prevalence ·
#           top AMRFinder virulence genes.
# Salmonella: single panel — top abricate VFDB hits (fallback: AMRFinder).

# Kleborate pathovar markers and their display labels
PATHOVAR_MARKERS = [
    ("kleborate_Stx1",     "Stx1 (Shiga toxin 1)"),
    ("kleborate_Stx2",     "Stx2 (Shiga toxin 2)"),
    ("kleborate_eae",      "eae (intimin)"),
    ("kleborate_ipaH",     "ipaH (invasion)"),
    ("kleborate_LT",       "LT (heat-labile toxin)"),
    ("kleborate_ST_toxin", "ST (heat-stable toxin)"),
]

# Pathotype display order and colours (most severe first)
PATHOVAR_ORDER = ["EHEC", "STEC", "EPEC", "EIEC", "ETEC", "-"]
PATHOVAR_LABEL = {
    "EHEC": "EHEC",
    "STEC": "STEC",
    "EPEC": "EPEC",
    "EIEC": "EIEC",
    "ETEC": "ETEC",
    "-":    "ExPEC",
}
PATHOVAR_COLOR = {
    "EHEC": "#d62728",
    "STEC": "#e15759",
    "EPEC": "#f28e2b",
    "EIEC": "#9467bd",
    "ETEC": "#1f77b4",
    "-":    "#bab0ac",
}


def _marker_present(val) -> bool:
    """Kleborate pathovar marker: '+' or gene hit string = present; '-' = absent."""
    s = str(val).strip()
    return s not in ("-", "NA", "nan", "None", "")


def fig_virulence(df: pd.DataFrame, outdir: Path, prefix: str, top_n: int = 25) -> None:
    is_ecoli = "ecoli" in prefix.lower() or "kleborate_pathovar" in df.columns
    if is_ecoli:
        _fig_virulence_ecoli(df, outdir, prefix, top_n)
    else:
        _fig_virulence_salmonella(df, outdir, prefix, top_n)


def _fig_virulence_ecoli(df: pd.DataFrame, outdir: Path, prefix: str, top_n: int) -> None:
    has_pathovar  = "kleborate_pathovar" in df.columns
    has_markers   = any(col in df.columns for col, _ in PATHOVAR_MARKERS)
    has_amr_vir   = "amrfinder_virulence_genes" in df.columns

    if not has_pathovar and not has_amr_vir:
        print("WARNING: no virulence data — skipping virulence plot", file=sys.stderr)
        return

    n_panels = sum([has_pathovar, has_markers, has_amr_vir])
    fig_w    = 5.5 * n_panels
    fig_h    = max(5.0, top_n * 0.28 + 2.0)
    fig, axes = plt.subplots(1, n_panels, figsize=(fig_w, fig_h))
    if n_panels == 1:
        axes = [axes]
    ax_idx = 0
    n = len(df)

    # ── Panel A: Pathotype distribution ──────────────────────────────────────
    if has_pathovar:
        ax = axes[ax_idx]; ax_idx += 1
        raw = df["kleborate_pathovar"].fillna("-").replace({"NA": "-", "nan": "-", "": "-"})
        ctr: Counter = Counter()
        for v in raw:
            # Kleborate can produce compound pathotypes like "STEC/EPEC"; count each
            for pt in str(v).split("/"):
                ctr[pt.strip()] += 1

        # Build ordered list: known pathotypes first, then any extras, then '-'
        known   = [p for p in PATHOVAR_ORDER if p in ctr]
        extras  = [p for p in sorted(ctr) if p not in PATHOVAR_ORDER]
        ordered = known + extras

        lbls   = [PATHOVAR_LABEL.get(p, p) for p in ordered]
        vals   = [ctr[p] for p in ordered]
        colors = [PATHOVAR_COLOR.get(p, "#76b7b2") for p in ordered]

        y = np.arange(len(lbls))
        ax.barh(y, vals, color=colors, edgecolor="white", linewidth=0.3, height=0.72)
        ax.set_yticks(y); ax.set_yticklabels(lbls, fontsize=8.5)
        ax.invert_yaxis()
        ax.set_xlabel("Number of isolates")
        ax.set_title("A. Pathotype (Kleborate)", fontweight="bold")
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=5))

    # ── Panel B: Pathovar marker prevalence ──────────────────────────────────
    if has_markers:
        ax = axes[ax_idx]; ax_idx += 1
        lbls, pcts, cols = [], [], []
        for col, label in PATHOVAR_MARKERS:
            if col not in df.columns:
                continue
            pct = 100 * df[col].apply(_marker_present).sum() / n
            lbls.append(label)
            pcts.append(pct)
            # Match colour to the pathotype each marker defines
            if "Stx" in col:  c = PATHOVAR_COLOR["STEC"]
            elif col == "kleborate_eae": c = PATHOVAR_COLOR["EPEC"]
            elif col == "kleborate_ipaH": c = PATHOVAR_COLOR["EIEC"]
            else:              c = PATHOVAR_COLOR["ETEC"]
            cols.append(c)

        y = np.arange(len(lbls))
        ax.hlines(y, 0, pcts, colors="#d0d0d0", linewidth=1.5, zorder=1)
        ax.scatter(pcts, y, color=cols, s=55, zorder=2)
        ax.set_yticks(y); ax.set_yticklabels(lbls, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Prevalence (% of isolates)")
        ax.set_title("B. Virulence gene markers (InPEC)", fontweight="bold")
        ax.set_xlim(0, 108)
        ax.xaxis.set_major_locator(plt.MultipleLocator(20))

    # ── Panel C: Top AMRFinder virulence genes ────────────────────────────────
    if has_amr_vir:
        ax = axes[ax_idx]
        gene_ctr: Counter = Counter()
        for val in df["amrfinder_virulence_genes"]:
            for g in parse_list(val):
                gene_ctr[g] += 1

        if not gene_ctr:
            ax.text(0.5, 0.5, "No virulence genes\ndetected (AMRFinder)",
                    ha="center", va="center", transform=ax.transAxes,
                    color="grey", fontsize=9)
            ax.axis("off")
        else:
            top  = gene_ctr.most_common(top_n)
            glbls = [g for g, _ in top]
            gpcts = [100 * v / n for _, v in top]
            gy = np.arange(len(glbls))
            ax.barh(gy, gpcts, color="#e07b54", edgecolor="white",
                    linewidth=0.3, height=0.76)
            ax.set_yticks(gy)
            ax.set_yticklabels(glbls, fontsize=7.5, fontstyle="italic")
            ax.invert_yaxis()
            ax.set_xlabel("Prevalence (% of isolates)")
            ax.set_title("C. Virulence genes (AMRFinder)", fontweight="bold")
            ax.set_xlim(0, 108)
            ax.xaxis.set_major_locator(plt.MultipleLocator(20))

    fig.suptitle("E. coli virulence profiling", fontsize=11, fontweight="bold", y=1.01)
    plt.tight_layout()
    _save(fig, outdir, f"{prefix}_fig5_virulence")


def _fig_virulence_salmonella(df: pd.DataFrame, outdir: Path, prefix: str, top_n: int) -> None:
    # Prefer abricate VFDB; fall back to AMRFinder virulence
    col = next((c for c in ["abricate_vfdb_genes", "amrfinder_virulence_genes"]
                if c in df.columns), None)
    if col is None:
        print("WARNING: no virulence gene data — skipping virulence plot", file=sys.stderr)
        return

    n = len(df)
    gene_ctr: Counter = Counter()
    for val in df[col]:
        for g in parse_list(val):
            gene_ctr[g] += 1

    if not gene_ctr:
        fig, ax = plt.subplots(figsize=(6.5, 2))
        ax.text(0.5, 0.5, "No virulence genes detected",
                ha="center", va="center", transform=ax.transAxes, color="grey", fontsize=10)
        ax.set_title("Virulence genes (VFDB / AMRFinder)", fontweight="bold")
        ax.axis("off")
        _save(fig, outdir, f"{prefix}_fig5_virulence")
        return

    top    = gene_ctr.most_common(top_n)
    labels = [g for g, _ in top]
    pcts   = [100 * v / n for _, v in top]

    fig, ax = plt.subplots(figsize=(6.5, len(labels) * 0.28 + 2.0))
    y = np.arange(len(labels))
    ax.barh(y, pcts, color="#e07b54", edgecolor="white", linewidth=0.3, height=0.76)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5, fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlabel("Prevalence (% of isolates)")
    src = "VFDB (abricate)" if col == "abricate_vfdb_genes" else "AMRFinder"
    ax.set_title(f"Virulence genes — {src}", fontweight="bold")
    ax.set_xlim(0, 108)
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    plt.tight_layout()
    _save(fig, outdir, f"{prefix}_fig5_virulence")


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
    fig_virulence(df, outdir, args.prefix)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
