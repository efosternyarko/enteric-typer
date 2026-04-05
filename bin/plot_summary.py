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
  {prefix}_fig5_virulence            — virulence gene prevalence
  {prefix}_fig7_amr_by_st            — AMRnet-style tile heatmap: drug class % by MLST ST
  {prefix}_fig8_amr_by_group         — AMRnet-style tile heatmap: drug class % by serovar/phylogroup
  {prefix}_fig9_shigella_serotypes   — Shigella species + serotype composition (stacked bars)
  {prefix}_fig10_shigella_features   — Shigella virulence & invasion feature panel (binary heatmap)
  {prefix}_fig11_shigella_is         — Shigella IS element landscape (copy-number heatmap)
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
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
    # Strip float artifact written by pandas (e.g. "68.0" → "68")
    if "." in s:
        try:
            s = str(int(float(s)))
        except ValueError:
            pass
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
    is_shigella   = "shigeifinder_serotype" in df.columns
    is_salmonella = "sistr_serovar" in df.columns or "salmonella" in prefix.lower()
    sp_label = ("Shigella spp." if is_shigella
                else "Salmonella enterica" if is_salmonella
                else "Escherichia coli")

    fig = plt.figure(figsize=(18, 11))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.50, wspace=0.55)
    _panel_st(df, fig.add_subplot(gs[0, 0]))

    if is_shigella:
        # Panel B: IS element landscape (replaces K-locus/serotype panel for Shigella)
        _panel_shigella_is(df, fig.add_subplot(gs[0, 1]))
    else:
        _panel_sero(df, fig.add_subplot(gs[0, 1]))

    _panel_amr(df, fig.add_subplot(gs[1, 0]))
    _panel_mdr(df, fig.add_subplot(gs[1, 1]))
    fig.suptitle(f"{sp_label} genomic surveillance summary  (n = {len(df)} isolates)",
                 fontsize=12, fontweight="bold", y=1.01)
    _save(fig, outdir, f"{prefix}_fig1_population_summary")


def _panel_st(df: pd.DataFrame, ax: plt.Axes, top_n: int = 15) -> None:
    sts = df["mlst_st"].apply(clean_st) if "mlst_st" in df.columns else pd.Series(
        ["Unknown"] * len(df), index=df.index)

    # Detect species: Shigella > Salmonella > E. coli
    is_shigella   = "shigeifinder_serotype" in df.columns
    is_salmonella = (not is_shigella and
                     ("sistr_serovar" in df.columns or "sistr_cgmlst_ST" in df.columns))

    ctr = Counter(sts)
    unk_n = ctr.pop("Unknown", 0)
    top   = ctr.most_common(top_n)
    other_n = sum(v for k, v in ctr.items() if k not in dict(top))

    labels, values = zip(*top) if top else ([], [])
    labels, values = list(labels), list(values)
    if other_n:  labels.append("Other");    values.append(other_n)
    if unk_n:    labels.append("Unknown");  values.append(unk_n)

    y = np.arange(len(labels))

    if is_shigella:
        # Stacked bars by Shigella species using _SHIGELLA_SPECIES_PALETTE
        sp_series = _infer_shigella_species(df)
        top_st_set = set(dict(top).keys())
        label_sp_counts: dict[str, Counter] = {}
        for st_val, sp_val in zip(sts, sp_series):
            if st_val == "Unknown":
                lbl = "Unknown"
            elif st_val in top_st_set:
                lbl = st_val
            else:
                lbl = "Other"
            label_sp_counts.setdefault(lbl, Counter())[sp_val] += 1

        sp_order = list(_SHIGELLA_SPECIES_PALETTE.keys())
        seen_sp: dict[str, str] = {}
        for y_idx, lbl in enumerate(labels):
            x = 0
            sp_ctr = label_sp_counts.get(lbl, Counter())
            for sp in sorted(sp_ctr.keys(),
                             key=lambda s: sp_order.index(s) if s in sp_order else 99):
                cnt = sp_ctr[sp]
                color = _SHIGELLA_SPECIES_PALETTE.get(sp, "#adb5bd")
                ax.barh(y_idx, cnt, left=x, height=0.78,
                        color=color, edgecolor="white", linewidth=0.6)
                seen_sp.setdefault(sp, color)
                x += cnt

        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Number of isolates", fontsize=8.5)
        ax.set_title("A. MLST sequence types", fontweight="bold")
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=5))
        ax.tick_params(axis="both", length=3)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_linewidth(0.8)
        ax.yaxis.set_ticks_position("none")
        patches = [mpatches.Patch(facecolor=c, label=sp, linewidth=0)
                   for sp, c in seen_sp.items()]
        ax.legend(handles=patches, title="Species", fontsize=7, title_fontsize=7.5,
                  bbox_to_anchor=(1.02, 1), loc="upper left",
                  frameon=False, handlelength=1.2, handleheight=1.2)
        return

    if is_salmonella:
        # Single-hue gradient: most frequent ST = deepest, least frequent = lightest
        n_real = sum(1 for l in labels if l not in ("Other", "Unknown"))
        cmap = plt.cm.Blues
        colors = []
        for i, l in enumerate(labels):
            if l in ("Other", "Unknown"):
                colors.append("#bab0ac")
            else:
                frac = 0.85 - (i / max(n_real - 1, 1)) * 0.50
                colors.append(cmap(frac))
        ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.6, height=0.78)
        title_suffix = ""
        show_pg_legend = False

    else:
        # E. coli: stacked bars by phylogroup — every bar (including "Other") is
        # split by the actual phylogroup composition of isolates in that bin.
        has_real_pg = (
            "clermont_phylogroup" in df.columns and
            df["clermont_phylogroup"].notna().any() and
            df["clermont_phylogroup"].astype(str).str.strip().ne("NA").any()
        )
        if has_real_pg:
            pg_series = df["clermont_phylogroup"].astype(str).str.strip()
            pg_series = pg_series.where(
                ~pg_series.isin({"NA", "nan", "None", "-", ""}), "Unknown")
            title_suffix = " (Clermont)"
        else:
            pg_series = sts.apply(get_phylogroup)
            title_suffix = " (ST-inferred)"

        # Accumulate phylogroup counts per display label (individual ST, "Other", "Unknown")
        top_st_set = set(dict(top).keys())
        label_pg_counts: dict[str, Counter] = {}
        for st_val, pg_val in zip(sts, pg_series):
            if st_val == "Unknown":
                lbl = "Unknown"
            elif st_val in top_st_set:
                lbl = st_val
            else:
                lbl = "Other"
            label_pg_counts.setdefault(lbl, Counter())[pg_val] += 1

        # Draw stacked horizontal bars — each segment = one phylogroup
        pg_order = list(PHYLOGROUP_COLORS.keys())
        seen_pg: dict[str, str] = {}
        for y_idx, lbl in enumerate(labels):
            x = 0
            pg_ctr = label_pg_counts.get(lbl, Counter())
            # Sort segments by canonical phylogroup order for visual consistency
            for pg in sorted(pg_ctr.keys(),
                             key=lambda p: pg_order.index(p) if p in pg_order else 99):
                cnt = pg_ctr[pg]
                color = PHYLOGROUP_COLORS.get(pg, PHYLOGROUP_COLORS["Unknown"])
                ax.barh(y_idx, cnt, left=x, height=0.78,
                        color=color, edgecolor="white", linewidth=0.6)
                seen_pg.setdefault(pg, color)
                x += cnt

        show_pg_legend = True

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

    if show_pg_legend:
        patches = [mpatches.Patch(facecolor=c, label=pg, linewidth=0)
                   for pg, c in seen_pg.items()]
        ax.legend(handles=patches, title="Phylogroup", fontsize=7, title_fontsize=7.5,
                  bbox_to_anchor=(1.02, 1), loc="upper left",
                  frameon=False, handlelength=1.2, handleheight=1.2)


# K-locus group colours (G1–G4 + Unknown) — consistent across all figures
KGROUP_COLORS = {
    "G1": "#4e79a7", "G2": "#f28e2b",
    "G3": "#59a14f", "G4": "#e15759",
    # Kaptive reports database groups as compound strings
    "G1/G4": "#9b59b6", "G2/G3": "#e67e22",
    "Unknown": "#bab0ac",
}
KGROUP_ORDER = ["G1", "G2", "G3", "G4", "G1/G4", "G2/G3", "Unknown"]


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
    use_klocus    = (species == "ecoli"      and "k_locus"   in df.columns)
    use_mlst_st   = (species == "salmonella" and "mlst_st"   in df.columns)
    use_stcomplex = (species == "salmonella" and not use_mlst_st
                     and "mlst_st_complex" in df.columns)

    # ── Build per-serotype fill counts ────────────────────────────────────────
    if use_klocus:
        raw_klocus = df["k_locus"].fillna("Unknown").replace(
            {"NA": "Unknown", "": "Unknown", "-": "Unknown"})
        # Build locus → group mapping (used for legend and colour lookup)
        locus_grp: dict[str, str] = {}
        if "k_group" in df.columns:
            for locus, grp in zip(raw_klocus,
                                  df["k_group"].fillna("Unknown").replace(
                                      {"NA": "Unknown", "": "Unknown"})):
                locus_grp.setdefault(str(locus), str(grp).strip())
        locus_freq = Counter(raw_klocus)
        top_loci   = [l for l, _ in locus_freq.most_common(12) if l != "Unknown"]
        has_other  = any(l not in top_loci and l != "Unknown" for l in locus_freq)

        # Fill bars by K-group (not individual locus) so that rare K-loci
        # still receive their proper group colour rather than collapsing to
        # a grey "Other" segment.
        if "k_group" in df.columns:
            raw_kgroup = df["k_group"].fillna("Unknown").replace(
                {"NA": "Unknown", "": "Unknown", "-": "Unknown"})
        else:
            raw_kgroup = raw_klocus.apply(lambda l: locus_grp.get(str(l), "Unknown"))
        grp_freq  = Counter(raw_kgroup)
        fill_vals = [g for g in KGROUP_ORDER if g in grp_freq and g != "Unknown"]
        if "Unknown" in grp_freq:
            fill_vals = fill_vals + ["Unknown"]
        fill_col   = raw_kgroup
        fill_color = {g: KGROUP_COLORS.get(g, "#bab0ac") for g in fill_vals}
        fill_label = {g: g for g in fill_vals}
        fill_title = f"{title}  ·  K-locus type (coloured by group)"

    elif use_mlst_st:
        raw_st     = df["mlst_st"].apply(clean_st)
        fill_col   = raw_st
        top_st     = [s for s, _ in Counter(raw_st).most_common(12) if s != "Unknown"]
        has_other  = len(Counter(raw_st)) > len(top_st) + (1 if "Unknown" in Counter(raw_st) else 0)
        fill_vals  = top_st + (["Other"] if has_other else []) + (["Unknown"] if "Unknown" in Counter(raw_st) else [])
        st_pal     = sns.color_palette("husl", len(top_st))
        fill_color = {s: st_pal[i] for i, s in enumerate(top_st)}
        fill_color.update({"Other": "#bab0ac", "Unknown": "#d3d3d3"})
        fill_label = {s: s for s in fill_vals}
        fill_title = f"{title}  ·  filled by sequence type (MLST)"

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
            fv = "Other" if "Other" in valid else "Unknown"

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
    if use_klocus:
        leg_title = "K-locus type"
    elif use_mlst_st:
        leg_title = "Sequence type"
    else:
        leg_title = "ST complex"

    locus_leg = ax.legend(handles=patches, fontsize=7, ncol=1,
                          bbox_to_anchor=(1.02, 1), loc="upper left",
                          frameon=False, handlelength=1.2, handleheight=1.2,
                          title=leg_title, title_fontsize=7.5)

    # K-locus group colour key — second legend below the locus legend
    # Only include groups actually observed in the top loci; add "Other" entry if
    # rare loci were collapsed into the "Other" bar segment.
    if use_klocus:
        ax.add_artist(locus_leg)   # keep the first legend when adding a second
        # Show all K-groups present in fill_vals (bars are stacked by group,
        # so every group colour actually appears in the chart)
        group_patches = [
            mpatches.Patch(facecolor=KGROUP_COLORS[g], label=g, linewidth=0)
            for g in ["G1", "G2", "G3", "G4", "G1/G4", "G2/G3"]
            if g in fill_vals
        ]
        if group_patches:
            ax.legend(handles=group_patches, fontsize=7, ncol=1,
                      bbox_to_anchor=(1.02, 0), loc="lower left",
                      frameon=False, handlelength=1.2, handleheight=1.2,
                      title="K-locus group", title_fontsize=7.5)


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


def _panel_shigella_is(df: pd.DataFrame, ax: plt.Axes) -> None:
    """Compact IS element copy-number heatmap for fig1 Panel B (Shigella)."""
    import re as _re

    if "is_elements" not in df.columns:
        ax.set_title("B   IS element landscape"); return

    df2 = df.copy()
    df2["_species"]  = _infer_shigella_species(df2)
    df2["_serotype"] = df2.get(
        "shigeifinder_serotype", pd.Series("NA", index=df2.index)
    ).fillna("NA").astype(str)
    sp_rank = {s: i for i, s in enumerate(_SHIGELLA_SPECIES_PALETTE)}
    df2["_sp_rank"] = df2["_species"].map(sp_rank).fillna(99)
    df2 = df2.sort_values(["_sp_rank", "_serotype", "sample"]).reset_index(drop=True)

    def _parse_is(val) -> dict:
        counts: dict = {}
        if pd.isna(val) or str(val).strip() in {"", "NA", "nan", "-"}:
            return counts
        for tok in str(val).split(";"):
            m = _re.match(r'([A-Za-z0-9_]+)\((\d+)\)', tok.strip())
            if m:
                counts[m.group(1)] = int(m.group(2))
        return counts

    parsed = df2["is_elements"].apply(_parse_is)
    found  = set()
    for d in parsed:
        found.update(d.keys())
    elem_order = [e for e in _IS_ELEMENTS if e in found] + sorted(found - set(_IS_ELEMENTS))
    if not elem_order:
        ax.set_title("B   IS element landscape"); return

    mat = pd.DataFrame(
        [{e: d.get(e, 0) for e in elem_order} for d in parsed],
        columns=elem_order, index=df2.index,
    )
    n_samp, n_elem = mat.shape
    vmax = max(mat.values.max(), 1)
    cmap = LinearSegmentedColormap.from_list("is_cmap", ["#f8f9fa", "#1d3557"])
    ax.imshow(mat.values.astype(float), aspect="auto", cmap=cmap,
              vmin=0, vmax=vmax, interpolation="none")

    # Copy numbers in cells
    for ri in range(n_samp):
        for ci in range(n_elem):
            val = int(mat.iloc[ri, ci])
            if val > 0:
                tc = "white" if val > vmax * 0.6 else "#1d3557"
                ax.text(ci, ri, str(val), ha="center", va="center",
                        fontsize=5.5, color=tc, fontweight="bold")

    ax.set_xticks(range(n_elem))
    ax.set_xticklabels(elem_order, rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(n_samp))
    ax.set_yticklabels(df2["sample"].tolist(), fontsize=5.5)
    ax.tick_params(length=0)

    # Colour ytick labels by species
    for lbl, sp in zip(ax.get_yticklabels(), df2["_species"]):
        lbl.set_color(_SHIGELLA_SPECIES_PALETTE.get(sp, "#333333"))

    ax.set_title("B.  IS element landscape", fontweight="bold")


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
    if "salmonella" in prefix.lower():
        # Salmonella VFDB virulence genes are uninformative (core genes present in
        # all isolates). Omitted pending systematic review of useful virulence markers.
        return
    if is_ecoli:
        _fig_virulence_ecoli(df, outdir, prefix, top_n)
    else:
        _fig_virulence_salmonella(df, outdir, prefix, top_n)


def _fig_virulence_ecoli(df: pd.DataFrame, outdir: Path, prefix: str, top_n: int) -> None:
    has_pathovar  = "kleborate_pathovar" in df.columns
    has_markers   = any(col in df.columns for col, _ in PATHOVAR_MARKERS)
    has_amr_vir   = "amrfinder_virulence_genes" in df.columns

    if not has_pathovar and not has_markers and not has_amr_vir:
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
            top   = gene_ctr.most_common(top_n)
            glbls = [g for g, _ in top]
            gpcts = [100 * v / n for _, v in top]
            gy    = np.arange(len(glbls))
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


# ── Figures 7 & 8: AMRnet-style tile heatmaps ────────────────────────────────

# Colour map: 0 % = mid-grey; >0 % → cream → orange → dark-red → purple
_AMRNET_GREY  = np.array([0.502, 0.502, 0.502])
_AMRNET_CMAP  = LinearSegmentedColormap("amrnet_pos", {
    "red":   [(0.0, 1.000, 1.000), (0.5, 1.000, 1.000),
              (0.8, 0.741, 0.741), (1.0, 0.290, 0.290)],
    "green": [(0.0, 0.973, 0.973), (0.5, 0.600, 0.600),
              (0.8, 0.031, 0.031), (1.0, 0.004, 0.004)],
    "blue":  [(0.0, 0.882, 0.882), (0.5, 0.302, 0.302),
              (0.8, 0.051, 0.051), (1.0, 0.549, 0.549)],
})

# Canonical short abbreviations for drug classes
_DC_ABBREV: dict[str, str] = {
    "TETRACYCLINE":        "TET",
    "AMINOGLYCOSIDE":      "AMG",
    "SULFONAMIDE":         "SUL",
    "BETA-LACTAM":         "BLA",
    "TRIMETHOPRIM":        "TMP",
    "FOSFOMYCIN":          "FOS",
    "MACROLIDE":           "MAC",
    "PHENICOL":            "PHE",
    "QUINOLONE":           "QNL",
    "QUINOLONE/TRICLOSAN": "QNL",   # merge with QUINOLONE
    "MULTIDRUG":           "MDR",
    "FOSMIDOMYCIN":        "FOSM",
    "COLISTIN":            "COL",
    "NITROFURAN":          "NIT",
    "STREPTOTHRICIN":      "STR",
}
_EXCLUDE_DC = {"EFFLUX"}   # near-universal in E. coli — excluded from heatmaps


def _amrnet_matrix(
    df: pd.DataFrame,
    row_col: str,
    top_n: int = 15,
) -> pd.DataFrame:
    """
    Build a (rows × drug-class) matrix of % isolates carrying each class.
    Rows = top_n most common values in row_col (Unknown pushed to end).
    Columns = drug classes sorted by overall prevalence (EFFLUX excluded).
    """
    def _parse(cell) -> set[str]:
        if pd.isna(cell) or str(cell).strip() in ("", "-"):
            return set()
        result: set[str] = set()
        for tok in str(cell).split(";"):
            tok = tok.strip()
            if not tok or tok in _EXCLUDE_DC:
                continue
            # Whole-token lookup first (e.g. QUINOLONE/TRICLOSAN → QNL)
            if tok in _DC_ABBREV:
                result.add(_DC_ABBREV[tok])
            else:
                # Expand compound classes (e.g. AMINOGLYCOSIDE/QUINOLONE → AMG + QNL)
                for part in tok.split("/"):
                    part = part.strip()
                    if part and part not in _EXCLUDE_DC:
                        result.add(_DC_ABBREV.get(part, part))
        return result

    # Values that are QC flags rather than real phylogroup/serovar labels
    _NON_GROUP = {"EC_control_fail", "Unknown", "NA", "nan", "-", ""}

    df = df.copy()
    df["_row"] = df[row_col].fillna("Unknown").apply(
        lambda v: clean_st(v) if "st" in row_col.lower() else str(v)
    )
    # Normalise EzClermont control-fail to "Unknown" so it doesn't appear as a row
    df["_row"] = df["_row"].apply(lambda v: "Unknown" if v in _NON_GROUP else v)
    dc_col = next((c for c in ["amrfinder_drug_classes", "amr_classes",
                                "drug_classes", "resistance_classes"]
                   if c in df.columns), None)
    if dc_col is None:
        return pd.DataFrame()
    df["_cls"] = df[dc_col].apply(_parse)

    counts  = Counter(df["_row"])
    ordered = [r for r, _ in counts.most_common() if r not in _NON_GROUP][:top_n]
    if "Unknown" in counts:
        ordered.append("Unknown")

    cls_counts = Counter(c for s in df["_cls"] for c in s)
    all_cls    = [c for c, _ in cls_counts.most_common()]

    rows, labels = [], []
    for cat in ordered:
        sub = df[df["_row"] == cat]
        n   = len(sub)
        rows.append([100 * sum(1 for s in sub["_cls"] if dc in s) / n
                     for dc in all_cls])
        labels.append(f"{cat}  (n={n})")

    return pd.DataFrame(rows, index=labels, columns=all_cls)


def _draw_amrnet_ax(ax: plt.Axes, mat: pd.DataFrame, title: str) -> None:
    """Draw one AMRnet tile heatmap onto *ax*."""
    nrows, ncols = mat.shape
    ax.set_xlim(-0.5, ncols - 0.5)
    ax.set_ylim(-0.5, nrows - 0.5)
    ax.invert_yaxis()
    TILE = 0.88

    for ri in range(nrows):
        for ci in range(ncols):
            pct  = mat.iloc[ri, ci]
            face = _AMRNET_GREY if pct == 0 else np.array(_AMRNET_CMAP(pct / 100))[:3]
            ax.add_patch(mpatches.FancyBboxPatch(
                (ci - TILE / 2, ri - TILE / 2), TILE, TILE,
                boxstyle="square,pad=0",
                facecolor=face, edgecolor="white", linewidth=0.5,
            ))
            if pct > 0:
                lum     = 0.299 * face[0] + 0.587 * face[1] + 0.114 * face[2]
                txt_col = "white" if lum < 0.45 else "#333333"
                label   = f"{pct:.0f}" if pct >= 1 else f"{pct:.1f}"
                ax.text(ci, ri, label, ha="center", va="center",
                        fontsize=7, color=txt_col)

    ax.set_xticks(range(ncols))
    ax.set_xticklabels(mat.columns, rotation=45, ha="left", fontsize=8)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_yticks(range(nrows))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.tick_params(axis="y", length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=9, fontweight="bold", pad=14)

    legend_items = [
        mpatches.Patch(facecolor=_AMRNET_GREY,                 label="0 %"),
        mpatches.Patch(facecolor=_AMRNET_CMAP(0.01)[:3],       label="1–25 %"),
        mpatches.Patch(facecolor=_AMRNET_CMAP(0.50)[:3],       label="50 %"),
        mpatches.Patch(facecolor=_AMRNET_CMAP(0.80)[:3],       label="75–80 %"),
        mpatches.Patch(facecolor=_AMRNET_CMAP(1.00)[:3],       label="100 %"),
    ]
    ax.legend(handles=legend_items, loc="lower right",
              bbox_to_anchor=(1.02, -0.12), ncol=len(legend_items),
              frameon=False, fontsize=7, handlelength=1.2,
              handletextpad=0.4, columnspacing=0.8)


def fig_amrnet_by_st(df: pd.DataFrame, outdir: Path, prefix: str) -> None:
    """Fig 7 — AMR drug class % by MLST ST."""
    mat = _amrnet_matrix(df, "mlst_st", top_n=12)
    if mat.empty or mat.shape[1] == 0:
        return
    nrows, ncols = mat.shape
    fig, ax = plt.subplots(figsize=(max(6, ncols * 0.85 + 3),
                                    max(4, nrows * 0.55 + 2)))
    sp = ("Salmonella enterica" if "sistr_serovar" in df.columns else "Escherichia coli")
    _draw_amrnet_ax(ax, mat,
                    f"{sp} — AMR drug class prevalence by sequence type (ST)")
    plt.tight_layout()
    _save(fig, outdir, f"{prefix}_fig7_amr_by_st")


def fig_amrnet_by_group(df: pd.DataFrame, outdir: Path, prefix: str) -> None:
    """Fig 8 — AMR drug class % by serovar (Salmonella) or phylogroup (E. coli)."""
    is_sal = "sistr_serovar" in df.columns
    if is_sal:
        row_col  = "sistr_serovar"
        sp       = "Salmonella enterica"
        grp_name = "serovar"
        top_n    = 15
    elif "clermont_phylogroup" in df.columns:
        row_col  = "clermont_phylogroup"
        sp       = "Escherichia coli"
        grp_name = "Clermont phylogroup"
        top_n    = 10
    elif "shigeifinder_serotype" in df.columns:
        row_col  = "shigeifinder_serotype"
        sp       = "Shigella spp."
        grp_name = "serotype"
        top_n    = 15
    else:
        return

    mat = _amrnet_matrix(df, row_col, top_n=top_n)
    if mat.empty or mat.shape[1] == 0:
        return
    nrows, ncols = mat.shape
    fig, ax = plt.subplots(figsize=(max(6, ncols * 0.85 + 3),
                                    max(4, nrows * 0.55 + 2)))
    _draw_amrnet_ax(ax, mat,
                    f"{sp} — AMR drug class prevalence by {grp_name}")
    plt.tight_layout()
    _save(fig, outdir, f"{prefix}_fig8_amr_by_group")


# ── Shigella-specific figures ─────────────────────────────────────────────────

_SHIGELLA_SPECIES_PALETTE = {
    "S. sonnei":      "#e63946",
    "S. flexneri":    "#457b9d",
    "S. dysenteriae": "#f4a261",
    "S. boydii":      "#2a9d8f",
    "Unknown":        "#adb5bd",
}

# PINV markers screened by pinv_screen module (gene names as stored in pinv_genes column)
_PINV_GENES = ["icsA_virG", "virF", "virB", "ipaB", "ipaC", "ipaD"]
# IS elements screened by is_screen module
_IS_ELEMENTS = ["IS1", "IS1A", "IS30", "IS186", "IS600", "IS629"]


def _infer_shigella_species(df: pd.DataFrame) -> pd.Series:
    """Return a Series of 'S. sonnei' / 'S. flexneri' / etc. from ShigEiFinder columns."""
    def _classify(row):
        for col in ("shigeifinder_cluster", "shigeifinder_serotype"):
            val = str(row.get(col, "")).lower()
            if "sonnei"      in val: return "S. sonnei"
            if "flexneri"    in val: return "S. flexneri"
            if "dysenteriae" in val or "dysenteri" in val: return "S. dysenteriae"
            if "boydii"      in val: return "S. boydii"
        return "Unknown"
    return df.apply(_classify, axis=1)


def fig_shigella_serotypes(df: pd.DataFrame, outdir: Path, prefix: str) -> None:
    """
    Fig 9 — Shigella species + serotype composition.

    Horizontal stacked bar chart: one bar per species, stacked by serotype.
    Serotype labels from shigeifinder_serotype; species inferred from cluster/serotype.
    """
    if "shigeifinder_serotype" not in df.columns:
        return

    df2 = df.copy()
    df2["_species"]  = _infer_shigella_species(df2)
    df2["_serotype"] = df2["shigeifinder_serotype"].fillna("Unknown").astype(str).str.strip()
    df2["_serotype"] = df2["_serotype"].where(
        ~df2["_serotype"].isin({"", "NA", "nan", "None", "-"}), "Unknown")

    species_order = [s for s in _SHIGELLA_SPECIES_PALETTE if s != "Unknown"]
    # Only species that have samples
    species_order = [s for s in species_order if (df2["_species"] == s).any()]
    if not species_order:
        return

    # Per-species serotype counts
    all_serotypes = sorted(df2["_serotype"].unique())
    # Build a colour palette for serotypes using tab20 / paired
    sero_colors = {s: c for s, c in zip(
        all_serotypes,
        plt.cm.tab20.colors[:len(all_serotypes)]   # type: ignore[attr-defined]
    )}

    fig_h = max(3, len(species_order) * 0.9 + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    lefts = {sp: 0 for sp in species_order}
    for sero in all_serotypes:
        counts = [int((df2[df2["_species"] == sp]["_serotype"] == sero).sum())
                  for sp in species_order]
        bars = ax.barh(species_order, counts, left=[lefts[sp] for sp in species_order],
                       color=sero_colors[sero], edgecolor="white", linewidth=0.5,
                       label=sero if any(c > 0 for c in counts) else "_nolegend_")
        # Annotate count inside bar if wide enough
        for bar, n in zip(bars, counts):
            if n > 0 and bar.get_width() >= 1:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_y() + bar.get_height() / 2,
                        str(n), ha="center", va="center", fontsize=7, color="white",
                        fontweight="bold")
        for sp, cnt in zip(species_order, counts):
            lefts[sp] += cnt

    # Total n per species on right
    for i, sp in enumerate(species_order):
        total = int((df2["_species"] == sp).sum())
        ax.text(lefts[sp] + 0.15, i, f"n={total}", va="center", fontsize=7.5, color="#444")

    ax.set_yticks(range(len(species_order)))
    ax.set_yticklabels([f"$\\it{{{s.replace(' ', '~')}}}$" for s in species_order], fontsize=9)
    ax.set_xlabel("Number of isolates", fontsize=9)
    ax.set_title("Shigella species and serotype composition", fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend — only real serotypes, outside plot
    handles, labels = ax.get_legend_handles_labels()
    seen, h2, l2 = set(), [], []
    for h, l in zip(handles, labels):
        if l not in seen and not l.startswith("_"):
            seen.add(l); h2.append(h); l2.append(l)
    ncol = max(1, len(l2) // 8 + 1)
    ax.legend(h2, l2, title="Serotype", fontsize=6.5, title_fontsize=7,
              loc="lower right", bbox_to_anchor=(1, 1), ncol=ncol, frameon=False)

    plt.tight_layout()
    _save(fig, outdir, f"{prefix}_fig9_shigella_serotypes")


def fig_shigella_features(df: pd.DataFrame, outdir: Path, prefix: str) -> None:
    """
    Fig 10 — Shigella virulence & invasion feature panel.

    Binary presence/absence heatmap modelled on S. flexneri eLife (2015) and
    S. sonnei Nature Communications (2021) comparative genomics figures.
    Columns: ipaH | virulence plasmid | pINV invasion genes | IS elements.
    Rows: samples sorted by species then serotype, with species colour strip.
    """
    needed = {"shigeifinder_ipaH", "shigeifinder_virulence_plasmid", "pinv_genes", "is_elements"}
    if not needed.intersection(df.columns):
        return

    df2 = df.copy()
    df2["_species"]  = _infer_shigella_species(df2)
    df2["_serotype"] = df2.get("shigeifinder_serotype", pd.Series("NA", index=df2.index)).fillna("NA").astype(str)
    # Sort: species order → serotype → sample name
    sp_rank = {s: i for i, s in enumerate(_SHIGELLA_SPECIES_PALETTE)}
    df2["_sp_rank"] = df2["_species"].map(sp_rank).fillna(99)
    df2 = df2.sort_values(["_sp_rank", "_serotype", "sample"]).reset_index(drop=True)

    def _present(val) -> int:
        v = str(val).strip().lower()
        return 0 if v in {"na", "nan", "none", "-", "", "0", "no", "n", "absent"} else 1

    # ── Build binary matrix ───────────────────────────────────────────────────
    records = []
    for _, row in df2.iterrows():
        r: dict = {}
        r["ipaH"]         = _present(row.get("shigeifinder_ipaH", ""))
        r["Vir.\nplasmid"] = _present(row.get("shigeifinder_virulence_plasmid", ""))

        pinv_raw = str(row.get("pinv_genes", "")).lower()
        for g in _PINV_GENES:
            # Match gene name allowing _/- variation and partial tokens
            r[g] = 1 if re.search(r'\b' + re.escape(g.lower().replace("_", ".?")) + r'\b',
                                   pinv_raw.replace("_", ".")) else 0

        is_raw = str(row.get("is_elements", ""))
        for elem in _IS_ELEMENTS:
            r[elem] = 1 if re.search(r'\b' + re.escape(elem) + r'\b', is_raw, re.I) else 0
        records.append(r)

    feat_cols = (["ipaH", "Vir.\nplasmid"] +
                 [g.replace("_", "\n") for g in _PINV_GENES] +
                 _IS_ELEMENTS)
    # Rename keys to match feat_cols
    records2 = []
    for rec in records:
        r2 = {}
        r2["ipaH"]          = rec["ipaH"]
        r2["Vir.\nplasmid"] = rec["Vir.\nplasmid"]
        for g in _PINV_GENES:
            r2[g.replace("_", "\n")] = rec[g]
        for e in _IS_ELEMENTS:
            r2[e] = rec[e]
        records2.append(r2)

    mat = pd.DataFrame(records2, columns=feat_cols)

    n_samp, n_feat = mat.shape
    cell_w, cell_h = 0.55, 0.22
    fig_w = max(9,  n_feat * cell_w + 4.5)
    fig_h = max(4,  n_samp * cell_h + 2.5)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # ── Heatmap ───────────────────────────────────────────────────────────────
    dark_red = "#9b2226"
    light    = "#f8f9fa"
    cmap_bin = LinearSegmentedColormap.from_list("shig_feat", [light, dark_red])
    ax.imshow(mat.values.astype(float), aspect="auto", cmap=cmap_bin,
              vmin=0, vmax=1, interpolation="none")

    # Column dividers between feature groups (after ipaH+plasmid, after PINV genes)
    boundaries = [1.5, 1.5 + len(_PINV_GENES)]
    for xb in boundaries:
        ax.axvline(xb, color="white", lw=2.5)

    # Column group labels above the heatmap
    group_info = [
        (0,   1,                       "ShigEiFinder"),
        (2,   2 + len(_PINV_GENES) - 1, "pINV genes"),
        (2 + len(_PINV_GENES), n_feat - 1, "IS elements"),
    ]
    for x0, x1, label in group_info:
        mid = (x0 + x1) / 2
        ax.annotate(label, xy=(mid, -0.8), xycoords=("data", "axes fraction"),
                    ha="center", va="bottom", fontsize=7.5, fontweight="bold",
                    annotation_clip=False)
        ax.annotate("", xy=(x0 - 0.4, -0.6), xytext=(x1 + 0.4, -0.6),
                    xycoords=("data", "axes fraction"),
                    arrowprops=dict(arrowstyle="-", color="#555", lw=1),
                    annotation_clip=False)

    # Axes ticks
    ax.set_xticks(range(n_feat))
    ax.set_xticklabels(feat_cols, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_samp))
    ax.set_yticklabels(df2["sample"].tolist(), fontsize=6)
    ax.tick_params(length=0)

    # Colour ytick labels by species (replaces colour strip)
    for lbl, sp in zip(ax.get_yticklabels(), df2["_species"]):
        lbl.set_color(_SHIGELLA_SPECIES_PALETTE.get(sp, "#333333"))

    # Species legend
    sp_handles = [mpatches.Patch(color=v, label=k)
                  for k, v in _SHIGELLA_SPECIES_PALETTE.items() if k != "Unknown"]
    leg = ax.legend(handles=sp_handles, title="Species", fontsize=7, title_fontsize=7.5,
                    loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    ax.add_artist(leg)

    # Presence legend
    pres_handles = [
        mpatches.Patch(color=dark_red, label="Present"),
        mpatches.Patch(color=light,    label="Absent", edgecolor="#adb5bd", lw=0.5),
    ]
    ax.legend(handles=pres_handles, fontsize=7, loc="upper left",
              bbox_to_anchor=(1.01, 0.55), frameon=False)

    ax.set_title("Shigella virulence & invasion feature panel", fontsize=10, fontweight="bold")
    plt.tight_layout()
    _save(fig, outdir, f"{prefix}_fig10_shigella_features")


def fig_shigella_is_elements(df: pd.DataFrame, outdir: Path, prefix: str) -> None:
    """
    Fig 11 — Shigella IS element copy-number heatmap.

    Rows = samples sorted by species/serotype.
    Columns = IS element type.
    Cell colour encodes copy number (white = 0, dark = high).
    Literature basis: IS elements account for >50 % of Shigella virulence
    plasmid ORFs; copy-number variation tracks pathoadaptive evolution
    (Parkhill et al.; S. flexneri eLife 2015 pangenome figure).
    """
    if "is_elements" not in df.columns:
        return

    df2 = df.copy()
    df2["_species"]  = _infer_shigella_species(df2)
    df2["_serotype"] = df2.get("shigeifinder_serotype", pd.Series("NA", index=df2.index)).fillna("NA").astype(str)
    sp_rank = {s: i for i, s in enumerate(_SHIGELLA_SPECIES_PALETTE)}
    df2["_sp_rank"] = df2["_species"].map(sp_rank).fillna(99)
    df2 = df2.sort_values(["_sp_rank", "_serotype", "sample"]).reset_index(drop=True)

    # Parse is_elements: "IS600(3);IS629(1)" → {IS600: 3, IS629: 1}
    def _parse_is(val) -> dict:
        counts: dict = {}
        if pd.isna(val) or str(val).strip() in {"", "NA", "nan", "-"}:
            return counts
        for tok in str(val).split(";"):
            m = re.match(r'([A-Za-z0-9_]+)\((\d+)\)', tok.strip())
            if m:
                counts[m.group(1)] = int(m.group(2))
        return counts

    parsed = df2["is_elements"].apply(_parse_is)

    # Determine column order: screened elements first, then any others found
    found_elements = set()
    for d in parsed:
        found_elements.update(d.keys())
    elem_order = [e for e in _IS_ELEMENTS if e in found_elements] + \
                 sorted(found_elements - set(_IS_ELEMENTS))
    if not elem_order:
        return

    mat = pd.DataFrame(
        [{e: d.get(e, 0) for e in elem_order} for d in parsed],
        columns=elem_order,
        index=df2.index,
    )

    n_samp, n_elem = mat.shape
    fig_w = max(5, n_elem * 0.7 + 3.5)
    fig_h = max(3, n_samp * 0.22 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vmax = max(mat.values.max(), 1)
    cmap = LinearSegmentedColormap.from_list("is_cmap", ["#f8f9fa", "#1d3557"])
    im = ax.imshow(mat.values.astype(float), aspect="auto", cmap=cmap,
                   vmin=0, vmax=vmax, interpolation="none")

    # Annotate cells with copy number where > 0
    for row_i in range(n_samp):
        for col_i in range(n_elem):
            val = int(mat.iloc[row_i, col_i])
            if val > 0:
                text_col = "white" if val > vmax * 0.6 else "#1d3557"
                ax.text(col_i, row_i, str(val), ha="center", va="center",
                        fontsize=6.5, color=text_col, fontweight="bold")

    ax.set_xticks(range(n_elem))
    ax.set_xticklabels(elem_order, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(n_samp))
    ax.set_yticklabels(df2["sample"].tolist(), fontsize=6)
    ax.tick_params(length=0)

    # Colour ytick labels by species (replaces colour strip)
    for lbl, sp in zip(ax.get_yticklabels(), df2["_species"]):
        lbl.set_color(_SHIGELLA_SPECIES_PALETTE.get(sp, "#333333"))

    # Species legend
    sp_handles = [mpatches.Patch(color=v, label=k)
                  for k, v in _SHIGELLA_SPECIES_PALETTE.items() if k != "Unknown"]
    ax.legend(handles=sp_handles, title="Species", fontsize=7, title_fontsize=7.5,
              loc="upper left", bbox_to_anchor=(1.12, 1), frameon=False)

    # Colour bar
    cbar = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("Copy number", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)

    ax.set_title("Shigella IS element landscape", fontsize=10, fontweight="bold")
    plt.tight_layout()
    _save(fig, outdir, f"{prefix}_fig11_shigella_is")


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
    fig_amrnet_by_st(df, outdir, args.prefix)
    fig_amrnet_by_group(df, outdir, args.prefix)
    # Shigella-specific figures (silently skipped for non-Shigella datasets)
    fig_shigella_serotypes(df, outdir, args.prefix)
    fig_shigella_features(df, outdir, args.prefix)
    fig_shigella_is_elements(df, outdir, args.prefix)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
