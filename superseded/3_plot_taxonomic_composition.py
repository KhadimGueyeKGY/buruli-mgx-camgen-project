#!/usr/bin/env python3
"""
Figures of the bacterial composition of the wound swabs, both sites.

Figure 1 shows where the reads go at domain level and how the bacterial fraction
breaks down by phylum. Figure 2 shows the leading genera per sample. Only Bacteria
enter the composition panels, as requested; the domain panel is kept so the size of
the bacterial fraction stays visible.

Run 2b_build_abundance_tables.py first.

Usage: python3 3_plot_taxonomic_composition.py

Auteur : Khadim Gueye
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

ROOT = Path.home() / "khadim/Francis-job"
OUT_ROOT = ROOT / "data/16s_metagenomics/analysis_2026"
FIG_DIR = OUT_ROOT / "figures"
SITES = {"ayos": "A", "akonolinga": "B"}
LABEL = {"ayos": "Ayos", "akonolinga": "Akonolinga"}

CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#6250d6", "#e34948"]
OTHER = "#9a9890"
COLOUR_SITE = {"ayos": "#2a78d6", "akonolinga": "#eb6834"}
SEQUENTIAL = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec",
              "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK = "#0b0b0b"
MUTED = "#52514e"
TOP_PHYLA = 8
TOP_GENERA = 20
GAP = 2.0

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.6,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def load(level):
    frames = []
    for site, letter in SITES.items():
        path = OUT_ROOT / f"{site}_set{letter}/abundance_tables/bacteria_{level}_relative.tsv"
        frame = pd.read_csv(path, sep="\t", index_col=0)
        # Sort first, then label: assigning sorted names to unsorted columns
        # would attach every sample's name to another sample's data.
        frame = frame[sorted(frame.columns)]
        frame.columns = pd.MultiIndex.from_arrays([[site] * frame.shape[1], frame.columns])
        frames.append(frame)
    return pd.concat(frames, axis=1).fillna(0)


def load_domains():
    frames = []
    for site, letter in SITES.items():
        path = OUT_ROOT / f"{site}_set{letter}/abundance_tables/domain_composition.tsv"
        frame = pd.read_csv(path, sep="\t", index_col=0)
        frame["site"] = site
        frames.append(frame)
    return pd.concat(frames)


def positions(columns):
    """One x position per sample, with a gap between the two sites."""
    xs, ticks, offset, prev = [], [], 0.0, None
    for site, _ in columns:
        if prev is not None and site != prev:
            offset += GAP
        xs.append(len(xs) + offset)
        prev = site
    for site in SITES:
        idx = [i for i, (s, _) in enumerate(columns) if s == site]
        ticks.append((np.mean([xs[i] for i in idx]), LABEL[site]))
    return np.array(xs), ticks


def collapse(frame, keep):
    top = frame.sum(axis=1).sort_values(ascending=False).head(keep).index.tolist()
    out = frame.loc[top].copy()
    rest = frame.drop(index=top).sum(axis=0)
    if (rest > 0).any():
        out.loc["Other"] = rest
    return out


def stacked(ax, frame, xs, colours, title):
    bottom = np.zeros(frame.shape[1])
    for name, colour in zip(frame.index, colours):
        values = frame.loc[name].to_numpy(dtype=float)
        ax.bar(xs, values, bottom=bottom, width=0.82, color=colour,
               edgecolor="white", linewidth=0.5, label=name)
        bottom += values
    ax.set_xlim(xs[0] - 0.8, xs[-1] + 0.8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Relative abundance", color=INK)
    ax.set_title(title, loc="left", fontsize=9, fontweight="bold", color=INK)
    ax.tick_params(length=2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def figure1(phyla):
    columns = list(phyla.columns)
    xs, ticks = positions(columns)
    samples = [s for _, s in columns]
    top = collapse(phyla, TOP_PHYLA)
    colours = (CATEGORICAL[:len(top) - 1] + [OTHER]) if "Other" in top.index \
        else CATEGORICAL[:len(top)]

    fig = plt.figure(figsize=(13, 6.9))
    grid = fig.add_gridspec(2, 1, height_ratios=[2.0, 0.62], hspace=0.38)

    ax = fig.add_subplot(grid[0])
    stacked(ax, top, xs, colours, "A  Composition by phylum")
    ax.title.set_position((0, 1.09))
    ax.set_xticks(xs)
    ax.set_xticklabels(samples, rotation=90, fontsize=5.5)
    ax.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), frameon=False,
              fontsize=7.5, labelcolor=INK)
    for x, label in ticks:
        ax.text(x, 1.025, label, ha="center", va="bottom", fontsize=8.5,
                color=MUTED, transform=ax.get_xaxis_transform())

    ax = fig.add_subplot(grid[1])
    sites = list(SITES)
    left = np.zeros(len(sites))
    for name, colour in zip(top.index, colours):
        values = np.array([top.loc[name, [c for c in top.columns if c[0] == s]].mean() * 100
                           for s in sites])
        ax.barh(np.arange(len(sites)), values, left=left, height=0.5, color=colour,
                edgecolor="white", linewidth=0.6)
        for i, (v, l) in enumerate(zip(values, left)):
            if v >= 3:
                ax.text(l + v / 2, i, f"{v:.0f}", ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
        left += values
    ax.set_yticks(np.arange(len(sites)))
    ax.set_yticklabels([LABEL[s] for s in sites])
    ax.set_xlim(0, 100)
    ax.set_xlabel("Mean percent of bacterial reads", color=INK)
    ax.set_title("B  Mean composition per site, same colours as panel A",
                 loc="left", fontsize=9, fontweight="bold", color=INK)
    ax.tick_params(length=2)
    ax.invert_yaxis()
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    fig.tight_layout(rect=[0, 0, 0.87, 0.97])
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"Figure1_bacterial_composition.{ext}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)


def figure2(genera):
    columns = list(genera.columns)
    xs, ticks = positions(columns)
    samples = [s for _, s in columns]
    top = genera.sum(axis=1).sort_values(ascending=False).head(TOP_GENERA).index.tolist()
    data = genera.loc[top] * 100

    cmap = LinearSegmentedColormap.from_list("blues", SEQUENTIAL)
    fig, ax = plt.subplots(figsize=(13, 6.0))
    mesh = ax.pcolormesh(np.arange(len(samples) + 1), np.arange(len(top) + 1),
                         data.to_numpy()[::-1], cmap=cmap,
                         norm=PowerNorm(0.5, vmin=0, vmax=100),
                         edgecolors="white", linewidth=0.4)
    ax.set_yticks(np.arange(len(top)) + 0.5)
    ax.set_yticklabels(top[::-1], fontsize=7)
    ax.set_xticks(np.arange(len(samples)) + 0.5)
    ax.set_xticklabels(samples, rotation=90, fontsize=5.5)
    split = sum(1 for s, _ in columns if s == "ayos")
    ax.axvline(split, color=INK, linewidth=1.2)
    for x, label in ticks:
        pos = x if label == LABEL["ayos"] else x - GAP
        ax.text(pos + 0.5, len(top) + 0.4, label, ha="center", va="bottom",
                fontsize=8.5, color=MUTED)
    ax.set_title("Twenty most abundant bacterial genera", loc="left", fontsize=9,
                 fontweight="bold", color=INK, pad=22)
    bar = fig.colorbar(mesh, ax=ax, pad=0.012, fraction=0.025,
                       ticks=[0, 1, 5, 10, 25, 50, 100])
    bar.set_label("Percent of bacterial reads, square root scale", color=INK, fontsize=8)
    bar.outline.set_edgecolor(MUTED)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"Figure2_top_genera.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    domains = load_domains()
    phyla = load("phylum")
    genera = load("genus")
    figure1(phyla)
    figure2(genera)

    def per_site(frame, site):
        return frame[[c for c in frame.columns if c[0] == site]].mean(axis=1) * 100

    top_phyla = collapse(phyla, TOP_PHYLA)
    table1 = top_phyla * 100
    table1.columns = [f"{LABEL[site]}_{sample}" for site, sample in table1.columns]
    for site in reversed(list(SITES)):
        table1.insert(0, f"mean_percent_{LABEL[site]}", per_site(top_phyla, site))
    table1.to_csv(FIG_DIR / "Table1_phylum_composition.tsv", sep="\t",
                  float_format="%.4f")

    summary = pd.DataFrame({
        "mean_percent_all": genera.mean(axis=1) * 100,
        "mean_percent_Ayos": per_site(genera, "ayos"),
        "mean_percent_Akonolinga": per_site(genera, "akonolinga"),
        "median_percent_all": genera.median(axis=1) * 100,
        "samples_detected": (genera > 0).sum(axis=1),
        "samples_above_1pct": (genera > 0.01).sum(axis=1),
    }).sort_values("mean_percent_all", ascending=False)
    summary.head(40).to_csv(FIG_DIR / "Table2_top_genera.tsv", sep="\t",
                            float_format="%.3f")
    domains.to_csv(FIG_DIR / "Table3_domain_composition.tsv", sep="\t",
                   float_format="%.3f")
    print(f"figures and tables written to {FIG_DIR}")


if __name__ == "__main__":
    main()
