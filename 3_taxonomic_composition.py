#!/usr/bin/env python3
"""
Bacterial composition, point 3 of the analysis plan.

Samples are ordered by site and then by Buruli status, so the two comparisons the
plan asks for can be read off the same axis. Only Bacteria enter these panels.
SILVA resolves to genus, so there is no species panel.

Run 2b_build_abundance_tables.py and 1c_prepare_metadata.py first.

Usage: python3 3_taxonomic_composition.py

Auteur : Khadim Gueye
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

import common as c

TOP_PHYLA = 8
TOP_GENERA = 20
SITE_GAP = 2.4
GROUP_GAP = 1.0


def ordered_samples():
    """One row per sample: site, group, and an x position with gaps between blocks."""
    rows, x = [], 0.0
    for i, site in enumerate(c.SITES):
        if i:
            x += SITE_GAP
        meta = c.load_metadata(site)
        for j, group in enumerate(c.GROUPS):
            if j:
                x += GROUP_GAP
            for sample in sorted(meta.index[meta["group"] == group]):
                rows.append({"site": site, "group": group, "sample_id": sample, "x": x})
                x += 1
    return pd.DataFrame(rows)


def block_labels(order, key):
    return [(sub["x"].mean(), value) for value, sub in order.groupby(key, sort=False)]


def wide(level, order, relative=True):
    """Taxa by sample over both sites, columns in the plotting order."""
    frames = [c.load_table(site, level, relative=relative) for site in c.SITES]
    table = pd.concat(frames, axis=1).fillna(0)
    return table[order["sample_id"].tolist()]


def collapse(frame, keep):
    top = frame.sum(axis=1).sort_values(ascending=False).head(keep).index.tolist()
    out = frame.loc[top].copy()
    rest = frame.drop(index=top).sum(axis=0)
    if (rest > 0).any():
        out.loc["Other"] = rest
    return out


def colours_for(frame):
    return (c.CATEGORICAL[:len(frame) - 1] + [c.OTHER]) if "Other" in frame.index \
        else c.CATEGORICAL[:len(frame)]


def group_bands(ax, order, y=1.012, height=0.026, label_sites=True):
    """A coloured band marking Buruli status above the bars, with the block names."""
    for (site, group), sub in order.groupby(["site", "group"], sort=False):
        left, right = sub["x"].min() - 0.42, sub["x"].max() + 0.42
        ax.add_patch(plt.Rectangle((left, y), right - left, height,
                                   transform=ax.get_xaxis_transform(), clip_on=False,
                                   facecolor=c.GROUP_COLOUR[group], linewidth=0))
        ax.text((left + right) / 2, y + height + 0.008, group, ha="center", va="bottom",
                fontsize=7, color=c.MUTED, transform=ax.get_xaxis_transform())
    if label_sites:
        for x, site in block_labels(order, "site"):
            ax.text(x, y + height + 0.062, c.LABEL[site], ha="center", va="bottom",
                    fontsize=9, color=c.MUTED, transform=ax.get_xaxis_transform())


def stacked(ax, frame, xs, colours, title, pad=16):
    bottom = np.zeros(frame.shape[1])
    for name, colour in zip(frame.index, colours):
        values = frame.loc[name].to_numpy(dtype=float)
        ax.bar(xs, values, bottom=bottom, width=0.82, color=colour,
               edgecolor="white", linewidth=0.5, label=name)
        bottom += values
    ax.set_xlim(xs[0] - 0.8, xs[-1] + 0.8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Relative abundance", color=c.INK)
    c.panel_title(ax, title, pad=pad)
    c.tidy_axes(ax)


def figure1(order, phyla):
    xs = order["x"].to_numpy()
    top = collapse(phyla, TOP_PHYLA)
    colours = colours_for(top)

    fig = plt.figure(figsize=(13.5, 7.6))
    grid = fig.add_gridspec(2, 1, height_ratios=[2.0, 0.85], hspace=0.62)

    ax = fig.add_subplot(grid[0])
    stacked(ax, top, xs, colours, "A  Composition by phylum", pad=50)
    ax.set_xticks(xs)
    ax.set_xticklabels(order["sample_id"], rotation=90, fontsize=5.5)
    ax.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), frameon=False,
              fontsize=7.5, labelcolor=c.INK)
    group_bands(ax, order)

    ax = fig.add_subplot(grid[1])
    blocks = [(site, group) for site in c.SITES for group in c.GROUPS]
    left = np.zeros(len(blocks))
    for name, colour in zip(top.index, colours):
        values = []
        for site, group in blocks:
            picked = order.loc[(order["site"] == site) & (order["group"] == group),
                               "sample_id"]
            values.append(top.loc[name, picked].mean() * 100 if len(picked) else 0.0)
        values = np.array(values)
        ax.barh(np.arange(len(blocks)), values, left=left, height=0.58, color=colour,
                edgecolor="white", linewidth=0.6)
        for i, (v, l) in enumerate(zip(values, left)):
            if v >= 4:
                ax.text(l + v / 2, i, f"{v:.0f}", ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
        left += values
    ax.set_yticks(np.arange(len(blocks)))
    ax.set_yticklabels([f"{c.LABEL[s]}, {g}  (n={int(((order['site']==s)&(order['group']==g)).sum())})"
                        for s, g in blocks], fontsize=7.5)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Mean percent of bacterial reads", color=c.INK)
    c.panel_title(ax, "B  Mean composition per site and group, same colours as panel A")
    ax.invert_yaxis()
    c.tidy_axes(ax, hide=("top", "right", "left"))

    c.save(fig, "Figure1_bacterial_composition")
    return top


def figure2(order, genera):
    top = genera.sum(axis=1).sort_values(ascending=False).head(TOP_GENERA).index.tolist()
    data = genera.loc[top] * 100
    samples = order["sample_id"].tolist()

    cmap = LinearSegmentedColormap.from_list("blues", c.SEQUENTIAL)
    fig, ax = plt.subplots(figsize=(13.5, 6.2))
    mesh = ax.pcolormesh(np.arange(len(samples) + 1), np.arange(len(top) + 1),
                         data.to_numpy()[::-1], cmap=cmap,
                         norm=PowerNorm(0.5, vmin=0, vmax=100),
                         edgecolors="white", linewidth=0.4)
    ax.set_yticks(np.arange(len(top)) + 0.5)
    ax.set_yticklabels(top[::-1], fontsize=7)
    ax.set_xticks(np.arange(len(samples)) + 0.5)
    ax.set_xticklabels(samples, rotation=90, fontsize=5.5)

    compact = order.copy()
    compact["x"] = np.arange(len(order)) + 0.5
    edges = np.cumsum([0] + [int(((order["site"] == s) & (order["group"] == g)).sum())
                             for s in c.SITES for g in c.GROUPS])
    for edge in edges[1:-1]:
        ax.axvline(edge, color=c.MUTED, linewidth=0.7, linestyle=(0, (3, 2)))
    ax.axvline(edges[2], color=c.INK, linewidth=1.3, linestyle="-")
    group_bands(ax, compact, y=1.010, height=0.022)
    c.panel_title(ax, "Twenty most abundant bacterial genera", pad=58)
    bar = fig.colorbar(mesh, ax=ax, pad=0.012, fraction=0.025,
                       ticks=[0, 1, 5, 10, 25, 50, 100])
    bar.set_label("Percent of bacterial reads, square root scale", color=c.INK, fontsize=8)
    bar.outline.set_edgecolor(c.MUTED)
    c.save(fig, "Figure2_top_genera")


def tables(order, phyla, genera, top_phyla):
    def mean_for(frame, site=None, group=None):
        picked = order
        if site:
            picked = picked[picked["site"] == site]
        if group:
            picked = picked[picked["group"] == group]
        return frame[picked["sample_id"].tolist()].mean(axis=1) * 100

    table1 = top_phyla * 100
    table1.columns = [f"{c.LABEL[r['site']]}_{r['group']}_{r['sample_id']}"
                      for _, r in order.iterrows()]
    for site in reversed(list(c.SITES)):
        for group in reversed(c.GROUPS):
            table1.insert(0, f"mean_{c.LABEL[site]}_{group}",
                          mean_for(top_phyla, site, group))
    table1.to_csv(c.FIG_DIR / "Table1_phylum_composition.tsv", sep="\t",
                  float_format="%.4f")

    summary = pd.DataFrame({"mean_percent_all": mean_for(genera)})
    for site in c.SITES:
        for group in c.GROUPS:
            summary[f"mean_{c.LABEL[site]}_{group}"] = mean_for(genera, site, group)
    summary["median_percent_all"] = genera.median(axis=1) * 100
    summary["samples_detected"] = (genera > 0).sum(axis=1)
    summary["samples_above_1pct"] = (genera > 0.01).sum(axis=1)
    summary.sort_values("mean_percent_all", ascending=False).head(40).to_csv(
        c.FIG_DIR / "Table2_top_genera.tsv", sep="\t", float_format="%.3f")

    domains = pd.concat([
        pd.read_csv(c.site_dir(site) / "abundance_tables/domain_composition.tsv",
                    sep="\t", index_col=0).assign(site=c.LABEL[site])
        for site in c.SITES])
    domains.to_csv(c.FIG_DIR / "Table3_domain_composition.tsv", sep="\t",
                   float_format="%.3f")


def main():
    order = ordered_samples()
    phyla = wide("phylum", order)
    genera = wide("genus", order)
    top_phyla = figure1(order, phyla)
    figure2(order, genera)
    tables(order, phyla, genera, top_phyla)
    print(order.groupby(["site", "group"]).size().to_string())
    print(f"figures and tables written to {c.FIG_DIR}")


if __name__ == "__main__":
    main()
