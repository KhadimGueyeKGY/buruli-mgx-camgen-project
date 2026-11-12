#!/usr/bin/env python3
"""
Beta diversity, point 5 of the analysis plan.

Bray-Curtis on the bacterial genus relative abundances, with a principal coordinates
ordination. PERMANOVA tests site, Buruli status within each site, and the clinical
variables, each over 9,999 free permutations. UniFrac is not available because SILVA
gives no phylogeny with these reports.

Usage: python3 5_beta_diversity.py

Auteur : Khadim Gueye

"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import common as c

FACTORS = ("group", "lesion_category", "lesion_size", "treatment")


def combined():
    frames, metas = [], []
    for site in c.SITES:
        frame, meta = c.aligned(site, "genus", relative=True)
        frames.append(frame)
        metas.append(meta.assign(site=site))
    table = pd.concat(frames, axis=1).fillna(0)
    meta = pd.concat(metas)
    order = sorted(table.columns)
    return table[order], meta.loc[order]


def ordination(ax, coords, explained, meta, colour_by, title,
               legend_title=None):
    values = meta[colour_by].fillna("unknown")
    palette = c.GROUP_COLOUR if colour_by == "group" else None
    categories = [g for g in c.GROUPS if g in set(values)] if colour_by == "group" \
        else sorted(set(values))
    for i, value in enumerate(categories):
        mask = (values == value).to_numpy()
        colour = palette[value] if palette else c.CATEGORICAL[i % len(c.CATEGORICAL)]
        ax.scatter(coords[mask, 0], coords[mask, 1], s=46, color=colour,
                   edgecolor="white", linewidth=0.8, label=str(value))
    ax.axhline(0, color=c.MUTED, linewidth=0.6, zorder=0)
    ax.axvline(0, color=c.MUTED, linewidth=0.6, zorder=0)
    ax.set_xlabel(f"PCo1 ({100 * explained[0]:.1f} percent of the variation)",
                  color=c.INK)
    ax.set_ylabel(f"PCo2 ({100 * explained[1]:.1f} percent)", color=c.INK)
    c.panel_title(ax, title, pad=12)
    if legend_title:
        c.side_legend(ax, title=legend_title)
    c.tidy_axes(ax)


def main():
    rng = np.random.default_rng(c.SEED)
    table, meta = combined()
    distance = c.bray_curtis(table)
    coords, explained = c.pcoa(distance)

    rows = []
    f_site, p_site = c.permanova(distance, meta["site"].to_numpy(), rng)
    rows.append({"scope": "both sites", "factor": "site", "n": len(meta),
                 "levels": 2, "F": f_site, "p_value": p_site})

    per_site = {}
    for site in c.SITES:
        sub_table, sub_meta = c.aligned(site, "genus", relative=True)
        sub_distance = c.bray_curtis(sub_table)
        per_site[site] = (sub_table, sub_meta, sub_distance, c.pcoa(sub_distance))
        for factor in FACTORS:
            labels = sub_meta[factor].fillna("unknown").to_numpy()
            levels = [v for v in set(labels) if (labels == v).sum() >= 3]
            keep = np.isin(labels, levels)
            if len(levels) < 2 or keep.sum() <= len(levels):
                continue
            f, p = c.permanova(sub_distance[np.ix_(keep, keep)], labels[keep], rng)
            rows.append({"scope": c.LABEL[site], "factor": factor, "n": int(keep.sum()),
                         "levels": len(levels), "F": f, "p_value": p})

    stats = pd.DataFrame(rows)
    stats["p_adjusted_bh"] = c.benjamini_hochberg(stats["p_value"].to_numpy())
    stats.to_csv(c.FIG_DIR / "Table6_beta_diversity_permanova.tsv", sep="\t",
                 index=False, float_format="%.5f")

    fig = plt.figure(figsize=(6.9, 7.2))
    grid = fig.add_gridspec(2, 2, hspace=0.62, wspace=1.15)

    ax = fig.add_subplot(grid[0, 0])
    for site in c.SITES:
        mask = (meta["site"] == site).to_numpy()
        ax.scatter(coords[mask, 0], coords[mask, 1], s=46,
                   color=c.SITE_COLOUR[site], edgecolor="white", linewidth=0.8,
                   label=c.LABEL[site])
    ax.axhline(0, color=c.MUTED, linewidth=0.6, zorder=0)
    ax.axvline(0, color=c.MUTED, linewidth=0.6, zorder=0)
    ax.set_xlabel(f"PCo1 ({100 * explained[0]:.1f} percent of the variation)",
                  color=c.INK)
    ax.set_ylabel(f"PCo2 ({100 * explained[1]:.1f} percent)", color=c.INK)
    c.panel_title(ax, f"A  Both sites, F = {f_site:.2f}, p = {p_site:.4f}", pad=12)
    c.side_legend(ax, title="Site")
    c.tidy_axes(ax)

    for position, letter, site in zip(((0, 1), (1, 0)), "BC", c.SITES):
        ax = fig.add_subplot(grid[position])
        _, sub_meta, _, (sub_coords, sub_explained) = per_site[site]
        row = stats[(stats["scope"] == c.LABEL[site]) & (stats["factor"] == "group")]
        title = f"{letter}  {c.LABEL[site]}"
        if not row.empty:
            title += (f", F = {row['F'].iloc[0]:.2f}, "
                      f"p = {row['p_value'].iloc[0]:.4f}")
        ordination(ax, sub_coords, sub_explained, sub_meta, "group", title,
                   legend_title="Buruli status" if letter == "C" else None)

    ax = fig.add_subplot(grid[1, 1])
    shown = stats.copy()
    shown["name"] = [f"{row['scope']}: {row['factor'].replace('_', ' ')}"
                     for _, row in shown.iterrows()]
    shown = shown.iloc[::-1]
    positions = np.arange(len(shown))
    colours = [c.CATEGORICAL[0] if sig else c.OTHER
               for sig in shown["p_adjusted_bh"] < 0.05]
    ax.barh(positions, shown["F"], height=0.66, color=colours, edgecolor="white",
            linewidth=0.8)
    for y, (f, p) in enumerate(zip(shown["F"], shown["p_adjusted_bh"])):
        ax.text(f + 0.08, y, f"BH p = {p:.3f}", va="center", fontsize=8.5,
                color=c.MUTED)
    ax.set_yticks(positions)
    ax.set_yticklabels(shown["name"], fontsize=9)
    ax.set_xlim(0, float(shown["F"].max()) * 1.75)
    ax.set_xlabel("PERMANOVA pseudo-F", color=c.INK)
    ax.set_ylabel("Site and factor tested", color=c.INK, labelpad=6)
    c.panel_title(ax, "D  Every factor tested", pad=12)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c.CATEGORICAL[0],
                             label="BH p below 0.05"),
               plt.Rectangle((0, 0), 1, 1, color=c.OTHER, label="not significant")]
    c.side_legend(ax, handles=handles, labels=[h.get_label() for h in handles])
    c.tidy_axes(ax, hide=("top", "right", "left"))

    c.save(fig, "Figure5_beta_diversity")

    names = table.columns
    pd.DataFrame(distance, index=names, columns=names).to_csv(
        c.FIG_DIR / "Table7_bray_curtis.tsv", sep="\t", float_format="%.4f")
    print(stats.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
