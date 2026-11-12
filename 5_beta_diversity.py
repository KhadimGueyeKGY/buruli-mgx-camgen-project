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


def ordination(ax, coords, explained, meta, colour_by, title, legend_loc="best"):
    values = meta[colour_by].fillna("unknown")
    palette = c.GROUP_COLOUR if colour_by == "group" else None
    categories = [g for g in c.GROUPS if g in set(values)] if colour_by == "group" \
        else sorted(set(values))
    for i, value in enumerate(categories):
        mask = (values == value).to_numpy()
        colour = palette[value] if palette else c.CATEGORICAL[i % len(c.CATEGORICAL)]
        ax.scatter(coords[mask, 0], coords[mask, 1], s=26, color=colour,
                   edgecolor="white", linewidth=0.6, label=str(value))
    ax.axhline(0, color=c.MUTED, linewidth=0.4, zorder=0)
    ax.axvline(0, color=c.MUTED, linewidth=0.4, zorder=0)
    ax.set_xlabel(f"PCo1 ({100 * explained[0]:.1f} percent)", color=c.INK)
    ax.set_ylabel(f"PCo2 ({100 * explained[1]:.1f} percent)", color=c.INK)
    c.panel_title(ax, title)
    ax.legend(frameon=False, fontsize=7, labelcolor=c.INK, loc=legend_loc)
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

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for site in c.SITES:
        mask = (meta["site"] == site).to_numpy()
        axes[0].scatter(coords[mask, 0], coords[mask, 1], s=26,
                        color=c.SITE_COLOUR[site], edgecolor="white", linewidth=0.6,
                        label=c.LABEL[site])
    axes[0].axhline(0, color=c.MUTED, linewidth=0.4, zorder=0)
    axes[0].axvline(0, color=c.MUTED, linewidth=0.4, zorder=0)
    axes[0].set_xlabel(f"PCo1 ({100 * explained[0]:.1f} percent)", color=c.INK)
    axes[0].set_ylabel(f"PCo2 ({100 * explained[1]:.1f} percent)", color=c.INK)
    c.panel_title(axes[0], f"A  Both sites, PERMANOVA F = {f_site:.2f}, p = {p_site:.4f}")
    axes[0].legend(frameon=False, fontsize=7.5, labelcolor=c.INK)
    c.tidy_axes(axes[0])

    for ax, letter, site in zip(axes[1:], "BC", c.SITES):
        _, sub_meta, _, (sub_coords, sub_explained) = per_site[site]
        row = stats[(stats["scope"] == c.LABEL[site]) & (stats["factor"] == "group")]
        caption = f"{letter}  {c.LABEL[site]}, Buruli status"
        if not row.empty:
            caption += f", F = {row['F'].iloc[0]:.2f}, p = {row['p_value'].iloc[0]:.4f}"
        ordination(ax, sub_coords, sub_explained, sub_meta, "group", caption)

    fig.tight_layout()
    c.save(fig, "Figure4_beta_diversity")

    names = table.columns
    pd.DataFrame(distance, index=names, columns=names).to_csv(
        c.FIG_DIR / "Table7_bray_curtis.tsv", sep="\t", float_format="%.4f")
    print(stats.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
