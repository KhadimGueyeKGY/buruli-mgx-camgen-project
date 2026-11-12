#!/usr/bin/env python3
"""
Alpha diversity, point 4 of the analysis plan.

Observed genera, Shannon, Simpson and Chao1 on the bacterial genus counts. Every
sample is rarefied to the smallest bacterial library so depth cannot drive the
comparison, and the mean over 100 draws is reported. Buruli against non Buruli is
tested within each site with a two sided Mann-Whitney test, and the two sites are
tested against each other on the pooled samples.

Chao1 is computed on the rarefied counts, so it estimates the richness reachable at
that depth rather than the richness of the full library.

Usage: python3 4_alpha_diversity.py

Auteur : Khadim Gueye
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import common as c

METRICS = {"richness": "Observed genera", "shannon": "Shannon index",
           "simpson": "Simpson index", "chao1": "Chao1"}
REPEATS = 100


def table():
    rows, counts = [], {}
    for site in c.SITES:
        frame, meta = c.aligned(site, "genus")
        for sample in frame.columns:
            counts[(site, sample)] = frame[sample]
            rows.append({"site": site, "sample_id": sample,
                         "group": meta.loc[sample, "group"],
                         "lesion_category": meta.loc[sample, "lesion_category"],
                         "treatment": meta.loc[sample, "treatment"],
                         "bacterial_reads": int(frame[sample].sum())})
    frame = pd.DataFrame(rows)
    depth = int(frame["bacterial_reads"].min())
    rng = np.random.default_rng(c.SEED)
    values = [c.rarefy(counts[(r["site"], r["sample_id"])], depth, rng, REPEATS)
              for _, r in frame.iterrows()]
    frame[list(METRICS)] = pd.DataFrame(values, index=frame.index)
    return frame, depth


def tests(frame):
    rows = []
    for metric in METRICS:
        for site in c.SITES:
            sub = frame[frame["site"] == site]
            a = sub.loc[sub["group"] == "BU", metric].dropna()
            b = sub.loc[sub["group"] == "non-BU", metric].dropna()
            if len(a) >= 3 and len(b) >= 3:
                p = mannwhitneyu(a, b, alternative="two-sided").pvalue
                rows.append({"metric": metric, "comparison": f"{c.LABEL[site]} BU vs non-BU",
                             "n1": len(a), "n2": len(b), "median1": a.median(),
                             "median2": b.median(), "p_value": p})
        a = frame.loc[frame["site"] == "ayos", metric].dropna()
        b = frame.loc[frame["site"] == "akonolinga", metric].dropna()
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        rows.append({"metric": metric, "comparison": "Ayos vs Akonolinga",
                     "n1": len(a), "n2": len(b), "median1": a.median(),
                     "median2": b.median(), "p_value": p})
    out = pd.DataFrame(rows)
    out["p_adjusted_bh"] = c.benjamini_hochberg(out["p_value"].to_numpy())
    return out


def figure(frame, stats, depth):
    blocks = [(site, group) for site in c.SITES for group in c.GROUPS]
    rng = np.random.default_rng(c.SEED)
    fig = plt.figure(figsize=(6.3, 6.2))
    grid = fig.add_gridspec(2, 2, hspace=0.92, wspace=0.38)
    axes = [fig.add_subplot(grid[k // 2, k % 2]) for k in range(len(METRICS))]

    for letter, ax, (metric, name) in zip("ABCD", axes, METRICS.items()):
        data, positions = [], []
        for i, (site, group) in enumerate(blocks):
            values = frame.loc[(frame["site"] == site) & (frame["group"] == group),
                               metric].dropna().to_numpy()
            data.append(values)
            positions.append(i + (0.45 if i >= 2 else 0))
        for pos, values, (site, group) in zip(positions, data, blocks):
            if not len(values):
                continue
            parts = ax.boxplot(values, positions=[pos], widths=0.56,
                               patch_artist=True, showfliers=False)
            parts["boxes"][0].set(facecolor=c.GROUP_COLOUR[group], alpha=0.32,
                                  edgecolor=c.MUTED, linewidth=1.0)
            for element in ("whiskers", "caps", "medians"):
                for item in parts[element]:
                    item.set(color=c.MUTED, linewidth=1.0)
            ax.scatter(rng.normal(pos, 0.06, len(values)), values, s=24,
                       color=c.GROUP_COLOUR[group], edgecolor="white",
                       linewidth=0.6, zorder=3)
        ax.set_xticks(positions)
        ax.set_xticklabels([g for _, g in blocks], fontsize=9.5)
        ax.set_ylabel(name, color=c.INK)
        ax.set_xlabel("Site and Buruli status", color=c.INK, labelpad=30)
        c.panel_title(ax, f"{letter}  {name}", pad=12)
        c.tidy_axes(ax)

        top = max((v.max() for v in data if len(v)), default=1.0)
        bottom = min((v.min() for v in data if len(v)), default=0.0)
        span = top - bottom or 1.0
        ax.set_ylim(bottom - 0.08 * span, top + 0.30 * span)
        for k, site in enumerate(c.SITES):
            left, right = positions[2 * k], positions[2 * k + 1]
            middle = (left + right) / 2
            ax.text(middle, -0.17, c.LABEL[site], ha="center", va="top", fontsize=10,
                    fontweight="bold", color=c.MUTED,
                    transform=ax.get_xaxis_transform())
            row = stats[(stats["metric"] == metric) &
                        (stats["comparison"] == f"{c.LABEL[site]} BU vs non-BU")]
            if row.empty:
                continue
            y = top + 0.10 * span
            ax.plot([left, left, right, right],
                    [y, y + 0.035 * span, y + 0.035 * span, y],
                    color=c.MUTED, linewidth=1.0)
            ax.text(middle, y + 0.06 * span,
                    f"p = {row['p_value'].iloc[0]:.3f}", ha="center", va="bottom",
                    fontsize=9, color=c.MUTED)

    handles = [plt.Line2D([], [], marker="o", linestyle="", color=c.GROUP_COLOUR[g],
                          markersize=7, label=g) for g in c.GROUPS]
    c.figure_legend(fig, handles, [h.get_label() for h in handles],
                    title="Buruli status")
    c.save(fig, "Figure4_alpha_diversity")


def main():
    frame, depth = table()
    stats = tests(frame)
    frame_out = frame.copy()
    frame_out["site"] = frame_out["site"].map(c.LABEL)
    frame_out.to_csv(c.FIG_DIR / "Table4_alpha_diversity.tsv", sep="\t",
                     index=False, float_format="%.4f")
    stats.to_csv(c.FIG_DIR / "Table5_alpha_diversity_tests.tsv", sep="\t",
                 index=False, float_format="%.5f")
    figure(frame, stats, depth)
    print(f"rarefaction depth {depth:,}")
    print(stats.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
