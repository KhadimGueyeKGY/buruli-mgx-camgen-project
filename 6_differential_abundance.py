#!/usr/bin/env python3
"""
Differential abundance, point 6 of the analysis plan.

Buruli against non Buruli, separately for each site. Counts are turned into centred
log ratios on a half count pseudocount, which is the compositional treatment used by
ALDEx2 and ANCOM, then each genus is tested with a two sided Mann-Whitney and the p
values are corrected with Benjamini-Hochberg within the site. DESeq2 assumes a
negative binomial on counts that are not compositional, which these are, and the
group sizes here are small, so a rank test on log ratios is the safer choice.

The effect size is the difference of mean centred log ratios, a log2 fold change on
the log ratio scale. It can disagree in sign with the difference of the mean
percentages when one sample carries most of a genus, because an arithmetic mean of
proportions follows that sample while a mean of log ratios follows the typical one.
Both summaries are reported, with the medians, so the disagreement is visible.

Gram status is assigned from the lineage: Actinomycetota, Bacillota and Deinococcota
are Gram positive, except the class Negativicutes which is Gram negative; the other
phyla present are Gram negative. Phyla whose wall structure is variable or absent are
left unassigned rather than forced into a class.

Usage: python3 6_differential_abundance.py

Auteur : Khadim Gueye
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu

import common as c

MIN_PREVALENCE = 0.25
MIN_MEAN_ABUNDANCE = 0.001
ALPHA = 0.05
TOP_HEATMAP = 20

def clr(counts):
    """Centred log ratio in base 2, on a half count pseudocount."""
    values = counts.to_numpy(dtype=float) + 0.5
    values = values / values.sum(axis=0, keepdims=True)
    log = np.log2(values)
    return pd.DataFrame(log - log.mean(axis=0, keepdims=True),
                        index=counts.index, columns=counts.columns)


def test_site(site):
    counts, meta = c.aligned(site, "genus")
    relative = counts.div(counts.sum(axis=0), axis=1)
    bu = meta.index[meta["group"] == "BU"]
    non = meta.index[meta["group"] == "non-BU"]
    if len(bu) < 3 or len(non) < 3:
        return pd.DataFrame()

    prevalence = (counts > 0).mean(axis=1)
    abundance = pd.concat([relative[bu].mean(axis=1), relative[non].mean(axis=1)],
                          axis=1).max(axis=1)
    keep = counts.index[(prevalence >= MIN_PREVALENCE) & (abundance >= MIN_MEAN_ABUNDANCE)]
    transformed = clr(counts).loc[keep]

    rows = []
    for genus in keep:
        a = transformed.loc[genus, bu].to_numpy()
        b = transformed.loc[genus, non].to_numpy()
        if np.allclose(a.min(), a.max()) and np.allclose(b.min(), b.max()) \
                and np.allclose(a.mean(), b.mean()):
            p = 1.0
        else:
            p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        rows.append({"site": c.LABEL[site], "genus": genus,
                     "n_BU": len(bu), "n_nonBU": len(non),
                     "prevalence": prevalence[genus],
                     "mean_percent_BU": relative.loc[genus, bu].mean() * 100,
                     "mean_percent_nonBU": relative.loc[genus, non].mean() * 100,
                     "median_percent_BU": relative.loc[genus, bu].median() * 100,
                     "median_percent_nonBU": relative.loc[genus, non].median() * 100,
                     "log2_fold_change": a.mean() - b.mean(), "p_value": p})
    frame = pd.DataFrame(rows)
    frame["p_adjusted_bh"] = c.benjamini_hochberg(frame["p_value"].to_numpy())
    frame["significant"] = frame["p_adjusted_bh"] < ALPHA
    return frame.sort_values("p_value")


def gram_table(gram):
    rows = []
    for site in c.SITES:
        counts, meta = c.aligned(site, "genus")
        relative = counts.div(counts.sum(axis=0), axis=1)
        status = gram.reindex(relative.index).fillna("unassigned")
        shares = relative.groupby(status).sum() * 100
        for sample in shares.columns:
            for label, value in shares[sample].items():
                rows.append({"site": c.LABEL[site], "sample_id": sample,
                             "group": meta.loc[sample, "group"],
                             "gram": label, "percent": value})
    return pd.DataFrame(rows)


def gram_tests(frame):
    rows = []
    for site in frame["site"].unique():
        for label in ("Gram positive", "Gram negative"):
            sub = frame[(frame["site"] == site) & (frame["gram"] == label)]
            a = sub.loc[sub["group"] == "BU", "percent"]
            b = sub.loc[sub["group"] == "non-BU", "percent"]
            if len(a) >= 3 and len(b) >= 3:
                rows.append({"site": site, "gram": label, "n_BU": len(a),
                             "n_nonBU": len(b), "median_BU": a.median(),
                             "median_nonBU": b.median(),
                             "p_value": mannwhitneyu(a, b, alternative="two-sided").pvalue})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_adjusted_bh"] = c.benjamini_hochberg(out["p_value"].to_numpy())
    return out


def volcano(ax, frame, title):
    if frame.empty:
        ax.text(0.5, 0.5, "too few samples", ha="center", va="center",
                transform=ax.transAxes, color=c.MUTED)
        c.panel_title(ax, title)
        c.tidy_axes(ax)
        return
    y = -np.log10(frame["p_value"].clip(lower=1e-300))
    for mask, colour, size, order, label in (
            (~frame["significant"], c.OTHER, 26, 2, "not significant"),
            (frame["significant"], c.CATEGORICAL[0], 54, 3, "BH p below 0.05")):
        ax.scatter(frame.loc[mask, "log2_fold_change"], y[mask], s=size, color=colour,
                   edgecolor="white", linewidth=0.7, zorder=order, label=label)
    if frame["significant"].any():
        cut = frame.loc[frame["significant"], "p_value"].max()
        ax.axhline(-np.log10(cut), color=c.CATEGORICAL[1], linewidth=1.2,
                   linestyle=(0, (4, 2)), zorder=1, label="Benjamini-Hochberg 0.05")
    ax.axvline(0, color=c.MUTED, linewidth=0.6, zorder=0)
    for _, row in frame[frame["significant"]].head(8).iterrows():
        ax.annotate(row["genus"],
                    (row["log2_fold_change"], -np.log10(max(row["p_value"], 1e-300))),
                    textcoords="offset points", xytext=(8, 2), fontsize=9.5,
                    fontstyle="italic", color=c.INK)
    ax.set_xlabel("log2 fold change, BU over non-BU", color=c.INK)
    ax.set_ylabel("-log10 p value", color=c.INK)
    c.panel_title(ax, title, pad=12)
    c.tidy_axes(ax)


def heatmap_data(site, frame):
    """Standardised centred log ratios for the lowest p values, patients in order."""
    counts, meta = c.aligned(site, "genus")
    chosen = frame.sort_values("p_value")["genus"].head(TOP_HEATMAP).tolist()
    order = [s for g in c.GROUPS for s in sorted(meta.index[meta["group"] == g])]
    data = clr(counts).loc[chosen, order]
    # Standardised within each genus: the centred log ratios of an abundant genus are
    # all positive, so an unscaled diverging map would read as one colour throughout.
    spread = data.std(axis=1).replace(0, np.nan)
    data = data.sub(data.mean(axis=1), axis=0).div(spread, axis=0).fillna(0)
    return data, chosen, order, meta


def heatmap(ax, fig, site, frame, limit, show_bar):
    if frame.empty:
        ax.axis("off")
        return
    data, chosen, order, meta = heatmap_data(site, frame)
    mesh = ax.pcolormesh(np.arange(len(order) + 1), np.arange(len(chosen) + 1),
                         data.to_numpy()[::-1],
                         cmap=LinearSegmentedColormap.from_list("div", c.DIVERGING),
                         vmin=-limit, vmax=limit, edgecolors="white", linewidth=0.3)
    ax.set_yticks(np.arange(len(chosen)) + 0.5)
    ax.set_yticklabels([g[:18] + "..." if len(g) > 21 else g for g in chosen[::-1]],
                       fontsize=8)
    ax.set_xticks(np.arange(len(order)) + 0.5)
    ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_xlabel("Patient, grouped by Buruli status", color=c.INK, labelpad=6)
    ax.set_ylabel("Genus", color=c.INK)
    split = int((meta["group"] == "BU").sum())
    ax.axvline(split, color=c.INK, linewidth=1.6)
    ax.text(split / 2, len(chosen) + 0.35, "BU", ha="center", va="bottom",
            fontsize=10.5, fontweight="bold", color=c.GROUP_COLOUR["BU"])
    ax.text((split + len(order)) / 2, len(chosen) + 0.35, "non-BU", ha="center",
            va="bottom", fontsize=10.5, fontweight="bold",
            color=c.GROUP_COLOUR["non-BU"])
    if show_bar:
        bar = fig.colorbar(mesh, ax=ax, pad=0.020, fraction=0.034)
        bar.set_label("Centred log ratio, standardised per genus", color=c.INK,
                      fontsize=9.5, fontweight="bold")
        bar.ax.tick_params(labelsize=9)
        bar.outline.set_edgecolor(c.MUTED)


def figure_gram(gram_long, gram_stats):
    """Point 6: Gram composition, averaged per group and per patient."""
    fig, axes = plt.subplots(1, 2, figsize=(5.6, 3.4),
                             gridspec_kw={"wspace": 0.62})
    blocks = [(c.LABEL[s], g) for s in c.SITES for g in c.GROUPS]

    ax = axes[0]
    shares = gram_long.groupby(["site", "group", "gram"])["percent"].mean()
    left = np.zeros(len(blocks))
    for name in ("Gram negative", "Gram positive", "unassigned"):
        values = np.array([shares.get((site, group, name), 0.0)
                           for site, group in blocks])
        ax.barh(np.arange(len(blocks)), values, left=left, height=0.62,
                color=c.GRAM_COLOUR[name], edgecolor="white", linewidth=0.8,
                label=name)
        for i, (v, l) in enumerate(zip(values, left)):
            if v >= 6:
                ax.text(l + v / 2, i, f"{v:.0f}", ha="center", va="center",
                        fontsize=8.5, color="white", fontweight="bold")
        left += values
    ax.set_yticks(np.arange(len(blocks)))
    ax.set_yticklabels([f"{site}, {group}" for site, group in blocks], fontsize=8.5)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Mean percent of bacterial reads", color=c.INK)
    ax.set_ylabel("Site and group", color=c.INK)
    ax.invert_yaxis()
    c.panel_title(ax, "A  Mean Gram composition", pad=12)
    gram_keys = ax.get_legend_handles_labels()
    c.tidy_axes(ax, hide=("top", "right", "left"))

    group_keys = gram_panel(axes[1], gram_long, gram_stats, "B")
    # Figure legends, not axes legends: a key beside the left hand panel would land
    # in the gutter, and an artist added to an axes is cropped at save time.
    for (entries, title), anchor in zip(
            ((gram_keys, "Gram status"), (group_keys, "Buruli status")),
            (0.95, 0.45)):
        legend = fig.legend(entries[0], entries[1], loc="upper left",
                            bbox_to_anchor=(1.0, anchor), frameon=False,
                            labelcolor=c.INK, title=title, handlelength=1.4,
                            borderaxespad=0.0)
        legend.get_title().set_fontsize(9)
        legend.get_title().set_fontweight("bold")
    c.save(fig, "Figure7_gram_composition")


def gram_panel(ax, gram_long, gram_stats, letter):
    """Gram positive share per patient, Buruli against non Buruli at each site."""
    blocks = [(c.LABEL[s], g) for s in c.SITES for g in c.GROUPS]
    rng = np.random.default_rng(c.SEED)
    positions = []
    for i, (site, group) in enumerate(blocks):
        pos = i + (0.45 if i >= 2 else 0)
        positions.append(pos)
        values = gram_long.loc[(gram_long["site"] == site) &
                               (gram_long["group"] == group) &
                               (gram_long["gram"] == "Gram positive"), "percent"]
        if values.empty:
            continue
        parts = ax.boxplot(values, positions=[pos], widths=0.56, patch_artist=True,
                           showfliers=False)
        parts["boxes"][0].set(facecolor=c.GROUP_COLOUR[group], alpha=0.32,
                              edgecolor=c.MUTED, linewidth=1.0)
        for element in ("whiskers", "caps", "medians"):
            for item in parts[element]:
                item.set(color=c.MUTED, linewidth=1.0)
        ax.scatter(rng.normal(pos, 0.06, len(values)), values, s=24,
                   color=c.GROUP_COLOUR[group], edgecolor="white", linewidth=0.6,
                   zorder=3, label=group)
    ax.set_xticks(positions)
    ax.set_xticklabels([g for _, g in blocks], fontsize=9.5)
    ax.set_ylabel("Percent Gram positive", color=c.INK)
    ax.set_xlabel("Site and Buruli status", color=c.INK, labelpad=40)
    ax.set_ylim(bottom=-3)
    c.panel_title(ax, f"{letter}  Gram positive share", pad=12)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.18)
    for k, site in enumerate(c.SITES):
        middle = (positions[2 * k] + positions[2 * k + 1]) / 2
        ax.text(middle, -0.17, c.LABEL[site], ha="center", va="top", fontsize=10,
                fontweight="bold", color=c.MUTED, transform=ax.get_xaxis_transform())
        row = gram_stats[(gram_stats["site"] == c.LABEL[site]) &
                         (gram_stats["gram"] == "Gram positive")]
        if not row.empty:
            ax.text(middle, 0.90, f"p = {row['p_value'].iloc[0]:.3f}", ha="center",
                    va="bottom", fontsize=9, color=c.MUTED,
                    transform=ax.get_xaxis_transform())
    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=8,
                          color=c.GROUP_COLOUR[g], label=g) for g in c.GROUPS]
    c.tidy_axes(ax)
    return handles, [h.get_label() for h in handles]


def figure_differential(results):
    fig = plt.figure(figsize=(5.4, 11.4))
    grid = fig.add_gridspec(3, 2, height_ratios=[0.90, 1.0, 1.0], hspace=0.78,
                            wspace=0.46)
    # One colour scale over both heatmaps, so a cell means the same thing in each.
    limit = max((float(np.abs(heatmap_data(site,
                results[results["site"] == c.LABEL[site]])[0].to_numpy()).max())
                 for site in c.SITES
                 if not results[results["site"] == c.LABEL[site]].empty), default=1.0)
    for k, site in enumerate(c.SITES):
        frame = results[results["site"] == c.LABEL[site]]
        n_sig = int(frame["significant"].sum()) if not frame.empty else 0
        top = fig.add_subplot(grid[0, k])
        volcano(top, frame,
                f"{'AB'[k]}  {c.LABEL[site]}, {n_sig} of {len(frame)} significant")
        if k == len(c.SITES) - 1:
            c.side_legend(top)
        ax = fig.add_subplot(grid[1 + k, :])
        heatmap(ax, fig, site, frame, limit,
                show_bar=(k == len(c.SITES) - 1))
        c.panel_title(ax, f"{'CD'[k]}  {c.LABEL[site]}", pad=30)
    c.save(fig, "Figure6_differential_abundance")


def main():
    lineage = pd.read_csv(c.OUT_ROOT / "bacteria_lineage.tsv", sep="\t").fillna("")
    gram = c.gram_status(lineage)

    results = pd.concat([test_site(site) for site in c.SITES], ignore_index=True)
    results = results.merge(gram.rename("gram"), left_on="genus", right_index=True,
                            how="left")
    results.to_csv(c.FIG_DIR / "Table8_differential_abundance.tsv", sep="\t",
                   index=False, float_format="%.5f")

    gram_long = gram_table(gram)
    gram_long.to_csv(c.FIG_DIR / "Table9_gram_composition.tsv", sep="\t",
                     index=False, float_format="%.3f")
    gram_stats = gram_tests(gram_long)
    gram_stats.to_csv(c.FIG_DIR / "Table10_gram_tests.tsv", sep="\t", index=False,
                      float_format="%.5f")

    figure_differential(results)
    figure_gram(gram_long, gram_stats)
    for site in results["site"].unique():
        sub = results[results["site"] == site]
        print(f"{site}: {len(sub)} genera tested, "
              f"{int(sub['significant'].sum())} significant after Benjamini-Hochberg")
        print(sub.head(8)[["genus", "mean_percent_BU", "mean_percent_nonBU",
                           "median_percent_BU", "median_percent_nonBU",
                           "log2_fold_change", "p_value", "p_adjusted_bh"]]
              .to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    if not gram_stats.empty:
        print(gram_stats.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
