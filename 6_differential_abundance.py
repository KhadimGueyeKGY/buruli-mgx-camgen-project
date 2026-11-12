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
TOP_HEATMAP = 25

GRAM_POSITIVE_PHYLA = {"Actinomycetota", "Bacillota", "Deinococcota"}
GRAM_NEGATIVE_CLASSES = {"Negativicutes"}
UNASSIGNED_PHYLA = {"Mycoplasmatota", "Chloroflexota", "Thermomicrobiota"}


def gram_status(lineage):
    """One row per genus: Gram positive, Gram negative or unassigned."""
    out = {}
    for _, row in lineage.iterrows():
        phylum, klass = row["phylum"], row["class"]
        if phylum in UNASSIGNED_PHYLA or not phylum:
            status = "unassigned"
        elif phylum in GRAM_POSITIVE_PHYLA:
            status = "Gram negative" if klass in GRAM_NEGATIVE_CLASSES else "Gram positive"
        else:
            status = "Gram negative"
        out[row["genus"]] = status
    return pd.Series(out, name="gram")


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
    for mask, colour, size, order in ((~frame["significant"], c.OTHER, 15, 2),
                                      (frame["significant"], c.CATEGORICAL[0], 30, 3)):
        ax.scatter(frame.loc[mask, "log2_fold_change"], y[mask], s=size, color=colour,
                   edgecolor="white", linewidth=0.5, zorder=order)
    if frame["significant"].any():
        cut = frame.loc[frame["significant"], "p_value"].max()
        ax.axhline(-np.log10(cut), color=c.CATEGORICAL[1], linewidth=0.9,
                   linestyle=(0, (4, 2)), zorder=1)
        ax.text(ax.get_xlim()[1], -np.log10(cut), "  BH 0.05", va="center",
                fontsize=6.5, color=c.CATEGORICAL[1])
    ax.axvline(0, color=c.MUTED, linewidth=0.4, zorder=0)
    for _, row in frame[frame["significant"]].head(8).iterrows():
        ax.annotate(row["genus"],
                    (row["log2_fold_change"], -np.log10(max(row["p_value"], 1e-300))),
                    textcoords="offset points", xytext=(5, 3), fontsize=6, color=c.INK)
    ax.set_xlabel("log2 fold change of centred log ratios, BU over non-BU", color=c.INK)
    ax.set_ylabel("-log10 p", color=c.INK)
    c.panel_title(ax, title)
    c.tidy_axes(ax)


def heatmap(ax, fig, site, frame):
    if frame.empty:
        ax.axis("off")
        return
    counts, meta = c.aligned(site, "genus")
    chosen = frame.sort_values("p_value")["genus"].head(TOP_HEATMAP).tolist()
    order = [s for g in c.GROUPS for s in sorted(meta.index[meta["group"] == g])]
    data = clr(counts).loc[chosen, order]
    # Standardised within each genus: the centred log ratios of an abundant genus are
    # all positive, so an unscaled diverging map would read as one colour throughout.
    spread = data.std(axis=1).replace(0, np.nan)
    data = data.sub(data.mean(axis=1), axis=0).div(spread, axis=0).fillna(0)
    limit = float(np.abs(data.to_numpy()).max())
    mesh = ax.pcolormesh(np.arange(len(order) + 1), np.arange(len(chosen) + 1),
                         data.to_numpy()[::-1],
                         cmap=LinearSegmentedColormap.from_list("div", c.DIVERGING),
                         vmin=-limit, vmax=limit, edgecolors="white", linewidth=0.3)
    ax.set_yticks(np.arange(len(chosen)) + 0.5)
    ax.set_yticklabels(chosen[::-1], fontsize=6)
    ax.set_xticks(np.arange(len(order)) + 0.5)
    ax.set_xticklabels(order, rotation=90, fontsize=5)
    split = int((meta["group"] == "BU").sum())
    ax.axvline(split, color=c.INK, linewidth=1.2)
    ax.text(split / 2, len(chosen) + 0.3, "BU", ha="center", va="bottom",
            fontsize=7.5, color=c.MUTED)
    ax.text((split + len(order)) / 2, len(chosen) + 0.3, "non-BU", ha="center",
            va="bottom", fontsize=7.5, color=c.MUTED)
    bar = fig.colorbar(mesh, ax=ax, pad=0.014, fraction=0.030)
    bar.set_label("Centred log ratio, standardised per genus", color=c.INK,
                  fontsize=7)
    bar.outline.set_edgecolor(c.MUTED)


def figure_differential(results):
    fig, axes = plt.subplots(2, 2, figsize=(13, 10.4),
                             gridspec_kw={"height_ratios": [1.0, 1.55], "hspace": 0.40,
                                          "wspace": 0.26})
    for k, site in enumerate(c.SITES):
        frame = results[results["site"] == c.LABEL[site]]
        n_sig = int(frame["significant"].sum()) if not frame.empty else 0
        volcano(axes[0, k], frame,
                f"{'AB'[k]}  {c.LABEL[site]}, {len(frame)} genera tested, "
                f"{n_sig} significant")
        heatmap(axes[1, k], fig, site, frame)
        c.panel_title(axes[1, k],
                      f"{'CD'[k]}  {c.LABEL[site]}, {TOP_HEATMAP} lowest p values",
                      pad=16)
    c.save(fig, "Figure5_differential_abundance")


def figure_gram(gram_long, gram_stats):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw={"wspace": 0.26})
    blocks = [(c.LABEL[s], g) for s in c.SITES for g in c.GROUPS]
    rng = np.random.default_rng(c.SEED)

    ax = axes[0]
    positions = []
    for i, (site, group) in enumerate(blocks):
        pos = i + (0.35 if i >= 2 else 0)
        positions.append(pos)
        values = gram_long.loc[(gram_long["site"] == site) & (gram_long["group"] == group)
                               & (gram_long["gram"] == "Gram positive"), "percent"]
        if values.empty:
            continue
        parts = ax.boxplot(values, positions=[pos], widths=0.52, patch_artist=True,
                           showfliers=False)
        parts["boxes"][0].set(facecolor=c.GROUP_COLOUR[group], alpha=0.32,
                              edgecolor=c.MUTED, linewidth=0.8)
        for element in ("whiskers", "caps", "medians"):
            for item in parts[element]:
                item.set(color=c.MUTED, linewidth=0.8)
        ax.scatter(rng.normal(pos, 0.055, len(values)), values, s=14,
                   color=c.GROUP_COLOUR[group], edgecolor="white", linewidth=0.5,
                   zorder=3)
    ax.set_xticks(positions)
    ax.set_xticklabels([g for _, g in blocks], fontsize=7.5)
    ax.set_ylabel("Percent Gram positive", color=c.INK)
    ax.set_ylim(bottom=-2)
    c.panel_title(ax, "A  Gram positive share per sample", pad=14)
    for k, site in enumerate(c.SITES):
        middle = (positions[2 * k] + positions[2 * k + 1]) / 2
        ax.text(middle, -0.11, c.LABEL[site], ha="center", va="top", fontsize=8,
                color=c.MUTED, transform=ax.get_xaxis_transform())
        row = gram_stats[(gram_stats["site"] == c.LABEL[site]) &
                         (gram_stats["gram"] == "Gram positive")]
        if not row.empty:
            ax.text(middle, 1.01, f"p = {row['p_value'].iloc[0]:.3f}", ha="center",
                    va="bottom", fontsize=7, color=c.MUTED,
                    transform=ax.get_xaxis_transform())
    c.tidy_axes(ax)

    ax = axes[1]
    shares = gram_long.groupby(["site", "group", "gram"])["percent"].mean().reset_index()
    labels = [f"{s}\n{g}" for s, g in blocks]
    bottom = np.zeros(len(labels))
    for colour, label in zip((c.CATEGORICAL[0], c.CATEGORICAL[1], c.OTHER),
                             ("Gram negative", "Gram positive", "unassigned")):
        values = np.array([
            shares.loc[(shares["site"] == s) & (shares["group"] == g) &
                       (shares["gram"] == label), "percent"].sum() for s, g in blocks])
        ax.bar(np.arange(len(labels)), values, bottom=bottom, width=0.62, color=colour,
               edgecolor="white", linewidth=0.8, label=label)
        for i, (v, b) in enumerate(zip(values, bottom)):
            if v >= 6:
                ax.text(i, b + v / 2, f"{v:.0f}", ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
        bottom += values
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Mean percent of bacterial reads", color=c.INK)
    c.panel_title(ax, "B  Mean Gram composition", pad=14)
    ax.legend(frameon=False, fontsize=7.5, labelcolor=c.INK,
              loc="upper left", bbox_to_anchor=(1.01, 1.0))
    c.tidy_axes(ax)

    c.save(fig, "Figure6_gram_composition")


def main():
    lineage = pd.read_csv(c.OUT_ROOT / "bacteria_lineage.tsv", sep="\t").fillna("")
    gram = gram_status(lineage)

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
