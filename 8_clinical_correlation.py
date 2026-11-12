#!/usr/bin/env python3
"""
Clinical correlation, point 8 of the analysis plan.

Lesion size, lesion category and treatment status against the bacterial community.
Lesion size and category are ordinal and are encoded as 1, 2, 3 from the smallest or
mildest upward, so a positive Spearman coefficient means the taxon rises with
severity. Treatment is binary.

Genus abundances are centred log ratios, so the correlations are not forced negative
by the constant sum constraint of relative abundances.

Usage: python3 8_clinical_correlation.py

Auteur : Khadim Gueye
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu, spearmanr

import common as c

SIZE_ORDER = {"<= 5 cm": 1, "5-15 cm": 2, "> 15 cm": 3}
CATEGORY_ORDER = {"I": 1, "II": 2, "III": 3}
TREATMENT_ORDER = {"No": 0, "Yes": 1}
VARIABLES = {"lesion_size": SIZE_ORDER, "lesion_category": CATEGORY_ORDER,
             "treatment": TREATMENT_ORDER}
NICE = {"lesion_size": "Lesion size", "lesion_category": "Lesion category",
        "treatment": "Treatment"}
TOP_GENERA = 20
MIN_PAIRS = 6
MIN_PREVALENCE = 0.25


def clr(counts):
    values = counts.to_numpy(dtype=float) + 0.5
    values = values / values.sum(axis=0, keepdims=True)
    log = np.log2(values)
    return pd.DataFrame(log - log.mean(axis=0, keepdims=True),
                        index=counts.index, columns=counts.columns)


def encoded(meta):
    out = pd.DataFrame(index=meta.index)
    for name, mapping in VARIABLES.items():
        out[name] = meta[name].map(mapping)
    return out


def correlations(site):
    counts, meta = c.aligned(site, "genus")
    keep = counts.index[(counts > 0).mean(axis=1) >= MIN_PREVALENCE]
    transformed = clr(counts).loc[keep]
    chosen = counts.loc[keep].sum(axis=1).sort_values(ascending=False) \
        .head(TOP_GENERA).index.tolist()
    scores = encoded(meta)

    rows = []
    for variable in VARIABLES:
        values = scores[variable]
        usable = values.dropna().index
        if len(usable) < MIN_PAIRS or values.loc[usable].nunique() < 2:
            continue
        for genus in chosen:
            rho, p = spearmanr(transformed.loc[genus, usable].to_numpy(),
                               values.loc[usable].to_numpy())
            rows.append({"site": c.LABEL[site], "variable": variable, "genus": genus,
                         "n": len(usable), "spearman_rho": rho, "p_value": p})
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["p_adjusted_bh"] = c.benjamini_hochberg(frame["p_value"].to_numpy())
    return frame, chosen


def diversity_against_clinical(alpha):
    rows = []
    for site in c.SITES:
        meta = c.load_metadata(site)
        sub = alpha[alpha["site"] == c.LABEL[site]].set_index("sample_id")
        shared = sorted(set(sub.index) & set(meta.index))
        scores = encoded(meta.loc[shared])
        for metric in ("shannon", "richness", "simpson", "chao1"):
            for variable in VARIABLES:
                values = scores[variable].dropna()
                usable = sorted(set(values.index) & set(shared))
                if len(usable) < MIN_PAIRS or values.loc[usable].nunique() < 2:
                    continue
                rho, p = spearmanr(sub.loc[usable, metric].to_numpy(),
                                   values.loc[usable].to_numpy())
                rows.append({"site": c.LABEL[site], "metric": metric,
                             "variable": variable, "n": len(usable),
                             "spearman_rho": rho, "p_value": p})
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["p_adjusted_bh"] = c.benjamini_hochberg(frame["p_value"].to_numpy())
    return frame


def figure(correlation, alpha, diversity_stats, chosen):
    fig = plt.figure(figsize=(5.8, 10.4))
    grid = fig.add_gridspec(2, 2, height_ratios=[2.25, 1.0], hspace=0.28,
                            wspace=0.46)
    cmap = LinearSegmentedColormap.from_list("div", c.DIVERGING)

    ax = fig.add_subplot(grid[0, :])
    if correlation.empty:
        ax.axis("off")
    else:
        columns = [(s, v) for s in correlation["site"].unique() for v in VARIABLES
                   if not correlation[(correlation["site"] == s) &
                                      (correlation["variable"] == v)].empty]
        genera = sorted({g for gs in chosen.values() for g in gs})
        data = pd.DataFrame(index=genera,
                            columns=pd.MultiIndex.from_tuples(columns), dtype=float)
        for site, variable in columns:
            sub = correlation[(correlation["site"] == site) &
                              (correlation["variable"] == variable)]
            data[(site, variable)] = sub.set_index("genus")["spearman_rho"]
        mesh = ax.pcolormesh(np.arange(len(columns) + 1), np.arange(len(genera) + 1),
                             data.to_numpy(dtype=float)[::-1], cmap=cmap,
                             vmin=-1, vmax=1, edgecolors="white", linewidth=0.5)
        for i, (site, variable) in enumerate(columns):
            sub = correlation[(correlation["site"] == site) &
                              (correlation["variable"] == variable)].set_index("genus")
            for j, genus in enumerate(genera[::-1]):
                if genus in sub.index and sub.loc[genus, "p_value"] < 0.05:
                    marker = "**" if sub.loc[genus, "p_adjusted_bh"] < 0.05 else "*"
                    ax.text(i + 0.5, j + 0.5, marker, ha="center", va="center",
                            fontsize=11, fontweight="bold", color=c.INK)
        ax.set_yticks(np.arange(len(genera)) + 0.5)
        ax.set_yticklabels([g[:18] + "..." if len(g) > 21 else g
                            for g in genera[::-1]], fontsize=8)
        ax.set_ylabel("Genus", color=c.INK)
        ax.set_xticks(np.arange(len(columns)) + 0.5)
        ax.set_xticklabels([NICE[v] for _, v in columns], rotation=28, ha="right",
                           fontsize=9.5)
        ax.set_xlabel("Clinical variable, by site", color=c.INK, labelpad=6)
        for edge in np.cumsum([sum(1 for s, _ in columns if s == site)
                               for site in dict.fromkeys(s for s, _ in columns)])[:-1]:
            ax.axvline(edge, color=c.INK, linewidth=1.2)
        for site in dict.fromkeys(s for s, _ in columns):
            idx = [i for i, (s, _) in enumerate(columns) if s == site]
            ax.text(np.mean(idx) + 0.5, len(genera) + 0.4, site, ha="center",
                    va="bottom", fontsize=11, fontweight="bold", color=c.MUTED)
        c.panel_title(ax, "A  Genus abundance against the clinical variables",
                      pad=22)
        bar = fig.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.30,
                           fraction=0.038, aspect=46)
        bar.set_label("Spearman rho    * p below 0.05    "
                      "** Benjamini-Hochberg below 0.05",
                      color=c.INK, fontsize=8.5)
        bar.ax.tick_params(labelsize=8.5)
        bar.outline.set_edgecolor(c.MUTED)

    rng = np.random.default_rng(c.SEED)
    for k, variable in enumerate(("lesion_category", "treatment")):
        ax = fig.add_subplot(grid[1, k])
        levels = [l for l, _ in sorted(VARIABLES[variable].items(), key=lambda kv: kv[1])]
        positions, ticks = [], []
        pos = 0.0
        for site in c.SITES:
            meta = c.load_metadata(site)
            sub = alpha[alpha["site"] == c.LABEL[site]].set_index("sample_id")
            for level in levels:
                members = [s for s in sub.index
                           if s in meta.index and meta.loc[s, variable] == level]
                values = sub.loc[members, "shannon"].dropna() if members else pd.Series(dtype=float)
                if len(values) >= 2:
                    parts = ax.boxplot(values, positions=[pos], widths=0.56,
                                       patch_artist=True, showfliers=False)
                    parts["boxes"][0].set(facecolor=c.SITE_COLOUR[site], alpha=0.30,
                                          edgecolor=c.MUTED, linewidth=0.8)
                    for element in ("whiskers", "caps", "medians"):
                        for item in parts[element]:
                            item.set(color=c.MUTED, linewidth=0.8)
                if len(values):
                    ax.scatter(rng.normal(pos, 0.06, len(values)), values, s=24,
                               color=c.SITE_COLOUR[site], edgecolor="white",
                               linewidth=0.6, zorder=3)
                positions.append(pos)
                ticks.append(level)
                pos += 1
            pos += 1.0
        ax.set_xticks(positions)
        ax.set_xticklabels(ticks, fontsize=9)
        ax.set_ylabel("Shannon index", color=c.INK)
        ax.set_xlabel(f"{NICE[variable]}, by site", color=c.INK, labelpad=40)
        row = diversity_stats[(diversity_stats["metric"] == "shannon") &
                              (diversity_stats["variable"] == variable)]
        c.panel_title(ax, f"{'BC'[k]}  Shannon by {NICE[variable].lower()}", pad=24)
        if not row.empty:
            ax.text(0.0, 1.015, "   ".join(
                f"{r['site']} rho {r['spearman_rho']:.2f}"
                for _, r in row.iterrows()), transform=ax.transAxes, fontsize=8.5,
                color=c.MUTED, va="bottom")
        c.tidy_axes(ax)
        for site, idx in zip(c.SITES, (range(0, len(levels)),
                                       range(len(levels), 2 * len(levels)))):
            chosen_pos = [positions[i] for i in idx if i < len(positions)]
            if chosen_pos:
                ax.text(np.mean(chosen_pos), -0.23, c.LABEL[site], ha="center",
                        va="top", fontsize=10.5, fontweight="bold", color=c.MUTED,
                        transform=ax.get_xaxis_transform())

    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=7,
                          color=c.SITE_COLOUR[site], label=c.LABEL[site])
               for site in c.SITES]
    c.figure_legend(fig, handles, [h.get_label() for h in handles], title="Site")
    c.save(fig, "Figure9_clinical_correlation")


def main():
    alpha = pd.read_csv(c.FIG_DIR / "Table4_alpha_diversity.tsv", sep="\t")
    frames, chosen = [], {}
    for site in c.SITES:
        frame, picked = correlations(site)
        frames.append(frame)
        chosen[site] = picked
    correlation = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    correlation.to_csv(c.FIG_DIR / "Table11_clinical_genus_correlation.tsv", sep="\t",
                       index=False, float_format="%.5f")

    diversity_stats = diversity_against_clinical(alpha)
    diversity_stats.to_csv(c.FIG_DIR / "Table12_clinical_diversity_correlation.tsv",
                           sep="\t", index=False, float_format="%.5f")

    treatment = []
    for site in c.SITES:
        meta = c.load_metadata(site)
        sub = alpha[alpha["site"] == c.LABEL[site]].set_index("sample_id")
        for metric in ("shannon", "richness"):
            groups = [sub.loc[[s for s in sub.index
                               if s in meta.index and meta.loc[s, "treatment"] == level],
                              metric].dropna() for level in ("Yes", "No")]
            if min(len(g) for g in groups) >= 3:
                treatment.append({"site": c.LABEL[site], "metric": metric,
                                  "n_treated": len(groups[0]), "n_untreated": len(groups[1]),
                                  "p_value": mannwhitneyu(groups[0], groups[1],
                                                          alternative="two-sided").pvalue})
    pd.DataFrame(treatment).to_csv(c.FIG_DIR / "Table13_treatment_tests.tsv", sep="\t",
                                   index=False, float_format="%.5f")

    figure(correlation, alpha, diversity_stats, chosen)
    if not correlation.empty:
        top = correlation.reindex(correlation["p_value"].sort_values().index).head(10)
        print(top.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    if not diversity_stats.empty:
        print(diversity_stats.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    if not treatment:
        print("treatment: too few untreated patients for a test at either site")


if __name__ == "__main__":
    main()
