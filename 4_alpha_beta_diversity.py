#!/usr/bin/env python3
"""
Alpha and beta diversity of the bacterial communities, both sites.

Richness, Shannon and Simpson are computed on the bacterial genus counts after
rarefying every sample to the smallest bacterial library, so the comparison is not
driven by sequencing depth. Beta diversity uses Bray-Curtis on relative abundance
with a principal coordinates ordination. Samples are labelled by site only; the
clinical grouping needs the metadata sheets, which are not on this machine.

Run 2b_build_abundance_tables.py first.

Usage: python3 4_alpha_beta_diversity.py

Auteur : Khadim Gueye
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path.home() / "khadim/Francis-job"
OUT_ROOT = ROOT / "data/16s_metagenomics/analysis_2026"
FIG_DIR = OUT_ROOT / "figures"
SITES = {"ayos": "A", "akonolinga": "B"}
LABEL = {"ayos": "Ayos", "akonolinga": "Akonolinga"}
COLOUR = {"ayos": "#2a78d6", "akonolinga": "#eb6834"}
INK = "#0b0b0b"
MUTED = "#52514e"
RAREFY_REPEATS = 100
SEED = 20261111

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.6,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def load_counts():
    frames = []
    for site, letter in SITES.items():
        path = OUT_ROOT / f"{site}_set{letter}/abundance_tables/bacteria_genus_counts.tsv"
        frame = pd.read_csv(path, sep="\t", index_col=0)
        # Sort first, then label: assigning sorted names to unsorted columns
        # would attach every sample's name to another sample's data.
        frame = frame[sorted(frame.columns)]
        frame.columns = pd.MultiIndex.from_arrays([[site] * frame.shape[1], frame.columns])
        frames.append(frame)
    return pd.concat(frames, axis=1).fillna(0).astype(int)


def alpha(counts, depth, rng):
    """Mean richness, Shannon and Simpson over repeated rarefactions."""
    values = counts.to_numpy(dtype=np.int64)
    total = values.sum()
    if total < depth:
        return np.nan, np.nan, np.nan
    probabilities = values / total
    out = []
    for _ in range(RAREFY_REPEATS):
        drawn = rng.multinomial(depth, probabilities)
        drawn = drawn[drawn > 0]
        p = drawn / drawn.sum()
        out.append((len(drawn), float(-(p * np.log(p)).sum()), float(1 - (p ** 2).sum())))
    return tuple(np.mean(out, axis=0))


def bray_curtis(relative):
    values = relative.to_numpy().T
    n = values.shape[0]
    distance = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = values[i], values[j]
            d = np.abs(a - b).sum() / (a.sum() + b.sum())
            distance[i, j] = distance[j, i] = d
    return distance


def pcoa(distance):
    n = distance.shape[0]
    centring = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centring @ (distance ** 2) @ centring
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    positive = values > 0
    coords = vectors[:, positive] * np.sqrt(values[positive])
    explained = values[positive] / values[positive].sum()
    return coords, explained


def permanova(distance, groups, rng, permutations=9999):
    groups = np.asarray(groups)
    n = len(groups)
    labels = np.unique(groups)

    def statistic(assignment):
        total = (distance ** 2).sum() / (2 * n)
        within = 0.0
        for label in labels:
            idx = np.where(assignment == label)[0]
            if len(idx) > 1:
                within += (distance[np.ix_(idx, idx)] ** 2).sum() / (2 * len(idx))
        between = total - within
        return (between / (len(labels) - 1)) / (within / (n - len(labels)))

    observed = statistic(groups)
    null = np.array([statistic(rng.permutation(groups)) for _ in range(permutations)])
    return observed, (1 + (null >= observed).sum()) / (permutations + 1)


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    counts = load_counts()
    depth = int(counts.sum(axis=0).min())

    rows = []
    for site, sample in counts.columns:
        richness, shannon, simpson = alpha(counts[(site, sample)], depth, rng)
        rows.append({"sample_id": sample, "site": LABEL[site],
                     "bacterial_reads": int(counts[(site, sample)].sum()),
                     "richness": richness, "shannon": shannon, "simpson": simpson})
    table = pd.DataFrame(rows)
    table.to_csv(FIG_DIR / "Table4_alpha_diversity.tsv", sep="\t",
                 index=False, float_format="%.4f")

    relative = counts.div(counts.sum(axis=0), axis=1)
    distance = bray_curtis(relative)
    coords, explained = pcoa(distance)
    sites = [LABEL[s] for s, _ in counts.columns]
    f_stat, p_value = permanova(distance, sites, rng)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    for ax, metric, name in zip(axes[:2], ("shannon", "richness"),
                                ("Shannon index", "Observed genera")):
        data = [table.loc[table["site"] == LABEL[s], metric].to_numpy() for s in SITES]
        parts = ax.boxplot(data, widths=0.5, patch_artist=True, showfliers=False)
        for patch, site in zip(parts["boxes"], SITES):
            patch.set(facecolor=COLOUR[site], alpha=0.35, edgecolor=MUTED, linewidth=0.8)
        for element in ("whiskers", "caps", "medians"):
            for item in parts[element]:
                item.set(color=MUTED, linewidth=0.8)
        for i, site in enumerate(SITES, start=1):
            y = table.loc[table["site"] == LABEL[site], metric].to_numpy()
            ax.scatter(rng.normal(i, 0.05, len(y)), y, s=14, color=COLOUR[site],
                       edgecolor="white", linewidth=0.5, zorder=3)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([LABEL[s] for s in SITES])
        ax.set_ylabel(name, color=INK)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    tests = {}
    for metric in ("shannon", "richness"):
        a, b = (table.loc[table["site"] == LABEL[s], metric].to_numpy() for s in SITES)
        tests[metric] = mannwhitneyu(a, b, alternative="two-sided").pvalue
    axes[0].set_title(f"A  Shannon, rarefied to {depth:,} reads, "
                      f"Mann-Whitney p = {tests['shannon']:.3f}",
                      loc="left", fontsize=9, fontweight="bold", color=INK)
    axes[1].set_title(f"B  Observed genera, Mann-Whitney p = {tests['richness']:.3f}",
                      loc="left", fontsize=9, fontweight="bold", color=INK)

    ax = axes[2]
    for site in SITES:
        idx = [i for i, (s, _) in enumerate(counts.columns) if s == site]
        ax.scatter(coords[idx, 0], coords[idx, 1], s=22, color=COLOUR[site],
                   edgecolor="white", linewidth=0.6, label=LABEL[site])
    ax.axhline(0, color=MUTED, linewidth=0.4, zorder=0)
    ax.axvline(0, color=MUTED, linewidth=0.4, zorder=0)
    ax.set_xlabel(f"PCo1 ({100 * explained[0]:.1f} percent)", color=INK)
    ax.set_ylabel(f"PCo2 ({100 * explained[1]:.1f} percent)", color=INK)
    ax.set_title(f"C  Bray-Curtis, PERMANOVA F = {f_stat:.2f}, p = {p_value:.4f}",
                 loc="left", fontsize=9, fontweight="bold", color=INK)
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"Figure3_diversity.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(distance, index=[s for _, s in counts.columns],
                 columns=[s for _, s in counts.columns]).to_csv(
        FIG_DIR / "Table5_bray_curtis.tsv", sep="\t", float_format="%.4f")
    print(f"rarefaction depth {depth}, PERMANOVA F = {f_stat:.3f}, p = {p_value:.4f}")
    for metric, value in tests.items():
        print(f"{metric}: Mann-Whitney p = {value:.4f}")
    print(f"figures and tables written to {FIG_DIR}")


if __name__ == "__main__":
    main()
