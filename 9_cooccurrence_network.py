#!/usr/bin/env python3
"""
Genus co-occurrence, point 9 of the analysis plan.

SparCC and CoNet are not installed here, so the same idea is implemented directly:
Spearman correlations between genera on centred log ratio abundances, which removes
the spurious negative correlation that relative abundances impose, with significance
from a permutation null and a Benjamini-Hochberg correction. One network per site.

The surviving correlations are drawn as a matrix rather than as a force directed
graph. With 35 genera and 125 edges per site the graph collapses into a core that
no label can reach without covering a node, while the matrix shows every pair at
once and stays readable at the size a page gives it. Every pair, with rho and the
adjusted p value, is in Table14.

Genera are ordered once, over both sites, by the leading eigenvector of the mean
absolute correlation, so genera that move together sit together and the two panels
can be compared row by row.

Usage: python3 9_cooccurrence_network.py

Auteur : Khadim Gueye
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import common as c

TOP_GENERA = 35
MIN_PREVALENCE = 0.30
PERMUTATIONS = 999
ALPHA = 0.05
MIN_ABS_RHO = 0.5


def clr(counts):
    values = counts.to_numpy(dtype=float) + 0.5
    values = values / values.sum(axis=0, keepdims=True)
    log = np.log2(values)
    return pd.DataFrame(log - log.mean(axis=0, keepdims=True),
                        index=counts.index, columns=counts.columns)


def network(site, rng):
    counts, meta = c.aligned(site, "genus")
    keep = counts.index[(counts > 0).mean(axis=1) >= MIN_PREVALENCE]
    chosen = counts.loc[keep].sum(axis=1).sort_values(ascending=False) \
        .head(TOP_GENERA).index.tolist()
    data = clr(counts).loc[chosen]
    values = data.to_numpy()
    n, samples = len(chosen), values.shape[1]
    if n < 3 or samples < 5:
        return pd.DataFrame(), chosen, counts

    observed = spearmanr(values, axis=1).statistic
    observed = np.atleast_2d(observed)
    null = np.empty((PERMUTATIONS, n, n))
    for k in range(PERMUTATIONS):
        shuffled = np.array([row[rng.permutation(samples)] for row in values])
        null[k] = np.atleast_2d(spearmanr(shuffled, axis=1).statistic)

    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            rho = observed[i, j]
            p = (1 + int((np.abs(null[:, i, j]) >= abs(rho)).sum())) / (PERMUTATIONS + 1)
            rows.append({"site": c.LABEL[site], "genus_a": chosen[i],
                         "genus_b": chosen[j], "spearman_rho": rho, "p_value": p})
    frame = pd.DataFrame(rows)
    frame["p_adjusted_bh"] = c.benjamini_hochberg(frame["p_value"].to_numpy())
    frame["edge"] = (frame["p_adjusted_bh"] < ALPHA) & \
                    (frame["spearman_rho"].abs() >= MIN_ABS_RHO)
    return frame, chosen, counts


def shorten(name, limit=20):
    """SILVA joins ambiguous genera with hyphens, which no axis can hold."""
    if len(name) <= limit:
        return name
    head = name.split("-")[0]
    return f"{head}-..." if len(head) + 4 <= limit else name[:limit - 3] + "..."


def seriate(matrix):
    """An order that puts correlated genera next to each other.

    The leading eigenvector of the absolute correlation gives each genus a position
    on one axis; sorting by it is the cheapest spectral seriation and needs no
    clustering library.
    """
    values = np.nan_to_num(np.abs(matrix.to_numpy(dtype=float)))
    np.fill_diagonal(values, 0.0)
    eigenvalues, eigenvectors = np.linalg.eigh(values)
    leading = eigenvectors[:, int(np.argmax(eigenvalues))]
    return [matrix.index[i] for i in np.argsort(leading)]


def matrix_for(frame, order):
    """Square matrix of the surviving correlations, blank where none survived."""
    out = pd.DataFrame(np.nan, index=order, columns=order, dtype=float)
    for _, row in frame[frame["edge"]].iterrows():
        a, b = row["genus_a"], row["genus_b"]
        if a in out.index and b in out.index:
            out.loc[a, b] = out.loc[b, a] = row["spearman_rho"]
    return out


def draw(ax, fig, site, matrix, letter, show_bar):
    shown = matrix.to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list("div", c.DIVERGING)
    cmap.set_bad("#f7f7f6")
    mesh = ax.pcolormesh(np.arange(len(matrix) + 1), np.arange(len(matrix) + 1),
                         np.ma.masked_invalid(shown[::-1]), cmap=cmap,
                         vmin=-1, vmax=1, edgecolors="white", linewidth=0.4)
    # The matrix is symmetric, so the columns carry the same names as the rows and
    # printing them twice only costs room. The names are given once, on panel A.
    ax.set_xticks([])
    ax.set_yticks(np.arange(len(matrix)) + 0.5 if letter == "A" else [])
    if letter == "A":
        ax.set_yticklabels([shorten(name) for name in list(matrix.index)[::-1]],
                           fontsize=7)
        ax.set_ylabel("Genus", color=c.INK)
    ax.set_xlabel("Genus, in the order of the rows", color=c.INK)
    edges = int(np.isfinite(shown).sum() // 2)
    c.panel_title(ax, f"{letter}  {c.LABEL[site]}, {edges} pairs", pad=12)
    # Square cells would tie the row height to the panel width, which leaves the
    # genus names less room than they need. The cells are allowed to be tall.
    ax.tick_params(length=2, width=0.7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(left=(letter == "A"))
    if show_bar:
        bar = fig.colorbar(mesh, ax=ax, pad=0.03, fraction=0.045)
        bar.set_label("Spearman rho", color=c.INK, fontsize=9)
        bar.ax.tick_params(labelsize=8.5)
        bar.outline.set_edgecolor(c.MUTED)


def main():
    rng = np.random.default_rng(c.SEED)
    frames, chosen_by_site = {}, {}
    for site in c.SITES:
        frame, chosen, _ = network(site, rng)
        frames[site] = frame
        chosen_by_site[site] = chosen

    # Only the genera that carry a surviving correlation somewhere: a row of blanks
    # names a genus the figure has nothing to say about.
    linked = set()
    for frame in frames.values():
        if not frame.empty:
            kept = frame[frame["edge"]]
            linked |= set(kept["genus_a"]) | set(kept["genus_b"])
    shared = sorted(linked)
    mean = pd.DataFrame(0.0, index=shared, columns=shared)
    for frame in frames.values():
        if frame.empty:
            continue
        for _, row in frame.iterrows():
            a, b = row["genus_a"], row["genus_b"]
            if a in mean.index and b in mean.index:
                mean.loc[a, b] += abs(row["spearman_rho"]) / len(frames)
                mean.loc[b, a] = mean.loc[a, b]
    order = seriate(mean)

    fig, axes = plt.subplots(1, 2, figsize=(6.6, 5.6),
                             gridspec_kw={"wspace": 0.10})
    for ax, letter, site in zip(axes, "AB", c.SITES):
        draw(ax, fig, site, matrix_for(frames[site], order), letter,
             show_bar=(letter == "B"))
    c.save(fig, "Figure11_cooccurrence_matrix")

    table = pd.concat([f for f in frames.values() if not f.empty], ignore_index=True) \
        if any(not f.empty for f in frames.values()) else pd.DataFrame()
    table.to_csv(c.FIG_DIR / "Table14_cooccurrence_edges.tsv", sep="\t", index=False,
                 float_format="%.5f")
    for site, frame in frames.items():
        if not frame.empty:
            print(f"{c.LABEL[site]}: {len(frame)} pairs tested, "
                  f"{int(frame['edge'].sum())} edges retained")


if __name__ == "__main__":
    main()
