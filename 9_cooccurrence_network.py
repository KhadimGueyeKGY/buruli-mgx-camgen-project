#!/usr/bin/env python3
"""
Co-occurrence networks, point 9 of the analysis plan.

SparCC and CoNet are not installed here, so the same idea is implemented directly:
Spearman correlations between genera on centred log ratio abundances, which removes
the spurious negative correlation that relative abundances impose, with significance
from a permutation null and a Benjamini-Hochberg correction. One network per site.

Nodes are the leading genera, sized by mean abundance and coloured by whether the
genus is higher in Buruli or in non Buruli lesions. Edges are the correlations that
survive the correction, blue for positive and red for negative.

Layout is a deterministic Fruchterman-Reingold placement seeded from the shared
random seed, so the picture is reproducible.

Run 6_differential_abundance.py first, its table supplies the node colours.

Usage: python3 9_cooccurrence_network.py

Auteur : Khadim Gueye
"""

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
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


LABEL_SIZE = 6.5
LABEL_OFFSETS = ((0, 9), (0, -17), (38, 0), (-38, 0), (0, 20), (0, -28),
                 (38, 16), (-38, 16), (38, -16), (-38, -16), (72, 0), (-72, 0))


def layout(nodes, edges, rng, iterations=400):
    """Fruchterman-Reingold on a unit square, written out so the result does not
    depend on networkx. Repulsion is k squared over the distance and attraction is
    the distance squared over k, the forces of the original paper, and positions are
    kept inside the frame so a few isolated nodes cannot squash the rest."""
    n = len(nodes)
    index = {name: i for i, name in enumerate(nodes)}
    pos = rng.random((n, 2))
    adjacency = np.zeros((n, n))
    for _, row in edges.iterrows():
        i, j = index[row["genus_a"]], index[row["genus_b"]]
        adjacency[i, j] = adjacency[j, i] = 1.0
    k = np.sqrt(1.0 / n)
    temperature = 0.1
    for _ in range(iterations):
        delta = pos[:, None, :] - pos[None, :, :]
        distance = np.linalg.norm(delta, axis=-1)
        # The diagonal is zero, so it is floored rather than set to infinity: an
        # infinite diagonal multiplied by a zero adjacency gives nan, which then
        # spreads to every coordinate.
        safe = np.maximum(distance, 1e-6)
        direction = delta / safe[..., None]
        repulsion = (direction * (k ** 2 / safe)[..., None]).sum(axis=1)
        attraction = (direction * (adjacency * safe ** 2 / k)[..., None]).sum(axis=1)
        move = repulsion - attraction
        length = np.linalg.norm(move, axis=1, keepdims=True)
        length[length == 0] = 1.0
        pos += move / length * np.minimum(length, temperature)
        pos = np.clip(pos, 0.0, 1.0)
        temperature *= 0.985
    span = np.ptp(pos, axis=0)
    span[span < 1e-9] = 1.0
    return (pos - pos.min(axis=0)) / span


def draw(ax, site, frame, chosen, counts, enrichment, letter):
    if frame.empty:
        ax.text(0.5, 0.5, "too few samples", ha="center", va="center",
                transform=ax.transAxes, color=c.MUTED)
        c.panel_title(ax, f"{letter}  {c.LABEL[site]}")
        ax.axis("off")
        return 0
    edges = frame[frame["edge"]]
    rng = np.random.default_rng(c.SEED)
    pos = layout(chosen, edges, rng)
    index = {name: i for i, name in enumerate(chosen)}
    relative = counts.div(counts.sum(axis=0), axis=1)
    sizes = relative.loc[chosen].mean(axis=1).to_numpy()
    # Square root so one dominant genus does not swallow the rest of the picture.
    sizes = 30 + 420 * np.sqrt(sizes / sizes.max() if sizes.max() else sizes)

    for _, row in edges.iterrows():
        i, j = index[row["genus_a"]], index[row["genus_b"]]
        colour = c.CATEGORICAL[0] if row["spearman_rho"] > 0 else c.CATEGORICAL[7]
        ax.plot(*zip(pos[i], pos[j]), color=colour, linewidth=0.5 + 1.8 * abs(row["spearman_rho"]),
                alpha=0.55, zorder=1, solid_capstyle="round")

    lookup = enrichment[enrichment["site"] == c.LABEL[site]].set_index("genus") \
        if not enrichment.empty else pd.DataFrame()
    colours = []
    for genus in chosen:
        if genus in getattr(lookup, "index", []):
            fold = lookup.loc[genus, "log2_fold_change"]
            colours.append(c.GROUP_COLOUR["BU"] if fold > 0 else c.GROUP_COLOUR["non-BU"])
        else:
            colours.append(c.OTHER)
    ax.scatter(pos[:, 0], pos[:, 1], s=sizes, c=colours, edgecolor="white",
               linewidth=0.9, zorder=2)
    degree = {name: 0 for name in chosen}
    for _, row in edges.iterrows():
        degree[row["genus_a"]] += 1
        degree[row["genus_b"]] += 1
    ax.set_xlim(-0.10, 1.10)
    ax.set_ylim(-0.10, 1.10)
    ax.set_aspect("equal")
    ax.axis("off")

    labelled = sorted(degree, key=degree.get, reverse=True)[:10]
    # A label goes to the first offset whose box clears the boxes already placed.
    # The box is measured with the renderer rather than guessed from the name
    # length, because a guess puts a long genus underneath a short neighbour.
    font = FontProperties(family=plt.rcParams["font.family"], size=LABEL_SIZE)
    renderer = ax.figure.canvas.get_renderer()
    per_point = ax.figure.dpi / 72.0
    pad = 0.12 * LABEL_SIZE * per_point
    placed = []
    for genus in labelled:
        width, height, descent = renderer.get_text_width_height_descent(
            genus, font, False)
        anchor = ax.transData.transform(pos[index[genus]])
        chosen_offset, chosen_box = LABEL_OFFSETS[0], None
        for candidate in LABEL_OFFSETS:
            x = anchor[0] + candidate[0] * per_point
            y = anchor[1] + candidate[1] * per_point
            box = (x - width / 2 - pad, x + width / 2 + pad,
                   y - descent - pad, y + height - descent + pad)
            if all(box[1] < q[0] or box[0] > q[1] or box[3] < q[2] or box[2] > q[3]
                   for q in placed):
                chosen_offset, chosen_box = candidate, box
                break
        if chosen_box is None:
            x = anchor[0] + chosen_offset[0] * per_point
            y = anchor[1] + chosen_offset[1] * per_point
            chosen_box = (x - width / 2 - pad, x + width / 2 + pad,
                          y - descent - pad, y + height - descent + pad)
        placed.append(chosen_box)
        ax.annotate(genus, pos[index[genus]], textcoords="offset points",
                    xytext=chosen_offset, ha="center", fontsize=LABEL_SIZE,
                    color=c.INK, zorder=4,
                    bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                              edgecolor="none", alpha=0.8))
    c.panel_title(ax, f"{letter}  {c.LABEL[site]}, {len(edges)} edges")
    return len(edges)


def main():
    rng = np.random.default_rng(c.SEED)
    enrichment = pd.read_csv(c.FIG_DIR / "Table8_differential_abundance.tsv", sep="\t") \
        if (c.FIG_DIR / "Table8_differential_abundance.tsv").exists() else pd.DataFrame()

    frames = {}
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.8), gridspec_kw={"wspace": 0.04})
    for ax, letter, site in zip(axes, "AB", c.SITES):
        frame, chosen, counts = network(site, rng)
        frames[site] = frame
        draw(ax, site, frame, chosen, counts, enrichment, letter)

    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=c.GROUP_COLOUR["BU"],
                   markersize=7, label="higher in BU"),
        plt.Line2D([], [], marker="o", linestyle="", color=c.GROUP_COLOUR["non-BU"],
                   markersize=7, label="higher in non-BU"),
        plt.Line2D([], [], marker="o", linestyle="", color=c.OTHER, markersize=7,
                   label="not tested"),
        plt.Line2D([], [], color=c.CATEGORICAL[0], linewidth=2, label="positive edge"),
        plt.Line2D([], [], color=c.CATEGORICAL[7], linewidth=2, label="negative edge"),
    ]
    fig.legend(handles=handles, frameon=False, fontsize=7.5, labelcolor=c.INK,
               ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Genus co-occurrence, Spearman on centred log ratios, "
                 "Benjamini-Hochberg below 0.05 and absolute rho at least 0.5",
                 fontsize=9, fontweight="bold", color=c.INK, x=0.012, ha="left", y=0.99)
    c.save(fig, "Figure8_cooccurrence_network")

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
