#!/usr/bin/env python3
"""
Shared paths, loaders, statistics and plot settings for the Buruli 16S analysis.

Every downstream script reads its data through here, so a sample is filtered,
labelled and coloured the same way everywhere.

Auteur : Khadim Gueye
"""

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path.home() / "khadim/Francis-job"
OUT_ROOT = ROOT / "data/16s_metagenomics/analysis_2026"
FIG_DIR = OUT_ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
SITES = {"ayos": "A", "akonolinga": "B"}
LABEL = {"ayos": "Ayos", "akonolinga": "Akonolinga"}
LEVELS = ("phylum", "family", "genus")
GROUPS = ("BU", "non-BU")

CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#6250d6", "#e34948"]
OTHER = "#9a9890"
SEQUENTIAL = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec",
              "#3987e5", "#256abf", "#184f95", "#0d366b"]
DIVERGING = ["#184f95", "#6da7ec", "#f0efec", "#f0a080", "#b93d1c"]
# One palette per variable, the same in every figure. Site and Buruli status share
# a figure more than once, so they cannot share colours: the beta diversity figure
# coloured a panel by site and the next by group, with the same blue for both.
GROUP_COLOUR = {"BU": "#2a78d6", "non-BU": "#eb6834"}
SITE_COLOUR = {"ayos": "#1baf7a", "akonolinga": "#e87ba4"}
INK = "#0b0b0b"
MUTED = "#52514e"
SEED = 20261111

# Figures are sized for a Word page, 6.5 inches of text width. A figure wider than
# that is scaled down when it is placed, and the text shrinks with it, so the widths
# below stay near the page and the type is set large enough to survive the scaling.
PAGE = 6.9
WIDE = 9.2

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.labelsize": 10, "axes.labelweight": "normal",
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 9,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.8,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "figure.dpi": 400,
})


def site_dir(site):
    return OUT_ROOT / f"{site}_set{SITES[site]}"


REPEAT_RUN = re.compile(r"_\d+$")
CLINICAL = ("group", "ub_status", "lesion_size", "lesion_category", "treatment")


def patient_id(sample):
    """The patient behind a sequencing run.

    The sample sheet names a repeat run of an already listed patient by adding a
    numeric suffix, AY21-F25_2 beside AY21-F25. The two share one lesion and one
    clinical record, so they are one patient measured twice, not two patients.
    """
    return REPEAT_RUN.sub("", sample)


def load_metadata(site, patients_only=True, merge_repeats=True):
    """Metadata indexed by sample, restricted to what was classified.

    With merge_repeats the index is the patient, which is what every statistical
    comparison needs; without it the index is the sequencing run, which is what
    the read level quality control needs.
    """
    frame = pd.read_csv(site_dir(site) / "metadata.tsv", sep="\t", dtype=str)
    frame = frame[frame["barcode"].notna()]
    names = pd.read_csv(site_dir(site) / "barcode_sample_map.tsv", sep="\t", dtype=str)
    frame = frame[frame["barcode"].isin(set(names["barcode"]))]
    if patients_only:
        frame = frame[frame["sample_type"] == "patient"]
    if not merge_repeats:
        return frame.set_index("sample_id").sort_index()
    frame = frame.copy()
    frame["patient"] = frame["sample_id"].map(patient_id)
    for patient, rows in frame.groupby("patient"):
        for column in CLINICAL:
            if rows[column].nunique(dropna=False) > 1:
                raise SystemExit(f"{site}: runs of {patient} disagree on {column}")
    merged = frame.groupby("patient", sort=True).agg(
        {**{column: "first" for column in frame.columns if column != "patient"},
         "barcode": lambda values: ",".join(sorted(values))})
    merged["n_runs"] = frame.groupby("patient").size()
    merged = merged.drop(columns=["sample_id"])
    merged.index.name = "sample_id"
    return merged


def load_table(site, level, relative=False, merge_repeats=True):
    """Abundance table, with the repeat runs of one patient added together.

    Relative abundance is recomputed after the merge, because adding two shares
    of two different run totals is not a share of anything.
    """
    path = site_dir(site) / f"abundance_tables/bacteria_{level}_counts.tsv"
    frame = pd.read_csv(path, sep="\t", index_col=0)
    if merge_repeats:
        frame = frame.T.groupby(patient_id, sort=True).sum().T
    frame = frame[sorted(frame.columns)]
    if relative:
        totals = frame.sum(axis=0)
        frame = frame.div(totals.where(totals > 0, 1), axis=1)
    return frame


def aligned(site, level, relative=False, patients_only=True):
    """Abundance table and metadata restricted to the samples present in both."""
    table = load_table(site, level, relative)
    meta = load_metadata(site, patients_only)
    shared = sorted(set(table.columns) & set(meta.index))
    return table[shared], meta.loc[shared]


GRAM_POSITIVE_PHYLA = {"Actinomycetota", "Bacillota", "Deinococcota"}
GRAM_NEGATIVE_CLASSES = {"Negativicutes"}
GRAM_UNASSIGNED_PHYLA = {"Mycoplasmatota", "Chloroflexota", "Thermomicrobiota"}
# Deliberately not the group colours: Gram status and Buruli status share a figure,
# and the same blue meaning two things in one figure is unreadable.
GRAM_COLOUR = {"Gram negative": "#6250d6", "Gram positive": "#eda100",
               "unassigned": "#9a9890"}


def load_lineage():
    return pd.read_csv(OUT_ROOT / "bacteria_lineage.tsv", sep="\t", dtype=str).fillna("")


def gram_status(lineage=None):
    """One row per genus: Gram positive, Gram negative or unassigned.

    Assigned from the lineage, because SILVA carries no Gram field. The phyla above
    stain positive, except the class Negativicutes, which has a second membrane.
    """
    lineage = load_lineage() if lineage is None else lineage
    out = {}
    for _, row in lineage.iterrows():
        phylum, klass = row["phylum"], row["class"]
        if phylum in GRAM_UNASSIGNED_PHYLA or not phylum:
            status = "unassigned"
        elif phylum in GRAM_POSITIVE_PHYLA:
            status = "Gram negative" if klass in GRAM_NEGATIVE_CLASSES else "Gram positive"
        else:
            status = "Gram negative"
        out[row["genus"]] = status
    return pd.Series(out, name="gram")


def shannon(counts):
    p = counts[counts > 0]
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def simpson(counts):
    p = counts[counts > 0]
    p = p / p.sum()
    return float(1 - (p ** 2).sum())


def chao1(counts):
    """Bias corrected Chao1, as in Chao 1987."""
    counts = np.asarray(counts)
    counts = counts[counts > 0]
    observed = len(counts)
    singletons = int((counts == 1).sum())
    doubletons = int((counts == 2).sum())
    return float(observed + singletons * (singletons - 1) / (2 * (doubletons + 1)))


def rarefy(counts, depth, rng, repeats=100):
    """Mean richness, Shannon, Simpson and Chao1 over repeated rarefactions."""
    values = np.asarray(counts, dtype=np.int64)
    total = values.sum()
    if total < depth:
        return (np.nan,) * 4
    probabilities = values / total
    out = []
    for _ in range(repeats):
        drawn = rng.multinomial(depth, probabilities)
        kept = drawn[drawn > 0]
        out.append((len(kept), shannon(kept), simpson(kept), chao1(drawn)))
    return tuple(np.mean(out, axis=0))


def bray_curtis(table):
    """Pairwise Bray-Curtis over the columns of a taxa by sample table."""
    values = table.to_numpy(dtype=float).T
    n = values.shape[0]
    distance = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = values[i], values[j]
            total = a.sum() + b.sum()
            d = np.abs(a - b).sum() / total if total else 0.0
            distance[i, j] = distance[j, i] = d
    return distance


def pcoa(distance):
    n = distance.shape[0]
    centring = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centring @ (distance ** 2) @ centring
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    positive = values > 1e-12
    coords = vectors[:, positive] * np.sqrt(values[positive])
    return coords, values[positive] / values[positive].sum()


def permanova(distance, groups, rng, permutations=9999):
    """One way PERMANOVA, Anderson 2001, with a free permutation null."""
    groups = np.asarray(groups)
    n = len(groups)
    labels = np.unique(groups)
    if len(labels) < 2 or n <= len(labels):
        return np.nan, np.nan

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
    return float(observed), float((1 + (null >= observed).sum()) / (permutations + 1))


def benjamini_hochberg(pvalues):
    """Adjusted p values, monotone, same order as the input."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    adjusted = np.empty(n)
    running = 1.0
    for rank, index in zip(range(n, 0, -1), order[::-1]):
        running = min(running, p[index] * n / rank)
        adjusted[index] = running
    return np.minimum(adjusted, 1.0)


def tidy_axes(ax, hide=("top", "right")):
    ax.tick_params(length=3, width=0.9)
    for side in hide:
        ax.spines[side].set_visible(False)


def panel_title(ax, text, size=11, **kwargs):
    ax.set_title(text, loc="left", fontsize=size, fontweight="bold", color=INK, **kwargs)


def side_legend(ax, handles=None, labels=None, title=None, **kwargs):
    """A key placed outside the axes, on the right, so it never covers the data."""
    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                       frameon=False, labelcolor=INK, title=title,
                       handlelength=1.4, borderaxespad=0.0, **kwargs)
    if title:
        legend.get_title().set_fontsize(9)
        legend.get_title().set_fontweight("bold")
    return legend


def figure_legend(fig, handles, labels, title=None, **kwargs):
    """One key for a whole figure, outside it on the right."""
    legend = fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
                        frameon=False, labelcolor=INK, title=title,
                        handlelength=1.4, borderaxespad=0.0, **kwargs)
    if title:
        legend.get_title().set_fontsize(9)
        legend.get_title().set_fontweight("bold")
    return legend


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=400, bbox_inches="tight",
                    pad_inches=0.08)
    plt.close(fig)
