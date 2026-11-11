#!/usr/bin/env python3
"""
Compare the previous classification with the SILVA one.

The database used in the first run held no bacterial sequences, so the reports
named only archaea and left most reads unassigned. This script reads the archived
reports, matches each block to a sample by its read total, and puts the two runs
side by side.

Usage: python3 2c_classification_comparison.py

Auteur : Khadim Gueye
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path.home() / "khadim/Francis-job"
ARCHIVE = ROOT / "data/16s_metagenomics/plots_updates"
OUT_ROOT = ROOT / "data/16s_metagenomics/analysis_2026"
FIG_DIR = OUT_ROOT / "figures"
SITES = {"ayos": ("A", "plot_ayos_setA"), "akonolinga": ("B", "plot_akonolinga_setB")}
LABEL = {"ayos": "Ayos", "akonolinga": "Akonolinga"}
FRACTIONS = ["Bacteria", "Archaea", "Eukaryota", "Unclassified"]
COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#9a9890"]
INK = "#0b0b0b"
MUTED = "#52514e"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.6,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def archived_blocks(path):
    """Split the concatenated report into one Kraken2 report per sample.

    Each sample contributes its Kraken2 report followed by two Bracken reports,
    which carry no unclassified line. Only the first is wanted, so the block is cut
    at the second root node."""
    blocks, current, roots = [], [], 0
    for line in open(path):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        if fields[3] == "U":
            if current:
                blocks.append(current)
            current, roots = [], 0
        if fields[3] == "R":
            roots += 1
        if roots < 2:
            current.append(fields)
    if current:
        blocks.append(current)
    return blocks


def block_totals(block):
    out = {name: 0 for name in FRACTIONS}
    for fields in block:
        name, reads = fields[5].strip(), int(fields[1])
        if fields[3] == "U":
            out["Unclassified"] += reads
        elif name in ("Bacteria", "Archaea", "Eukaryota"):
            out[name] += reads
    out["total"] = sum(out[name] for name in FRACTIONS)
    return out


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for site, (letter, folder) in SITES.items():
        out_dir = OUT_ROOT / f"{site}_set{letter}"
        mapping = pd.read_csv(out_dir / "barcode_sample_map.tsv", sep="\t")
        by_total = dict(zip(mapping["filtered_reads"], mapping["sample_id"]))

        previous = {}
        for block in archived_blocks(ARCHIVE / folder / "krona_plots/combined_kraken2_report.txt"):
            totals = block_totals(block)
            sample = by_total.get(totals["total"])
            if sample:
                previous[sample] = totals

        current = pd.read_csv(out_dir / "abundance_tables/domain_composition.tsv",
                              sep="\t", index_col=0)
        for sample, totals in previous.items():
            for run, source in (("Previous database", totals),
                                ("SILVA 138.2", {
                                    "Bacteria": current.loc[sample, "Bacteria"],
                                    "Archaea": current.loc[sample, "Archaea"],
                                    "Eukaryota": current.loc[sample, "Eukaryota"],
                                    "Unclassified": current.loc[sample, "unclassified"],
                                    "total": current.loc[sample, "total_reads"]})):
                row = {"site": LABEL[site], "sample_id": sample, "run": run}
                row.update({name: 100 * source[name] / source["total"] for name in FRACTIONS})
                rows.append(row)

    table = pd.DataFrame(rows)
    table.to_csv(FIG_DIR / "TableS1_classification_comparison.tsv", sep="\t",
                 index=False, float_format="%.3f")

    summary = table.groupby(["site", "run"])[FRACTIONS].mean()
    order = [(LABEL[s], r) for s in SITES for r in ("Previous database", "SILVA 138.2")]
    summary = summary.loc[order]

    fig, ax = plt.subplots(figsize=(8.6, 2.9))
    left = np.zeros(len(order))
    for name, colour in zip(FRACTIONS, COLOURS):
        values = summary[name].to_numpy()
        ax.barh(np.arange(len(order)), values, left=left, height=0.56, color=colour,
                edgecolor="white", linewidth=0.6, label=name)
        for i, (v, l) in enumerate(zip(values, left)):
            if v >= 5:
                ax.text(l + v / 2, i, f"{v:.1f}", ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
        left += values
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([f"{site}\n{run}" for site, run in order], fontsize=7.5)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Mean percent of filtered reads", color=INK)
    ax.set_title("Reads assigned before and after the database was rebuilt",
                 loc="left", fontsize=9, fontweight="bold", color=INK)
    ax.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), frameon=False,
              fontsize=7.5, labelcolor=INK)
    ax.invert_yaxis()
    ax.tick_params(length=2)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"FigureS1_classification_comparison.{ext}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print(summary.round(3).to_string())
    print(f"\n{len(table) // 2} samples compared, written to {FIG_DIR}")


if __name__ == "__main__":
    main()
