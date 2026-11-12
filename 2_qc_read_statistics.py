#!/usr/bin/env python3
"""
Read statistics, point 2 of the analysis plan.

Total reads per barcode before and after filtering, the share retained, and the read
length distribution measured by FastQC on the filtered files. Barcodes below the
inclusion cutoff are drawn, and marked, rather than hidden.

Usage: python3 2_qc_read_statistics.py

Auteur : Khadim Gueye
"""

import zipfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import common as c

CUTOFF_FRACTION = 0.001


def length_profile(site):
    """Sequence length distribution per barcode, read out of the FastQC archives."""
    out = {}
    for archive in sorted((c.site_dir(site) / "filtered_fastqc").glob("*_fastqc.zip")):
        barcode = archive.name.split(".filtered")[0]
        with zipfile.ZipFile(archive) as zf:
            name = next(n for n in zf.namelist() if n.endswith("fastqc_data.txt"))
            inside, rows = False, []
            for raw in zf.read(name).decode().splitlines():
                if raw.startswith(">>Sequence Length Distribution"):
                    inside = True
                    continue
                if inside and raw.startswith(">>END_MODULE"):
                    break
                if inside and not raw.startswith("#"):
                    length, count = raw.split("\t")
                    middle = np.mean([float(v) for v in length.split("-")])
                    rows.append((middle, float(count)))
        if rows:
            out[barcode] = pd.Series(dict(rows)).sort_index()
    return out


def main():
    frames, profiles = [], {}
    for site in c.SITES:
        counts = pd.read_csv(c.site_dir(site) / "read_counts.tsv", sep="\t")
        meta = pd.read_csv(c.site_dir(site) / "metadata.tsv", sep="\t", dtype=str)
        counts = counts.merge(meta[["barcode", "sample_id", "sample_type", "group"]],
                              on="barcode", how="left")
        counts["site"] = c.LABEL[site]
        counts["label"] = counts["sample_id"].fillna(counts["barcode"])
        counts["retained_percent"] = 100 * counts["filtered_reads"] / counts["raw_reads"]
        counts["cutoff"] = CUTOFF_FRACTION * counts["filtered_reads"].sum()
        frames.append(counts)
        profiles[site] = length_profile(site)
    table = pd.concat(frames, ignore_index=True)
    table.to_csv(c.FIG_DIR / "Table0_read_statistics.tsv", sep="\t", index=False,
                 float_format="%.3f")

    fig = plt.figure(figsize=(6.9, 6.6))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.55, 1.0], width_ratios=[1.0, 1.0],
                            hspace=0.58, wspace=1.05)

    # A barcode with no row in the sample sheet has no identity to show, so it stays
    # in Table0 and out of the figure.
    order = table[table["sample_id"].notna()].copy()
    unnamed = int(table["sample_id"].isna().sum())
    order["site_rank"] = order["site"].map({c.LABEL[s]: i for i, s in enumerate(c.SITES)})
    order = order.sort_values(["site_rank", "label"]).reset_index(drop=True)
    xs = np.arange(len(order))

    ax = fig.add_subplot(grid[0, :])
    ax.bar(xs, order["raw_reads"] / 1000, width=0.78, color="#cde2fb",
           edgecolor="white", linewidth=0.5, label="Before filtering")
    ax.bar(xs, order["filtered_reads"] / 1000, width=0.78, color=c.CATEGORICAL[0],
           edgecolor="white", linewidth=0.5, label="After filtering")
    ax.set_yscale("log")
    floor = max(order["filtered_reads"].replace(0, np.nan).min() / 1000 * 0.5, 5e-4)
    ax.set_ylim(floor, order["raw_reads"].max() / 1000 * 2.6)
    for site, sub in order.groupby("site", sort=False):
        ax.hlines(sub["cutoff"].iloc[0] / 1000, sub.index.min() - 0.5,
                  sub.index.max() + 0.5, color=c.CATEGORICAL[1], linewidth=1.3,
                  linestyle=(0, (4, 2)), zorder=4,
                  label="0.1 percent cutoff" if site == order["site"].iloc[0] else None)
    dropped = order.index[order["kept"] == "no"]
    ax.scatter(dropped, np.full(len(dropped), floor * 1.35), marker="x", s=40,
               linewidth=1.8, color=c.CATEGORICAL[1], zorder=5,
               label="Below the cutoff")
    ax.set_ylabel("Reads, thousands, log scale", color=c.INK)
    ax.set_xlabel("Sample", color=c.INK, labelpad=6)
    ax.set_xticks(xs)
    ax.set_xticklabels(order["label"], rotation=90, fontsize=6.5)
    ax.set_xlim(-0.8, len(order) - 0.2)
    c.panel_title(ax, "A  Reads per sample", pad=18)
    c.side_legend(ax)
    c.tidy_axes(ax)
    for name, sub in order.groupby("site", sort=False):
        ax.text(sub.index.to_numpy().mean(), 0.98, name, ha="center", va="top",
                fontsize=10.5, fontweight="bold", color=c.MUTED,
                transform=ax.get_xaxis_transform())


    ax = fig.add_subplot(grid[1, 0])
    rng = np.random.default_rng(c.SEED)
    for i, site in enumerate(c.SITES):
        values = table.loc[table["site"] == c.LABEL[site], "retained_percent"].dropna()
        parts = ax.boxplot(values, positions=[i], widths=0.45, patch_artist=True,
                           showfliers=False)
        parts["boxes"][0].set(facecolor=c.SITE_COLOUR[site], alpha=0.35,
                              edgecolor=c.MUTED, linewidth=0.8)
        for element in ("whiskers", "caps", "medians"):
            for item in parts[element]:
                item.set(color=c.MUTED, linewidth=0.8)
        ax.scatter(rng.normal(i, 0.055, len(values)), values, s=26,
                   color=c.SITE_COLOUR[site], edgecolor="white", linewidth=0.6, zorder=3)
    ax.set_xlim(-0.6, len(c.SITES) - 0.4)
    ax.set_xticks(range(len(c.SITES)))
    ax.set_xticklabels([c.LABEL[s] for s in c.SITES])
    ax.set_ylabel("Percent of reads retained", color=c.INK)
    ax.set_xlabel("Site", color=c.INK)
    c.panel_title(ax, "B  Reads retained by filtering", pad=12)
    handles = [plt.Line2D([], [], marker="o", linestyle="none", markersize=6,
                          color=c.SITE_COLOUR[site], label=c.LABEL[site])
               for site in c.SITES]
    c.side_legend(ax, handles=handles, labels=[h.get_label() for h in handles],
                  title="Site")
    c.tidy_axes(ax)

    ax = fig.add_subplot(grid[1, 1])
    for site in c.SITES:
        for series in profiles[site].values():
            ax.plot(series.index, series / series.sum(), color=c.SITE_COLOUR[site],
                    linewidth=0.9, alpha=0.35)
        ax.plot([], [], color=c.SITE_COLOUR[site], linewidth=2.4, label=c.LABEL[site])
    ax.set_xlim(1000, 2000)
    ax.set_xlabel("Read length, bp", color=c.INK)
    ax.set_ylabel("Share of reads", color=c.INK)
    c.panel_title(ax, "C  Read length after filtering", pad=12)
    c.side_legend(ax, title="Site")
    c.tidy_axes(ax)

    c.save(fig, "Figure1_read_statistics")
    kept = table[table["kept"] == "yes"]
    print(f"{len(table)} barcodes, {len(kept)} above the cutoff, "
          f"{unnamed} with no sample sheet row left out of the figure")
    print(f"median retained {table['retained_percent'].median():.1f} percent")
    print(kept.groupby(["site", "group"]).size().to_string())


if __name__ == "__main__":
    main()
