#!/usr/bin/env python3
"""
Functional profiling, point 7 of the analysis plan.

PICRUSt2 predicts gene content from 16S, but it needs study sequences and a count
table, which Kraken2 and Bracken do not produce. Clustering the reads themselves
does not work here: after quality filtering a Nanopore 16S read still carries a few
percent of error, so full length reads almost never repeat and a de novo operational
taxonomic unit is built out of noise. On this run 282158 reads gave 282119 distinct
sequences, and clustering kept 35 units covering 7 percent of the data.

This script takes the other route. Bracken has already placed every read on a genus
against SILVA, so each detected genus is represented by its own SILVA reference
sequence and weighted by the Bracken count. The sequences PICRUSt2 places are then
reference quality rather than raw reads, and the abundances are the ones used
everywhere else in this analysis.

The predictions remain inferences from marker gene content, not measured gene
abundances, so the output supports comparisons between groups rather than absolute
statements about pathways.

Usage: python3 7_functional_profiling.py
       [--rebuild-reference] [--skip-picrust] [--analyse-only]

Auteur : Khadim Gueye
"""

import argparse
import os
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import common as c

WORK = c.OUT_ROOT / "functional"
SILVA = Path.home() / "DATABASE/kraken2_silva138.2_16S/library/silva.fna"
REFERENCE = WORK / "genus_reference.fasta"
PROVENANCE = WORK / "genus_reference_provenance.tsv"
TABLE = WORK / "genus_table.tsv"
THREADS = 48
MIN_LENGTH = 900
FULL_LENGTH = (1200, 2000)


def run(command, **kwargs):
    print("  " + " ".join(str(part) for part in command))
    subprocess.run(command, check=True, **kwargs)


def identifier(genus, taken):
    """A genus name PICRUSt2 can carry through a FASTA, an alignment and a tree."""
    base = re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]", "_", genus)).strip("_") or "genus"
    name, suffix = base, 2
    while name in taken:
        name, suffix = f"{base}_{suffix}", suffix + 1
    taken.add(name)
    return name


def wanted_genera():
    """Every genus in the abundance tables, keyed the way SILVA writes it.

    A genus label was disambiguated by its family where SILVA reuses a name, so both
    the plain label and the label without the family prefix are accepted. Some SILVA
    genus names carry their family themselves, Lachnospiraceae UCG-004 for one.
    """
    lineage = pd.read_csv(c.OUT_ROOT / "bacteria_lineage.tsv", sep="\t",
                          dtype=str).fillna("")
    keys = {}
    for _, row in lineage.iterrows():
        family, label = row["family"], row["genus"]
        forms = [label]
        if family and label.startswith(family + " "):
            forms.append(label[len(family) + 1:])
        for form in forms:
            held = keys.setdefault((family, form), label)
            if held != label:
                raise SystemExit(f"{family} {form} is claimed by {held} and {label}")
    return keys


def rank(sequence):
    """Prefer an unambiguous sequence, then a full length one, then a long one."""
    unambiguous = not set(sequence) - set("ACGT")
    full = FULL_LENGTH[0] <= len(sequence) <= FULL_LENGTH[1]
    return (unambiguous + full, len(sequence))


def build_reference():
    """One SILVA sequence per detected genus, chosen once and cached."""
    keys = wanted_genera()
    best = {}

    def consider(header, parts):
        path = header.split(" ", 1)[1].split(";") if " " in header else []
        if len(path) < 7 or path[0] != "Bacteria":
            return
        key = (path[4], path[5])
        if key not in keys:
            return
        sequence = "".join(parts)
        if len(sequence) < MIN_LENGTH:
            return
        genus = keys[key]
        current = best.get(genus)
        if current is None or (rank(sequence), header) > (rank(current[1]), current[0]):
            best[genus] = (header, sequence)

    header, parts = None, []
    with open(SILVA) as handle:
        for line in handle:
            if line.startswith(">"):
                if header:
                    consider(header, parts)
                header, parts = line[1:].rstrip(), []
            else:
                parts.append(line.strip())
    if header:
        consider(header, parts)

    taken, rows = set(), []
    with open(REFERENCE, "w") as out:
        for genus in sorted(best):
            name, sequence = best[genus]
            key = identifier(genus, taken)
            out.write(f">{key}\n")
            for start in range(0, len(sequence), 70):
                out.write(sequence[start:start + 70] + "\n")
            rows.append({"genus": genus, "reference_id": key,
                         "silva_accession": name.split(" ", 1)[0],
                         "length": len(sequence),
                         "silva_taxonomy": name.split(" ", 1)[1]})
    pd.DataFrame(rows).to_csv(PROVENANCE, sep="\t", index=False)
    print(f"{len(rows)} genera represented, {len(keys)} lookup keys")
    return pd.DataFrame(rows)


def build_table(provenance):
    """Bracken genus counts for both sites, one column per patient."""
    frames, named = [], set()
    for site in c.SITES:
        counts, _ = c.aligned(site, "genus")
        clash = named & set(counts.columns)
        if clash:
            raise SystemExit(f"the two sites share a sample name: {sorted(clash)}")
        named |= set(counts.columns)
        frames.append(counts)
    table = pd.concat(frames, axis=1).fillna(0)
    names = dict(zip(provenance["genus"], provenance["reference_id"]))
    covered = table.loc[[genus for genus in table.index if genus in names]]
    share = 100 * covered.to_numpy().sum() / max(table.to_numpy().sum(), 1)
    covered.index = [names[genus] for genus in covered.index]
    covered = covered.loc[covered.sum(axis=1) > 0]
    covered.index.name = "#OTU ID"
    covered.astype(int).to_csv(TABLE, sep="\t")
    print(f"{covered.shape[0]} genera over {covered.shape[1]} patients, "
          f"{share:.3f} percent of the bacterial reads represented")
    return TABLE


def picrust():
    out = WORK / "picrust2_out"
    if out.exists():
        shutil.rmtree(out)
    binary = shutil.which("picrust2_pipeline.py")
    if binary is None:
        env = Path.home() / "miniconda3/envs/picrust2/bin/picrust2_pipeline.py"
        binary = str(env) if env.exists() else None
    if binary is None:
        print("picrust2_pipeline.py not found, the table and references are in", WORK)
        return None
    # The pipeline calls place_seqs.py and its other steps by bare name, so its own
    # bin directory has to be on PATH even when the pipeline itself is called by path.
    env = dict(os.environ)
    env["PATH"] = str(Path(binary).parent) + os.pathsep + env.get("PATH", "")
    run([binary, "-s", str(REFERENCE), "-i", str(TABLE), "-o", str(out),
         "-p", str(THREADS), "--in_traits", "EC,KO"], env=env)
    return out


def clr(frame):
    values = frame.to_numpy(dtype=float) + 0.5
    values = values / values.sum(axis=0, keepdims=True)
    log = np.log2(values)
    return pd.DataFrame(log - log.mean(axis=0, keepdims=True),
                        index=frame.index, columns=frame.columns)


def pathway_names():
    """MetaCyc identifier to its description, from the PICRUSt2 installation."""
    binary = shutil.which("picrust2_pipeline.py")
    if binary is None:
        env = Path.home() / "miniconda3/envs/picrust2/bin/picrust2_pipeline.py"
        binary = str(env) if env.exists() else None
    if binary is None:
        return {}
    root = Path(binary).resolve().parent.parent
    matches = list(root.glob("lib/*/site-packages/picrust2/default_files/"
                             "description_mapfiles/metacyc_pathways_info.txt.gz"))
    if not matches:
        return {}
    frame = pd.read_csv(matches[0], sep="\t", header=None, names=["pathway", "description"],
                        compression="gzip", dtype=str)
    return dict(zip(frame["pathway"], frame["description"]))


def shorten(text, limit=30):
    return text if len(text) <= limit else text[:limit - 1].rstrip() + "\u2026"


def pathway_table(out):
    """MetaCyc pathway abundances per sample, whichever way PICRUSt2 named the file."""
    candidates = sorted((out / "pathways_out").glob("path_abun_unstrat*"))
    if not candidates:
        return pd.DataFrame()
    path = candidates[0]
    return pd.read_csv(path, sep="\t", index_col=0,
                       compression="gzip" if path.suffix == ".gz" else None)


def compare_pathways(abundance, names):
    """Buruli against non Buruli within each site, on centred log ratios."""
    rows, per_site = [], {}
    for site in c.SITES:
        meta = c.load_metadata(site)
        shared = sorted(set(abundance.columns) & set(meta.index))
        if len(shared) < 8:
            continue
        frame = abundance[shared]
        frame = frame.loc[frame.sum(axis=1) > 0]
        transformed = clr(frame)
        bu = [s for s in shared if meta.loc[s, "group"] == "BU"]
        non = [s for s in shared if meta.loc[s, "group"] == "non-BU"]
        if len(bu) < 3 or len(non) < 3:
            continue
        per_site[site] = (frame, transformed, bu, non)
        for pathway in transformed.index:
            a = transformed.loc[pathway, bu].to_numpy()
            b = transformed.loc[pathway, non].to_numpy()
            p = mannwhitneyu(a, b, alternative="two-sided").pvalue
            rows.append({"site": c.LABEL[site], "pathway": pathway,
                         "description": names.get(pathway, ""),
                         "n_BU": len(bu), "n_nonBU": len(non),
                         "mean_BU": frame.loc[pathway, bu].mean(),
                         "mean_nonBU": frame.loc[pathway, non].mean(),
                         "log2_fold_change": a.mean() - b.mean(), "p_value": p})
    results = pd.DataFrame(rows)
    if not results.empty:
        adjusted = []
        for site, sub in results.groupby("site", sort=False):
            adjusted.append(pd.Series(c.benjamini_hochberg(sub["p_value"].to_numpy()),
                                      index=sub.index))
        results["p_adjusted_bh"] = pd.concat(adjusted).sort_index()
        results["significant"] = results["p_adjusted_bh"] < 0.05
    return results, per_site


def figure_pathways(results):
    fig = plt.figure(figsize=(5.4, 9.4))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 2.0], hspace=0.52, wspace=0.46)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]),
            fig.add_subplot(grid[1, :])]
    for k, site in enumerate(c.SITES):
        ax = axes[k]
        frame = results[results["site"] == c.LABEL[site]] if not results.empty \
            else pd.DataFrame()
        if frame.empty:
            ax.text(0.5, 0.5, "too few samples", ha="center", va="center",
                    transform=ax.transAxes, color=c.MUTED)
            c.panel_title(ax, f"{'AB'[k]}  {c.LABEL[site]}")
            c.tidy_axes(ax)
            continue
        y = -np.log10(frame["p_value"].clip(lower=1e-300))
        for mask, colour, size, order, label in (
                (~frame["significant"], c.OTHER, 24, 2, "not significant"),
                (frame["significant"], c.CATEGORICAL[0], 50, 3, "BH p below 0.05")):
            ax.scatter(frame.loc[mask, "log2_fold_change"], y[mask], s=size,
                       color=colour, edgecolor="white", linewidth=0.6, zorder=order,
                       label=label)
        if frame["significant"].any():
            cut = frame.loc[frame["significant"], "p_value"].max()
            ax.axhline(-np.log10(cut), color=c.CATEGORICAL[1], linewidth=0.9,
                       linestyle=(0, (4, 2)), zorder=1)
        ax.axvline(0, color=c.MUTED, linewidth=0.4, zorder=0)
        ax.set_xlabel("log2 fold change, BU over non-BU", color=c.INK)
        ax.set_ylabel("-log10 p value", color=c.INK)
        c.panel_title(ax, f"{'AB'[k]}  {c.LABEL[site]}, "
                          f"{int(frame['significant'].sum())} of {len(frame)} "
                          f"significant", pad=12)
        if k == len(c.SITES) - 1:
            c.side_legend(ax)
        c.tidy_axes(ax)

    ax = axes[2]
    if results.empty:
        ax.axis("off")
    else:
        ranked = results.reindex(results["p_value"].sort_values().index).head(14)
        ranked = ranked.iloc[::-1]
        positions = np.arange(len(ranked))
        colours = [c.GROUP_COLOUR["BU"] if v > 0 else c.GROUP_COLOUR["non-BU"]
                   for v in ranked["log2_fold_change"]]
        ax.barh(positions, ranked["log2_fold_change"], height=0.68, color=colours,
                edgecolor="white", linewidth=0.8)
        ax.set_yticks(positions)
        # Site on one line and the pathway on the next, so the names do not push
        # the figure wider than a page.
        ax.set_yticklabels([f"{s}\n{shorten(d or p)}" for p, d, s in
                            zip(ranked["pathway"], ranked["description"],
                                ranked["site"])], fontsize=8.0,
                           linespacing=1.3)
        ax.set_ylim(-0.7, len(ranked) - 0.3)
        ax.axvline(0, color=c.MUTED, linewidth=0.6)
        ax.set_xlabel("log2 fold change, BU over non-BU", color=c.INK)
        ax.set_ylabel("Site and MetaCyc pathway", color=c.INK, labelpad=6)
        c.panel_title(ax, "C  The fourteen lowest p values", pad=12)
        c.tidy_axes(ax, hide=("top", "right", "left"))
        handles = [plt.Rectangle((0, 0), 1, 1, color=c.GROUP_COLOUR[g],
                                 label=g) for g in c.GROUPS]
        c.side_legend(ax, handles=handles,
                      labels=[h.get_label() for h in handles], title="Higher in")
    c.save(fig, "Figure8_functional_profiling")


def analyse(out):
    abundance = pathway_table(out)
    if abundance.empty:
        print("no pathway table produced, nothing to compare")
        return
    abundance.to_csv(c.FIG_DIR / "Table15_pathway_abundance.tsv", sep="\t",
                     float_format="%.3f")
    names = pathway_names()
    results, _ = compare_pathways(abundance, names)
    results.to_csv(c.FIG_DIR / "Table16_pathway_differential.tsv", sep="\t",
                   index=False, float_format="%.5f")
    figure_pathways(results)
    if not results.empty:
        for site, sub in results.groupby("site"):
            print(f"{site}: {len(sub)} pathways, "
                  f"{int(sub['significant'].sum())} significant after Benjamini-Hochberg")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-reference", action="store_true")
    parser.add_argument("--skip-picrust", action="store_true")
    parser.add_argument("--analyse-only", action="store_true",
                        help="redraw the comparison from a PICRUSt2 run already on disk")
    args = parser.parse_args()

    WORK.mkdir(parents=True, exist_ok=True)
    if args.analyse_only:
        analyse(WORK / "picrust2_out")
        return
    if args.rebuild_reference or not PROVENANCE.exists():
        if not SILVA.exists():
            raise SystemExit(f"SILVA sequences not found at {SILVA}")
        provenance = build_reference()
    else:
        provenance = pd.read_csv(PROVENANCE, sep="\t")
        print(f"{len(provenance)} genus references read from {PROVENANCE.name}")
    build_table(provenance)
    if args.skip_picrust:
        return
    out = picrust()
    if out is not None:
        analyse(out)


if __name__ == "__main__":
    main()
