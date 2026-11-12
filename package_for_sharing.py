#!/usr/bin/env python3
"""
Pack the figures and tables a co-author receives.

Collects everything in the figures folder, writes a short note describing what is in
the archive and what the results say, and zips the lot.

Run the analysis first.

Usage: python3 package_for_sharing.py

Auteur : Khadim Gueye
"""

import textwrap
import zipfile
from datetime import date

import pandas as pd

import common as c

ARCHIVE = c.OUT_ROOT / f"buruli_16s_analysis_{date.today():%Y-%m-%d}.zip"

NOTE = """Buruli ulcer wound microbiome, full length 16S rRNA, Ayos and Akonolinga
Buruli against non Buruli, following the project analysis plan
{today}
Khadim Gueye

WHY THE EARLIER FIGURES WERE ARCHAEAL

The Kraken2 database used in the first run contained no bacterial sequences. Its
bacterial download had failed and the rerun skipped index construction because a hash
table was already on disk, so the reports could only name archaea and 96 percent of
reads stayed unassigned. The samples were never the problem.

The database has been rebuilt from SILVA 138.2 SSU Ref NR99, the reference designed
for full length 16S amplicons. Classification goes from 0 percent bacterial to 100
percent bacterial, with no archaeal read left. The archaeal taxa in the current draft
should be removed, not corrected. See FigureS1 and TableS1.

SILVA resolves to genus, not species, so everything here stops at genus. Species level
statements from the first run cannot be carried over.

SAMPLES

The analysis plan sets the inclusion cutoff at 0.1 percent of the run read count, so
that is what was used. It is more permissive than the 10,000 read cutoff of the first
run and admits 8 more samples. 60 of 64 barcodes pass.

{samples}

The sample sheet lists five Ayos patients twice, on two barcodes, AY06-J25 beside
AY06-J25_2 and so on, with the same clinical record both times. Those are one patient
sequenced twice, so their counts are added together and they count once. Keeping them
apart would have inflated the Ayos non Buruli group from 15 to 20 and broken the
independence the tests assume.

Both water controls and the positive control fall below the cutoff, which is what a
clean run should look like.

Two barcodes pass the cutoff but have no row in the sample sheet, barcode32 at Ayos
and barcode40 at Akonolinga. They are left out. Please tell us what they are.

WHAT THE RESULTS SAY

{results}

FILES

figures/   the {n_figures} figures, PNG and PDF
tables/    every table behind them, tab separated
README.txt this note

Figure0  reads per barcode, share retained, read length
Figure1  phylum composition per sample and per group
Figure2  the twenty most abundant genera
Figure3  alpha diversity, observed genera, Shannon, Simpson, Chao1
Figure4  beta diversity, principal coordinates
Figure5  differential abundance, Buruli against non Buruli
Figure6  Gram composition
Figure7  the clinical variables against the community
Figure8  genus co-occurrence networks
Figure9  predicted pathways, Buruli against non Buruli
FigureS1 classification before and after the database was rebuilt
"""


def summary():
    lines = []
    beta = pd.read_csv(c.FIG_DIR / "Table6_beta_diversity_permanova.tsv", sep="\t")
    site = beta[beta["factor"] == "site"]
    if not site.empty:
        lines.append(
            f"The two sites carry different communities, PERMANOVA F = "
            f"{site['F'].iloc[0]:.2f}, p = {site['p_value'].iloc[0]:.4f}. Buruli status "
            f"does not separate them at either site, and neither does lesion size, "
            f"lesion category or treatment status.")
    alpha = pd.read_csv(c.FIG_DIR / "Table5_alpha_diversity_tests.tsv", sep="\t")
    group = alpha[alpha["comparison"].str.contains("BU vs non-BU")]
    if not group.empty:
        lines.append(
            f"No diversity index separates Buruli from non Buruli at either site; the "
            f"smallest adjusted p value is {group['p_adjusted_bh'].min():.2f}.")
    diff = pd.read_csv(c.FIG_DIR / "Table8_differential_abundance.tsv", sep="\t")
    for site_name, sub in diff.groupby("site"):
        hits = sub[sub["significant"]]
        if hits.empty:
            lines.append(f"At {site_name} no genus of the {len(sub)} tested survives "
                         f"correction.")
        else:
            named = ", ".join(f"{r['genus']} (log2 fold change "
                              f"{r['log2_fold_change']:.2f}, adjusted p "
                              f"{r['p_adjusted_bh']:.4f})" for _, r in hits.iterrows())
            word = "genus" if len(hits) == 1 else "genera"
            verb = "survives" if len(hits) == 1 else "survive"
            lines.append(f"At {site_name} {len(hits)} {word} of the {len(sub)} tested "
                         f"{verb} correction: {named}.")
    gram = pd.read_csv(c.FIG_DIR / "Table9_gram_composition.tsv", sep="\t")
    negative = gram[gram["gram"] == "Gram negative"]["percent"].mean()
    lines.append(f"The communities are overwhelmingly Gram negative, "
                 f"{negative:.0f} percent of bacterial reads on average, with no "
                 f"difference between the groups.")
    pathways = c.FIG_DIR / "Table16_pathway_differential.tsv"
    if pathways.exists():
        table = pd.read_csv(pathways, sep="\t")
        tested = ", ".join(f"{len(sub)} at {name}"
                           for name, sub in table.groupby("site", sort=False))
        best = table.loc[table["p_adjusted_bh"].idxmin()]
        lines.append(
            f"No predicted MetaCyc pathway survives correction either ({tested} "
            f"tested). The lowest adjusted p value, {best['p_adjusted_bh']:.3f}, is "
            f"{best['description'] or best['pathway']} at {best['site']}, lower in "
            f"Buruli lesions. That pathway and the three behind it are Corynebacterium "
            f"cell wall and mycothiol chemistry, tracking the Corynebacterium share at "
            f"Spearman rho 0.99, so they are the same finding as the genus above "
            f"rather than a second one.")
    edges = c.FIG_DIR / "Table14_cooccurrence_edges.tsv"
    if edges.exists():
        table = pd.read_csv(edges, sep="\t")
        kept = ", ".join(f"{int(sub['edge'].sum())} at {name}"
                         for name, sub in table.groupby("site", sort=False))
        lines.append(
            f"The co-occurrence networks are built around a dense Enterobacteriaceae "
            f"cluster that correlates negatively with Alcaligenes and the anaerobes "
            f"({kept} of {len(table) // table['site'].nunique()} genus pairs tested "
            f"per site).")
    lines.append("Mycobacterium is effectively absent, a single read across all "
                 "samples. Full length 16S does not recover M. ulcerans from these "
                 "swabs, so the pathogen cannot be quantified from this data.")
    lines.append("The main caveat is depth. The 0.1 percent cutoff admits samples with "
                 "around 3,000 reads, so rarefaction runs at 3,300 bacterial reads. "
                 "That is enough for genus level composition but it limits the power "
                 "of the diversity comparisons.")
    return "\n\n".join(textwrap.fill(line, width=84) for line in lines)


def main():
    counts = {}
    for site in c.SITES:
        meta = c.load_metadata(site)
        counts[site] = (len(meta), int((meta["group"] == "BU").sum()),
                        int((meta["group"] == "non-BU").sum()))

    samples = textwrap.fill(
        f"After removing controls, 58 patient runs remain, which are "
        f"{counts['ayos'][0] + counts['akonolinga'][0]} patients: "
        f"{counts['ayos'][0]} at Ayos ({counts['ayos'][1]} Buruli, "
        f"{counts['ayos'][2]} non Buruli) and {counts['akonolinga'][0]} at Akonolinga "
        f"({counts['akonolinga'][1]} Buruli, {counts['akonolinga'][2]} non Buruli).",
        width=79, subsequent_indent="")

    figures = sorted(c.FIG_DIR.glob("Figure*.png")) + sorted(c.FIG_DIR.glob("Figure*.pdf"))

    note = NOTE.format(
        today=f"{date.today():%d %B %Y}", samples=samples, results=summary(),
        n_figures=len({path.stem for path in figures}))
    tables = sorted(c.FIG_DIR.glob("Table*.tsv"))
    if not figures or not tables:
        raise SystemExit("no figures or tables found, run the analysis first")

    with zipfile.ZipFile(ARCHIVE, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("README.txt", note)
        for path in figures:
            zf.write(path, f"figures/{path.name}")
        for path in tables:
            zf.write(path, f"tables/{path.name}")
        lineage = c.OUT_ROOT / "bacteria_lineage.tsv"
        if lineage.exists():
            zf.write(lineage, f"tables/{lineage.name}")

    print(f"{ARCHIVE}  ({ARCHIVE.stat().st_size / 1024:.0f} KB)")
    print(f"  {len(figures)} figure files, {len(tables)} tables")


if __name__ == "__main__":
    main()
