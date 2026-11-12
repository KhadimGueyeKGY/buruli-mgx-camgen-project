#!/usr/bin/env python3
"""
Cross-check the barcode to sample_id mapping against the archived read counts.

The sample sheet is the source of the mapping. Before it was supplied, the mapping had
been reconstructed by matching filtered read counts against the archived QC table,
which was possible because the filter settings here reproduce the first run exactly.
That reconstruction is kept as a check: for every sample the first run carried, the
barcode whose read count matches must be the barcode the sample sheet names.

A disagreement means either the filter settings drifted or the sample sheet is wrong,
and both are worth knowing before any result is read.

Run 1c_prepare_metadata.py first.

Usage: python3 1b_recover_sample_names.py <site> <set_letter>

Auteur : Khadim Gueye
"""

import sys
from collections import Counter
from pathlib import Path

ROOT = Path.home() / "khadim/Francis-job"
ARCHIVE = ROOT / "data/16s_metagenomics/plots_updates"
OUT_ROOT = ROOT / "data/16s_metagenomics/analysis_2026"
ARCHIVE_DIR = {"ayos": "plot_ayos_setA", "akonolinga": "plot_akonolinga_setB"}


def archived(site):
    """sample_id to the FastQC read count of the first run."""
    path = ARCHIVE / ARCHIVE_DIR[site] / "qc_filtered_pass_cutoff.tsv"
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        column = header.index("FastQC_mqc.generalstats.fastqc.total_sequences")
        return {f[0].replace(".filtered", ""): int(float(f[column]))
                for f in (line.rstrip("\n").split("\t") for line in fh)}


def observed(out_dir):
    """barcode to the filtered read count measured here."""
    with open(out_dir / "read_counts.tsv") as fh:
        fh.readline()
        return {f[0]: int(f[2]) for f in (line.rstrip("\n").split("\t") for line in fh)}


def official(out_dir):
    """barcode to sample_id, as the sample sheet states it."""
    with open(out_dir / "metadata.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        b, s = header.index("barcode"), header.index("sample_id")
        rows = [line.rstrip("\n").split("\t") for line in fh]
    return {r[b]: r[s] for r in rows if r[b]}


def main():
    site, letter = sys.argv[1], sys.argv[2]
    out_dir = OUT_ROOT / f"{site}_set{letter}"
    want = archived(site)
    have = observed(out_dir)
    sheet = official(out_dir)

    repeated = [n for n, count in Counter(have.values()).items() if count > 1]
    by_count = {n: b for b, n in have.items() if n not in repeated}

    agree, disagree, unmatched = [], [], []
    for sample, count in want.items():
        barcode = by_count.get(count)
        if barcode is None:
            unmatched.append(sample)
        elif sheet.get(barcode) == sample:
            agree.append(sample)
        else:
            disagree.append((sample, barcode, sheet.get(barcode)))

    print(f"{site}: {len(agree)} of {len(want)} archived samples confirm the sample sheet")
    for sample, barcode, named in disagree:
        print(f"  DISAGREEMENT {sample} reproduces the read count of {barcode}, "
              f"which the sheet calls {named}")
    for sample in unmatched:
        print(f"  no barcode reproduces the archived read count of {sample}")
    if repeated:
        print(f"  read counts shared by more than one barcode, not usable as a key: {repeated}")
    return 1 if disagree else 0


if __name__ == "__main__":
    sys.exit(main())
