#!/usr/bin/env python3
"""
Recover the barcode to sample_id mapping lost with the metadata sheets.

The archived QC table holds the FastQC read count of every sample that passed the
10,000 read cutoff, measured on the filtered FASTQ. Rerunning the same NanoFilt
settings reproduces those counts exactly, so each sample maps back to one barcode.
The mapping is written only when every count is unique on both sides.

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
    path = ARCHIVE / ARCHIVE_DIR[site] / "qc_filtered_pass_cutoff.tsv"
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        column = header.index("FastQC_mqc.generalstats.fastqc.total_sequences")
        return {f[0].replace(".filtered", ""): int(float(f[column]))
                for f in (line.rstrip("\n").split("\t") for line in fh)}


def observed(out_dir):
    with open(out_dir / "read_counts.tsv") as fh:
        fh.readline()
        return {f[0]: int(f[2])
                for f in (line.rstrip("\n").split("\t") for line in fh) if f[3] == "yes"}


def main():
    site, letter = sys.argv[1], sys.argv[2]
    out_dir = OUT_ROOT / f"{site}_set{letter}"
    want = archived(site)
    have = observed(out_dir)

    for label, table in (("archived", want), ("observed", have)):
        repeated = [n for n, c in Counter(table.values()).items() if c > 1]
        if repeated:
            print(f"{site}: {label} read counts are not unique: {repeated}")
            return 1

    by_count = {n: b for b, n in have.items()}
    mapping, missing = {}, []
    for sample, count in want.items():
        if count in by_count:
            mapping[by_count[count]] = sample
        else:
            missing.append(sample)

    path = out_dir / "barcode_sample_map.tsv"
    with open(path, "w") as fh:
        fh.write("barcode\tsample_id\tfiltered_reads\n")
        for barcode in sorted(mapping):
            fh.write(f"{barcode}\t{mapping[barcode]}\t{have[barcode]}\n")

    print(f"{site}: {len(mapping)} of {len(want)} archived samples matched, written to {path}")
    for sample in missing:
        print(f"  no barcode reproduces the archived count of {sample}")
    extra = sorted(set(have) - set(mapping))
    if extra:
        print(f"  barcodes kept here but absent from the archived set: {', '.join(extra)}")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
