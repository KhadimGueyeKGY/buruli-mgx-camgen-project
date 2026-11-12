#!/usr/bin/env python3
"""
Turn the sample sheet into one tidy metadata table per site.

The sheet carries patients and controls together, a few rows whose barcode reads
"fail", and the lesion size written with a non standard greater than sign. This
script normalises all of it and marks the water and positive controls, which are
kept out of the group comparisons but classified so contamination stays visible.

It also writes barcode_sample_map.tsv, the barcode to sample_id mapping every later
step reads, restricted to the barcodes that passed the inclusion cutoff.

UB positive is Buruli ulcer, UB negative is non Buruli, as stated in the analysis
plan.

Usage: python3 1c_prepare_metadata.py

Auteur : Khadim Gueye
"""

from pathlib import Path

import pandas as pd

ROOT = Path.home() / "khadim/Francis-job"
SHEET = ROOT / "data/16s_metagenomics/samplesheet_16srna_positive_ub_setAB_final.xlsx"
OUT_ROOT = ROOT / "data/16s_metagenomics/analysis_2026"
SITES = {"fastq_pass_ayos": ("ayos", "A"), "fastq_pass_akonolinga": ("akonolinga", "B")}
COLUMNS = ["site", "sample_id", "barcode", "ub_status",
           "lesion_size", "lesion_category", "treatment"]
CONTROLS = {"NFW": "water control", "F-64": "positive control"}
GROUP = {"positive": "BU", "negative": "non-BU"}
SIZE = {"≤ 5cm": "<= 5 cm", "5-15cm": "5-15 cm", "˃ 15cm": "> 15 cm"}


def tidy(frame, site):
    frame = frame.copy()
    frame.columns = COLUMNS
    for column in ("sample_id", "barcode", "ub_status", "lesion_size",
                   "lesion_category", "treatment"):
        frame[column] = frame[column].astype(str).str.strip().replace({"nan": pd.NA})
    frame["site"] = site
    frame = frame.dropna(subset=["sample_id"]).drop_duplicates(subset=["sample_id", "barcode"])

    frame["sample_type"] = frame["sample_id"].map(CONTROLS).fillna("patient")
    frame["group"] = frame["ub_status"].map(GROUP)
    frame.loc[frame["sample_type"] != "patient", "group"] = pd.NA
    frame["lesion_size"] = frame["lesion_size"].map(SIZE).astype("object")
    frame["sequenced"] = frame["barcode"].ne("fail")
    frame.loc[~frame["sequenced"], "barcode"] = pd.NA
    return frame[["site", "sample_id", "barcode", "sample_type", "group", "ub_status",
                  "lesion_size", "lesion_category", "treatment", "sequenced"]]


def main():
    book = pd.ExcelFile(SHEET)
    for sheet, (site, letter) in SITES.items():
        frame = tidy(book.parse(sheet), site)
        out = OUT_ROOT / f"{site}_set{letter}"
        out.mkdir(parents=True, exist_ok=True)
        path = out / "metadata.tsv"
        frame.to_csv(path, sep="\t", index=False)

        counts = pd.read_csv(out / "read_counts.tsv", sep="\t")
        kept = set(counts.loc[counts["kept"] == "yes", "barcode"])
        analysed = frame[frame["barcode"].isin(kept) & (frame["sample_type"] == "patient")]
        controls = frame[frame["barcode"].isin(kept) & (frame["sample_type"] != "patient")]
        unused = frame[frame["sequenced"] & ~frame["barcode"].isin(kept)]

        passing = frame[frame["barcode"].isin(kept)]
        duplicated = passing["sample_id"][passing["sample_id"].duplicated()].tolist()
        if duplicated:
            raise SystemExit(f"{site}: sample_id used by more than one barcode: {duplicated}")
        mapping = passing[["barcode", "sample_id", "sample_type"]].copy()
        mapping = mapping.merge(counts[["barcode", "filtered_reads"]], on="barcode")
        mapping.sort_values("barcode").to_csv(out / "barcode_sample_map.tsv",
                                              sep="\t", index=False)
        unnamed = sorted(kept - set(frame["barcode"].dropna()))
        if unnamed:
            print(f"{site}: barcodes above the cutoff with no row in the sheet: "
                  f"{', '.join(unnamed)}")

        print(f"{site}: {len(frame)} rows, {int(frame['sequenced'].sum())} sequenced")
        print(f"  analysed patients {len(analysed)}: "
              f"{analysed['group'].value_counts().to_dict()}")
        print(f"  controls passing the cutoff: {list(controls['sample_id'])}")
        print(f"  sequenced but below the cutoff ({len(unused)}): "
              f"{', '.join(unused['sample_id'])}")
        print(f"  written to {path} and {out / 'barcode_sample_map.tsv'}")


if __name__ == "__main__":
    main()
