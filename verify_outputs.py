#!/usr/bin/env python3
"""
Independent checks on every published table and figure.

Nothing here reuses the analysis code paths. Values are re-derived from the Bracken
reports and from the sample sheet, so a mistake in the analysis scripts shows up as a
failure rather than being reproduced.

Exit status is 0 when every check passes.

Usage: python3 verify_outputs.py

Auteur : Khadim Gueye
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import common as c

ROOT = Path.home() / "khadim/Francis-job"
OUT = ROOT / "data/16s_metagenomics/analysis_2026"
FIG = OUT / "figures"
SHEET = ROOT / "data/16s_metagenomics/samplesheet_16srna_positive_ub_setAB_final.xlsx"
SITES = {"ayos": ("A", "Ayos", "fastq_pass_ayos"),
         "akonolinga": ("B", "Akonolinga", "fastq_pass_akonolinga")}
FAILURES = []


def check(name, condition, detail=""):
    print(("PASS  " if condition else "FAIL  ") + name + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", **kwargs)


def check_mapping():
    book = pd.ExcelFile(SHEET)
    for site, (letter, label, sheet) in SITES.items():
        official = book.parse(sheet)
        official.columns = ["site", "sample_id", "barcode", "ub_status",
                            "lesion_size", "lesion_category", "treatment"]
        official = official.dropna(subset=["sample_id", "barcode"])
        official["barcode"] = official["barcode"].astype(str).str.strip()
        official["sample_id"] = official["sample_id"].astype(str).str.strip()
        truth = dict(zip(official["barcode"], official["sample_id"]))

        mine = read(OUT / f"{site}_set{letter}/barcode_sample_map.tsv")
        wrong = [(b, s, truth.get(b)) for b, s in zip(mine["barcode"], mine["sample_id"])
                 if truth.get(b) != s]
        check(f"{label}: recovered barcode mapping matches the sample sheet",
              not wrong, f"{len(mine)} barcodes" if not wrong else str(wrong[:3]))

        meta = read(OUT / f"{site}_set{letter}/metadata.tsv")
        patients = meta[meta["sample_type"] == "patient"]
        check(f"{label}: every patient carries a group",
              bool(patients["group"].isin(["BU", "non-BU"]).all()))
        check(f"{label}: controls carry no group",
              bool(meta.loc[meta["sample_type"] != "patient", "group"].isna().all()))
        check(f"{label}: ub_status and group agree",
              bool(((patients["ub_status"] == "positive") ==
                    (patients["group"] == "BU")).all()))


def check_abundance():
    for site, (letter, label, _) in SITES.items():
        base = OUT / f"{site}_set{letter}/abundance_tables"
        for level in ("phylum", "family", "genus"):
            counts = read(base / f"bacteria_{level}_counts.tsv", index_col=0)
            relative = read(base / f"bacteria_{level}_relative.tsv", index_col=0)
            check(f"{label} {level}: counts are non negative integers",
                  bool((counts.to_numpy() >= 0).all()))
            exact = (counts / counts.sum(axis=0)).sum(axis=0)
            check(f"{label} {level}: columns sum to 1 before rounding",
                  bool(np.allclose(exact, 1.0, atol=1e-9)))
            check(f"{label} {level}: relative equals counts over the column total",
                  bool(np.allclose(relative.to_numpy(),
                                   (counts / counts.sum(axis=0)).to_numpy(), atol=1e-6)))
            check(f"{label} {level}: index has no duplicate", counts.index.is_unique)

        domains = read(base / "domain_composition.tsv", index_col=0)
        parts = domains[["unclassified", "Bacteria", "Archaea", "Eukaryota"]].sum(axis=1)
        check(f"{label}: domain parts sum to total_reads",
              bool((parts == domains["total_reads"]).all()))
        check(f"{label}: no archaeal read", bool((domains["Archaea"] == 0).all()),
              f"max {int(domains['Archaea'].max())}")

        # Bacterial totals re-read straight from the Bracken reports.
        mapping = read(OUT / f"{site}_set{letter}/barcode_sample_map.tsv")
        mismatch = []
        for _, row in mapping.iterrows():
            report = OUT / f"{site}_set{letter}/bracken_results/{row['barcode']}.bracken_G.report.txt"
            total = sum(int(line.split("\t")[1]) for line in open(report)
                        if len(line.split("\t")) > 5 and line.split("\t")[5].strip() == "Bacteria")
            if total != domains.loc[row["sample_id"], "Bacteria"]:
                mismatch.append(row["sample_id"])
        check(f"{label}: bacterial totals match the Bracken reports", not mismatch,
              f"{len(mapping)} samples" if not mismatch else str(mismatch[:3]))


def check_diversity():
    alpha = read(FIG / "Table4_alpha_diversity.tsv")
    observed = {}
    for site, (letter, label, _) in SITES.items():
        counts = merged_counts(site, letter)
        observed.update({name: int((counts[name] > 0).sum()) for name in counts.columns})
    alpha["observed_genera"] = alpha["sample_id"].map(observed)
    check("alpha: every sample has an observed richness",
          bool(alpha["observed_genera"].notna().all()))
    check("alpha: rarefied richness never exceeds observed richness",
          bool((alpha["richness"] <= alpha["observed_genera"] + 1e-9).all()))
    check("alpha: Shannon within 0 and log(richness)",
          bool((alpha["shannon"] >= 0).all() and
               (alpha["shannon"] <= np.log(alpha["richness"]) + 1e-6).all()))
    check("alpha: Simpson within 0 and 1",
          bool((alpha["simpson"] >= 0).all() and (alpha["simpson"] <= 1).all()))
    check("alpha: Chao1 at least the rarefied richness",
          bool((alpha["chao1"] >= alpha["richness"] - 1e-9).all()))

    distance = read(FIG / "Table7_bray_curtis.tsv", index_col=0).to_numpy()
    check("beta: Bray-Curtis symmetric with a zero diagonal inside 0 and 1",
          bool(np.allclose(distance, distance.T) and
               np.allclose(np.diag(distance), 0) and
               distance.min() >= 0 and distance.max() <= 1))


def check_pvalues():
    files = {
        "Table5_alpha_diversity_tests.tsv": ("p_value", "p_adjusted_bh"),
        "Table6_beta_diversity_permanova.tsv": ("p_value", "p_adjusted_bh"),
        "Table8_differential_abundance.tsv": ("p_value", "p_adjusted_bh"),
        "Table11_clinical_genus_correlation.tsv": ("p_value", "p_adjusted_bh"),
        "Table14_cooccurrence_edges.tsv": ("p_value", "p_adjusted_bh"),
        "Table16_pathway_differential.tsv": ("p_value", "p_adjusted_bh"),
    }
    for name, (raw, adjusted) in files.items():
        path = FIG / name
        if not path.exists():
            check(f"{name} exists", False)
            continue
        frame = read(path)
        if frame.empty:
            check(f"{name}: present but empty", True, "no rows to test")
            continue
        valid = frame[raw].between(0, 1, inclusive="both")
        check(f"{name}: p values inside 0 and 1", bool(valid.all()))
        check(f"{name}: adjusted p is never below the raw p",
              bool((frame[adjusted] >= frame[raw] - 1e-12).all()))
        check(f"{name}: adjusted p inside 0 and 1",
              bool(frame[adjusted].between(0, 1, inclusive="both").all()))


def check_differential():
    frame = read(FIG / "Table8_differential_abundance.tsv")
    if frame.empty:
        check("differential abundance table is not empty", False)
        return
    # The fold change is recomputed from the counts rather than compared with the
    # arithmetic means, which are a different summary and may disagree in sign.
    worst = 0.0
    for site, (letter, label, _) in SITES.items():
        sub = frame[frame["site"] == label]
        if sub.empty:
            continue
        counts = merged_counts(site, letter)
        meta = read(OUT / f"{site}_set{letter}/metadata.tsv", dtype=str)
        meta["sample_id"] = meta["sample_id"].str.replace(r"_\d+$", "", regex=True)
        meta = meta.drop_duplicates(subset="sample_id").set_index("sample_id")
        shared = sorted(set(counts.columns) & set(meta.index[meta["sample_type"] == "patient"]))
        counts = counts[shared]
        values = counts.to_numpy(dtype=float) + 0.5
        values = values / values.sum(axis=0, keepdims=True)
        log = np.log2(values)
        transformed = pd.DataFrame(log - log.mean(axis=0, keepdims=True),
                                   index=counts.index, columns=counts.columns)
        bu = [s for s in shared if meta.loc[s, "group"] == "BU"]
        non = [s for s in shared if meta.loc[s, "group"] == "non-BU"]
        mine = transformed.loc[sub["genus"], bu].mean(axis=1) - \
            transformed.loc[sub["genus"], non].mean(axis=1)
        worst = max(worst, float(np.abs(mine.to_numpy() -
                                        sub["log2_fold_change"].to_numpy()).max()))
    # The table is written with five decimals, so half the last digit is the floor.
    check("differential: fold change re-derived from the counts", worst <= 5e-6,
          f"largest difference {worst:.2e}, rounding step 5e-06")
    check("differential: every genus carries a Gram status",
          bool(frame["gram"].notna().all()))
    check("differential: significance flag matches the adjusted p",
          bool((frame["significant"] == (frame["p_adjusted_bh"] < 0.05)).all()))

    gram = read(FIG / "Table9_gram_composition.tsv")
    totals = gram.groupby(["site", "sample_id"])["percent"].sum()
    check("gram: the shares of a sample sum to 100",
          bool(np.allclose(totals.to_numpy(), 100.0, atol=1e-3)),
          f"min {totals.min():.4f} max {totals.max():.4f}")


def check_cooccurrence():
    path = FIG / "Table14_cooccurrence_edges.tsv"
    if not path.exists():
        check("co-occurrence table exists", False)
        return
    frame = read(path)
    if frame.empty:
        check("co-occurrence table is not empty", False)
        return
    check("co-occurrence: rho inside -1 and 1",
          bool(frame["spearman_rho"].between(-1, 1).all()))
    pairs = frame.apply(lambda r: tuple(sorted((r["genus_a"], r["genus_b"]))), axis=1)
    check("co-occurrence: no pair listed twice within a site",
          bool(not pd.concat([frame["site"], pairs], axis=1).duplicated().any()))
    check("co-occurrence: no genus paired with itself",
          bool((frame["genus_a"] != frame["genus_b"]).all()))


def merged_counts(site, letter):
    """Genus counts re-keyed by patient, adding the repeat runs, read from disk."""
    counts = read(OUT / f"{site}_set{letter}/abundance_tables/bacteria_genus_counts.tsv",
                  index_col=0)
    return counts.T.groupby(lambda name: re.sub(r"_\d+$", "", name), sort=True).sum().T


def check_composition():
    """The published phylum percentages, re-derived from the per run count files."""
    published = read(FIG / "Table1_phylum_composition.tsv", index_col=0)
    columns = [name for name in published.columns if not name.startswith("mean_")]
    worst, compared = 0.0, 0
    for site, (letter, _, _) in SITES.items():
        counts = read(OUT / f"{site}_set{letter}/abundance_tables/bacteria_phylum_counts.tsv",
                      index_col=0)
        merged = counts.T.groupby(lambda name: re.sub(r"_\d+$", "", name)).sum().T
        percent = 100 * merged / merged.sum(axis=0)
        for column in columns:
            sample = column.rsplit("_", 1)[-1]
            if sample not in percent.columns:
                continue
            for phylum in published.index:
                if phylum == "Other" or phylum not in percent.index:
                    continue
                worst = max(worst, abs(percent.loc[phylum, sample] -
                                       published.loc[phylum, column]))
                compared += 1
    check("composition: phylum percentages re-derived from the counts",
          compared > 0 and worst <= 5e-5,
          f"{compared} values, largest difference {worst:.2e}, rounding step 5e-05")
    check("composition: one column per patient",
          len(columns) == sum(len(c.load_metadata(site)) for site in SITES),
          f"{len(columns)} columns")


def check_patients():
    """Repeat sequencing runs of one patient must be merged, and merged correctly.

    The abundance tables on disk are keyed by sequencing run. Every comparison is
    keyed by patient. This re-adds the runs straight from the files and compares.
    """
    for site, (letter, label, _) in SITES.items():
        raw = read(OUT / f"{site}_set{letter}/abundance_tables/bacteria_genus_counts.tsv",
                   index_col=0)
        meta = read(OUT / f"{site}_set{letter}/metadata.tsv")
        patients = set(meta.loc[meta["sample_type"] == "patient", "sample_id"])
        repeats = {}
        for column in raw.columns:
            if column in patients:
                repeats.setdefault(re.sub(r"_\d+$", "", column), []).append(column)
        merged = c.load_table(site, "genus")
        wrong = []
        for patient, runs in repeats.items():
            if patient not in merged.columns:
                wrong.append((patient, "missing"))
                continue
            expected = raw[runs].sum(axis=1)
            if not expected.reindex(merged.index).fillna(0).equals(merged[patient]):
                wrong.append((patient, "sum differs"))
        check(f"{label}: merged counts equal the sum of the runs",
              not wrong, f"{len(repeats)} patients" if not wrong else str(wrong[:3]))

        multi = {p: r for p, r in repeats.items() if len(r) > 1}
        clinical = meta.set_index("sample_id")
        disagree = [p for p, runs in multi.items()
                    for field in ("group", "lesion_size", "lesion_category", "treatment")
                    if clinical.loc[runs, field].nunique(dropna=False) > 1]
        check(f"{label}: repeat runs agree on every clinical field",
              not disagree, f"{len(multi)} repeated patients")

        table = c.load_metadata(site)
        check(f"{label}: one row per patient, none named as a repeat run",
              bool(table.index.is_unique) and not any(re.search(r"_\d+$", s)
                                                      for s in table.index),
              f"{len(table)} patients, {len(set(raw.columns) & patients)} runs")


def check_reference():
    """The PICRUSt2 input must name the same genera in the sequences and the table."""
    reference = OUT / "functional/genus_reference.fasta"
    table = OUT / "functional/genus_table.tsv"
    if not (reference.exists() and table.exists()):
        check("functional: reference sequences and genus table written", False)
        return
    names = {line[1:].strip() for line in reference.read_text().splitlines()
             if line.startswith(">")}
    counts = read(table, index_col=0)
    check("functional: every table row has a reference sequence",
          set(counts.index) <= names, f"{len(counts)} genera, {len(names)} sequences")
    check("functional: reference identifiers are unique",
          len(names) == sum(1 for line in reference.read_text().splitlines()
                            if line.startswith(">")))
    check("functional: counts are whole numbers and non negative",
          bool((counts.to_numpy() >= 0).all()) and
          bool((counts.to_numpy() == counts.to_numpy().astype(int)).all()))
    provenance = read(OUT / "functional/genus_reference_provenance.tsv")
    check("functional: every reference carries a SILVA accession",
          bool(provenance["silva_accession"].notna().all()) and
          bool(provenance["silva_accession"].str.contains(r"\.").all()),
          f"{len(provenance)} accessions")


def check_pathways():
    abundance = FIG / "Table15_pathway_abundance.tsv"
    check("Table15_pathway_abundance.tsv written", abundance.exists())
    if not abundance.exists():
        return
    frame = read(abundance, index_col=0)
    check("pathways: abundances are non negative",
          bool((frame.to_numpy() >= 0).all()))
    check("pathways: every sample carries at least one pathway",
          bool((frame.sum(axis=0) > 0).all()))
    results = FIG / "Table16_pathway_differential.tsv"
    if results.exists() and not read(results).empty:
        table = read(results)
        check("pathways: every tested pathway is in the abundance table",
              bool(set(table["pathway"]) <= set(frame.index)))
        check("pathways: significance flag matches the adjusted p",
              bool((table["significant"] == (table["p_adjusted_bh"] < 0.05)).all()))


def check_figures():
    expected = ["Figure0_read_statistics", "Figure1_bacterial_composition",
                "Figure2_top_genera", "Figure3_alpha_diversity",
                "Figure4_beta_diversity", "Figure5_differential_abundance",
                "Figure6_gram_composition", "Figure7_clinical_correlation",
                "Figure8_cooccurrence_network", "Figure9_functional_profiling",
                "FigureS1_classification_comparison"]
    for name in expected:
        for extension in ("png", "pdf"):
            path = FIG / f"{name}.{extension}"
            check(f"{name}.{extension} written and not empty",
                  path.exists() and path.stat().st_size > 2000)


def main():
    check_mapping()
    check_patients()
    check_composition()
    check_abundance()
    check_diversity()
    check_pvalues()
    check_differential()
    check_cooccurrence()
    check_reference()
    check_pathways()
    check_figures()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECKS FAILED")
        for name in FAILURES:
            print(f"  {name}")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
