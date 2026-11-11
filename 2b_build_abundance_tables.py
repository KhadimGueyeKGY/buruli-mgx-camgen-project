#!/usr/bin/env python3
"""
Build per sample abundance tables restricted to Bacteria.

Francis asked for cellular organisms, and bacteria in particular. Bracken reports
carry the whole taxonomic tree, so each report is cut at the Bacteria node and every
descendant is kept. Relative abundance is expressed within the bacterial fraction, so
the columns of a sample sum to one. A domain summary is written alongside, so the
share of reads that is bacterial, archaeal, eukaryotic or unclassified stays visible.

Taxa are collected on taxid, never on name, because SILVA reuses labels such as
Incertae Sedis under unrelated families. Labels are resolved once across both sites,
so a taxon carries the same name in every table.

Usage: python3 2b_build_abundance_tables.py

Auteur : Khadim Gueye
"""

from pathlib import Path

import pandas as pd

ROOT = Path.home() / "khadim/Francis-job"
OUT_ROOT = ROOT / "data/16s_metagenomics/analysis_2026"
SITES = {"ayos": "A", "akonolinga": "B"}
LEVELS = {"P": "phylum", "F": "family", "G": "genus"}


def parse_report(path):
    """Read a Kraken2 style report and attach each row's parent name."""
    rows, stack = [], []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        name = f[5]
        indent = len(name) - len(name.lstrip(" "))
        while stack and stack[-1][0] >= indent:
            stack.pop()
        rows.append({
            "clade_reads": int(f[1]),
            "rank": f[3],
            "taxid": f[4],
            "name": name.strip(),
            "parent": stack[-1][1] if stack else "",
            "indent": indent,
        })
        stack.append((indent, name.strip()))
    return rows


def subtree(rows, domain):
    """Rows under a named node. Matched on name and indentation, because the rank
    code for a domain differs between a SILVA and an NCBI built database."""
    out, inside, base = [], False, None
    for row in rows:
        if inside and row["indent"] <= base:
            break
        if inside:
            out.append(row)
        elif row["name"] == domain:
            inside, base = True, row["indent"]
            out.append(row)
    return out


def domain_totals(rows):
    out = {"unclassified": sum(r["clade_reads"] for r in rows if r["rank"] == "U")}
    for name in ("Bacteria", "Archaea", "Eukaryota"):
        hit = subtree(rows, name)
        out[name] = hit[0]["clade_reads"] if hit else 0
    return out


def collect():
    """counts[(site, level)][sample][taxid], the labels, and the domain totals."""
    counts, labels, domains, names = {}, {}, {}, {}
    for site, letter in SITES.items():
        out = OUT_ROOT / f"{site}_set{letter}"
        mapping = {}
        if (out / "barcode_sample_map.tsv").exists():
            frame = pd.read_csv(out / "barcode_sample_map.tsv", sep="\t")
            mapping = dict(zip(frame["barcode"], frame["sample_id"]))
        names[site] = mapping
        for code, level in LEVELS.items():
            table = {}
            for report in sorted((out / "bracken_results").glob(f"*.bracken_{code}.report.txt")):
                barcode = report.name.split(".bracken_")[0]
                if mapping and barcode not in mapping:
                    continue  # kept by the read cutoff but not in the published set
                sample = mapping.get(barcode, barcode)
                rows = parse_report(report)
                if code == "G":
                    domains.setdefault(site, {})[sample] = domain_totals(rows)
                table[sample] = {}
                for row in subtree(rows, "Bacteria"):
                    if row["rank"] == code:
                        table[sample][row["taxid"]] = row["clade_reads"]
                        labels[row["taxid"]] = (row["name"], row["parent"])
            counts[(site, level)] = table
    return counts, labels, domains


def resolve(labels):
    """One label per taxid, decided once over both sites."""
    seen = {}
    for taxid, (name, _) in labels.items():
        seen.setdefault(name, set()).add(taxid)
    return {taxid: name if len(seen[name]) == 1 else f"{parent} {name}"
            for taxid, (name, parent) in labels.items()}


def main():
    counts, labels, domains = collect()
    label_of = resolve(labels)

    for site, letter in SITES.items():
        tables = OUT_ROOT / f"{site}_set{letter}/abundance_tables"
        tables.mkdir(exist_ok=True)
        for level in LEVELS.values():
            table = counts[(site, level)]
            if not table:
                continue
            frame = pd.DataFrame(table).fillna(0).astype(int)
            frame.index = [label_of[t] for t in frame.index]
            frame = frame.groupby(level=0).sum()
            frame = frame[sorted(frame.columns)]
            frame = frame.loc[frame.sum(axis=1).sort_values(ascending=False).index]
            frame.index.name = level
            frame.to_csv(tables / f"bacteria_{level}_counts.tsv", sep="\t")
            rel = frame.div(frame.sum(axis=0).replace(0, pd.NA), axis=1).fillna(0)
            rel.to_csv(tables / f"bacteria_{level}_relative.tsv", sep="\t", float_format="%.8f")
            print(f"{site} {level}: {frame.shape[0]} taxa, {frame.shape[1]} samples")

        dom = pd.DataFrame(domains[site]).T.sort_index()
        dom.index.name = "sample_id"
        dom["total_reads"] = dom.sum(axis=1)
        for col in ("Bacteria", "Archaea", "Eukaryota", "unclassified"):
            dom[f"pct_{col}"] = 100 * dom[col] / dom["total_reads"]
        dom.to_csv(tables / "domain_composition.tsv", sep="\t", float_format="%.3f")
        print(f"{site} domains: median bacterial share {dom['pct_Bacteria'].median():.1f} percent")


if __name__ == "__main__":
    main()
