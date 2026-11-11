# Buruli ulcer wound microbiome, full length 16S on Oxford Nanopore

Wound swabs from two Cameroonian treatment centres, Ayos and Akonolinga, sequenced
as full length 16S rRNA amplicons on Oxford Nanopore. The pipeline covers read
filtering, taxonomic classification, bacterial composition, diversity, differential
abundance and clinical correlation.

Reads are 1,572 bp on average and 1,589 bp at the median, which is the full 16S gene
from V1 to V9.

## The database had to be rebuilt

The first run of this project reported communities made almost entirely of archaea,
with 96 percent of reads unassigned. The cause was the Kraken2 database, not the
samples. Inspecting it showed Eukaryota 78.1 percent and Archaea 17.4 percent and no
bacterial sequence at all. Its build log shows the bacterial download failed:

```
-> Downloading bacteria
rsync_from_ncbi.pl: FTP binary mode error: [Net::FTP] Connection closed
Sequence ID to taxonomy ID map already present, skipping map creation.
Hash table already present, skipping database file build.
```

The download failed, and the rerun then skipped index construction because a hash
table was already on disk. The 252 GB of bacterial sequence sat in `library/bacteria`
without ever entering the index.

`0_build_silva_16s_db.sh` replaces it with SILVA 138.2 SSU Ref NR99, which is the
right reference for 16S amplicons and is 147 MB rather than 109 GB. The rebuilt index
holds Bacteria 72.9 percent, Eukaryota 22.4 percent and Archaea 4.7 percent.

Result across the 50 samples:

| | Bacteria | Archaea | Unclassified |
|---|---|---|---|
| Previous database, Ayos | 0.0 | 1.7 | 98.3 |
| SILVA 138.2, Ayos | 100.0 | 0.0 | 0.0 |
| Previous database, Akonolinga | 0.0 | 4.7 | 95.3 |
| SILVA 138.2, Akonolinga | 100.0 | 0.0 | 0.0 |

Mean percent of filtered reads. No archaeal read survives the rebuild, so the archaeal
taxa in the earlier draft should be removed rather than corrected.

SILVA resolves to genus, not species, so Bracken is run at phylum, family and genus.
Species level claims from the earlier run cannot be carried over.

## Scripts

Run in order. Steps 0 and 1 are the long ones.

| Script | What it does |
|---|---|
| `0_build_silva_16s_db.sh` | Builds the SILVA 138.2 Kraken2 database, then the Bracken database at 1,500 bp |
| `1_run_mgx_qc_taxonomy.sh` | Merges barcodes, filters at Q10 and 1,000 bp, drops samples below 10,000 reads, runs FastQC, Kraken2 and Bracken |
| `1b_recover_sample_names.py` | Maps barcode to sample_id by matching filtered read counts against the archived QC table |
| `2_qc_summary_plots.R` | Read length and read count summaries |
| `2b_build_abundance_tables.py` | Cuts each report at the Bacteria node and writes count and relative abundance tables per level |
| `2c_classification_comparison.py` | Compares the previous classification with the SILVA one |
| `3_plot_taxonomic_composition.py` | Figure 1 and Figure 2, Tables 1 to 3 |
| `4_alpha_beta_diversity.py` | Figure 3 and Tables 4 and 5 |
| `5_1`, `5_2` | Differential abundance, needs the clinical metadata |
| `6_1`, `6_2` | Co-occurrence networks, needs the clinical metadata |
| `7_clinical_correlation.py` | Lesion size and category against taxa, needs the clinical metadata |
| `8_1`, `8_2` | HUMAnN3 functional profiling, shotgun only, not applicable to amplicons |

```bash
bash 0_build_silva_16s_db.sh
bash 1_run_mgx_qc_taxonomy.sh
python3 1b_recover_sample_names.py ayos A
python3 1b_recover_sample_names.py akonolinga B
python3 2b_build_abundance_tables.py
python3 2c_classification_comparison.py
python3 3_plot_taxonomic_composition.py
python3 4_alpha_beta_diversity.py
```

Step 1 takes about half an hour on 48 threads for 64 barcodes. Run it under `tmux`.

## Paths

| | |
|---|---|
| Raw reads | `../../data/raw_16s_buruli/16srna_metagnomics_buruli_2025/` |
| Output | `../../data/16s_metagenomics/analysis_2026/` |
| Kraken2 database | `~/DATABASE/kraken2_silva138.2_16S` |
| Bracken | `~/DATABASE/Bracken/bracken` |
| Earlier run, kept for comparison | `../../data/16s_metagenomics/plots_updates/` |

## Method

Reads of each barcode are merged, then filtered with NanoFilt at mean quality 10 and
minimum length 1,000 bp. A sample is carried forward when at least 10,000 reads
survive. These are the settings of the first run, and they reproduce its read counts
exactly, which is how the barcode to sample mapping was recovered after the metadata
sheets were lost.

Classification uses Kraken2 against SILVA 138.2 SSU Ref NR99 at default confidence,
followed by Bracken at a read length of 1,500 bp, at phylum, family and genus. Each
report is cut at the Bacteria node and every descendant kept, so the composition
panels contain bacteria only. Taxa are collected on taxid rather than name, because
SILVA reuses labels such as `Incertae Sedis` under unrelated families; a label used by
more than one taxid is prefixed with its parent. Relative abundance is expressed
within the bacterial fraction of a sample.

Alpha diversity is computed on the bacterial genus counts after rarefying every sample
to the smallest bacterial library, 22,563 reads, averaged over 100 draws. Sites are
compared with a two sided Mann-Whitney test. Beta diversity uses Bray-Curtis on
relative abundance with a principal coordinates ordination and a PERMANOVA over 9,999
permutations.

## Results

50 samples, 22 from Ayos and 28 from Akonolinga.

Pseudomonadota dominates both sites, 97.2 percent of bacterial reads at Ayos and 87.1
percent at Akonolinga. Bacillota comes second at 2.5 and 9.9 percent. The leading
genera over both sites are *Pseudomonas* 28.7 percent, *Alcaligenes* 24.0 percent and
*Escherichia-Shigella* 15.0 percent. *Pseudomonas* is the stronger of the two at Ayos,
51.8 percent against 10.5 percent at Akonolinga.

Shannon diversity is higher at Akonolinga, Mann-Whitney p = 0.045, while observed
genus richness does not differ, p = 0.305. Community composition differs between the
sites, PERMANOVA F = 5.13, p = 0.0012, with PCo1 and PCo2 carrying 37.6 and 30.2
percent.

*Mycobacterium* is effectively absent, a single read across the 50 samples. Full
length 16S does not recover *M. ulcerans* from these swabs, so the pathogen cannot be
quantified from this data.

## Output

Figures and tables are written to `analysis_2026/figures/`.

| File | Content |
|---|---|
| `Figure1_bacterial_composition` | Phylum composition per sample, and the mean per site |
| `Figure2_top_genera` | Twenty most abundant genera per sample |
| `Figure3_diversity` | Shannon, observed genera, Bray-Curtis ordination |
| `FigureS1_classification_comparison` | Reads assigned before and after the rebuild |
| `Table1_phylum_composition.tsv` | Percent per phylum per sample, with the per site means |
| `Table2_top_genera.tsv` | Forty leading genera, mean and median percent, detection counts |
| `Table3_domain_composition.tsv` | Reads per domain per sample |
| `Table4_alpha_diversity.tsv` | Richness, Shannon and Simpson per sample |
| `Table5_bray_curtis.tsv` | Bray-Curtis distance matrix |
| `TableS1_classification_comparison.tsv` | Per sample, per run, percent assigned to each domain |

Per site count and relative abundance tables are in
`analysis_2026/<site>_set<letter>/abundance_tables/`.

## What is still blocked

`fastq_pass_ayos.tsv` and `fastq_pass_akonolinga.tsv` carried the barcode, sample_id,
Buruli status, lesion size, lesion category and treatment. They are not on this
machine. Sample names were recovered from the read counts, but Buruli against non
Buruli, lesion size and treatment cannot be. Scripts 5, 6 and 7 stay blocked until
those sheets are supplied. Empty templates carrying the real barcodes are in
`analysis_2026/metadata_template/`.

Scripts 8_1 and 8_2 run HUMAnN3, which profiles gene families from shotgun reads.
They do not apply to 16S amplicons and are kept only for reference.

## Requirements

Kraken2 2.0.8-beta, Bracken, NanoFilt, FastQC, seqkit, Python 3 with pandas, numpy,
scipy and matplotlib, R with tidyverse, vegan, phyloseq and igraph.

## Citation

Please cite this repository and acknowledge the Buruli MGX Camgen Project.
