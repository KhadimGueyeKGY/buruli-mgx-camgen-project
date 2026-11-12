# Buruli ulcer wound microbiome, full length 16S on Oxford Nanopore

Wound swabs from two Cameroonian treatment centres, Ayos and Akonolinga, sequenced as
full length 16S rRNA amplicons on Oxford Nanopore. Each site holds Buruli ulcer
patients and non Buruli patients, and the analysis compares the two groups within each
site, following the project analysis plan.

Reads are 1,572 bp on average and 1,589 bp at the median, which is the full 16S gene
from V1 to V9.

## The database had to be rebuilt first

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

`0_build_silva_16s_db.sh` replaces it with SILVA 138.2 SSU Ref NR99, the reference
designed for 16S amplicons, 147 MB rather than 109 GB. The rebuilt index holds
Bacteria 72.9 percent, Eukaryota 22.4 percent and Archaea 4.7 percent.

After the rebuild no archaeal read survives in any sample, so the archaeal taxa in the
earlier draft should be removed rather than corrected. `FigureS1` and `TableS1` record
the change.

SILVA resolves to genus, not species, so every analysis here stops at genus. Species
level claims from the earlier run cannot be carried over.

## Scripts

Numbered to match the analysis plan. `run_all.sh` runs them in order.

| Script | Plan item | What it does |
|---|---|---|
| `0_build_silva_16s_db.sh` | | Builds the SILVA 138.2 Kraken2 database, then the Bracken database at 1,500 bp |
| `1_run_mgx_qc_taxonomy.sh` | | Merges barcodes, filters at Q10 and 1,000 bp, applies the 0.1 percent inclusion cutoff, runs FastQC, Kraken2 and Bracken |
| `1c_prepare_metadata.py` | | Turns the sample sheet into one tidy metadata table per site, and writes the barcode to sample_id map |
| `1b_recover_sample_names.py` | | Cross-checks that map against the read counts of the first run |
| `2_qc_read_statistics.py` | 2 | Reads per barcode, share retained, read length distribution |
| `2b_build_abundance_tables.py` | | Cuts each report at the Bacteria node, writes count, relative abundance and lineage tables |
| `2c_classification_comparison.py` | | Compares the previous classification with the SILVA one |
| `3_taxonomic_composition.py` | 3 | Phylum and genus composition, by site and by Buruli status |
| `4_alpha_diversity.py` | 4 | Observed genera, Shannon, Simpson, Chao1, with group tests |
| `5_beta_diversity.py` | 5 | Bray-Curtis, principal coordinates, PERMANOVA |
| `6_differential_abundance.py` | 6 | Genus level Buruli against non Buruli, and the Gram composition |
| `7_functional_profiling.py` | 7 | SILVA reference sequence per detected genus, PICRUSt2 prediction, pathway comparison |
| `8_clinical_correlation.py` | 8 | Lesion size, lesion category and treatment against the community |
| `9_cooccurrence_network.py` | 9 | Genus co-occurrence networks per site |
| `verify_outputs.py` | | Re-derives every published number from the Bracken reports and the sample sheet |
| `package_for_sharing.py` | | Zips the figures and tables with a short note for a co-author |
| `common.py` | | Shared paths, loaders, statistics and plot settings |

```bash
bash run_all.sh            # everything except the database build
bash run_all.sh 3 4 5      # only those steps
bash run_all.sh 0          # rebuild the database, hours
```

Step 1 takes about forty minutes on 48 threads for 64 barcodes, step 7 about twenty
minutes. Run them under `tmux`.

## Paths

| | |
|---|---|
| Raw reads | `../../data/raw_16s_buruli/16srna_metagnomics_buruli_2025/` |
| Sample sheet | `../../data/16s_metagenomics/samplesheet_16srna_positive_ub_setAB_final.xlsx` |
| Output | `../../data/16s_metagenomics/analysis_2026/` |
| Kraken2 database | `~/DATABASE/kraken2_silva138.2_16S` |
| Bracken | `~/DATABASE/Bracken/bracken` |
| PICRUSt2 | `~/miniconda3/envs/picrust2` |
| Earlier run, kept for comparison | `../../data/16s_metagenomics/plots_updates/` |

## Method

**Filtering.** Reads of each barcode are merged, then filtered with NanoFilt at mean
quality 10 and minimum length 1,000 bp. A sample is carried forward when it holds at
least 0.1 percent of the filtered reads of its run, the cutoff the analysis plan
specifies. The first run used a flat 10,000 read cutoff instead, which admitted fewer
samples.

These filter settings are those of the first run and reproduce its read counts
exactly. Before the sample sheet was supplied, that fact was used to recover the
barcode to sample mapping from the archived QC table; the sample sheet later confirmed
every one of those 50 assignments. `1b_recover_sample_names.py` keeps the check, and
`verify_outputs.py` re-checks the published map against the sheet.

**Classification.** Kraken2 against SILVA 138.2 SSU Ref NR99 at default confidence,
then Bracken at a read length of 1,500 bp, at phylum, family and genus. Each report is
cut at the Bacteria node and every descendant kept, so the composition panels contain
bacteria only. Taxa are collected on taxid rather than name, because SILVA reuses
labels such as `Incertae Sedis` under unrelated families; a label used by more than one
taxid is prefixed with its parent, decided once over both sites so a taxon carries the
same name everywhere.

**Groups.** UB positive is Buruli ulcer, UB negative is non Buruli, as the plan states.
Water and positive controls are classified but kept out of the group comparisons.

**Alpha diversity.** Observed genera, Shannon, Simpson and Chao1 on the bacterial
genus counts, each sample rarefied to the smallest bacterial library and averaged over
100 draws. Chao1 is the bias corrected form of Chao 1987, computed on the rarefied
counts. Buruli against non Buruli is tested within each site with a two sided
Mann-Whitney test, and the sites against each other on the pooled samples, with
Benjamini-Hochberg correction across the family of tests.

**Beta diversity.** Bray-Curtis on genus relative abundance with a principal
coordinates ordination, and PERMANOVA over 9,999 free permutations for site, Buruli
status and each clinical variable. UniFrac is not available because these reports
carry no phylogeny.

**Differential abundance.** Counts are turned into centred log ratios on a half count
pseudocount, the compositional treatment used by ALDEx2 and ANCOM, and each genus is
tested with a two sided Mann-Whitney corrected with Benjamini-Hochberg within the
site. DESeq2 assumes a negative binomial on counts that are not compositional, which
these are, and the groups here are small, so a rank test on log ratios is the safer
choice. Genera are kept when present in at least 25 percent of the samples of the site
and reaching 0.1 percent mean abundance in one group.

**Gram status** is assigned from the lineage: Actinomycetota, Bacillota and
Deinococcota are Gram positive, except the class Negativicutes which is Gram negative;
the other phyla present are Gram negative. Phyla whose wall is variable or absent are
left unassigned rather than forced into a class.

**Clinical variables.** Lesion size and lesion category are ordinal and encoded 1, 2, 3
from the smallest or mildest upward, so a positive Spearman coefficient means the taxon
rises with severity. Treatment is binary. Correlations are computed on centred log
ratios, so they are not forced negative by the constant sum constraint.

**Co-occurrence.** SparCC and CoNet are not installed, so the same idea is implemented
directly: Spearman correlations between genera on centred log ratios, significance from
a 999 permutation null, Benjamini-Hochberg correction, and edges kept at absolute rho
of at least 0.5. Layout is a seeded Fruchterman-Reingold placement, so the picture is
reproducible.

**Repeat runs.** The sample sheet lists five Ayos patients twice, on two barcodes,
named `AY21-F25` and `AY21-F25_2` and so on, with one identical clinical record each
time. These are one patient measured twice, not two patients, so `common.py` adds their
counts together and every comparison runs on 24 Ayos patients rather than 29 runs.
Treating them as independent would have inflated the Ayos non Buruli group from 15 to
20 and broken the independence the tests assume. Read level quality control still reads
the runs separately, which is what it is about.

**Functional profiling.** Kraken2 and Bracken produce no study sequences, and
clustering the reads does not supply them either: after filtering, a Nanopore 16S read
still carries a few percent of error, so full length reads almost never repeat. On this
run 282,158 reads gave 282,119 distinct sequences and de novo clustering kept 35 units
covering 7 percent of the data, which is noise, not a community. Bracken has already
placed every read on a genus against SILVA, so each detected genus is instead
represented by its own SILVA reference sequence and weighted by its Bracken count:
1,194 genera, 100.000 percent of the bacterial reads. The sequences PICRUSt2 places are
reference quality and the abundances are the ones used everywhere else. These remain
predictions from marker gene content, not measured gene abundances, so the output
supports comparison between groups rather than absolute statements about pathways.

## Verification

`verify_outputs.py` re-derives every published number without reusing the analysis
code: the recovered barcode mapping against the sample sheet, bacterial totals against
the Bracken reports, column sums of every abundance table, the bounds of every
diversity index, the monotonicity of every adjusted p value, the sign of every fold
change against the mean abundance difference, and the presence and size of every
figure. It exits non zero on any failure.

The statistics in `common.py` were each checked against an independent implementation:
Shannon against `scipy.stats.entropy`, Bray-Curtis against
`scipy.spatial.distance.pdist`, Benjamini-Hochberg against
`scipy.stats.false_discovery_control`, Chao1 against a hand worked example, the
PERMANOVA pseudo-F against its definition, and the principal coordinates against a
reconstruction of the distance matrix.

## Results

60 of the 64 barcodes passed the 0.1 percent inclusion cutoff. After removing the
controls, 58 patient runs remain, which are 53 patients once the five Ayos repeat runs
are merged: 24 at Ayos, 9 Buruli and 15 non Buruli, and 29 at Akonolinga, 14 Buruli and
15 non Buruli. Both water controls and the positive control fell below the cutoff,
which is what a clean run should look like.

**Composition.** Pseudomonadota dominates every group, 95 and 96 percent of bacterial
reads in the Ayos Buruli and non Buruli groups and 90 and 85 percent at Akonolinga,
with Bacillota second at 4 to 11 percent. The leading genera over all samples are
*Pseudomonas* 28.1 percent, *Alcaligenes* 23.0 percent and *Escherichia-Shigella* 15.3
percent. *Pseudomonas* is much stronger at Ayos, around 45 to 50 percent of reads
against 6 to 14 percent at Akonolinga.

**Gram.** The communities are overwhelmingly Gram negative, 89 to 96 percent of
bacterial reads in every group, with no difference between Buruli and non Buruli at
either site.

**Alpha diversity.** No index separates Buruli from non Buruli at either site, and the
two sites do not differ either, once Benjamini-Hochberg is applied across the family of
tests. The smallest adjusted p value is 0.49.

**Beta diversity.** The two sites carry different communities, PERMANOVA F = 5.33,
p = 0.0006, adjusted 0.005. Buruli status does not separate the communities at either
site, Ayos F = 0.85 p = 0.46 and Akonolinga F = 0.41 p = 0.87, and neither does lesion
size, lesion category or treatment status.

**Differential abundance.** At Ayos no genus survives correction among the 35 tested.
At Akonolinga one of 49 does: *Corynebacterium* is lower in Buruli lesions, log2 fold
change -3.89, adjusted p = 0.011. The next lowest adjusted p value at that site is
0.56, so this is a single isolated signal and should be treated as a lead rather than
a result.

**Clinical variables.** No genus and no diversity index correlates with lesion size,
lesion category or treatment status after correction. Note that lesion size and lesion
category are the same information: at Ayos they map one to one, category I being 5 cm
or less, II being 5 to 15 cm and III being over 15 cm, and at Akonolinga they differ
for 2 patients out of 29. Testing both does not add evidence. Treatment is also close
to constant, 3 untreated patients at Ayos and 1 at Akonolinga, so it carries almost no
power.

**Predicted function.** All 1,194 genus reference sequences placed below the PICRUSt2
NSTI cutoff of 2.0, median 0.044, so the placement is sound. 594 MetaCyc pathways were
tested at Ayos and 618 at Akonolinga, and none survives Benjamini-Hochberg. The lowest
adjusted p value, 0.083, is the superpathway of mycolyl-arabinogalactan-peptidoglycan
complex biosynthesis at Akonolinga, lower in Buruli lesions, and the next three are
the same cell wall and mycothiol chemistry. That is the gene content of
*Corynebacterium*, the one genus that does survive correction at that site and in the
same direction: across the Akonolinga patients the predicted abundance of those two
cell wall pathways tracks the *Corynebacterium* share at Spearman rho 0.993, and the
two mycothiol pathways at 0.86. The pathway signal and the genus signal are one
finding seen twice, not two.

**Co-occurrence.** Of 595 genus pairs tested per site, 125 edges survive at Ayos and
125 at Akonolinga. Both networks are built around a dense Enterobacteriaceae cluster,
*Escherichia-Shigella*, *Klebsiella*, *Enterobacter*, *Citrobacter*, *Raoultella* and
*Salmonella*, which correlate positively with one another and negatively with
*Alcaligenes* and the anaerobes.

**What the data does not support.** *Mycobacterium* is effectively absent, a single
read across all samples. Full length 16S does not recover *M. ulcerans* from these
swabs, so the pathogen cannot be quantified from this data and the Buruli status used
here comes from the clinical sheet, not from the sequencing.

The main caveat is depth. The 0.1 percent cutoff admits samples with around 3,000
reads, so rarefaction runs at 3,300 bacterial reads. That is enough for genus level
composition in a community this uneven, but it limits the power of the diversity
comparisons. A flat 10,000 read cutoff would rarefy at 22,563 reads and drop 8 samples.

## Output

Figures and tables are written to `analysis_2026/figures/`, as PNG at 400 dpi and as
PDF. Figures follow the order of the analysis plan. Panels share a figure only where
they read together, and no figure is wider than 7.5 inches, so placing one in a 6.5
inch Word column scales it by at least 0.82 and the body type stays above 7 points.
Keys sit outside the axes on the right. Methods belong in this file, not on the
figures.

| File | Content |
|---|---|
| `Figure1_read_statistics` | Reads per sample, share retained, read length |
| `Figure2_phylum_composition` | Phylum composition, per patient and per group |
| `Figure3_top_genera` | The twenty most abundant genera |
| `Figure4_alpha_diversity` | Observed genera, Shannon, Simpson, Chao1 |
| `Figure5_beta_diversity` | Principal coordinates, both sites and each site |
| `Figure6_differential_abundance` | Volcano and heatmap per site |
| `Figure7_gram_composition` | Gram shares, mean and per patient |
| `Figure8_functional_profiling` | Predicted MetaCyc pathways, BU against non-BU |
| `Figure9_clinical_genus_correlation` | Genus abundance against the clinical variables |
| `Figure10_clinical_diversity` | Shannon against the clinical variables |
| `Figure11_cooccurrence_matrix` | Surviving genus correlations per site |
| `FigureS1_classification_comparison` | Classification before and after the rebuild |
| `Table0_read_statistics.tsv` | Reads per barcode with the inclusion decision |
| `Table1_phylum_composition.tsv` | Percent per phylum per sample, with group means |
| `Table2_top_genera.tsv` | Forty leading genera, per site and per group |
| `Table3_domain_composition.tsv` | Reads per domain per sample |
| `Table4_alpha_diversity.tsv` | Diversity indices per sample |
| `Table5_alpha_diversity_tests.tsv` | Mann-Whitney results for the indices |
| `Table6_beta_diversity_permanova.tsv` | PERMANOVA for site, group and clinical variables |
| `Table7_bray_curtis.tsv` | Bray-Curtis distance matrix |
| `Table8_differential_abundance.tsv` | Genus level results with fold change and adjusted p |
| `Table9_gram_composition.tsv` | Gram shares per sample |
| `Table10_gram_tests.tsv` | Gram shares, BU against non-BU |
| `Table11_clinical_genus_correlation.tsv` | Spearman, genus against clinical variables |
| `Table12_clinical_diversity_correlation.tsv` | Spearman, diversity against clinical variables |
| `Table13_treatment_tests.tsv` | Diversity by treatment status |
| `Table14_cooccurrence_edges.tsv` | Every tested genus pair with rho and adjusted p |
| `Table15_pathway_abundance.tsv` | Predicted pathway abundance per sample |
| `Table16_pathway_differential.tsv` | Pathways, BU against non-BU |
| `TableS1_classification_comparison.tsv` | Percent per domain, previous run and SILVA |
| `bacteria_lineage.tsv` | Genus to family, order, class and phylum |
| `functional/genus_reference_provenance.tsv` | SILVA accession behind each genus given to PICRUSt2 |

Per site count and relative abundance tables are in
`analysis_2026/<site>_set<letter>/abundance_tables/`. Those are keyed by sequencing
run, so the five repeated Ayos patients still appear twice there. Everything under
`figures/` is keyed by patient, which is what the comparisons use.

## Requirements

Kraken2 2.0.8-beta, Bracken, NanoFilt, FastQC, seqkit, PICRUSt2 in its own conda
environment with Python 3.10 and pandas 1.5, and Python 3 with pandas, numpy, scipy,
matplotlib and openpyxl. Step 7 also needs the SILVA sequences used to build the
database, at `~/DATABASE/kraken2_silva138.2_16S/library/silva.fna`, and takes about
half an hour, most of it inside the PICRUSt2 check of its own trait tables.

PICRUSt2 installed from the current bioconda package pulls a pandas that removed
`pandas.util.testing`, which PICRUSt2 still imports. The environment therefore pins
`python=3.10`, `pandas<2` and `numpy<2`.

## Citation

Please cite this repository and acknowledge the Buruli MGX Camgen Project.
