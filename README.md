# 16S rRNA Metagenomics Analysis Pipeline – Buruli MGX Camgen Project

This repository contains a modular pipeline for the analysis of 16S rRNA metagenomic sequencing data.  
The workflow is organised into sequential steps, covering quality control, taxonomic profiling, diversity analyses, differential abundance testing, co-occurrence network construction, and clinical correlations.

---

## Workflow Overview

The pipeline is composed of the following steps:

1. **Quality Control and Taxonomic Assignment**  
   - Script: `1_run_mgx_qc_taxonomy.sh`  
   - Performs raw read quality control, filtering, and taxonomic classification using **Kraken2** and **Bracken**.  
   - Outputs include cleaned FASTQ files, taxonomic profiles, and summary tables.  

2. **QC Summary Plots**  
   - Script: `2_qc_summary_plots.R`  
   - Generates visual summaries of read quality and taxonomic distributions across samples.  
   - Outputs include barplots and PDF figures.  

3. **Taxonomic Composition Visualization**  
   - Script: `3_plot_taxonomic_composition.sh` (calls `plot_taxonomic_composition.R`)  
   - Creates stacked barplots at different taxonomic levels (phylum, genus, species).  
   - Also generates Krona-compatible input files for interactive taxonomic exploration.  
   - Outputs are stored in `output/plots_combined_sites/`.  

4. **Alpha and Beta Diversity Analyses**  
   - Scripts:  
     - `4_1_alpha_beta_diversity.R` (all samples)  
     - `4_2_alpha_beta_diversity_two_sites.R` (comparisons between sites)  
   - Computes:  
     - Alpha diversity (Shannon, Simpson, etc.)  
     - Beta diversity (ordination, PCoA, PERMANOVA).  
   - Outputs include statistical summaries and ordination plots.  

5. **Differential Abundance Analysis**  
   - Scripts:  
     - `5_1_differential_abundance_analysis.R` (per site)  
     - `5_2_differential_abundance_analysis_both_sites.R` (site comparisons)  
   - Performs statistical testing to identify taxa significantly enriched or depleted across conditions or sites.  
   - Generates volcano plots, heatmaps, and result tables.  

6. **Co-occurrence Network Analysis**  
   - Scripts:  
     - `6_1_cooccurrence_network.R` (per site)  
     - `6_2_cooccurrence_network_both_sites.R` (combined analysis)  
   - Constructs microbial co-occurrence networks based on correlation measures.  
   - Outputs include network graphs and adjacency matrices.  

7. **Clinical Correlation Analysis**  
   - Script: `7_clinical_correlation.py`  
   - Links microbiome features to available clinical metadata.  
   - Performs correlation analyses (e.g., Spearman/Pearson tests) and generates summary tables and plots.  

---

## Requirements

- **R (≥ 4.0.0)** with the following packages:  
  `tidyverse`, `vegan`, `phyloseq`, `igraph`, `ggplot2`  
- **Python (≥ 3.8)** with:  
  `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`  
- **Kraken2**, **Bracken**, **KronaTools**  

---

## Usage

1. Clone the repository:  
   
   ```
   https://github.com/KhadimGueyeKGY/buruli-mgx-camgen-project.git
   cd buruli-mgx-camgen-project
   ```
  

2. Run each step sequentially, adapting inputs and paths as required:

   ```bash
   bash 1_run_mgx_qc_taxonomy.sh
   Rscript 2_qc_summary_plots.R
   bash 3_plot_taxonomic_composition.sh
   Rscript 4_1_alpha_beta_diversity.R
   Rscript 5_1_differential_abundance_analysis.R
   Rscript 6_1_cooccurrence_network.R
   python3 7_clinical_correlation.py
   ```

3. Outputs are organised into dedicated subdirectories (`output/`).

---

## Repository Structure

```
.
├── 1_run_mgx_qc_taxonomy.sh
├── 2_qc_summary_plots.R
├── 3_plot_taxonomic_composition.sh
├── 4_1_alpha_beta_diversity.R
├── 4_2_alpha_beta_diversity_two_sites.R
├── 5_1_differential_abundance_analysis.R
├── 5_2_differential_abundance_analysis_both_sites.R
├── 6_1_cooccurrence_network.R
├── 6_2_cooccurrence_network_both_sites.R
├── 7_clinical_correlation.py
├── Bracken/
├── Krona/
└── output/
```

---

## Citation

If you use this pipeline in your research, please cite this repository and acknowledge the **Buruli MGX Camgen Project**.


