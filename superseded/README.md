# Superseded scripts

Kept for reference only, do not run. Each was written against the previous Kraken2
database, which held no bacterial sequences, or against species level output that a
SILVA built database cannot support.

| Superseded | Replaced by |
|---|---|
| `1_run_mgx_qc_taxonomy_1.sh` | `../1_run_mgx_qc_taxonomy.sh`, the old one points at `~/kraken2_db`, which does not exist |
| `2_qc_summary_plots.R` | `../2_qc_read_statistics.py` |
| `3_plot_taxonomic_composition.sh`, `plot_taxonomic_composition.R`, `3_plot_taxonomic_composition.py` | `../3_taxonomic_composition.py` |
| `4_1_alpha_beta_diversity.R`, `4_2_alpha_beta_diversity_two_sites.R`, `4_alpha_beta_diversity.py` | `../4_alpha_diversity.py` and `../5_beta_diversity.py` |
| `5_1_differential_abundance_analysis.R`, `5_2_differential_abundance_analysis_both_sites.R` | `../6_differential_abundance.py` |
| `6_1_cooccurrence_network.R`, `6_2_cooccurrence_network_both_sites.R` | `../9_cooccurrence_network.py` |
| `7_clinical_correlation.py` | `../8_clinical_correlation.py` |
| `8_1_run_humann3.sh`, `8_2_plot_functional_profile.R` | `../7_functional_profiling.py`, HUMAnN3 profiles shotgun reads and does not apply to 16S amplicons |
