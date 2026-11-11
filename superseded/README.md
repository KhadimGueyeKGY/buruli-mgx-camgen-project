# Superseded scripts

Kept for reference only, do not run.

`1_run_mgx_qc_taxonomy_1.sh` points at `~/kraken2_db`, which does not exist, and is
replaced by `../1_run_mgx_qc_taxonomy.sh`.

`3_plot_taxonomic_composition.sh` and `plot_taxonomic_composition.R` plot species
level composition. SILVA resolves to genus, so species level output from this data is
not supportable. They are replaced by `../3_plot_taxonomic_composition.py`.
