#!/usr/bin/env bash
# Run the analysis end to end, in the order of the analysis plan.
#
#   bash run_all.sh              every step except the database build
#   bash run_all.sh 3 4 5        only those steps
#   bash run_all.sh 0            rebuild the SILVA database, hours
#
# Step 1 reads every FASTQ and takes about forty minutes on 48 threads. Step 7 runs
# PICRUSt2 and takes about two hours, most of it inside the PICRUSt2 check of its own
# trait tables, which it reads in full whatever --in_traits asks for. Run both under
# tmux.
#
# Auteur : Khadim Gueye

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"

STEPS=("$@")
[[ ${#STEPS[@]} -eq 0 ]] && STEPS=(1 1c 1b 2 2b 2c 3 4 5 6 7 8 9 verify pack)

say() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
wanted() { [[ " ${STEPS[*]} " == *" $1 "* ]]; }

cd "$HERE"

wanted 0  && { say "0  build the SILVA Kraken2 and Bracken databases"
               bash 0_build_silva_16s_db.sh; }
wanted 1  && { say "1  QC and taxonomic classification"; bash 1_run_mgx_qc_taxonomy.sh; }
wanted 1c && { say "1c prepare the metadata and the barcode map"; "$PYTHON" 1c_prepare_metadata.py; }
wanted 1b && { say "1b cross-check the mapping against the archived read counts"
               "$PYTHON" 1b_recover_sample_names.py ayos A
               "$PYTHON" 1b_recover_sample_names.py akonolinga B; }
wanted 2  && { say "2  read statistics"; "$PYTHON" 2_qc_read_statistics.py; }
wanted 2b && { say "2b bacterial abundance tables"; "$PYTHON" 2b_build_abundance_tables.py; }
wanted 2c && { say "2c classification before and after the rebuild"
               "$PYTHON" 2c_classification_comparison.py; }
wanted 3  && { say "3  taxonomic composition"; "$PYTHON" 3_taxonomic_composition.py; }
wanted 4  && { say "4  alpha diversity"; "$PYTHON" 4_alpha_diversity.py; }
wanted 5  && { say "5  beta diversity"; "$PYTHON" 5_beta_diversity.py; }
wanted 6  && { say "6  differential abundance"; "$PYTHON" 6_differential_abundance.py; }
wanted 7  && { say "7  functional profiling"; "$PYTHON" 7_functional_profiling.py; }
wanted 8  && { say "8  clinical correlation"; "$PYTHON" 8_clinical_correlation.py; }
wanted 9  && { say "9  co-occurrence networks"; "$PYTHON" 9_cooccurrence_network.py; }
wanted verify && { say "verify the published tables and figures"
                   "$PYTHON" verify_outputs.py; }
wanted pack   && { say "pack the figures and tables for a co-author"
                   "$PYTHON" package_for_sharing.py; }

say "done"
