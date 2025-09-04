#!/bin/bash

# -------------------------------
# 0. Parameters
# -------------------------------
SITE='akonolinga'      #  "ayos" -> A or "akonolinga" -> B 
SITESET='B'   

FASTQ_DIR="../output/16srna_buruli_${SITE}_set${SITESET}/filtered_fastq"
OUT_DIR="../output/16srna_buruli_${SITE}_set${SITESET}/functional_profile_${SITE}_set${SITESET}"
THREADS=8

mkdir -p $OUT_DIR

# -------------------------------
# 1. Loop over all fastq files
# -------------------------------
for fq in $FASTQ_DIR/*.fastq.gz; do
    sample=$(basename $fq .fastq.gz)
    echo "Processing sample $sample"
    
    humann --input $fq \
           --output $OUT_DIR/$sample \
           --threads $THREADS \
           --output-format tsv
done

# -------------------------------
# 2. Combine tables
# -------------------------------
humann_join_tables --input $OUT_DIR --output $OUT_DIR/combined_pathabundance.tsv --file_name pathabundance
humann_join_tables --input $OUT_DIR --output $OUT_DIR/combined_genefamilies.tsv --file_name genefamilies

# -------------------------------
# 3. Normalize to relative abundance
# -------------------------------
humann_renorm_table --input $OUT_DIR/combined_pathabundance.tsv \
                    --output $OUT_DIR/combined_pathabundance_relab.tsv \
                    --units relab
