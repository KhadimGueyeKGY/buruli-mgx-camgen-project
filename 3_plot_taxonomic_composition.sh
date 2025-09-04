#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
SITE="ayos"   # "akonolinga" ou "ayos"
SITESET="A"         # "A" ou "B"
BRACKEN_DIR="../output/16srna_buruli_${SITE}_set${SITESET}/bracken_results"
METADATA="../fastq_pass_${SITE}.tsv"
OUTPUT_DIR="../output/plot_${SITE}_set${SITESET}"
KRONA_DIR="$OUTPUT_DIR/krona_species"
KRAKEN_DIR="../output/16srna_buruli_${SITE}_set${SITESET}/kraken2_reports"
KRAKEN_KRONA_DIR="../output/krona_plots"
NCORES=4

# === Create directories ===
mkdir -p "$OUTPUT_DIR"
mkdir -p "$KRONA_DIR"
mkdir -p "$KRAKEN_KRONA_DIR"

echo "=== Step 1: Running R script to generate barplot + Krona inputs ==="
Rscript plot_taxonomic_composition.R \
    --bracken_dir "$BRACKEN_DIR" \
    --metadata "$METADATA" \
    --output_dir "$OUTPUT_DIR" \
    --krona_dir "$KRONA_DIR"

echo "=== Step 2: Generating Krona plots per-sample from Bracken ==="
for f in "$KRONA_DIR"/*.txt; do
    base=$(basename "$f" .txt)
    ktImportText -o "$KRONA_DIR/${base}.html" "$f"
done

echo "=== Step 3: Generating combined Krona plot from all Kraken2 reports ==="

KRONA_TAX_PATH="Krona/KronaTools/taxonomy"

cat "$KRAKEN_DIR"/*.txt > "$KRAKEN_KRONA_DIR/combined_kraken2_report.txt"
ktImportTaxonomy -tax "$KRONA_TAX_PATH" "$KRAKEN_KRONA_DIR/combined_kraken2_report.txt" -o "$OUTPUT_DIR/combined_krona.html"

echo "Finished: Taxonomic composition barplot + Krona HTMLs in $OUTPUT_DIR and $KRAKEN_KRONA_DIR"
