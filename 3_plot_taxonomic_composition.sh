#!/bin/bash
set -euo pipefail

# === CONFIGURATION ===
SITE="akonolinga"         # "akonolinga" ou "ayos"
SITESET="B"         # "A" ou "B"
BRACKEN_DIR="../output/16srna_buruli_${SITE}_set${SITESET}/bracken_results"
METADATA="../fastq_pass_${SITE}.tsv"
OUTPUT_DIR="../output/plot_${SITE}_set${SITESET}"
KRONA_DIR="$OUTPUT_DIR/krona_species"
KRAKEN_DIR="../output/16srna_buruli_${SITE}_set${SITESET}/kraken2_reports"
KRAKEN_KRONA_DIR="$OUTPUT_DIR/krona_plots"
NCORES=4

# === Create directories ===
mkdir -p "$OUTPUT_DIR"
mkdir -p "$KRONA_DIR"
mkdir -p "$KRAKEN_KRONA_DIR"

echo "=== Step 1: Running R script to generate stacked barplots + Krona input files ==="
Rscript plot_taxonomic_composition.R \
    --bracken_dir "$BRACKEN_DIR" \
    --metadata "$METADATA" \
    --output_dir "$OUTPUT_DIR" \
    --krona_dir "$KRONA_DIR" \

echo "=== Step 2: Generating Krona plots per sample from Bracken ==="
for f in "$KRONA_DIR"/*.txt; do
    base=$(basename "$f" .txt)
    ktImportText -o "$KRONA_DIR/${base}.html" "$f"
done

echo "=== Step 3: Generating combined Krona plot from all Kraken2 reports ==="
KRONA_TAX_PATH="Krona/KronaTools/taxonomy"

# Concatenate all Kraken2 report files and generate combined Krona plot
cat "$KRAKEN_DIR"/*.txt > "$KRAKEN_KRONA_DIR/combined_kraken2_report.txt"
ktImportTaxonomy -tax "$KRONA_TAX_PATH" \
    "$KRAKEN_KRONA_DIR/combined_kraken2_report.txt" \
    -o "$OUTPUT_DIR/combined_krona_all.html"

echo "=== Step 4: Generating Krona plots from R-prepared inputs (all + 10% cutoff) ==="
# Full Krona (all organisms)
if [[ -f "$KRONA_DIR/Krona_${SITE}_all.txt" ]]; then
    ktImportText -o "$OUTPUT_DIR/Krona_${SITE}_all.html" "$KRONA_DIR/Krona_${SITE}_all.txt"
fi

# Krona filtered at 10% cutoff (publication-ready)
if [[ -f "$KRONA_DIR/Krona_${SITE}_10pct.txt" ]]; then
    ktImportText -o "$OUTPUT_DIR/Krona_${SITE}_10pct.html" "$KRONA_DIR/Krona_${SITE}_10pct.txt"
fi

echo "✅ Finished:"
echo "   - Stacked barplots at species/genus/phylum: $OUTPUT_DIR"
echo "   - Per-sample Krona plots: $KRONA_DIR"
echo "   - Combined Kraken2 Krona plot: $OUTPUT_DIR/combined_krona_all.html"
echo "   - Site Krona plots (all + 10% cutoff): $OUTPUT_DIR/Krona_${SITE}_*.html"
