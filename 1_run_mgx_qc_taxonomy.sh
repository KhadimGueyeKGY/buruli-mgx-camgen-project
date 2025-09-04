#!/bin/bash
set -e

# === CONFIGURATION ===
SITE="ayos"   # Change to "akonolinga" or "ayos"
site="A"            # ayos = A, akonolinga = B
KRAKEN_DB="$(realpath ~/kraken2_db)"
BRACKEN_DB="$KRAKEN_DB"   # Bracken database built from same Kraken2 DB
NCORES=4

# === Input and Output Directories ===
input_dir="../16srna_metagnomics_buruli_2025/16srna_buruli_${SITE}_set${site}/fastq_pass/all_fastq"
output_dir="../output/16srna_buruli_${SITE}_set${site}"
metadata_tsv="../fastq_pass_${SITE}.tsv"

# === Create Output Directories ===
mkdir -p "$output_dir"/{renamed_fastq,raw_fastqc,filtered_fastq,filtered_fastqc,multiqc_raw,multiqc_filtered,kraken2_results,kraken2_reports,bracken_results,krona_plots}

# === Step 0: Create barcode → sample_id mapping from metadata ===
declare -A sample_map
while IFS=$'\t' read -r site sample_id barcode ub_status lesion_size lesion_category treatment; do
  sample_map["$barcode"]="$sample_id"
done < <(tail -n +2 "$metadata_tsv")

# === Step 1–4 (Renaming, QC, Filtering, MultiQC) ===
echo "Renaming FASTQ files..."
for f in "$input_dir"/barcode*.fastq.gz; do
  barcode=$(basename "$f" .fastq.gz)
  sample_id="${sample_map[$barcode]}"
  if [[ -n "$sample_id" ]]; then
    cp "$f" "$output_dir/renamed_fastq/${sample_id}.fastq.gz"
    echo "Renamed $barcode.fastq.gz → ${sample_id}.fastq.gz"
  else
    echo "Warning: No mapping found for $barcode — file not renamed"
  fi
done

echo "Running FastQC on raw files..."
fastqc "$output_dir/renamed_fastq"/*.fastq.gz -o "$output_dir/raw_fastqc" --threads $NCORES
multiqc "$output_dir/raw_fastqc" -o "$output_dir/multiqc_raw"

echo "Filtering FASTQ files..."
for f in "$output_dir/renamed_fastq"/*.fastq.gz; do
  base=$(basename "$f" .fastq.gz)
  zcat "$f" | NanoFilt -q 10 -l 1000 | gzip > "$output_dir/filtered_fastq/${base}.filtered.fastq.gz"
done

echo "Running FastQC on filtered files..."
fastqc "$output_dir/filtered_fastq"/*.filtered.fastq.gz -o "$output_dir/filtered_fastqc" --threads $NCORES
multiqc "$output_dir/filtered_fastqc" -o "$output_dir/multiqc_filtered"

# === Step 5: Kraken2 Classification + Bracken ===
echo "Running Kraken2 + Bracken..."
for fq in "$output_dir/filtered_fastq"/*.filtered.fastq.gz; do
  base=$(basename "$fq" .filtered.fastq.gz)
  echo "Classifying $base..."

  kraken2 \
    --db "$KRAKEN_DB" \
    --threads $NCORES \
    --gzip-compressed \
    --report "$output_dir/kraken2_reports/${base}.report.txt" \
    --output "$output_dir/kraken2_results/${base}.kraken" \
    "$fq"

  bracken \
    -d "$BRACKEN_DB" \
    -i "$output_dir/kraken2_reports/${base}.report.txt" \
    -o "$output_dir/bracken_results/${base}.bracken_genus.txt" \
    -r 150 -l G -t $NCORES

  bracken \
    -d "$BRACKEN_DB" \
    -i "$output_dir/kraken2_reports/${base}.report.txt" \
    -o "$output_dir/bracken_results/${base}.bracken_species.txt" \
    -r 150 -l S -t $NCORES
done

# === Step 6: Merge Bracken Reports (phylum/genus/species) ===
echo "Merging Bracken outputs..."
python3 Bracken/analysis_scripts/combine_bracken_outputs.py --files "$output_dir"/bracken_results/*.bracken_genus.txt -o "$output_dir/bracken_results/all_samples_genus.txt"
python3 Bracken/analysis_scripts/combine_bracken_outputs.py --files "$output_dir"/bracken_results/*.bracken_species.txt -o "$output_dir/bracken_results/all_samples_species.txt"

# === Step 7: Krona Visualization ===
echo "Generating Krona plots..."
KRONA_TAX_PATH="Krona/KronaTools/taxonomy"

for kr in "$output_dir/kraken2_results"/*.kraken; do
  base=$(basename "$kr" .kraken)
  ktImportTaxonomy -tax "$KRONA_TAX_PATH" "$kr" -o "$output_dir/krona_plots/${base}_krona.html"
done


echo "✅ Taxonomic composition analysis complete for site: $SITE"
echo "✅ Bracken results: $output_dir/bracken_results"
echo "✅ Krona plots: $output_dir/krona_plots"
