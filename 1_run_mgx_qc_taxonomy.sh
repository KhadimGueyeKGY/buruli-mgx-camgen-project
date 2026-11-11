#!/bin/bash
# QC and taxonomic classification of the full length 16S amplicon runs, both sites.
#
# Reads are merged per barcode, filtered with NanoFilt at Q10 and 1000 bp, and
# classified against a SILVA SSU Ref NR99 Kraken2 database. The earlier database
# contained no bacterial sequences, so every report was archaeal and 96 percent of
# reads were unclassified. SILVA stops at genus, so Bracken is run at genus, family
# and phylum, not species.
#
# Usage: bash 1_run_mgx_qc_taxonomy.sh [ayos|akonolinga]
#
# Auteur : Khadim Gueye

set -euo pipefail

ROOT="$HOME/khadim/Francis-job"
RAW_ROOT="$ROOT/data/raw_16s_buruli/16srna_metagnomics_buruli_2025"
OUT_ROOT="$ROOT/data/16s_metagenomics/analysis_2026"
KRAKEN_DB="${KRAKEN_DB:-$HOME/DATABASE/kraken2_silva138.2_16S}"
BRACKEN="$HOME/DATABASE/Bracken/bracken"
NCORES="${NCORES:-48}"
MIN_QUAL=10
MIN_LEN=1000
MIN_READS=10000
BRACKEN_READ_LEN=1500
export PATH="$HOME/.local/bin:$PATH"

declare -A SET=( [ayos]=A [akonolinga]=B )
SITES=("${1:-ayos akonolinga}")
read -r -a SITES <<< "${SITES[*]}"

for SITE in "${SITES[@]}"; do
  s="${SET[$SITE]}"
  in_dir="$RAW_ROOT/16srna_buruli_${SITE}_set${s}/fastq_pass"
  out="$OUT_ROOT/${SITE}_set${s}"
  mkdir -p "$out"/{merged_fastq,filtered_fastq,filtered_fastqc,kraken2_reports,kraken2_results,bracken_results,krona_plots}

  echo "[$(date '+%F %T')] $SITE: merging barcodes"
  find "$in_dir" -mindepth 1 -maxdepth 1 -type d -name 'barcode*' -printf '%f\n' | sort | \
  xargs -P "$NCORES" -I{} bash -c 'cat "$1/$2"/*.fastq.gz > "$3/$2.fastq.gz"' _ "$in_dir" {} "$out/merged_fastq"

  echo "[$(date '+%F %T')] $SITE: filtering Q>=$MIN_QUAL L>=$MIN_LEN"
  export MIN_QUAL MIN_LEN
  find "$out/merged_fastq" -name '*.fastq.gz' -printf '%f\n' | sed 's/.fastq.gz//' | sort | \
  xargs -P "$NCORES" -I{} bash -c 'zcat "$1/$2.fastq.gz" | NanoFilt -q "$MIN_QUAL" -l "$MIN_LEN" | gzip > "$3/$2.filtered.fastq.gz"' \
        _ "$out/merged_fastq" {} "$out/filtered_fastq"

  echo -e "barcode\traw_reads\tfiltered_reads\tkept" > "$out/read_counts.tsv"
  for f in "$out/filtered_fastq"/*.filtered.fastq.gz; do
    b=$(basename "$f" .filtered.fastq.gz)
    raw=$(( $(zcat "$out/merged_fastq/$b.fastq.gz" | wc -l) / 4 ))
    n=$(( $(zcat "$f" | wc -l) / 4 ))
    if (( n < MIN_READS )); then rm -f "$f"; k=no; else k=yes; fi
    echo -e "$b\t$raw\t$n\t$k" >> "$out/read_counts.tsv"
  done

  echo "[$(date '+%F %T')] $SITE: FastQC"
  fastqc "$out/filtered_fastq"/*.filtered.fastq.gz -o "$out/filtered_fastqc" --threads "$NCORES" --quiet

  echo "[$(date '+%F %T')] $SITE: Kraken2"
  for f in "$out/filtered_fastq"/*.filtered.fastq.gz; do
    b=$(basename "$f" .filtered.fastq.gz)
    kraken2 --db "$KRAKEN_DB" --threads "$NCORES" --gzip-compressed \
            --report "$out/kraken2_reports/$b.report.txt" \
            --output "$out/kraken2_results/$b.kraken" "$f"
  done

  echo "[$(date '+%F %T')] $SITE: Bracken"
  for r in "$out/kraken2_reports"/*.report.txt; do
    b=$(basename "$r" .report.txt)
    for lvl in P F G; do
      "$BRACKEN" -d "$KRAKEN_DB" -i "$r" -o "$out/bracken_results/$b.bracken_$lvl.txt" \
                 -w "$out/bracken_results/$b.bracken_$lvl.report.txt" \
                 -r "$BRACKEN_READ_LEN" -l "$lvl" -t 1 || echo "bracken failed at $lvl for $b"
    done
  done

  echo "[$(date '+%F %T')] $SITE: done, output in $out"
done
