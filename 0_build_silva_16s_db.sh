#!/bin/bash
# Build a Kraken2 database from SILVA SSU Ref NR99 for full-length 16S amplicons.
# The previous database held no bacterial sequences, which is why every report was
# archaeal and 96 percent of reads were unclassified.
#
# Usage: bash 0_build_silva_16s_db.sh
#
# Auteur : Khadim Gueye

set -euo pipefail

VERSION="138.2"
DB="${KRAKEN_DB:-$HOME/DATABASE/kraken2_silva${VERSION}_16S}"
THREADS="${THREADS:-48}"
BASE="https://ftp.arb-silva.de/release_${VERSION//./_}/Exports"
FASTA="SILVA_${VERSION}_SSURef_NR99_tax_silva.fasta"
TAXO="tax_slv_ssu_${VERSION}"
export PATH="/usr/lib/kraken2:$PATH"

mkdir -p "$DB"/{data,taxonomy,library}
cd "$DB/data"

[ -s "$FASTA" ] || { wget -c "$BASE/${FASTA}.gz"; gunzip -f "${FASTA}.gz"; }
[ -s "${TAXO}.acc_taxid" ] || { wget -c "$BASE/taxonomy/${TAXO}.acc_taxid.gz"; gunzip -f "${TAXO}.acc_taxid.gz"; }
[ -s "${TAXO}.txt" ] || { wget -c "$BASE/taxonomy/${TAXO}.txt.gz"; gunzip -f "${TAXO}.txt.gz"; }

build_silva_taxonomy.pl "${TAXO}.txt"
cd "$DB"
mv -f data/names.dmp data/nodes.dmp taxonomy/
cp -f "data/${TAXO}.acc_taxid" seqid2taxid.map
# SILVA stores RNA, Kraken2 indexes DNA.
sed -e '/^>/!y/U/T/' "data/$FASTA" > library/silva.fna

kraken2-build --db "$DB" --build --threads "$THREADS"
kraken2-inspect --db "$DB" --threads "$THREADS" | head -40
