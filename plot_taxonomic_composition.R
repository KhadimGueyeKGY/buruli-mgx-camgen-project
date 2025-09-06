#!/usr/bin/env Rscript

# === Load packages ===
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(readr)
})

# === CLI arguments ===
option_list <- list(
  make_option("--bracken_dir", type="character", help="Directory containing Bracken results"),
  make_option("--metadata", type="character", help="Metadata file (TSV)"),
  make_option("--output_dir", type="character", help="Output directory"),
  make_option("--krona_dir", type="character", help="Directory to store Krona inputs")
)
opt <- parse_args(OptionParser(option_list=option_list))

# === Functions ===
read_bracken <- function(file, level) {
  df <- read.delim(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  df$sample <- gsub(paste0(".bracken_", tolower(level), ".txt"), "", basename(file))
  df$level <- level
  return(df)
}

make_barplot <- function(df, level, outdir) {
  p <- df %>%
    group_by(sample, name) %>%
    summarise(mean_abundance = sum(fraction_total_reads), .groups="drop") %>%
    ggplot(aes(x = sample, y = mean_abundance, fill = name)) +
    geom_bar(stat="identity") +
    labs(title=paste("Taxonomic composition at", level, "level"),
         y="Fraction of total reads", x="Sample") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  
  ggsave(file.path(outdir, paste0("barplot_", tolower(level), ".pdf")),
         plot=p, width=12, height=6)
}

write_krona_inputs <- function(df, level, krona_dir) {
  df %>%
    group_by(sample, name) %>%
    summarise(reads = sum(new_est_reads), .groups="drop") %>%
    group_split(sample) %>%
    walk(function(d) {
      sample <- unique(d$sample)
      out_file <- file.path(krona_dir, paste0(sample, "_", tolower(level), ".txt"))
      write.table(d %>% select(reads, name),
                  file=out_file, sep="\t", quote=FALSE,
                  col.names=FALSE, row.names=FALSE)
    })
}

write_krona_top10 <- function(df, level, krona_dir) {
  df %>%
    group_by(sample, name) %>%
    summarise(reads = sum(new_est_reads), .groups="drop") %>%
    group_by(sample) %>%
    mutate(fraction = reads / sum(reads)) %>%
    filter(fraction >= 0.10) %>%  # keep only taxa ≥10%
    group_split(sample) %>%
    walk(function(d) {
      sample <- unique(d$sample)
      out_file <- file.path(krona_dir, paste0(sample, "_", tolower(level), "_top10.txt"))
      write.table(d %>% select(reads, name),
                  file=out_file, sep="\t", quote=FALSE,
                  col.names=FALSE, row.names=FALSE)
    })
}

# === Main ===
message("=== Loading Bracken outputs ===")
species_files <- list.files(opt$bracken_dir, pattern="bracken_species.txt$", full.names=TRUE)
genus_files   <- list.files(opt$bracken_dir, pattern="bracken_genus.txt$", full.names=TRUE)

df_species <- map_df(species_files, ~read_bracken(.x, "Species"))
df_genus   <- map_df(genus_files,   ~read_bracken(.x, "Genus"))

# Summarise phylum level from genus/species (optional coarse grouping)
df_phylum <- df_species %>%
  separate(name, into=c("genus","species"), sep=" ", fill="right", extra="merge") %>%
  mutate(phylum = ifelse(is.na(genus), name, genus)) %>%
  group_by(sample, phylum) %>%
  summarise(fraction_total_reads = sum(fraction_total_reads),
            new_est_reads = sum(new_est_reads), .groups="drop") %>%
  rename(name = phylum) %>%
  mutate(level="Phylum")

# === Barplots ===
message("=== Generating stacked barplots ===")
make_barplot(df_species, "Species", opt$output_dir)
make_barplot(df_genus,   "Genus",   opt$output_dir)
make_barplot(df_phylum,  "Phylum",  opt$output_dir)

# === Krona input files ===
message("=== Writing Krona input files ===")
dir.create(opt$krona_dir, showWarnings=FALSE, recursive=TRUE)
write_krona_inputs(df_species, "Species", opt$krona_dir)
write_krona_inputs(df_genus,   "Genus",   opt$krona_dir)

# === Krona top 10% taxa ===
message("=== Writing Krona input files (≥10% taxa only) ===")
write_krona_top10(df_species, "Species", opt$krona_dir)
write_krona_top10(df_genus,   "Genus",   opt$krona_dir)

message("✅ Done: barplots + Krona inputs generated")
