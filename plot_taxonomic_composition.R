suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# === Options ===
option_list <- list(
  make_option("--bracken_dir", type="character", help="Directory with Bracken outputs"),
  make_option("--metadata", type="character", help="Sample metadata file (TSV)"),
  make_option("--output_dir", type="character", help="Output directory"),
  make_option("--krona_dir", type="character", help="Krona input directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

# === Load Bracken species files ===
files <- list.files(opt$bracken_dir, pattern="*.bracken_species.txt", full.names = TRUE)
sample_names <- gsub(".bracken_species.txt", "", basename(files))

df_list <- lapply(seq_along(files), function(i) {
  read.delim(files[i], sep="\t", header=TRUE) %>%
    select(name, taxonomy_id, fraction_total_reads) %>%
    mutate(Sample = sample_names[i])
})

df <- bind_rows(df_list)

# === Select top taxa across all samples ===
top_taxa <- df %>%
  group_by(name) %>%
  summarise(mean_abundance = mean(fraction_total_reads)) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 20) %>%
  pull(name)

df_plot <- df %>%
  mutate(Taxon = ifelse(name %in% top_taxa, name, "Other")) %>%
  group_by(Sample, Taxon) %>%
  summarise(RelAbundance = sum(fraction_total_reads), .groups="drop") %>%
  group_by(Sample) %>%
  mutate(RelAbundance = RelAbundance / sum(RelAbundance)) %>%
  ungroup()

# === Single combined stacked barplot ===
p <- ggplot(df_plot, aes(x = Sample, y = RelAbundance, fill = Taxon)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Samples", y="Relative abundance (%)", fill="Species") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8))

# Save outputs
ggsave(file.path(opt$output_dir, "Taxonomic_Composition_Barplot.pdf"), p, width=12, height=6)
ggsave(file.path(opt$output_dir, "Taxonomic_Composition_Barplot.png"), p, width=12, height=6, dpi=300)

# === Prepare Krona input per group (all samples concatenated) ===
for (f in files) {
  base <- gsub(".bracken_species.txt", "", basename(f))
  krona_out <- file.path(opt$krona_dir, paste0(base, ".txt"))
  
  df_krona <- read.delim(f, sep="\t", header=TRUE) %>%
    select(fraction_total_reads, taxonomy_id, name)
  
  # Krona input format: <count> <taxonomy>
  write.table(df_krona, file=krona_out, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
