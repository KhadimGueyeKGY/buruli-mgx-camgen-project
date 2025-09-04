# ==============================
# 4_alpha_beta_diversity.R
# 16S Metagenomics - Alpha and Beta Diversity
# ==============================

# === Load Required Libraries ===
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)

# === Set input/output paths ===
SITE <- "akonolinga"  # "akonolinga" or "ayos"
site <- "B"     # ayos = A, akonolinga = B

kraken_reports_dir <- paste0("../output/16srna_buruli_", SITE, "_set", site, "/kraken2_reports/")
combined_csv_file <- paste0("../output/16srna_buruli_", SITE, "_set", site, "/combined_kraken_cleaned.csv")
output_dir <- paste0("../output/16srna_buruli_", SITE, "_set", site, "/4_alpha_beta_diversity")
plot_dir <- paste0("../output/plot_", SITE, "_set", site)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================
# 1️⃣ Combine Kraken2 Bracken Species reports
# ==============================
message("Combining Kraken2 reports...")

# List all report files
report_files <- list.files(kraken_reports_dir, pattern = "*.report_bracken_species.txt", full.names = TRUE)

# Read and combine all reports
kraken_list <- lapply(report_files, function(f) {
  read_table(f, col_names = c("percentage", "clade_reads", "taxon_reads", "rank", "taxid", "name")) %>%
    rename(TaxonName = name, Reads = taxon_reads) %>%
    mutate(Reads = as.numeric(Reads),
           Sample = gsub(".report_bracken_species.txt", "", basename(f)))
})

kraken_combined <- bind_rows(kraken_list) %>%
  select(Sample, TaxonName, Reads)

# Save combined CSV
write_csv(kraken_combined, combined_csv_file)
message("✅ Combined Kraken2 CSV saved: ", combined_csv_file)

# ==============================
# 2️⃣ Create phyloseq object
# ==============================
message("Creating phyloseq object...")

# Abundance matrix
abundance_matrix <- kraken_combined %>%
  group_by(TaxonName, Sample) %>%
  summarise(Reads = sum(Reads), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Reads, values_fill = 0) %>%
  column_to_rownames("TaxonName") %>%
  as.matrix()

otu_table_ps <- otu_table(abundance_matrix, taxa_are_rows = TRUE)

# Sample metadata (update Group info if available)
sample_metadata <- data.frame(SampleID = colnames(abundance_matrix))
rownames(sample_metadata) <- sample_metadata$SampleID
sample_metadata$Group <- sample_metadata$SampleID  # Replace with actual group/lesion info

sample_data_ps <- sample_data(sample_metadata)

physeq_object <- phyloseq(otu_table_ps, sample_data_ps)

# ==============================
# 3️⃣ Alpha Diversity
# ==============================
message("Calculating alpha diversity...")

alpha_df <- estimate_richness(physeq_object, measures = c("Shannon", "Simpson", "Chao1")) %>%
  rownames_to_column("SampleID") %>%
  left_join(sample_metadata, by = "SampleID")

# Alpha diversity plots (boxplot + violin)
alpha_plot_shannon <- ggplot(alpha_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.15, size = 2) +
  theme_minimal(base_size = 16) +
  labs(title = paste("Alpha Diversity - Shannon Index (", SITE, ")", sep = ""),
       x = "Group", y = "Shannon Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(file.path(plot_dir, paste0("4_alpha_shannon_violin_boxplot_", SITE, "_set", site, ".png")),
       plot = alpha_plot_shannon, width = 10, height = 7, dpi = 300)

# ==============================
# 4️⃣ Beta Diversity
# ==============================
message("Calculating beta diversity (Bray-Curtis, PCoA)...")

# Bray-Curtis distance and PCoA ordination
bray_dist <- phyloseq::distance(physeq_object, method = "bray")
ordination_bray <- ordinate(physeq_object, method = "PCoA", distance = bray_dist)

# Beta diversity plot
beta_plot <- plot_ordination(physeq_object, ordination_bray, color = "Group") +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal(base_size = 16) +
  labs(title = paste("Beta Diversity - Bray Curtis PCoA (", SITE, ")", sep = ""),
       x = "PCoA1", y = "PCoA2") +
  theme(legend.position = "right")

ggsave(file.path(plot_dir, paste0("4_beta_bray_pcoa_", SITE, "_set", site, ".png")),
       plot = beta_plot, width = 10, height = 7, dpi = 300)

# ==============================
# Done
# ==============================
message("✅ Alpha and Beta diversity plots generated successfully in: ", output_dir)
