# ==============================
# 4_alpha_beta_diversity_combined.R
# ==============================

# === Load Required Libraries ===
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)

# === Load combined CSVs from Ayos and Akonolinga ===
siteA_csv <- "../output/16srna_buruli_ayos_setA/combined_kraken_cleaned.csv"
siteB_csv <- "../output/16srna_buruli_akonolinga_setB/combined_kraken_cleaned.csv"

output_dir <- "../output/plots_combined_sites"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# === Load and merge ===
siteA <- read_csv(siteA_csv) %>% mutate(Site = "Ayos")
siteB <- read_csv(siteB_csv) %>% mutate(Site = "Akonolinga")
combined <- bind_rows(siteA, siteB) %>%
  select(Sample, TaxonName, Reads, Site)

# === Build abundance matrix ===
abundance_matrix <- combined %>%
  group_by(TaxonName, Sample) %>%
  summarise(Reads = sum(Reads), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Reads, values_fill = 0) %>%
  column_to_rownames("TaxonName") %>%
  as.matrix()

otu_ps <- otu_table(abundance_matrix, taxa_are_rows = TRUE)

# Sample metadata
sample_metadata <- combined %>%
  select(Sample, Site) %>%
  distinct() %>%
  column_to_rownames("Sample")

# 🔹 Assign color groups based on sample name prefix
sample_metadata$ColorGroup <- ifelse(grepl("^AK", rownames(sample_metadata)), "AK",
                              ifelse(grepl("^AY", rownames(sample_metadata)), "AY", "Other"))

sample_ps <- sample_data(sample_metadata)

# Build phyloseq object
physeq <- phyloseq(otu_ps, sample_ps)

# Define custom colors
custom_colors <- c("AK" = "#1f78b4",   # Blue
                   "AY" = "#33a02c",   # Green
                   "Other" = "grey50")

# ==============================
# Alpha Diversity (Shannon)
# ==============================
alpha_df <- estimate_richness(physeq, measures = c("Shannon")) %>%
  rownames_to_column("Sample") %>%
  left_join(sample_metadata %>% rownames_to_column("Sample"), by = "Sample")

alpha_plot <- ggplot(alpha_df, aes(x = ColorGroup, y = Shannon, fill = ColorGroup)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.15, size = 2) +
  scale_fill_manual(values = custom_colors) +
  theme_pubr(base_size = 14) +
  labs(title = "Alpha Diversity (Shannon)", x = "Sample Group", y = "Shannon Index")

# ==============================
# Beta Diversity (Bray-Curtis PCoA)
# ==============================
bray_dist <- phyloseq::distance(physeq, method = "bray")
ordination <- ordinate(physeq, method = "PCoA", distance = bray_dist)

beta_plot <- plot_ordination(physeq, ordination, color = "ColorGroup") +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", linetype = "dashed", level = 0.95) +
  scale_color_manual(values = custom_colors) +
  theme_pubr(base_size = 14) +
  labs(title = "Beta Diversity (PCoA, Bray-Curtis)", color = "Sample Group")

# ==============================
# Combine Plots (side by side)
# ==============================
final_plot <- ggarrange(alpha_plot, beta_plot,
                        labels = c("A", "B"),
                        ncol = 2, nrow = 1)

ggsave(file.path(output_dir, "alpha_beta_diversity_combined.png"),
       final_plot, width = 14, height = 6, dpi = 300)

message("✅ Combined alpha+beta diversity plot saved in: ", output_dir)
