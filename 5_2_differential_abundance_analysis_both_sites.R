#!/usr/bin/env Rscript

# ----------------------------------------
# Differential abundance analysis - Combined Sites (Ayos + Akonolinga)
# With site adjustment, Buruli effect, and direct site comparison
# ----------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(ggpubr)
})

# ----------------------------------------
# 0. Input files
# ----------------------------------------
sites <- c("ayos", "akonolinga")
site_sets <- c("A", "B")

kraken_files <- paste0("../output/16srna_buruli_", sites, "_set", site_sets, "/combined_kraken_cleaned.csv")
metadata_files <- paste0("../fastq_pass_", sites, ".tsv")
ub_cols <- paste0("UB_status_", sites)
sample_cols <- paste0("16srna_pos_", sites)

OUT_DIR <- "../output/plots_combined_sites/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------
# 1. Load Kraken + metadata for both sites
# ----------------------------------------
all_kraken <- list()
all_metadata <- list()

for (i in seq_along(sites)) {
  site <- sites[i]
  kraken_file <- kraken_files[i]
  metadata_file <- metadata_files[i]
  ub_col <- ub_cols[i]
  sample_col <- sample_cols[i]
  
  # Kraken
  kraken <- read_csv(kraken_file, col_types = cols()) %>%
    mutate(Reads = as.numeric(Reads),
           Site = site)
  
  # Metadata
  metadata <- read_tsv(metadata_file, col_types = cols()) %>%
    mutate(Sample = .data[[sample_col]],
           condition = ifelse(.data[[ub_col]] == "positive", "Buruli", "nonBuruli"),
           site = site)
  
  all_kraken[[site]] <- kraken
  all_metadata[[site]] <- metadata
}

kraken <- bind_rows(all_kraken) %>%
  group_by(Sample, TaxonName) %>%
  summarise(Reads = sum(Reads), .groups = "drop")

metadata <- bind_rows(all_metadata) %>%
  distinct(Sample, .keep_all = TRUE)

# ----------------------------------------
# 2. Pivot to wide format
# ----------------------------------------
kraken_wide <- kraken %>%
  pivot_wider(
    names_from = Sample,
    values_from = Reads,
    values_fill = 0
  ) %>%
  column_to_rownames("TaxonName") %>%
  as.matrix()

count_mat <- round(kraken_wide)
cat("✅ Combined Kraken table dimensions:", dim(count_mat), "\n")

# ----------------------------------------
# 3. Align metadata
# ----------------------------------------
metadata_filtered <- metadata %>%
  filter(Sample %in% colnames(count_mat)) %>%
  column_to_rownames("Sample")

common_samples <- intersect(colnames(count_mat), rownames(metadata_filtered))
if (length(common_samples) == 0) stop("❌ No matching samples")

count_mat <- count_mat[, common_samples, drop = FALSE]
metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]

metadata_filtered$condition <- factor(metadata_filtered$condition,
                                      levels = c("nonBuruli", "Buruli"))
metadata_filtered$site <- factor(metadata_filtered$site)

cat("✅ Samples after alignment:", length(common_samples), "\n")

# ----------------------------------------
# 4. DESeq2 with site adjustment
# ----------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = metadata_filtered,
                              design = ~ site + condition)

# Filter taxa
keep_taxa <- rowSums(counts(dds) > 0) >= 2
dds <- dds[keep_taxa, ]

dds <- DESeq(dds, sfType = "poscounts")
cat("✅ DESeq2 completed\n")

# Results
res_buruli <- results(dds, contrast = c("condition", "Buruli", "nonBuruli")) %>%
  as.data.frame() %>%
  rownames_to_column("Taxon") %>%
  mutate(
    padj = ifelse(is.na(padj), 1, padj),
    Significant = padj < 0.05 & abs(log2FoldChange) > 1
  )

res_site <- results(dds, contrast = c("site", "akonolinga", "ayos")) %>%
  as.data.frame() %>%
  rownames_to_column("Taxon") %>%
  mutate(
    padj = ifelse(is.na(padj), 1, padj),
    Significant = padj < 0.05 & abs(log2FoldChange) > 1
  )

# ----------------------------------------
# 5. Volcano plots
# ----------------------------------------
volcano_buruli <- ggplot(res_buruli, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = "Volcano - Buruli vs nonBuruli (adjusted for site)") +
  theme_minimal()

volcano_site <- ggplot(res_site, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
  labs(title = "Volcano - Akonolinga vs Ayos") +
  theme_minimal()

ggsave(file.path(OUT_DIR, "volcano_buruli.png"), volcano_buruli, width = 7, height = 5)
ggsave(file.path(OUT_DIR, "volcano_site.png"), volcano_site, width = 7, height = 5)

# ----------------------------------------
# 6. Heatmaps
# ----------------------------------------
top_buruli <- res_buruli %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20) %>%
  pull(Taxon)

top_site <- res_site %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20) %>%
  pull(Taxon)

mat_buruli <- counts(dds, normalized = TRUE)[top_buruli, ]
mat_site <- counts(dds, normalized = TRUE)[top_site, ]

pheatmap(log1p(mat_buruli),
         annotation_col = metadata_filtered[, c("condition", "site")],
         main = "Top taxa - Buruli vs nonBuruli",
         filename = file.path(OUT_DIR, "heatmap_buruli.png"))

pheatmap(log1p(mat_site),
         annotation_col = metadata_filtered[, c("condition", "site")],
         main = "Top taxa - Akonolinga vs Ayos",
         filename = file.path(OUT_DIR, "heatmap_site.png"))

# ----------------------------------------
# 7. Gram classification
# ----------------------------------------
gram_dict <- tibble(
  Genus = c("Staphylococcus", "Bacillus", "Clostridium", "Lactobacillus", "Mycobacterium", 
            "Escherichia", "Pseudomonas", "Salmonella"),
  Gram = c("Positive", "Positive", "Positive", "Positive", "Positive", "Negative", "Negative", "Negative")
)

kraken_genus <- kraken %>%
  mutate(Genus = word(TaxonName, 1)) %>%
  group_by(Sample, Genus) %>%
  summarise(Reads = sum(Reads), .groups = "drop") %>%
  left_join(metadata_filtered %>% rownames_to_column("Sample"), by = "Sample") %>%
  left_join(gram_dict, by = "Genus") %>%
  mutate(Gram = ifelse(is.na(Gram), "Unknown", Gram))

gram_plot <- kraken_genus %>%
  group_by(condition, site, Gram) %>%
  summarise(TotalReads = sum(Reads), .groups = "drop") %>%
  ggplot(aes(x = condition, y = TotalReads, fill = Gram)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~site) +
  theme_minimal() +
  labs(title = "Gram-positive vs Gram-negative by site",
       x = "Condition", y = "Total Reads")

ggsave(file.path(OUT_DIR, "gram_combined.png"), gram_plot, width = 9, height = 5)

# ----------------------------------------
# 8. Summary statistics
# ----------------------------------------
summary_stats <- kraken %>%
  left_join(metadata_filtered %>% rownames_to_column("Sample"), by = "Sample") %>%
  group_by(site, condition) %>%
  summarise(
    TotalReads = sum(Reads),
    TotalSpecies = n_distinct(TaxonName),
    .groups = "drop"
  )

write_csv(summary_stats, file.path(OUT_DIR, "summary_species_reads_combined.csv"))

cat("✅ Combined site analysis completed. Outputs in:", OUT_DIR, "\n")
