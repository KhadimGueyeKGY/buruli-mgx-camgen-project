#!/usr/bin/env Rscript

# ----------------------------------------
# Differential abundance analysis - Buruli vs Non-Buruli
# Extended with heatmap improvements, Gram classification, and summary stats
# ----------------------------------------

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(ggpubr)
})

# ----------------------------------------
# 0. User parameters
# ----------------------------------------
SITE <- "ayos"   # "akonolinga" or "ayos"
SITESET <- "A"         # "A" for Ayos, "B" for Akonolinga

KRAKEN_FILE <- paste0("../output/16srna_buruli_", SITE, "_set", SITESET, "/combined_kraken_cleaned.csv")
METADATA_FILE <- paste0("../fastq_pass_", SITE, ".tsv")
UB_COL <- paste0("UB_status_", SITE)  # Buruli status column

OUT_DIR <- paste0("../output/plot_", SITE, "_set", SITESET, "/differential_abundance/")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------
# 1. Load data
# ----------------------------------------
cat("🔍 Processing site:", SITE, "\n")

kraken <- read_csv(KRAKEN_FILE, col_types = cols()) %>%
  mutate(Reads = as.numeric(Reads))

metadata <- read_tsv(METADATA_FILE, col_types = cols())

# ----------------------------------------
# 2. Handle duplicate taxa
# ----------------------------------------
kraken <- kraken %>%
  group_by(Sample, TaxonName) %>%
  summarise(Reads = sum(Reads), .groups = "drop")

# ----------------------------------------
# 3. Pivot Kraken to wide format
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
cat("✅ Kraken table dimensions:", dim(count_mat), "\n")

# ----------------------------------------
# 4. Align metadata with count matrix
# ----------------------------------------
sample_col <- paste0("16srna_pos_", SITE)
if (!sample_col %in% colnames(metadata)) stop(paste("❌ Column", sample_col, "not found"))

metadata <- metadata %>%
  mutate(Sample = .data[[sample_col]])

metadata_filtered <- metadata %>%
  filter(Sample %in% colnames(count_mat)) %>%
  column_to_rownames("Sample")

common_samples <- intersect(colnames(count_mat), rownames(metadata_filtered))
if (length(common_samples) == 0) stop("❌ No matching samples between metadata and Kraken table")

count_mat <- count_mat[, common_samples, drop = FALSE]
metadata_filtered <- metadata_filtered[common_samples, , drop = FALSE]
cat("✅ Number of matching samples:", length(common_samples), "\n")

# ----------------------------------------
# 5. Encode condition
# ----------------------------------------
metadata_filtered$condition <- factor(
  ifelse(metadata_filtered[[UB_COL]] == "positive", "Buruli", "nonBuruli"),
  levels = c("nonBuruli", "Buruli")
)
cat("✅ Condition encoded\n")

# ----------------------------------------
# 6. Create DESeq2 object and filter low-prevalence taxa
# ----------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = metadata_filtered,
                              design = ~ condition)

# Filter taxa present in at least 2 samples
keep_taxa <- rowSums(counts(dds) > 0) >= 2
dds <- dds[keep_taxa, ]

# Run DESeq2 using poscounts
dds <- DESeq(dds, sfType = "poscounts")
cat("✅ DESeq2 analysis completed\n")

# ----------------------------------------
# 7. Extract differential abundance results
# ----------------------------------------
res <- results(dds, contrast = c("condition", "Buruli", "nonBuruli")) %>%
  as.data.frame() %>%
  rownames_to_column("Taxon") %>%
  mutate(
    padj = ifelse(is.na(padj), 1, padj),
    Significant = padj < 0.05 & abs(log2FoldChange) > 1
  )

# ----------------------------------------
# 8. Volcano plot
# ----------------------------------------
volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(title = paste("Volcano plot -", SITE),
       x = "Log2 Fold Change", y = "-log10 adjusted p-value") +
  theme_minimal()

#ggsave(file.path(OUT_DIR, "volcano_plot.png"), volcano, width = 7, height = 5)

# ----------------------------------------
# 9. Heatmap of significant taxa
# ----------------------------------------
top_taxa <- res %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20) %>%   # keep top 20 for clearer view
  pull(Taxon)

mat <- counts(dds, normalized = TRUE)[top_taxa, ]

pheatmap(log1p(mat),
         annotation_col = metadata_filtered["condition"],
         main = paste("Top taxa -", SITE),
         labels_col = colnames(mat),  # sample codes on X
         filename = file.path(OUT_DIR, "heatmap_top_taxa.png"))

# ----------------------------------------
# 10. Separate taxa by condition (Buruli-specific, nonBuruli-specific, shared)
# ----------------------------------------
buruli_taxa <- res %>% filter(Significant & log2FoldChange > 1) %>% pull(Taxon)
nonburuli_taxa <- res %>% filter(Significant & log2FoldChange < -1) %>% pull(Taxon)
shared_taxa <- res %>% filter(!Significant) %>% pull(Taxon)

cat("Buruli-specific taxa:", length(buruli_taxa), "\n")
cat("nonBuruli-specific taxa:", length(nonburuli_taxa), "\n")
cat("Shared taxa:", length(shared_taxa), "\n")

# ----------------------------------------
# 11. Gram-positive vs Gram-negative analysis (at genus level)
# ----------------------------------------
# Rough classification dictionary (extend if needed)
gram_dict <- tibble(
  Genus = c("Staphylococcus", "Bacillus", "Clostridium", "Lactobacillus", "Mycobacterium", "Escherichia", "Pseudomonas", "Salmonella"),
  Gram = c("Positive", "Positive", "Positive", "Positive", "Positive", "Negative", "Negative", "Negative")
)

kraken_genus <- kraken %>%
  mutate(Genus = word(TaxonName, 1)) %>%   # first word = genus
  group_by(Sample, Genus) %>%
  summarise(Reads = sum(Reads), .groups = "drop") %>%
  left_join(metadata_filtered %>% rownames_to_column("Sample"), by = "Sample") %>%
  left_join(gram_dict, by = "Genus") %>%
  mutate(Gram = ifelse(is.na(Gram), "Unknown", Gram))

gram_plot <- kraken_genus %>%
  group_by(condition, Gram) %>%
  summarise(TotalReads = sum(Reads), .groups = "drop") %>%
  ggplot(aes(x = condition, y = TotalReads, fill = Gram)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = paste("Gram-positive vs Gram-negative bacteria -", SITE),
       x = "Condition", y = "Total Reads")

#ggsave(file.path(OUT_DIR, "gram_positive_negative.png"), gram_plot, width = 7, height = 5)

# ----------------------------------------
# 12. Summary stats: species and reads per condition
# ----------------------------------------
summary_stats <- kraken %>%
  left_join(metadata_filtered %>% rownames_to_column("Sample"), by = "Sample") %>%
  group_by(condition) %>%
  summarise(
    TotalReads = sum(Reads),
    TotalSpecies = n_distinct(TaxonName)
  )

write_csv(summary_stats, file.path(OUT_DIR, "summary_species_reads.csv"))
cat("✅ Summary stats saved\n")

# ----------------------------------------
# 13. Panel of main plots
# ----------------------------------------
panel <- ggarrange(volcano, gram_plot, ncol = 2, labels = c("A", "B"))
ggsave(file.path(OUT_DIR, "panel_volcano_gram.png"), panel, width = 12, height = 6, dpi = 300)

cat("✅ All plots generated in:", OUT_DIR, "\n")
