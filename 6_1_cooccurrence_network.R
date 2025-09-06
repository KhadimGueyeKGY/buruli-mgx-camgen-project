#!/usr/bin/env Rscript

# ==========================================================
# Microbial co-occurrence network (Spearman correlation)
# For 16S rRNA amplicon sequencing data
#
# Author: [Your Name]
# Date: [YYYY-MM-DD]
# ==========================================================

# -----------------------------
#  Load required packages
# -----------------------------
suppressMessages({
  library(tidyverse)   # data wrangling & visualization
  library(Hmisc)       # correlation with p-values
  library(igraph)      # network creation
  library(ggraph)      # network visualization
})

# -----------------------------
#  Parameters (adapt as needed)
# -----------------------------

SITE <- "akonolinga"  # "akonolinga" or "ayos"
SITESET <- "B"     # ayos = A, akonolinga = B


KRKN_FILE <- file.path("..", "output", paste0("16srna_buruli_", SITE, "_set", SITESET),
                       "combined_kraken_cleaned.csv")
METADATA_FILE <- file.path("..", paste0("fastq_pass_", SITE, ".tsv"))
OUT_DIR <- file.path("..", "output", paste0("plot_", SITE, "_set", SITESET))
OUT_IMG <- file.path(OUT_DIR, "cooccurrence_network.png")

# Create output directory if not exists
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# -----------------------------
#  1. Read Kraken counts
# -----------------------------
kraken <- read_csv(KRKN_FILE, show_col_types = FALSE) %>%
  mutate(Reads = as.numeric(Reads))   # ensure numeric

# -----------------------------
#  2. Aggregate duplicates
# -----------------------------
kraken_clean <- kraken %>%
  group_by(Sample, TaxonName) %>%
  summarise(Reads = sum(Reads), .groups = "drop")

# -----------------------------
#  3. Pivot to wide format
# -----------------------------
kraken_wide <- kraken_clean %>%
  pivot_wider(
    names_from = Sample,
    values_from = Reads,
    values_fill = list(Reads = 0)
  ) %>%
  column_to_rownames("TaxonName")

# -----------------------------
#  4. Filter taxa (optional)
# -----------------------------
# Keep taxa with >100 total reads (remove rare taxa)
kraken_wide <- kraken_wide[rowSums(kraken_wide) > 100, ]

# -----------------------------
#  5. Compute Spearman correlation
# -----------------------------
if (nrow(kraken_wide) > 4) {
  
  corr_res   <- rcorr(as.matrix(kraken_wide), type = "spearman")
  corr_matrix <- corr_res$r
  p_matrix    <- corr_res$P
  
  # -----------------------------
  #  6. Extract significant edges
  # -----------------------------
  sig_edges <- which(abs(corr_matrix) > 0.6 & p_matrix < 0.05, arr.ind = TRUE)
  
  edges <- tibble(
    from   = rownames(corr_matrix)[sig_edges[, 1]],
    to     = colnames(corr_matrix)[sig_edges[, 2]],
    weight = corr_matrix[sig_edges]
  ) %>%
    filter(from != to) %>%
    distinct(from, to, .keep_all = TRUE)
  
  # -----------------------------
  #  7. Build co-occurrence network
  # -----------------------------
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # -----------------------------
  #  8. Visualization
  # -----------------------------
  set.seed(42) # reproducibility of layout
  
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.7) +
    geom_node_point(size = 5, color = "steelblue") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, color = "black") +
    scale_edge_color_gradient2(
      low = "#d73027", mid = "#f0f0f0", high = "#1a9850",
      midpoint = 0, name = "Spearman ρ"
    ) +
    theme_void(base_size = 14) +
    ggtitle(paste0("Microbial co-occurrence network\n(", SITE, " set ", SITESET, ")")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  # -----------------------------
  #  9. Save publication-ready plot
  # -----------------------------
  ggsave(OUT_IMG, p, width = 10, height = 8, dpi = 300)
  
  message("✅ Network saved to: ", OUT_IMG)
  
} else {
  message("Not enough taxa (n=", nrow(kraken_wide),
          ") for correlation network. Need >4 taxa.")
}
