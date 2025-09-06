#!/usr/bin/env Rscript

# ==========================================================
# Microbial co-occurrence networks (Spearman correlation)
# Sites: Ayos (set A) and Akonolinga (set B)
# ==========================================================

suppressMessages({
  library(tidyverse)
  library(Hmisc)
  library(igraph)
  library(ggraph)
  library(gridExtra)   # instead of patchwork
})

# -----------------------------
#  Function to build network plot
# -----------------------------
build_network_plot <- function(site, set_label) {
  
  KRKN_FILE <- file.path("..", "output",
                         paste0("16srna_buruli_", site, "_set", set_label),
                         "combined_kraken_cleaned.csv")
  OUT_DIR <- file.path("..", "output", "plots_combined_sites")
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  
  # Read data
  kraken <- read_csv(KRKN_FILE, show_col_types = FALSE) %>%
    mutate(Reads = as.numeric(Reads)) %>%
    group_by(Sample, TaxonName) %>%
    summarise(Reads = sum(Reads), .groups = "drop") %>%
    pivot_wider(names_from = Sample,
                values_from = Reads,
                values_fill = list(Reads = 0)) %>%
    column_to_rownames("TaxonName")
  
  # Filter taxa
  kraken <- kraken[rowSums(kraken) > 100, ]
  
  if (nrow(kraken) <= 4) {
    return(ggplot() +
             theme_void() +
             ggtitle(paste0("Not enough taxa in ", site, " set ", set_label)))
  }
  
  # Spearman correlation
  corr_res <- rcorr(as.matrix(kraken), type = "spearman")
  corr_matrix <- corr_res$r
  p_matrix <- corr_res$P
  
  sig_edges <- which(abs(corr_matrix) > 0.6 & p_matrix < 0.05, arr.ind = TRUE)
  
  edges <- tibble(
    from   = rownames(corr_matrix)[sig_edges[, 1]],
    to     = colnames(corr_matrix)[sig_edges[, 2]],
    weight = corr_matrix[sig_edges]
  ) %>%
    filter(from != to) %>%
    distinct(from, to, .keep_all = TRUE)
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  set.seed(42)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.7) +
    geom_node_point(size = 5, color = "steelblue") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, color = "black") +
    scale_edge_color_gradient2(
      low = "#d73027", mid = "#f0f0f0", high = "#1a9850",
      midpoint = 0, name = "Spearman ρ"
    ) +
    theme_void(base_size = 14) +
    ggtitle(paste0(site, " (set ", set_label, ")")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
  
  return(p)
}

# -----------------------------
#  Build plots for both sites
# -----------------------------
p1 <- build_network_plot("ayos", "A")
p2 <- build_network_plot("akonolinga", "B")

# Arrange side by side with gridExtra
final_plot <- grid.arrange(p1, p2, ncol = 2,
                           top = "Microbial co-occurrence networks (Panels A & B)")

# Save
OUT_IMG <- file.path("..", "output", "plots_combined_sites",
                     "cooccurrence_networks_combined.png")
ggsave(OUT_IMG, final_plot, width = 14, height = 7, dpi = 300)

message("✅ Combined panel plot saved to: ", OUT_IMG)
