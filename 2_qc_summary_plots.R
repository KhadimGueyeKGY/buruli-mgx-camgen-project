# === Load Required Libraries ===
library(tidyverse)
library(readr)
library(cowplot)  # for theme improvements

# === Read Environment Variables from Shell ===
SITE <- "akonolinga"  # Change to "akonolinga" or "ayos"
site <- "B"     # ayos = A, akonolinga = B

# === Construct Paths Dynamically ===
multiqc_raw <- paste0("../output/16srna_buruli_", SITE, "_set", site, "/multiqc_raw/multiqc_data/multiqc_general_stats.txt")
multiqc_filtered <- paste0("../output/16srna_buruli_", SITE, "_set", site, "/multiqc_filtered/multiqc_data/multiqc_general_stats.txt")
output_dir <- paste0("../output/plot_", SITE, "_set", site)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# === Load Raw MultiQC Summary ===
message("Generating QC plots for: ", SITE, " Set ", site, " (Before Filtering)")
qc_raw <- read_tsv(multiqc_raw, comment = "#")
names(qc_raw) <- make.names(names(qc_raw))

# === Plot 1: Total Reads per Sample (Before Filtering) ===
p1 <- ggplot(qc_raw, aes(x = reorder(Sample, -FastQC_mqc.generalstats.fastqc.total_sequences), 
                         y = FastQC_mqc.generalstats.fastqc.total_sequences)) +
  geom_col(fill = "#0072B2", width = 0.7) +
  theme_classic(base_size = 16) +
  labs(title = paste("Total Reads per Sample (Before Filtering) -", SITE),
       x = "Sample", y = "Total Reads") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(output_dir, "total_reads_before_filter.pdf"), plot = p1, width = 10, height = 5)
ggsave(file.path(output_dir, "total_reads_before_filter.png"), plot = p1, width = 10, height = 5, dpi = 300)

# === Plot 2: Read Length Distribution (Before Filtering) - per sample ===
p2 <- ggplot(qc_raw, aes(x = Sample, y = FastQC_mqc.generalstats.fastqc.avg_sequence_length)) +
  geom_violin(fill = "#009E73", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  theme_classic(base_size = 16) +
  labs(title = paste("Read Length Distribution (Before Filtering) -", SITE),
       x = "Sample", y = "Average Read Length") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(output_dir, "read_length_before_filter.pdf"), plot = p2, width = 10, height = 6)
ggsave(file.path(output_dir, "read_length_before_filter.png"), plot = p2, width = 10, height = 6, dpi = 300)

# === Load Filtered MultiQC Summary ===
message("Generating QC plots for: ", SITE, " Set ", site, " (After Filtering)")
qc_filtered <- read_tsv(multiqc_filtered, comment = "#")
names(qc_filtered) <- make.names(names(qc_filtered))

# === Apply 0.1% Cutoff ===
max_reads <- max(qc_filtered$FastQC_mqc.generalstats.fastqc.total_sequences, na.rm = TRUE)
cutoff <- 0.001 * max_reads  # 0.1% of maximum reads
qc_filtered_cut <- qc_filtered %>%
  filter(FastQC_mqc.generalstats.fastqc.total_sequences >= cutoff)

message("Cutoff applied: ", cutoff, " reads (0.1% of max = ", max_reads, ")")
message("Samples retained after cutoff: ", nrow(qc_filtered_cut), "/", nrow(qc_filtered))

# === Plot 3: Total Reads per Sample (After Filtering & Cutoff) ===
p3 <- ggplot(qc_filtered_cut, aes(x = reorder(Sample, -FastQC_mqc.generalstats.fastqc.total_sequences), 
                                  y = FastQC_mqc.generalstats.fastqc.total_sequences)) +
  geom_col(fill = "#D55E00", width = 0.7) +
  theme_classic(base_size = 16) +
  labs(title = paste("Total Reads per Sample (After Filtering, ≥0.1% Cutoff) -", SITE),
       x = "Sample", y = "Total Reads (retained)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(output_dir, "total_reads_after_filter_cutoff.pdf"), plot = p3, width = 10, height = 5)
ggsave(file.path(output_dir, "total_reads_after_filter_cutoff.png"), plot = p3, width = 10, height = 5, dpi = 300)

# === Plot 4: Read Length Distribution (After Filtering & Cutoff) - per sample ===
p4 <- ggplot(qc_filtered_cut, aes(x = Sample, y = FastQC_mqc.generalstats.fastqc.avg_sequence_length)) +
  geom_violin(fill = "#F0E442", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  theme_classic(base_size = 16) +
  labs(title = paste("Read Length Distribution (After Filtering, ≥0.1% Cutoff) -", SITE),
       x = "Sample", y = "Average Read Length") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(output_dir, "read_length_after_filter_cutoff.pdf"), plot = p4, width = 10, height = 6)
ggsave(file.path(output_dir, "read_length_after_filter_cutoff.png"), plot = p4, width = 10, height = 6, dpi = 300)

# === Save filtered metadata for downstream analysis ===
write_tsv(qc_filtered_cut, file.path(output_dir, "qc_filtered_pass_cutoff.tsv"))

message("✅ All publication-ready QC plots generated successfully with cutoff applied.")
