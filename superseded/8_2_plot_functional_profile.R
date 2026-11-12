library(tidyverse)
library(reshape2)

# -------------------------------
# Parameters
# -------------------------------
input_file <- "../output/functional_profile_ayos_setA/combined_pathabundance_relab.tsv"
metadata_file <- "../fastq_pass_ayos.tsv"
outdir <- "../output/functional_profile_ayos_setA/plots/"
dir.create(outdir, showWarnings = FALSE)

# -------------------------------
# 1. Load data
# -------------------------------
func_df <- read_tsv(input_file)
metadata <- read_tsv(metadata_file)

# Clean column names
colnames(metadata) <- tolower(trimws(colnames(metadata)))
sample_col <- "16srna_pos_ayos"
metadata <- metadata %>% rename(sample_id = !!sample_col)

# -------------------------------
# 2. Filter top pathways
# -------------------------------
func_df_long <- func_df %>%
    pivot_longer(-Name, names_to = "sample_id", values_to = "abundance") %>%
    group_by(Name) %>%
    summarise(mean_abundance = mean(abundance, na.rm = TRUE)) %>%
    top_n(20, mean_abundance) %>%
    pull(Name)

func_top <- func_df %>% filter(Name %in% func_df_long)

# -------------------------------
# 3. Merge with metadata (group: Buruli vs Non-Buruli)
# -------------------------------
func_long <- func_top %>%
    pivot_longer(-Name, names_to = "sample_id", values_to = "abundance") %>%
    left_join(metadata %>% select(sample_id, ub_status_ayos), by = "sample_id")

# -------------------------------
# 4. Plot: Barplot (stacked by group)
# -------------------------------
ggplot(func_long, aes(x = sample_id, y = abundance, fill = Name)) +
    geom_bar(stat="identity") +
    facet_wrap(~ub_status_ayos, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
    labs(y="Relative Abundance", x="Sample", fill="Pathway") +
    ggsave(file.path(outdir, "top20_pathways_barplot.png"), width=12, height=6)

# -------------------------------
# 5. Plot: Bubble plot
# -------------------------------
ggplot(func_long, aes(x=Name, y=ub_status_ayos, size=abundance, color=ub_status_ayos)) +
    geom_point(alpha=0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
    labs(x="Pathway", y="Group", size="Rel. Abundance") +
    ggsave(file.path(outdir, "top20_pathways_bubbleplot.png"), width=12, height=6)
