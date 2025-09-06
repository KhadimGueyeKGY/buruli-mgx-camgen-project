#!/usr/bin/env python3

# -------------------------------
# Clinical correlation analysis
# -------------------------------

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math

# -------------------------------
# 0. User parameters
# -------------------------------

SITE = "akonolinga"  # "akonolinga" or "ayos"
SITESET = "B"     # ayos = A, akonolinga = B

KRAKEN_FILE = f"../output/16srna_buruli_{SITE}_set{SITESET}/combined_kraken_cleaned.csv"
METADATA_FILE = f"../fastq_pass_{SITE}.tsv"
OUTDIR = f"../output/plot_{SITE}_set{SITESET}/clinical_corr_pub/"
os.makedirs(OUTDIR, exist_ok=True)

# Clinical variables to analyze
CLINICAL_NUMERIC = ["lesion_size"]  
CLINICAL_CATEGORICAL = ["lesion_category", "treatment_status"]

print(f"Processing site: {SITE}")

# -------------------------------
# 1. Load data
# -------------------------------
kraken = pd.read_csv(KRAKEN_FILE)
metadata = pd.read_csv(METADATA_FILE, sep="\t")

# clean column names: strip spaces, lowercase
metadata.columns = [c.strip().lower() for c in metadata.columns]

print("Metadata columns after cleaning:", metadata.columns.tolist())

# Force Reads to numeric
kraken["Reads"] = pd.to_numeric(kraken["Reads"], errors="coerce").fillna(0)

# -------------------------------
# 2. Pivot kraken (samples x taxa)
# -------------------------------
kraken_wide = kraken.pivot_table(index="Sample", columns="TaxonName",
                                 values="Reads", fill_value=0)
count_mat = kraken_wide.loc[:, (kraken_wide.sum(axis=0) > 0)]

# -------------------------------
# 3. Align metadata with count matrix
# -------------------------------
sample_col = f"16srna_pos_{SITE}".lower()

if sample_col not in metadata.columns:
    raise KeyError(f"Sample ID column '{sample_col}' not found in metadata")

metadata = metadata.rename(columns={sample_col: "sample_id"})
metadata = metadata.set_index("sample_id")

common_samples = metadata.index.intersection(count_mat.index)

metadata_filtered = metadata.loc[common_samples]
count_mat_filtered = count_mat.loc[common_samples]

print(f"Matched {len(common_samples)} samples between kraken and metadata")

# -------------------------------
# 4. Correlation: numeric clinical vars (panel scatterplots)
# -------------------------------
for var in CLINICAL_NUMERIC:
    if var in metadata_filtered.columns:
        top_taxa = count_mat_filtered.sum(axis=0).sort_values(ascending=False).head(20).index
        mat_df = count_mat_filtered[top_taxa].copy()
        mat_df[var] = pd.to_numeric(metadata_filtered[var], errors="coerce")

        # Create bins for X axis

        # correlation matrix
        cor_mat = mat_df.corr(method="pearson")
        plt.figure(figsize=(12, 10))
        sns.heatmap(cor_mat, cmap="coolwarm", center=0, annot=False)
        plt.title(f"Correlation matrix: {var}", fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTDIR, f"{var}_correlation_matrix.png"))
        plt.close()

        # scatter plots in a panel
        bins = [0, 5, 15, float('inf')]
        labels = ["<=5", "5-15", ">15"]
        metadata_filtered[f"{var}_bin"] = pd.cut(pd.to_numeric(metadata_filtered[var], errors="coerce"),
                                         bins=bins,
                                         labels=labels,
                                         include_lowest=True)
        n_taxa = len(top_taxa)
        ncols = 4
        nrows = math.ceil(n_taxa / ncols)

        fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), sharex=True)
        axes = axes.flatten()

        for i, taxon in enumerate(top_taxa):
            ordered_categories = ["≤ 5cm", "5-15cm", "˃ 15cm"]
            metadata_filtered[var] = pd.Categorical(metadata_filtered[var],
                                                    categories=ordered_categories,
                                                    ordered=True)

            sorted_idx = metadata_filtered[var].sort_values().index

            # Vérification
            print(metadata_filtered.loc[sorted_idx, var])
            print(count_mat_filtered.loc[sorted_idx, taxon])


            sns.scatterplot(x=metadata_filtered[var], y=count_mat_filtered[taxon], ax=axes[i])
            sns.regplot(x=pd.to_numeric(metadata_filtered[var], errors="coerce"),
                        y=count_mat_filtered[taxon],
                        scatter=False, color="blue", ax=axes[i])
            axes[i].set_title(taxon, fontsize=10)
            axes[i].set_xlabel(var)
            axes[i].set_ylabel("Abundance")

        # hide empty subplots if any
        for j in range(i+1, len(axes)):
            fig.delaxes(axes[j])

        fig.suptitle(f"Scatter plots: Top taxa vs {var}", fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(os.path.join(OUTDIR, f"{var}_scatter_panel.png"))
        plt.close()
    else:
        print(f"⚠️ {var} not found in metadata")

# -------------------------------
# 5. Correlation: categorical vars (panel boxplots)
# -------------------------------
for var in CLINICAL_CATEGORICAL:
    if var in metadata_filtered.columns:
        top_taxa = count_mat_filtered.sum(axis=0).sort_values(ascending=False).head(10).index

        # Define desired order for categorical variables
        desired_order = ["I", "II", "III"]  # adjust according to your data
        metadata_filtered[var] = pd.Categorical(metadata_filtered[var],
                                                categories=desired_order,
                                                ordered=True)

        n_taxa = len(top_taxa)
        ncols = 4
        nrows = math.ceil(n_taxa / ncols)

        fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), sharey=False)
        axes = axes.flatten()

        for i, taxon in enumerate(top_taxa):
            sns.boxplot(x=metadata_filtered[var], y=count_mat_filtered[taxon], ax=axes[i])
            sns.stripplot(x=metadata_filtered[var], y=count_mat_filtered[taxon],
                          color="black", alpha=0.6, jitter=True, ax=axes[i])
            axes[i].set_title(taxon, fontsize=10)
            axes[i].set_xlabel(var)
            axes[i].set_ylabel("Abundance")

        for j in range(i+1, len(axes)):
            fig.delaxes(axes[j])

        fig.suptitle(f"Boxplots: Top taxa by {var}", fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(os.path.join(OUTDIR, f"{var}_boxplot_panel.png"))
        plt.close()
    else:
        print(f"⚠️ {var} not found in metadata")

print("✅ Clinical correlation analysis completed")
