# Declaration:

# I, John Alexander Beggs, declare that I have employed a Chat-GPT-3.5, to assist in the creation

# of this .py script. Specifically, I used it to proactively learn how to carry out a four component PCA plot grid for analysis of infection status.

#Part 3: Clustered heat map

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Load in and merge features data and probe data

# Load the gene expression data saved as features.tsv
features_all = pd.read_csv("features.tsv", sep="\t", index_col=0)

# Load probe metadata (GPL570-55999.txt) containing the gene annotations: gene symbol and GO terms
probe_metadata = pd.read_csv("GPL570-55999.txt", sep="\t", index_col=0, comment="#", low_memory=False)

# Extract only relevant columns: Gene Symbol and GO terms
probe_metadata = probe_metadata[["Gene Symbol", "Gene Ontology Biological Process",
                                  "Gene Ontology Cellular Component", "Gene Ontology Molecular Function"]]

# Merge the gene expression data with the probe metadata on the gene ID (index)
merged_data = features_all.merge(probe_metadata, left_index=True, right_index=True, how="left")

# 2. Filter data to get top 100 most significant genes for both RSV and Flu (lowest p-values)

# create variable to hold top 100 genes based on the lowest p-values for RSV and Flu
rsv_top100 = merged_data.nsmallest(100, "rsv_pvalue")
flu_top100 = merged_data.nsmallest(100, "flu_pvalue")

# Combine top RSV and Flu genes, removing any duplicates using drop_duplicates()
top_genes = pd.concat([rsv_top100, flu_top100]).drop_duplicates()

# 3. Prepare data for the heatmap by extracting numerical data and excluding non-numeric columns

# Extract numerical data for clustered heatmap
heatmap_data = top_genes.iloc[:, :features_all.shape[1] - 6]  # .iloc used to exclude non-numerical columns (last 6) from top genes

# Add GO Gene Symbols as labels for the rows
heatmap_data.index = top_genes["Gene Symbol"]

# Ensure column names in heatmap_data and metadata index match
heatmap_data.columns = heatmap_data.columns.str.strip().str.lower()
sample_metadata = pd.read_csv("meta_data.csv").set_index("Sample_geo_accession")
sample_metadata.index = sample_metadata.index.str.strip().str.lower()

# Reindex metadata to match the heatmap_data columns showing infection status (flu, rsv, control)
infection_labels = sample_metadata['infection_status'].reindex(heatmap_data.columns)

# Ensure no there are no unknown entries appearing
if infection_labels.isnull().any():
    print("WARNING: Some samples are unmatched. Check your metadata.")
infection_labels = infection_labels.fillna("none")  # Replace the unmatched genes with "none"

# Define the colours for infection statuses
infection_palette = {"rsv": "red", "influenza": "blue", "none": "green"}
infection_colours = [infection_palette[status] for status in infection_labels]

# Create heatmap
g = sns.clustermap(
    heatmap_data,
    cmap="turbo",
    figsize=(28, 15),  # Increase width for visibility
    yticklabels=True,
    xticklabels=True,
    dendrogram_ratio=(0.15, 0.15),
    cbar_pos=(0.02, 0.8, 0.03, 0.18),
    col_colors=infection_colours
)

# Rotate x-axis and y-axis labels
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=2)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=2)

# Add title to heatmap
plt.title("Clustered Heatmap of Top 100 RSV and Flu Genes", fontsize=16, pad=20)

# Add the x-axis Samples label
g.ax_heatmap.set_xlabel("Samples", fontsize=12, labelpad=15)

# Move the legend to the side of plot
legend_labels = [plt.Line2D([0], [0], color=color, lw=4) for color in infection_palette.values()]
plt.legend(
    legend_labels, infection_palette.keys(),
    title="Infection Status", loc='center left',
    bbox_to_anchor=(1.3, 0.5), fontsize=10
)

# Save the figure as a pdf file for zoom
g.savefig("clustered_heatmap_top_100.pdf", bbox_inches="tight")
plt.close()

print("Clustered heatmap saved as 'clustered_heatmap_top_100.pdf'.")


# 5. Generation of boxplots of top 100 genes for flu v RSV v control

melted_data = heatmap_data.reset_index().melt(id_vars="Gene Symbol", var_name="Sample", value_name="Log Ratio")
melted_data = melted_data.merge(sample_metadata, left_on="Sample", right_index=True)

# Create a grid of boxplots
plt.figure(figsize=(25, 17))
g = sns.catplot(
    data=melted_data,
    x="Gene Symbol",
    y="Log Ratio",
    hue="infection_status",
    kind="box",
    aspect=3,
    height=5,
    palette="Set2",
    linewidth=0.5,
    showfliers=False  # Remove outliers
)

# Customize boxplot to make longer with smaller thinner plots, making it less cluttered
g.set_xticklabels(rotation=90, fontsize=6)
g.fig.subplots_adjust(top=0.9, bottom=0.4)  # Increased spacing at bottom
g.fig.suptitle("Boxplots of Log Ratios for Top 100 Flu and RSV Genes Compared with Controls", fontsize=16)
g.set_axis_labels("Gene Symbol", "Log Ratio")  # Add axis labels, gene symbol on x and log ratio on y

plt.savefig("boxplots_top_100_genes.pdf") #save top 100 genes boxplot figure
plt.close()

print("Boxplots saved as 'boxplots_top_100_genes.pdf'.")

#6. Generation of boxplots of top 100 genes for flu v RSV v control

# Filter for top 20 Flu and RSV genes by smallest p value
top_20_flu_genes = flu_top100.nsmallest(20, "flu_pvalue")
top_20_rsv_genes = rsv_top100.nsmallest(20, "rsv_pvalue")

# Combine the two datasets using panda concat function
top_20_genes_combined = pd.concat([top_20_flu_genes, top_20_rsv_genes]).drop_duplicates() #drop_duplicates function removes duplicate genes

# Prepare data for boxplots
melted_data_top_20 = heatmap_data.loc[top_20_genes_combined["Gene Symbol"]].reset_index().melt(
    id_vars="Gene Symbol", var_name="Sample", value_name="Log Ratio"
)

# Merge with metadata to include infection status
melted_data_top_20 = melted_data_top_20.merge(sample_metadata, left_on="Sample", right_index=True)

# Create boxplots for top 20 Flu and RSV genes
plt.figure(figsize=(25, 17))
g = sns.catplot(
    data=melted_data_top_20,
    x="Gene Symbol",
    y="Log Ratio",
    hue="infection_status",
    kind="box",
    aspect=3,
    height=5,
    palette="Set2",
    linewidth=0.5,
    showfliers=False  # Remove outliers
)

# Customize boxplot appearance
g.fig.subplots_adjust(top=0.9, bottom=0.4)  # Increases bottom spacing
g.fig.suptitle("Boxplots of Log Ratios for Top 20 Flu and RSV Genes Compared with Controls", fontsize=16)
g.set_axis_labels("Gene Symbol", "Log Ratio")  # Add axis labels, gene symbol on x and log ratio on y

# x-axis labels are rotated
for ax in g.axes.flat:  # Iterate over axes in the grid
    for label in ax.get_xticklabels():
        label.set_rotation(90)  # Rotate x-axis labels for visibility
        label.set_horizontalalignment("center")
    ax.tick_params(axis='x', labelsize=6)  # Set font size for x-axis labels to make visible

# Save the boxplot figure
plt.savefig("boxplots_top_20_flu_rsv_genes.pdf")
plt.close()

# Save as a CSV for the selected top 20 genes
top_20_genes_combined.to_csv("top_20_flu_rsv_genes.csv", index=False)

print("Boxplots saved as 'boxplots_top_20_flu_rsv_genes.pdf'.")










