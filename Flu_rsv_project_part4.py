# Declaration:

# I, John Alexander Beggs, declare that I have employed a Chat-GPT-3.5, to assist in the creation

# of this .py script. Specifically, I used it to proactively learn how to find overlapping genes (eg intersection function) and also to create an extract
#GO terms function using index slicing. Furthermore, I used Ai to embed my enrichment factor equation within an esay to use function, enabling fast
#calculations of enrichment values for GO terms. I also used it to learn how to add boxes to venn diagrams showing gene symbols/accession codes for top genes


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np

#Part 3 stretch: Overlapping gene venn + GO term analysis

#Venn Diagram: Find overlappping genes for venn diagram

# 1. Load the features.tsv and probe metadata
features = pd.read_csv("features.tsv", sep="\t", index_col=0)
probe_metadata = pd.read_csv("GPL570-55999.txt", sep="\t", index_col=0, comment="#", low_memory=False)

# Keep Gene Symbol and GO terms, BP, CC, MF
probe_metadata = probe_metadata[["Gene Symbol", "Gene Ontology Biological Process", "Gene Ontology Cellular Component", "Gene Ontology Molecular Function"]]

# Merge features file with probe metadata
merged_data = features.merge(probe_metadata, left_index=True, right_index=True, how="left")

#2. Retain top 100 significant genes for RSV and Flu by lowest p value
rsv_top100 = merged_data.nsmallest(100, "rsv_pvalue")
flu_top100 = merged_data.nsmallest(100, "flu_pvalue")

# 3. Use Gene symbol or Accession Number in cases where gene symbol isn't available
rsv_genes = rsv_top100.apply(lambda row: row["Gene Symbol"] if pd.notna(row["Gene Symbol"]) else row.name, axis=1) #use row name (accession number) if gene symbol = NaN
flu_genes = flu_top100.apply(lambda row: row["Gene Symbol"] if pd.notna(row["Gene Symbol"]) else row.name, axis=1)#use row name (accession number) if gene symbol = NaN

#4. Convert into gene sets for creation of Venn diagram
rsv_gene_set = set(rsv_genes)
flu_gene_set = set(flu_genes)

#5. Find unique and overlapping genes (top 10 for each category)
rsv_only_genes = sorted(rsv_gene_set - flu_gene_set)[:10]  # Top 10 RSV genes
flu_only_genes = sorted(flu_gene_set - rsv_gene_set)[:10]  # Top 10 Flu genes
overlapping_genes = sorted(rsv_gene_set & flu_gene_set)[:10]    # Top 10 overlapping genes

# 6. Generate Venn diagram
plt.figure(figsize=(10, 8))
venn = venn2(
    [rsv_gene_set, flu_gene_set],
    set_labels=("Top 100 RSV Genes", "Top 100 Flu Genes"),
    set_colors=("blue", "orange"),
)

# Add gene symbols around the counts
if venn.get_label_by_id('10'):
    venn.get_label_by_id('10').set_text(f"{len(rsv_gene_set - flu_gene_set)}")
    plt.text(
        -0.7, 0.2,
        "\n".join(rsv_only_genes),
        fontsize=10, ha="center", va="center",
        bbox=dict(facecolor='blue', alpha=0.3, boxstyle="round,pad=0.4"),
    )

if venn.get_label_by_id('01'):
    venn.get_label_by_id('01').set_text(f"{len(flu_gene_set - rsv_gene_set)}")
    plt.text(
        0.7, 0.2,  # Adjust position
        "\n".join(flu_only_genes),
        fontsize=10, ha="center", va="center",
        bbox=dict(facecolor='orange', alpha=0.3, boxstyle="round,pad=0.4"),
    )

if venn.get_label_by_id('11'):
    venn.get_label_by_id('11').set_text(f"{len(rsv_gene_set & flu_gene_set)}")
    plt.text(
        0, -0.6,  # Adjust position
        "\n".join(overlapping_genes),
        fontsize=10, ha="center", va="center",
        bbox=dict(facecolor='grey', alpha=0.3, boxstyle="round,pad=0.4"),
    )

# Customize venn diagram apperance
plt.title("Venn Diagram of Top 100 Significant Genes (with Gene Symbol or Accession Numbers)", fontsize=16)
plt.tight_layout()

# Save the Venn diagram
plt.savefig("venn_diagram_overlapping.png")
plt.show()

print("Venn diagram saved as 'venn_diagram_overlapping.png'.")

#GO TERM ANALYSIS: compare GO terms between RSV and flu

#1. Comparison of cellular component GO terms between top 100 RSV and Flu genes

# Function to extract the middle part of the GO term
def extract_middle_go_term(go_term):
    if isinstance(go_term, str): #isinstance checks that the go_term is a python str
        parts = go_term.split(" // ")  # Split by //
        if len(parts) > 1:  # Ensure there is a middle part
            return parts[1]  # Extract the middle part after first //
    return None  # Return None if term is invalid

# Count GO CC terms for RSV and Flu
rsv_go_cc = rsv_top100["Gene Ontology Cellular Component"].apply(extract_middle_go_term) #using function to extract just GO term
flu_go_cc = flu_top100["Gene Ontology Cellular Component"].apply(extract_middle_go_term)

# Count occurrences of each GO term for RSV and Flu
rsv_go_cc_count = rsv_go_cc.value_counts()
flu_go_cc_count = flu_go_cc.value_counts()
#.value_counts counts occurances of genes within the panda series

# Combine the counts into a DataFrame for comparison
go_cc_comparison = pd.DataFrame({
    "RSV": rsv_go_cc_count,
    "Flu": flu_go_cc_count
}).fillna(0)  # Replace NaN with 0 for terms not present in one dataset

# Sort by total occurrences across both conditions
go_cc_comparison["Total"] = go_cc_comparison["RSV"] + go_cc_comparison["Flu"]
go_cc_comparison = go_cc_comparison.sort_values("Total", ascending=False) #sorted values in decending order

# Plot comparison
go_cc_comparison[["RSV", "Flu"]].plot(kind="bar", figsize=(12, 8), width=0.8)
plt.title("Comparison of Gene Ontology Cellular Component Terms")
plt.xlabel("GO Cellular Component")
plt.ylabel("Frequency")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Condition")
plt.tight_layout()
plt.show() #decided not to include this figure in report as I found the top 10 cellular components clearer for analysis

#2. Bar chart comparison of flu and RSV: top 10 cellular components

# Get the top 10 most frequent GO terms for Flu and RSV
top_flu_go_cc = flu_go_cc_count.head(10) #extract just top 10 using .head(10)
top_rsv_go_cc = rsv_go_cc_count.head(10)

# Combine the top 10 counts of Flu and RSV into one DataFrame for comparison
go_comparison_df = pd.DataFrame({
    "RSV": top_rsv_go_cc,
    "Flu": top_flu_go_cc
}).fillna(0)  # Replace NaN with 0 for terms not present in one dataset

# Sort the GO terms by total frequency
go_comparison_df["Total"] = go_comparison_df["RSV"] + go_comparison_df["Flu"]
go_comparison_df = go_comparison_df.sort_values("Total", ascending=False) #sorted values in decending order

# Plot the comparison of top 10 GO terms in Cellular Components
plt.figure(figsize=(12, 8))
go_comparison_df[["RSV", "Flu"]].plot(kind="bar", figsize=(12, 8), width=0.8)
plt.title("Comparison of Top 10 Gene Ontology Cellular Component Terms for RSV and Flu")
plt.xlabel("GO Cellular Component Terms")
plt.ylabel("Frequency")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Condition")
plt.tight_layout()
plt.show()

#3. Bar chart comparison of flu and RSV: Top 10 Biological processes
# Extract the middle GO term for Biological Process
rsv_go_bp = rsv_top100["Gene Ontology Biological Process"].apply(extract_middle_go_term) #using function to extract just GO term
flu_go_bp = flu_top100["Gene Ontology Biological Process"].apply(extract_middle_go_term)

# Count occurrences of each GO term for RSV and Flu
rsv_go_bp_count = rsv_go_bp.value_counts()
flu_go_bp_count = flu_go_bp.value_counts()
#.value_counts counts occurances of genes within the panda series

# Get the top 10 most frequent GO terms for Flu and RSV
top_flu_go_bp = flu_go_bp_count.head(10) #extract just top 10 using .head(10)
top_rsv_go_bp = rsv_go_bp_count.head(10)

# Combine the top 10 counts of Flu and RSV into one panda DataFrame for comparison
go_comparison_bp_df = pd.DataFrame({
    "RSV": top_rsv_go_bp,
    "Flu": top_flu_go_bp
}).fillna(0)  # Replace NaN with 0 for terms not present in dataset

# Sort the GO terms by total frequency
go_comparison_bp_df["Total"] = go_comparison_bp_df["RSV"] + go_comparison_bp_df["Flu"]
go_comparison_bp_df = go_comparison_bp_df.sort_values("Total", ascending=False) #sorted values in decending order

# Plot the comparison of top 10 GO terms in Biological Process
plt.figure(figsize=(12, 8))
go_comparison_bp_df[["RSV", "Flu"]].plot(kind="bar", figsize=(12, 8), width=0.8)
plt.title("Comparison of Top 10 Gene Ontology Biological Process Terms for RSV and Flu")
plt.xlabel("GO Biological Process Terms")
plt.ylabel("Frequency")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Condition")
plt.tight_layout()
plt.show()

#3. Bar chart comparison of flu and rsv: Top 10 molecular functions
# Extract the middle GO term for Molecular Function
rsv_go_mf = rsv_top100["Gene Ontology Molecular Function"].apply(extract_middle_go_term) #using function to extract just GO term
flu_go_mf = flu_top100["Gene Ontology Molecular Function"].apply(extract_middle_go_term)

# Count occurrences of each GO term for RSV and Flu
rsv_go_mf_count = rsv_go_mf.value_counts()
flu_go_mf_count = flu_go_mf.value_counts()
#.value_counts counts occurances of genes within the panda series

# Get the top 10 most frequent GO terms for Flu and RSV
top_flu_go_mf = flu_go_mf_count.head(10)  #extract just top 10 using .head(10)
top_rsv_go_mf = rsv_go_mf_count.head(10)

# Combine the top 10 counts of Flu and RSV into one DataFrame for comparison
go_comparison_mf_df = pd.DataFrame({
    "RSV": top_rsv_go_mf,
    "Flu": top_flu_go_mf
}).fillna(0)  # Replace NaN with 0 for terms not present in dataset

# Sort the GO terms by total frequency
go_comparison_mf_df["Total"] = go_comparison_mf_df["RSV"] + go_comparison_mf_df["Flu"]
go_comparison_mf_df = go_comparison_mf_df.sort_values("Total", ascending=False) #sorted values in decending order

# Plot the comparison of top 10 GO terms in Molecular Functions
plt.figure(figsize=(12, 8))
go_comparison_mf_df[["RSV", "Flu"]].plot(kind="bar", figsize=(12, 8), width=0.8)
plt.title("Comparison of Top 10 Gene Ontology Molecular Function Terms for RSV and Flu")
plt.xlabel("GO Molecular Function Terms")
plt.ylabel("Frequency")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Condition")
plt.tight_layout()
plt.show()


#Summary Table Generation: OVERLAPPING genes

# Filter show summary table of 8 overlapping genes from venn diagram

# genes on interest from the venn diagrams
overlapping_genes_interest = ['CLEC4C', 'COX20', 'GLG1', 'H1F0', 'IFI27', 'LDLR', 'SIGLEC1']
filtered_df = merged_data[merged_data['Gene Symbol'].isin(overlapping_genes_interest)].copy()
filtered_df = filtered_df.drop_duplicates(subset=['Gene Symbol'], keep='first')

GO_to_clean = ['Gene Ontology Biological Process', 'Gene Ontology Cellular Component', 'Gene Ontology Molecular Function']
# use for for loop to apply extract_middle_go_term function to columns

for col in GO_to_clean:
    filtered_df[col] = filtered_df[col].apply(extract_middle_go_term)

# Create a table with GO term information headings
selected_columns = [
    "Gene Symbol",
    "Gene Ontology Biological Process",
    "Gene Ontology Cellular Component",
    "Gene Ontology Molecular Function",
    "rsv_pvalue",
    "flu_pvalue"
]
# use these headings to extract relevant information
genes_of_interest = filtered_df[selected_columns]

# Save as csv then convert to html and then write + save as html
genes_of_interest.to_csv("overlapping_genes_summary.html.csv", index=False)
html_table = genes_of_interest.to_html(index=False)
with open("overlapping_genes_summary.html", "w") as file:
    file.write(html_table)


#Summary Table Generation: Top 10 Flu genes

# 1. Filter Flu data to exclude rows with NaN (using .dropna) in GO terms
flu_with_go_terms = flu_top100.dropna(subset=[
    "Gene Ontology Biological Process",
    "Gene Ontology Cellular Component",
    "Gene Ontology Molecular Function"
])

# 2. Select the top 10 significant Flu genes with GO terms
top_10_flu_genes = flu_with_go_terms.nsmallest(10, "flu_pvalue") #show top 10 most sig genes

# 3. Display the table with relevant information including gene symbol, GO terms and p values
table_columns = [
    "Gene Symbol",
    "Gene Ontology Biological Process",
    "Gene Ontology Cellular Component",
    "Gene Ontology Molecular Function",
    "flu_pvalue"
]

# 4. Extract relevant information
top_10_flu_info = top_10_flu_genes[table_columns]


# Apply the function to the Flu gene GO term columns to extract only middle using function created above
top_10_flu_info["Gene Ontology Biological Process"] = top_10_flu_info["Gene Ontology Biological Process"].apply(extract_middle_go_term)
top_10_flu_info["Gene Ontology Cellular Component"] = top_10_flu_info["Gene Ontology Cellular Component"].apply(extract_middle_go_term)
top_10_flu_info["Gene Ontology Molecular Function"] = top_10_flu_info["Gene Ontology Molecular Function"].apply(extract_middle_go_term)

# 5. Save table to a CSV file
top_10_flu_info.to_csv("top_10_flu_genes_info_cleaned.csv", index=False)

# 6. Convert table to HTML and save it as an HTML file
html_table_cleaned = top_10_flu_info.to_html(index=False)

with open("top_10_flu_genes_info_cleaned.html", "w") as f:
    f.write(html_table_cleaned)

# 7. Print the HTML table for review
print(html_table_cleaned)

#Summary Table Generation: Top 10 RSV genes

# 1. Filter Flu data to exclude rows with NaN in GO terms
rsv_with_go_terms = rsv_top100.dropna(subset=[
    "Gene Ontology Biological Process",
    "Gene Ontology Cellular Component",
    "Gene Ontology Molecular Function"
])

# 2. Select the top 10 significant (lowest pvalue) Flu genes with GO terms
top_10_rsv_genes = rsv_with_go_terms.nsmallest(10, "rsv_pvalue") #show top 10 most sig genes

# 3. Display the table with relevant information including gene symbol, GO terms and p values
table_columns = [
    "Gene Symbol",
    "Gene Ontology Biological Process",
    "Gene Ontology Cellular Component",
    "Gene Ontology Molecular Function",
    "rsv_pvalue"
]

# 4. Extract relevant information
top_10_rsv_info = top_10_rsv_genes[table_columns]

# Apply the function to the RSV gene GO term columns to extract only middle using function created above
top_10_rsv_info["Gene Ontology Biological Process"] = top_10_rsv_info["Gene Ontology Biological Process"].apply(extract_middle_go_term)
top_10_rsv_info["Gene Ontology Cellular Component"] = top_10_rsv_info["Gene Ontology Cellular Component"].apply(extract_middle_go_term)
top_10_rsv_info["Gene Ontology Molecular Function"] = top_10_rsv_info["Gene Ontology Molecular Function"].apply(extract_middle_go_term)

# 5. Save the table to CSV file
top_10_rsv_info.to_csv("top_10_rsv_genes_info_cleaned.csv", index=False)

# 6. Convert table to HTML and save it as an HTML file
html_table_rsv_cleaned = top_10_rsv_info.to_html(index=False)

with open("top_10_rsv_genes_info_cleaned.html", "w") as f:
    f.write(html_table_rsv_cleaned)

# 7. Print the HTML table string for review
print(html_table_rsv_cleaned)


#Enrichment GO terms analysis


# Function to calculate enrichment factors for each GO category
def calculate_enrichment_factors(go_terms, all_go_terms, de_genes):
    """
    Calculate Enrichment Factor (EF) for each GO term.
    go_terms: GO terms from the differential expressed genes (DE).
    all_go_terms: GO terms from the background gene set.
    de_genes: Differentially expressed genes.
    """
    background_go_count = all_go_terms.value_counts()  # Background gene count for each GO term
    de_go_count = go_terms.value_counts()  # DE gene count for each GO term

    b = len(de_genes)  # Total DE genes
    B = len(all_go_terms)  # Total genes in the background

    enrichment_factors = {}

    for go_term in de_go_count.index:
        a = de_go_count[go_term]  # Count of the differentially expressed genes mapped to this GO term
        A = background_go_count.get(go_term, 0)  # Count of background genes mapped to this GO term

        if A > 0:  # Avoid dividing by zero
            EF = (a / b) / (A / B)  # Enrichment Factor formula from report
            enrichment_factors[go_term] = EF

    return enrichment_factors


# 1. Calculate Enrichment Factors values for each GO category (Biological Process(BP), Cellular Component (CC), Molecular Function (MF)

# RSV Genes Enrichment Analysis

# Cellular Component (GO CC) - RSV
rsv_go_cc = rsv_top100["Gene Ontology Cellular Component"].apply(extract_middle_go_term)
all_genes_go_cc = merged_data["Gene Ontology Cellular Component"].dropna().apply(extract_middle_go_term)

rsv_go_cc_enrichment = calculate_enrichment_factors(rsv_go_cc, all_genes_go_cc, rsv_go_cc)

# Biological Process (GO BP) - RSV
rsv_go_bp = rsv_top100["Gene Ontology Biological Process"].apply(extract_middle_go_term)
all_genes_go_bp = merged_data["Gene Ontology Biological Process"].dropna().apply(extract_middle_go_term)

rsv_go_bp_enrichment = calculate_enrichment_factors(rsv_go_bp, all_genes_go_bp, rsv_go_bp)

# For Molecular Function (GO MF) - RSV
rsv_go_mf = rsv_top100["Gene Ontology Molecular Function"].apply(extract_middle_go_term)
all_genes_go_mf = merged_data["Gene Ontology Molecular Function"].dropna().apply(extract_middle_go_term)

rsv_go_mf_enrichment = calculate_enrichment_factors(rsv_go_mf, all_genes_go_mf, rsv_go_mf)

#GO terms extracted from top 100 genes for RSV. All gene go terms extracted dropping any NA values
# before feeding this into enrichment function to get enrichment factor value, repeated for CC, BP AND MF

# Flu Genes Enrichment Analysis

# Cellular Component (GO CC) - Flu
flu_go_cc = flu_top100["Gene Ontology Cellular Component"].apply(extract_middle_go_term)
all_genes_go_cc = merged_data["Gene Ontology Cellular Component"].dropna().apply(extract_middle_go_term)

flu_go_cc_enrichment = calculate_enrichment_factors(flu_go_cc, all_genes_go_cc, flu_go_cc)

# Biological Process (GO BP) - Flu
flu_go_bp = flu_top100["Gene Ontology Biological Process"].apply(extract_middle_go_term)
all_genes_go_bp = merged_data["Gene Ontology Biological Process"].dropna().apply(extract_middle_go_term)

flu_go_bp_enrichment = calculate_enrichment_factors(flu_go_bp, all_genes_go_bp, flu_go_bp)

# Molecular Function (GO MF) - Flu
flu_go_mf = flu_top100["Gene Ontology Molecular Function"].apply(extract_middle_go_term)
all_genes_go_mf = merged_data["Gene Ontology Molecular Function"].dropna().apply(extract_middle_go_term)

flu_go_mf_enrichment = calculate_enrichment_factors(flu_go_mf, all_genes_go_mf, flu_go_mf)

#GO terms extracted from top 100 genes for flu. All gene go terms extracted dropping any NA values
# before feeding this into enrichment function to get enrichment factor value, repeated for CC, BP AND MF


# For RSV
top_rsv_cc = max(rsv_go_cc_enrichment.items(), key=lambda x: x[1])
top_rsv_bp = max(rsv_go_bp_enrichment.items(), key=lambda x: x[1])
top_rsv_mf = max(rsv_go_mf_enrichment.items(), key=lambda x: x[1])

# For Flu
top_flu_cc = max(flu_go_cc_enrichment.items(), key=lambda x: x[1])
top_flu_bp = max(flu_go_bp_enrichment.items(), key=lambda x: x[1])
top_flu_mf = max(flu_go_mf_enrichment.items(), key=lambda x: x[1])

#rsv/flu_go_cc_enrichment.items dictionaires have GO terms as keys and enrichment scores as values
#max(..., key=lambda x: x[1]) finds the GO terms that have the highest enrichment scores

# DataFrames created for RSV and Flu Results

# DataFrame for RSV showing the top enriched GO terms and their EF for each category
top_enriched_rsv_df = pd.DataFrame({
    "GO Term (Cellular Component) - RSV": [top_rsv_cc[0]],
    "Enrichment Factor (CC) - RSV": [top_rsv_cc[1]],
    "GO Term (Biological Process) - RSV": [top_rsv_bp[0]],
    "Enrichment Factor (BP) - RSV": [top_rsv_bp[1]],
    "GO Term (Molecular Function) - RSV": [top_rsv_mf[0]],
    "Enrichment Factor (MF) - RSV": [top_rsv_mf[1]],
})

# DataFrame for Flu showing the top enriched GO terms and their EF for each category
top_enriched_flu_df = pd.DataFrame({
    "GO Term (Cellular Component) - Flu": [top_flu_cc[0]],
    "Enrichment Factor (CC) - Flu": [top_flu_cc[1]],
    "GO Term (Biological Process) - Flu": [top_flu_bp[0]],
    "Enrichment Factor (BP) - Flu": [top_flu_bp[1]],
    "GO Term (Molecular Function) - Flu": [top_flu_mf[0]],
    "Enrichment Factor (MF) - Flu": [top_flu_mf[1]],
})

# DataFrames converted to HTML Format

# Convert RSV DataFrame to HTML
html_table_rsv = top_enriched_rsv_df.to_html(index=False)

# Convert Flu DataFrame to HTML
html_table_flu = top_enriched_flu_df.to_html(index=False)

# Save the RSV HTML table to a file
with open("top_enriched_go_terms_rsv_final.html", "w") as f:
    f.write(html_table_rsv)

# Save the Flu HTML table to a file
with open("top_enriched_go_terms_flu_final.html", "w") as f:
    f.write(html_table_flu)

#print the HTML tables for review
print("RSV HTML Table:")
print(html_table_rsv)

print("Flu HTML Table:")
print(html_table_flu)


