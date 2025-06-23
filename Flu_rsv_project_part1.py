# Declaration:

# I, John Alexander Beggs, declare that I have employed a Chat-GPT-3.5, to assist in the creation

# of this .py script. Specifically, I used it to proactively learn how to make my t test and volcano plot creation embeded in functions.

#This made the code more concise and flow more smoothly. I was also taught how to create html tables for my report and also make my bar

#charts more visually appealing. eg gene count inside bar. I also had to ask it how to calculate FDR as i had issues with scipy, so i used statsmodels.stats.multitest instead.

import pandas as pd
import numpy as np


#Part 1: Differential Gene expression analysis + volcano plot

# 1. Load metadata
metadata = pd.read_csv("meta_data.csv")
#check headers
print(metadata.head())
#check infection status column
print(metadata['infection_status'].unique())
#rsv, none, influenza

# Extract sample IDs into lists based on infection_status
flu_samples = metadata.loc[metadata['infection_status'] == 'influenza']['Sample_geo_accession'].tolist()
rsv_samples = metadata.loc[metadata['infection_status'] == 'rsv']['Sample_geo_accession'].tolist()
control_samples = metadata.loc[metadata['infection_status'] == 'none']['Sample_geo_accession'].tolist()
#.loc[] filters to select the specific rows from the metadata dataframe, based on infection status
#.tolist() converts values in the sample_geo_accession column from the filtered rows into a list

print("Flu Samples:", flu_samples)
print("RSV Samples:", rsv_samples)
print("Control Samples:", control_samples)

# 2.Read in matrix of gene expression values
#num_rows = 1000  # Testing using first 1000 rows, total rows: 54675
df_matrix = pd.read_csv("GSE34205_series_matrix_clean.txt", sep="\t", index_col=0)

print(df_matrix.head())

# Reorder columns: flu samples, RSV samples, control samples
reordered_columns = flu_samples + rsv_samples + control_samples #combines the three lists into single list(reordered_columns)
df_matrix = df_matrix[reordered_columns] #columns of df_matrix are rearranged according to the order of reordered_columns list

print("\nReordered Gene Expression Matrix:")
print(df_matrix.head())  # Confirm the column order

# Convert all values to their natural log using np.log function
df_matrix = np.log(df_matrix) #numpy package has log conversion function built in

# Preview the transformed DataFrame
print("\nLog-transformed Gene Expression Matrix:")
print(df_matrix.head())

#3. Function to calculate a one-sample t-test (testing against zero)
def one_group_ttest(row, group_start, group_end):
    group = row.iloc[group_start:group_end]
    t_stat, p_value = ttest_1samp(group, popmean=0, nan_policy='omit') #popmean=0, checks if samples mean is sig different from 0
    return p_value

#row.iloc extracts a slice of the row from group start to group end saved as group
#group represents the data values for the t test

# Calculate p-values for flu and RSV samples using t test
df_matrix['flu_pvalue'] = df_matrix.apply(lambda row: one_group_ttest(row, 0, len(flu_samples)), axis=1)
df_matrix['rsv_pvalue'] = df_matrix.apply(lambda row: one_group_ttest(row, len(flu_samples), len(flu_samples) + len(rsv_samples)), axis=1)
# the rsv samples start after the flu in the list, so len(flu) samples is set as the group start (instead of 0) and len(flu_samples) + len(rsv_samples))
#is set as the group_end to capture the full list of rsv samples for the t test function to apply to.

#4. Apply false discovery rate correction to p values
df_matrix['flu_fdr'] = multipletests(df_matrix['flu_pvalue'], method='fdr_bh')[1]
df_matrix['rsv_fdr'] = multipletests(df_matrix['rsv_pvalue'], method='fdr_bh')[1]
#use of statsmodels multipletests function, using Benjamini-Hochberg procedure to FDR correct p values
#I was unable to use scipy.false_discovery_control due to an issue with updating scipy due to M1 macbook compatabaility issues (spoke to demonstrator about this)

# Calculate the mean of logged values for flu and RSV samples
df_matrix['flu_mean_log_ratio'] = df_matrix.iloc[:, 0:len(flu_samples)].mean(axis=1)
df_matrix['rsv_mean_log_ratio'] = df_matrix.iloc[:,len(flu_samples):len(flu_samples) + len(rsv_samples)].mean(axis=1)

# Preview the DataFrame with flu and RSV values
print("\nDataFrame with flu and rsv p-values, FDR values, and mean log ratios:")
print(df_matrix[['rsv_pvalue', 'rsv_fdr', 'rsv_mean_log_ratio']].head())

#5. Create a new column 'significant_feature' based on conditions for 'flu_fdr' and 'rsv_fdr'
df_matrix['significant_feature'] = (df_matrix['flu_fdr'] < 0.01) & (df_matrix['rsv_fdr'] < 0.01) #only flu and rsv samples with FDR <0.01

# Preview updated DataFrame with 'significant_feature' column
print(df_matrix[['flu_fdr', 'rsv_fdr', 'significant_feature']].head())

#6. Filter the DataFrame to include only rows with significant_feature == True
significant_df = df_matrix[df_matrix['significant_feature'] == True]

#7. Write the filtered DataFrame to a TSV file
significant_df.to_csv('features.tsv', sep='\t', index=True)

# Confirm that the file was written
print("Filtered data has been written to 'Features.tsv'.")

#8.Merge with probe data and create volcano plots with labels

# Load microarray probe data file
try:
    probe_df = pd.read_csv(
        "GPL570-55999.txt",
        sep="\t",  # Tab-separated
        index_col=0,  # first column as the index
        comment="#",  # Skip lines starting with #
        engine="python"  # Use Python engine
    )
    print(f"Probe data loaded successfully with shape: {probe_df.shape}")
except Exception as e:
    print(f"Error loading probe data: {e}") #troubleshooting method
    raise


# Merge the data matrix with probe data file
df_matrix = df_matrix.merge(probe_df, left_index=True, right_index=True, how="left")


def create_labeled_volcano_plot(df_matrix, mean_col, fdr_col, group_name, pval_col, gene_symbol_col, y_limit=None):
    sns.set_style("darkgrid")  # Set grey background in the plot area
    plt.figure(figsize=(10, 8))
    df_matrix['-log10_fdr'] = -np.log10(df_matrix[fdr_col]) #numoy used to get the 10 log of FDR values, - sign gives the negative
    #log10fdr now contains the negative -10 log of FDR values

    # Create the volcano plot
    sns.scatterplot(
        x=df_matrix[mean_col],
        y=df_matrix['-log10_fdr'],
        hue=df_matrix[fdr_col] < 0.01,  # Significant vs non significant
        palette={True: 'red', False: 'green'}, #set dot colours
        legend='full',  # Include legend
        s=6,  # smaller dots
        edgecolor=None  # Remove white outline
    )

    # Label the 25 most significant points
    top25 = df_matrix.nsmallest(25, pval_col) #retrieves 25 lowest p value points
    for i, row in top25.iterrows():
        plt.text(
            row[mean_col],
            row['-log10_fdr'],
            str(row[gene_symbol_col]),
            fontsize=3,  # Smaller font for labels for visability
            color='black',
            ha='right'
        )

    # Add titles, labels, and grid to volcano plot
    plt.title(f"Volcano Plot for {group_name}", fontsize=14)
    plt.xlabel("Mean Log Ratio", fontsize=12)
    plt.ylabel("-log10(FDR)", fontsize=12)
    plt.axhline(y=2, color='grey', linestyle='--', linewidth=0.8)  # Threshold line
    plt.axvline(x=0, color='grey', linestyle='--', linewidth=0.8)  # Zero mean log ratio line, left under expressed, right overexpressed
    plt.xlim(-5, 6)  # Set x-axis range to -5 to 5

    # Set y-axis limit based on the argument y_limit from function, ensures both plots are to scale
    if y_limit:
        plt.ylim(-0.5, y_limit)  # Use the provided y_limit: given below
    else:
        plt.ylim(-0.5, 25)  # Default y-axis range

    # Customize the legend
    handles, labels = plt.gca().get_legend_handles_labels()
    labels = ['Not Significant', 'Significant']  # Modification of legend labels
    plt.legend(handles, labels, title="Significance", fontsize=10, title_fontsize=12)

    # Save the plot as pdf file for zoom
    plt.savefig(f"volcano_plot_{group_name}.pdf", dpi=300)  # High-quality output
    plt.show()  # Display the plot
    plt.close()

# Create volcano plots using function for influenza and RSV with different y-axis limits
create_labeled_volcano_plot(df_matrix, 'flu_mean_log_ratio', 'flu_fdr', 'Influenza', 'flu_pvalue', 'Gene Symbol', y_limit=15)
create_labeled_volcano_plot(df_matrix, 'rsv_mean_log_ratio', 'rsv_fdr', 'RSV', 'rsv_pvalue', 'Gene Symbol', y_limit=25)
#given y limit to make volcano plots scaled

print("Analysis complete. Results saved to 'features.tsv'. Labeled volcano plots saved.")

#9. Create summary tables for flu and rsv by fdr

# Calculate the number of significant genes for each condition
flu_significant_genes = df_matrix[df_matrix['flu_fdr'] < 0.01].shape[0]
rsv_significant_genes = df_matrix[df_matrix['rsv_fdr'] < 0.01].shape[0]

# Create a summary DataFrame
summary_table = pd.DataFrame({
    'Condition': ['Influenza', 'RSV'],
    'Significant Genes': [flu_significant_genes, rsv_significant_genes]
})

# Print the summary table
print("\nSummary Table of Significant Genes:")
print(summary_table)

# Save the summary table to a CSV file
summary_table.to_csv("summary_table_significant_genes.csv", index=False)

#10. Create overexpressed / under expressed table for flu and RSV

# Categorise significant genes into overexpressed and underexpressed for Influenza
flu_overexpressed = df_matrix[(df_matrix['flu_fdr'] < 0.01) & (df_matrix['flu_mean_log_ratio'] > 0)].shape[0] #significant overexpressed
flu_underexpressed = df_matrix[(df_matrix['flu_fdr'] < 0.01) & (df_matrix['flu_mean_log_ratio'] < 0)].shape[0] #significant underexpressed

# Categorise significant genes into overexpressed and underexpressed for RSV
rsv_overexpressed = df_matrix[(df_matrix['rsv_fdr'] < 0.01) & (df_matrix['rsv_mean_log_ratio'] > 0)].shape[0] #significant overexpressed
rsv_underexpressed = df_matrix[(df_matrix['rsv_fdr'] < 0.01) & (df_matrix['rsv_mean_log_ratio'] < 0)].shape[0] #significant underexpressed

# Create summary DataFrame
over_underexpressed_summary = pd.DataFrame({
    'Condition': ['Influenza', 'Influenza', 'RSV', 'RSV'],
    'Category': ['Overexpressed', 'Underexpressed', 'Overexpressed', 'Underexpressed'],
    'Count': [flu_overexpressed, flu_underexpressed, rsv_overexpressed, rsv_underexpressed]
})

# Print the summary table
print("\nSummary of Overexpressed and Underexpressed Genes:")
print(over_underexpressed_summary)


#11. Create bar plot with separate bars for overexpressed and underexpressed genes
plt.figure(figsize=(10, 6))
bar_plot = sns.barplot(
    data=over_underexpressed_summary,
    x='Condition',
    y='Count',
    hue='Category',
    palette='Set2'
)

# Annotate each bar with the count value inside the bar
for p in bar_plot.patches:
    height = p.get_height()
    if height > 0:  # Only annotate bars with positive height
        bar_plot.text(
            x=p.get_x() + p.get_width() / 2,  # Center of the bar
            y=height / 2,  # Vertically centered
            s=f"{int(height)}",  # Annotation text
            ha="center",
            va="center",
            color="black",
            fontsize=12
        )

# Add titles and labels
plt.title('Significantly Overexpressed and Underexpressed Genes for RSV and Influenza', fontsize=14)
plt.ylabel('Number of Genes (fdr < 0.01)', fontsize=12)
plt.xlabel('Infection Status', fontsize=12)
plt.legend(title='Category', fontsize=10, title_fontsize=12)

# Save the plot
plt.savefig("over_under_expressed_genes_barplot.png", dpi=300)
plt.show()


