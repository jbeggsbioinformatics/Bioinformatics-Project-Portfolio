# Declaration:

# I, John Alexander Beggs, declare that I have employed a Chat-GPT-3.5, to assist in the creation

# of this .py script. Specifically, I used it to proactively learn how to carry out a four component PCA plot grid for analysis of infection status.

#Part 2: PCA plots

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load in features.tsv file
features = pd.read_csv("features.tsv", sep="\t", index_col=0)

# 2: Load in metadata file to merge
metadata = pd.read_csv("meta_data.csv")

# Ensure metadata matches the sample columns in features file
sample_metadata = metadata.set_index('Sample_geo_accession')

#3. Carry out principal components analysis on the data, using all samples and significant features
# Remove any non-numerical columns
numerical_features = features.select_dtypes(include=[np.number])

# Produce principal components analysis
pca = PCA(n_components=2) # 2 components PCA
principal_components = pca.fit_transform(numerical_features.T)  # Transpose to have rows as samples

# Creation of dataFrame for the principal components analysis results
pca_dataframe = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=numerical_features.columns)

# Merge PCA dataframe with metadata to add infection status, age group and sex
pca_dataframe = pca_dataframe.merge(sample_metadata, left_index=True, right_index=True)

# 4. Create principal components analysis plots

# By infection status:

plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_dataframe, x='PC1', y='PC2', hue='infection_status', palette='Set2')
plt.title('PCA Plot by Infection Status')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Infection Status')
plt.savefig('pca_by_infection_status.png')
plt.close()

#By sex:

palette = {'F': 'red', 'M': 'blue'}

plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_dataframe, x='PC1', y='PC2', hue='gender', palette=palette)
plt.title('PCA Plot by Sex')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Sex')
plt.savefig('pca_by_sex.png')
plt.close()

#By age group: (>6 months vs <=6 months)

pca_dataframe['age_group'] = np.where(pca_dataframe['age_months'] > 6, '>6 months', '<=6 months')
plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_dataframe, x='PC1', y='PC2', hue='age_group', palette='bright6')
plt.title('PCA Plot by Age Group')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Age Group')
plt.savefig('pca_by_age_group.png')
plt.close()

#5. 4 components principal components analysis
pca = PCA(n_components=4)
principal_components = pca.fit_transform(numerical_features.T)  # Transpose to rows as samples

# Create a dataframe for the principal components analysis results
pca_dataframe = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2', 'PC3', 'PC4'], index=numerical_features.columns)

# Merge with metadata to add infection status
pca_dataframe = pca_dataframe.merge(sample_metadata, left_index=True, right_index=True)

#Creation of PCA plots grid colored by infection status
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Pairs of PCs to plot
pc_pairs = [("PC1", "PC2"), ("PC1", "PC3"), ("PC2", "PC3"), ("PC2", "PC4")]

# Generation of scatter plots for each pair
# Coloured by infection status
for ax, (x_pc, y_pc) in zip(axes.flatten(), pc_pairs): #axes.flatten, flattens grid of subplot axes making it a 1D arrary to iterate over
    sns.scatterplot(data=pca_dataframe, x=x_pc, y=y_pc, hue='infection_status', palette='Set1', ax=ax, alpha=0.8)
    ax.set_title(f"{x_pc} vs {y_pc}")
    ax.set_xlabel(x_pc)
    ax.set_ylabel(y_pc)
    ax.legend(title='Infection Status', loc='upper right',  fontsize=7.5, markerscale=0.7, title_fontsize=8.5)

    ax.set_facecolor((0.92, 0.92, 0.92))  # Light gray background
    ax.grid(color="white", linestyle="-", linewidth=0.5, zorder=0)  # ensure lines of grid are behind points

# Adjustment of layout
plt.tight_layout()
plt.savefig('4_component_pca_4_grid.png')
plt.close()

print("4 components PCA grid plot colored by infection status saved as '4_component_pca_4_grid.png'.")


