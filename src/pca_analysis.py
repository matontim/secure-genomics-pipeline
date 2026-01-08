#Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path

# File paths
DATA_PATH = Path("data/raw/GSE260586_FPKM+1.txt")
fig_dir = Path("results/figures")
fig_dir.mkdir(parents=True, exist_ok=True)

# Load dataset
df = pd.read_csv(DATA_PATH, sep="\t")

# Set GeneID as index
gene_column = df.columns[7]
df = df.set_index(gene_column)

#Drop metadata columns
metadata_rows = ['Chromosome', 'Start', 'Stop', 'Length', 'Strand', 'Gene Symbol', 'gene_biotype', 'gene_name', 'gene_source', 'gene_version']
df = df.drop(columns=metadata_rows)

# QC and log2 transform
log2expr = np.log2(df + 1)

# Transpose data
log2expr_T = log2expr.T   

# Scale data
scaler = StandardScaler()
scaled_T = scaler.fit_transform(log2expr_T)

# PCA
pca = PCA(n_components=8)
pca_T = pca.fit_transform(scaled_T)

# Explained variance
explained_var = pca.explained_variance_ratio_
print("Explained variance ratios:", pca.explained_variance_ratio_)

# Create a Scree plot for explained variance
plt.figure(figsize=(6,4))
plt.plot(range(1, len(explained_var) + 1), explained_var * 100, marker='o')
plt.xlabel('Principal Component')
plt.ylabel('Explained Variance (%)')
plt.title('Scree Plot of PCA on Expression Data')
plt.grid(True)
plt.tight_layout()
# Elbow at PC2 suggests using first two PCs for visualization
plt.savefig(fig_dir/"scree_plot.png")
plt.close()

# Add cumulative variance
cum_var = np.cumsum(explained_var) * 100
plt.plot(range(1, len(cum_var) + 1), cum_var, marker='s', color='orange')
plt.axhline(80, linestyle='--', color='red', label='80% Variance')
plt.xlabel('Principal Component')
plt.ylabel('Cumulative Explained Variance (%)')
plt.title('Cumulative Explained Variance by PCA Components')
plt.tight_layout()
plt.savefig(fig_dir/"cumulative_variance.png")
plt.close()
#The first three componenets make up 80% of the variance

# PCA with 2 components
pca2 = PCA(n_components=2)
pcs2 = pca2.fit_transform(scaled_T)

# Create new DataFrame with PCs
pca_df = pd.DataFrame(data=pcs2, columns=['PC1', 'PC2'], index=log2expr.columns)

# Variance labels
pc1_var = pca2.explained_variance_ratio_[0] * 100
pc2_var = pca2.explained_variance_ratio_[1] * 100   

# Plot 
plt.figure(figsize=(6,5))
plt.scatter(pca_df['PC1'], pca_df['PC2'], s=50)
plt.xlabel(f'PC1 ({pc1_var:.1f}% variance)')
plt.ylabel(f'PC2 ({pc2_var:.1f}% variance)')
plt.title('PCA of GSE260586 Expression Data')
plt.tight_layout()
plt.savefig(fig_dir/"pca_scatter.png")
plt.show()