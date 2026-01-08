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