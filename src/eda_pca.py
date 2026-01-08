#Importing necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#Load dataset
DATA_PATH = "data/raw/GSE260586_FPKM+1.txt"

df = pd.read_csv(DATA_PATH, sep="\t")

print("Dataset shape:", df.shape)
print(df.head())

#Column 7 is GeneID
gene_column = df.columns[7]

#Set GeneID as index
df = df.set_index(gene_column)

print("After setting GeneID as index:")
print(df.head())

#Drop metadata rows
metadata_rows = ['Chromosome', 'Start', 'Stop', 'Length', 'Strand', 'Gene Symbol', 'gene_biotype', 'gene_name', 'gene_source', 'gene_version']
df = df.drop(columns=metadata_rows)

print(df.head())
print("Dataset shape after dropping metadata columns:", df.shape)

# Raw QC to catch parsing errors or missing values

# Ensure numeric
df = df.apply(pd.to_numeric, errors='coerce')

# Check for missing values 
print("Missing values:", df.isna().sum().sum()) 

# Check basic stats
print(df.describe())

# FPKM data must be log transformed for PCA
log2expr = np.log2(df + 1)  # log2(FPKM + 1) transformation

# Transpose data to have samples as rows and genes as columns
log2expr_T = log2expr.T     

print("Log2 transformed data shape:", log2expr_T.shape )

# Post transform QC making sure log transform behaves as expected
plt.hist(log2expr_T.values.flatten(), bins=50)
plt.title("Distribution of Log-transformed Expression Values")
plt.xlabel("Log2(FPKM + 1)")
plt.ylabel("Frequency")

# Scale the data
scaler = StandardScaler()
scaled_T = scaler.fit_transform(log2expr_T)

# PCA with multiple components
pca = PCA(n_components=8)
pca_T = pca.fit_transform(scaled_T)

# Explained variance
explained_var = pca.explained_variance_ratio_
print("Explained variance ratios:", pca.explained_variance_ratio_)

# Save PCA scores
pd.DataFrame(pca_T, columns=[f'PC{i+1}' for i in range(pca_T.shape[1])], index=log2expr_T.index).to_csv("results/pca_scores.csv")