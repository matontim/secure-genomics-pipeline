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

#Transpose so rows are samples and columns are genes
df_t = df.T

print("Transposed DataFrame shape:", df_t.shape)
print(df.T.head())

#Scale the data
scaler = StandardScaler()
df_scaled = scaler.fit_transform(df_t)
