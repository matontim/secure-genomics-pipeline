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