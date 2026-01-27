# Imports and processed Data
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from scipy.stats import ttest_rel
from pathlib import Path

# Loading data
expr = pd.read_csv('data/processed/log2_expression.csv', index_col=0)

# File paths
DATA_PATH = Path("data/raw/GSE260586_FPKM+1.txt")
fig_dir = Path("results/figures")
fig_dir.mkdir(parents=True, exist_ok=True)

print(expr.shape)
expr.head()

# Building metadata manually 
expr.columns
metadata = pd.DataFrame({
    'sample_id': expr.columns,
    'geo_accession': [
        'GSM8119803', 'GSM8119804', 'GSM8119805', 'GSM8119806',
        'GSM8119807', 'GSM8119808', 'GSM8119809', 'GSM8119810',
    ],
    'cell_line': [
        'HCC1937', 'HCC1937', 'SUM159-M1a-p44', 'SUM159-M1a-p44', 'MDA-468', 'MDA-468', '1833TR-p94', '1833TR-p94'
    ],
    'condition': ['DMSO control', 'MBE1.5 treated', 'DMSO control', 'MBE1.5 treated', 
                  'DMSO control', 'MBE1.5 treated', 'DMSO control', 'MBE1.5 treated'
    ]
}).set_index('sample_id')

print(expr.columns.tolist())
print()
print(metadata)

# Per-gene paired t-test across all cell lines
results = []

for gene in expr.index:
    treated = []
    control = []

    for cl in metadata['cell_line'].unique():
        sub = metadata[metadata['cell_line'] == cl]

        treated_sample = sub[sub['condition'] == 'MBE1.5 treated'].index[0]
        control_sample = sub[sub['condition'] == 'DMSO control'].index[0]

        treated.append(expr.loc[gene, treated_sample])
        control.append(expr.loc[gene, control_sample])

    stat, pval = ttest_rel(treated, control)
    log2fc = np.mean(treated) - np.mean(control)

    results.append([gene, log2fc, pval])

# Build DE results table
de_results = pd.DataFrame(
    results, 
    columns=['gene_id', 'log2FC', 'pval']
)

de_results['adj_pval'] = (de_results['pval'].rank(method='min') / len(de_results))

# Create significance column
log2FC_thresh = 1
padj_thresh = 0.05

de_results["significance"] = "Not Significant"

de_results.loc[
    (de_results["log2FC"] >= log2FC_thresh) & 
    (de_results["adj_pval"] < padj_thresh),
    "significance"
] = "Upregulated"

de_results.loc[
    (de_results["log2FC"] <= -log2FC_thresh) &
    (de_results["adj_pval"] < padj_thresh),
    "significance"
] = "Downregulated"


# Volcano plot
colors = {
    "Upregulated": "red",
    "Downregulated": "blue",
    "Not Significant": "lightgrey"
}

plt.figure(figsize=(8, 6))

for category, color in colors.items():
    subset = de_results[de_results['significance'] == category]
    plt.scatter(
        subset['log2FC'],
        -np.log10(subset['pval']),
        c=color,
        label=category,
        alpha=0.7
    )

# Threshold lines
plt.axvline(log2FC_thresh, color='black', linestyle='--', linewidth=1)
plt.axvline(-log2FC_thresh, color='black', linestyle='--', linewidth=1)
plt.axhline(-np.log10(padj_thresh), color='black', linestyle='--', linewidth=1)

plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted p-value')
plt.title('Differential Expression: MBE1.5 vs DMSO')
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(fig_dir/'volcano_plot.png')
plt.show()