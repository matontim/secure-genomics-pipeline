Secure Genomic Data Analysis 

OVERVIEW:
Beginner-friendly RNA-seq analysis workflow demonstrating QC and exploratory analysis, normalization PCA, basic differential exploration, and a short data integrity demo 

PROJECT GOALS:
Demonstrate QC, normalization PCA, and basic DE analysis on a small public dataset
Show simple data integrity and encryption demonstration for secure genoics workflows 

DATA: GSE260586 (FPKM+1)

REPO STRUCTURE
-'data/raw/' raw input files (GSE260586_FPKM+1.txt)
-'data/processed/' Log2 expression data
-'notebooks/' Jupyter notebooks
-'src/' scripts
-'results/figures/' figures and outputs

NOTEBOOKS
01_exploratory_analysis
02_stats_and_DE
    Differential expression was assess across four breast cancer cell lines, comparing MBE1.5 treatment to DMSO controls within each cell line. Log2 fold changes were computed as the mean paired difference, and significance was assessed using a paired t-test

TOOLS:
Python, pandas, numpy, matplotlib, scikit-learn, Git

NOTES & LIMITATIONS:

-Exploratory Analysis 
-Stats and DE
    Sample grouping was derived from GEO Series metadata describing treatment with both MBE1.5 and DMSO control across four breast cancer cell lines. Due to limited sample size, differential expression analysis serves as methodological demonstration for educational purposes. More robust methods, such as DESeq2, are used for definitive biological inference.