# Secure Genomic Data Analysis 

## OVERVIEW:
This project explores transcriptomic changes associated with ALDH1A3 inhibition in breast cancer cell lines using an RNA-seq analysis workflow built in Python and R. The pipeline includes quality control, normalization, principal component analysis (PCA), differential expression analysis, and a short data integrity demonstration for secure genomics workflows. ALDH1A3 is an enzyme involved in retinoic acid synthesis and has been linked to tumor progression and immune suppression in several cancers. Recent work suggests that ALDH1A3-derived retinoic acid may suppress anti-tumor immune responses through paracrine signaling rather than directly affecting tumor cell growth. Using publicly available RNA-seq data from ALDH1A3 inhibition experiments, this project investigates how blocking ALDH1A3 alters gene expression patterns and compares findings against published results from Esposito et al. (2025).

## PROJECT GOALS:
- Demonstrate QC, normalization, PCA, and basic DE analysis on a small public dataset
- Show simple data integrity and encryption demonstration for secure genomics workflows
- Compare results with published findings from Esposito et al. (2025).

## BIOLOGICAL CONTEXT: 
ALDH1A3 converts retinaldehyde into all-trans retinoic acid (atRA), a signaling molecule historically associated with tumor suppression. However, recent evidence suggests that ALDH1A3-expressing tumors may become resistant to retinoid signaling while still using atRA to suppress anti-tumor immunity in the tumor microenvironment. In the referenced study, researchers developed the ALDH1A3 inhibitor MBE1.5 to investigate how ALDH1A3 inhibition reshapes gene expression in cancer cells (Esposito et al., 2025, iScience 28, 113675).

## EXPERIMENTAL DESIGN: 
Four breast cancer cell lines - MDA-MB-468, SUM159-M1a, and HCC1937 (ALDH1A3 positive) and 1833TR-p94 (ALDH1A3 negative) - were treated with inhibitor MBE1.5 vs 0.1% DMSO control for 9 hours. Total RNA was harvested and transcriptome-wide expression was measured by RNA-seq.

## DATA: 
GSE260586 (FPKM+1)

## REPO STRUCTURE
- 'data/raw/' raw input files (GSE260586_FPKM+1.txt)
- 'data/processed/' Log2 expression data
- 'notebooks/' Jupyter notebooks
- 'src/' scripts
- 'results/figures/' figures and outputs

## NOTEBOOKS
1. 01_exploratory_analysis:
Performed initial inspection and normalization of the expression matrix and explored global variance structure. Goals were to identify large-scale variance patterns, verify normalization behavior, and provide initial assessment of treatment effects. 
2. 02_stats_and_DE:
Differential expression was assessed across four breast cancer cell lines, comparing MBE1.5 treatment to DMSO controls within each cell line. Log2 fold changes were computed as the mean paired difference, and significance was assessed using a paired t-test.
3. 03_pca_analysis:
Validated PCA results using R to demonstrate cross-language reproducibility. Confirmed that observed variance structure is not dependent on a specific software system.
4. 04_security_demo:
Demonstrated basic secure data handling concepts relevant to genomic workflows. Demonstrated hashing and encryption to show awareness of data integrity and secure computational practices.

## TOOLS:
Python, R, pandas, numpy, scipy, matplotlib, scikit-learn, Git

## COMPARISON TO PUBLISHED FINDINGS:
Calculated PCA scores (PC1: 32.9%, PC2: 25.3%) closely match published values (PC1: 32.46%, PC2: 25.63%) (Esposito et al., 2025, iScience 28, 113675). Consistent with Figure 4E of the paper, clustering observed in PCA is due to cell line grouping rather than treatment with MBE1.5 or DMSO. 

![PCA labeled by cell line and treatment](results/figures/pca_labeled.png)

CYP1A1 and CYP1B1 appear as top upregulated genes, consistent with experimental findings which identified these as off-target aryl hydrocarbon receptor responses to MBE1.5 (Esposito et al., 2025, iScience 28, 113675). The absence of other significant hits corroborates the paper's conclusion that ALDH1A3 inhibition produces no cell-intrinsic transcriptomic response in these retinoid-insensitive cell lines.

![Volcano plot MBE1.5 vs DMSO](results/figures/volcano_plot.png)

## Notes and Limitations
- Sample grouping was derived from GEO Series metadata describing treatment with both MBE1.5 and DMSO control across four breast cancer cell lines
- Due to limited sample size, differential expression analysis serves as methodological demonstration for educational purposes. With only 2 replicates per condition per cell line, statistical power is limited. This analysis serves as a methodological demonstration; definitive biological inference would require DESeq2 or edgeR with more replicates
- The security demonstration is conceptual and educational, not a production deployment

## Key Takeaways
 
- PCA scores (PC1: 32.9%, PC2: 25.3%) closely matched published values, confirming the pipeline was implemented correctly on publicly available data.
- Variance structure is driven by cell line identity, not treatment condition — MBE1.5 does not dramatically alter the transcriptome of retinoid-insensitive lines.
- CYP1A1 and CYP1B1 are the only significantly upregulated genes, consistent with an off-target aryl hydrocarbon receptor response to MBE1.5.
- Cross-language validation (Python + R) confirmed results are not dependent on a specific computational tool.
- The security notebook highlights data integrity and encryption as important considerations for any pipeline handling patient-adjacent genomic data.
- **Future direction:** Apply more robust methods (DESeq2, edgeR) with additional replicates, and expand to ALDH1A3-positive immune cell co-culture models to investigate the paracrine immunosuppression hypothesis directly.

## Citations
1. Esposito M, Fang C, Wei Y, Pozzan A, Beato C, Su X, Hutton JE III, Reed T, Hang X, Perini ED, Wang W, Cheng X, Pan Y, Yu J, Kane M, Manoharan M, Proudfoot J, Cristea IM, Kang Y. Development of retinoid nuclear receptor pathway antagonists through targeting aldehyde dehydrogenase 1A3. *iScience*. 2025;28:113675. https://doi.org/10.1016/j.isci.2025.113675