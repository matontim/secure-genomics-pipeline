# Pull processed data
data <- read.csv(
    "data/processed/log2_expression.csv", 
    row.names = 1, 
    check.names = FALSE
    )
# Transpose so rows are samples and columns are genes
data_t <- t(data)

# Perform PCA without scaling
pca <- prcomp(data_t, scale. = FALSE)
summary(pca)

# PCA visualization
var_explained <-- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

png("results/figures/pca_r.png", width=1200, height=1000, res=150)
par(
    cex.lab=1.2,
    cex.axis=1.0,
    cex.main=1.5
)
plot(
    pca$x[,1], pca$x[,2], col="blue", cex=2.0, pch=19,
    xlab = paste0("PC1 (", var_explained[1], "% variance)"),
    ylab = paste0("PC2 (", var_explained[2], "% variance)"),
    main = "PCA of GSE260586 Expression Data (R)")

dev.off()