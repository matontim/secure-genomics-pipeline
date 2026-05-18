# Pull processed data
data <- read.csv(
    "data/processed/log2_expression.csv", 
    row.names = 1, 
    check.names = FALSE
    )
# Transpose so rows are samples and columns are genes
data_t <- t(data)

# Perform PCA without scaling
pca <- prcomp(data_t, scale. = TRUE)
summary(pca)

# PCA visualization
var_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

# Sample metadata confirmed from GEO GSE260586
sample_names <- c('r542', 'r543', 'r544', 'r545', 'r546', 'r547', 'r548', 'r549')
cell_lines   <- c('HCC1937', 'HCC1937',
                  'SUM159-M1a', 'SUM159-M1a',
                  'MDA-468', 'MDA-468',
                  '1833TR-p94', '1833TR-p94')
treatments   <- c('DMSO', 'MBE1.5', 'DMSO', 'MBE1.5',
                  'DMSO', 'MBE1.5', 'DMSO', 'MBE1.5')

# Colors by cell line
color_map <- c(
    'HCC1937'    = 'steelblue',
    'SUM159-M1a' = 'darkorange',
    'MDA-468'    = 'green4',
    '1833TR-p94' = 'red'       
)

# Shapes by treatment
shape_map <- c('DMSO' = 16, 'MBE1.5' = 17)  

point_colors <- color_map[cell_lines]
point_shapes <- shape_map[treatments]

png("results/figures/pca_r.png", width=1400, height=1100, res=150)
par(cex.lab=1.2, cex.axis=1.0, cex.main=1.4,
    mar=c(5, 5, 4, 14), xpd=TRUE)   

plot(
    pca$x[,1], pca$x[,2],
    col = point_colors,
    pch = point_shapes,
    cex = 2.0,
    xlab = paste0("PC1 (", var_explained[1], "% variance)"),
    ylab = paste0("PC2 (", var_explained[2], "% variance)"),
    main = "PCA of GSE260586 Expression Data (R)"
)

# Cell line legend
legend("topright", inset=c(-0.3, 0),
       legend = names(color_map),
       col    = color_map,
       pch    = 16, pt.cex = 1.5,
       title  = "Cell Line", bty="n", cex=0.85)

# Treatment legend
legend("topright", inset=c(-0.25, 0.4),
       legend = names(shape_map),
       col    = "black",
       pch    = shape_map, pt.cex = 1.5,
       title  = "Treatment", bty="n", cex=0.85)

dev.off()