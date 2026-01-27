data <- read.csv("/data/processed/log2_expression.csv", row.names = 1)

pca <- prcomp(data, scale. = FALSE)
summary(pca)

