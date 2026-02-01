#!/usr/bin/env Rscript
library(DESeq2)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

counts_data <- read.table(
  snakemake@input[[1]], header = TRUE, row.names = 1, check.names = FALSE
)
counts_data <- counts_data[, order(names(counts_data))]

col_data <- read.table(
  snakemake@params[['samples']], header = TRUE, row.names = 1, check.names = FALSE, fill = TRUE
)
col_data <- col_data[order(rownames(col_data)), , drop = FALSE]
col_data <- col_data %>% select(starts_with('condition'))

for (i in 1:ncol(col_data)) {
  if (is.character(col_data[[i]]) || is.logical(col_data[[i]])) {
    col_data[[i]] <- as.factor(col_data[[i]])
  }
}

dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = as.formula(snakemake@params[['design']])
)
dds <- DESeq(dds)

rld <- rlog(dds, blind = TRUE)
sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- sprintf(
  '%s(%s.%s)', rownames(sample_dist_matrix), rld$condition_1, rld$condition_2
)
colnames(sample_dist_matrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, 'Oranges')))(255)
cluster_plot <- pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  col = colors,
  filename = snakemake@output[['cluster']],
  width = 8,
  height = 6
)

pca_data <- plotPCA(rld, intgroup = c('condition_1', 'condition_2'), returnData=TRUE)
percent_var <- round(100 * attr(pca_data, 'percentVar'))
pca_plot <- ggplot(
  pca_data,
  aes(PC1, PC2, color = condition_1, shape = condition_2)
) +
  theme_minimal() +
  geom_point(size = 3) +
  xlab(paste0('PC1: ', percent_var[1], '% variance')) +
  ylab(paste0('PC2: ', percent_var[2], '% variance')) +
  coord_fixed() +
  geom_text_repel(
    aes(label = rownames(pca_data)),
    size = 3,
    max.overlaps = Inf
  )

ggsave(snakemake@output[['pca']], plot = pca_plot)
