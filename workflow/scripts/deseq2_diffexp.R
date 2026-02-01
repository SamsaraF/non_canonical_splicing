#!/usr/bin/env Rscript
library(DESeq2)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)

counts_data <- read.table(
  snakemake@input[[1]], header = TRUE, row.names = 1, check.names = FALSE
)

col_data <- read.table(
  snakemake@params[['samples']], header = TRUE, row.names = 1, check.names = FALSE, fill = TRUE
)
col_data <- col_data %>% select(starts_with('condition'))
col_data <- col_data %>% filter(!!sym(snakemake@params[['filter_by']]) == snakemake@params[['filter_value']])

counts_data <- counts_data[, rownames(col_data)]
counts_data <- counts_data[, order(names(counts_data))]
col_data <- col_data[order(rownames(col_data)), , drop = FALSE]

for (i in 1:ncol(col_data)) {
  if (is.character(col_data[[i]]) || is.logical(col_data[[i]])) {
    col_data[[i]] <- as.factor(col_data[[i]])
  }
}

contrast <- snakemake@params[['contrast']]
design <- paste('~', contrast[[1]])
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = as.formula(design)
)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

res <- results(dds, contrast = contrast)
res <- lfcShrink(dds, contrast = contrast, res = res, type = 'ashr')
res <- res[order(res$padj), ]

# MA plot
pdf(snakemake@output[['ma_plot']])
plotMA(res, ylim = c(-2, 2))
dev.off()

# load ensembl id to symbol mapping
ensid_to_symbol <- read.table(
  snakemake@params[['id_to_symbol']], header = TRUE, sep = '\t', row.names = 1
)
symbol_map <- ensid_to_symbol[, 1]
names(symbol_map) <- rownames(ensid_to_symbol)

# save normalized counts with DESeq2 results
res_data <- merge(
  as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)),
  by = 'row.names', sort = FALSE
)
res_data$symbol <- symbol_map[res_data$Row.names]
write.table(
  res_data,
  file = snakemake@output[['norm']],
  row.names = FALSE,
  sep = '\t'
)

# volcano plot
volcano_pdf <- EnhancedVolcano(
  res, lab = Map(function(arr) symbol_map[[arr]], rownames(res)),
  x = 'log2FoldChange', y = 'padj', pCutoff = 0.05
)
ggsave(snakemake@output[['volcano']], plot = volcano_pdf, device = 'pdf')
