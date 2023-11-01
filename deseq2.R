#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 4) {
  stop("USAGE: Rscript deseq2.R sample_a sample_b counts.matrix dge_results.tsv", call.=FALSE)
}

if (! suppressMessages(require(DESeq2))) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
  library(DESeq2)
}

if (! suppressMessages(require(apeglm))) {
  install.packages("BiocManager")
  BiocManager::install("apeglm")
  library(apeglm)
}

if (! suppressMessages(require(edgeR))) {
  install.packages("BiocManager")
  BiocManager::install("edgeR")
  library(edgeR)
}

ctrl <- args[1]
treated <- args[2]
counts_file <- args[3]
results_file <- args[4]
# cut-off for the significance filtering
# alpha_thrshold <- as.numeric(args[5])
# fc_thrshold <- as.numeric(args[6])
alpha_thrshold <- 0.05
fc_thrshold <- 1

counts_data <- read.csv(file = counts_file, row.names = 1, header = TRUE, sep = "\t") # nolint
counts_data <- round(counts_data)
head(counts_data)
counts_data = counts_data[rowSums(cpm(counts_data) > 1) >= 2,]
head(counts_data)
conditions = data.frame(conditions=factor(c(rep(ctrl, 3), rep(treated, 3))))
rownames(conditions) = colnames(counts_data)
conditions
dds <- DESeqDataSetFromMatrix(countData=counts_data, 
                              colData=conditions, 
                              design=~ conditions)
# filter lowly expressed genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds <- DESeq(dds)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == ctrl])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == treated])
contrast <- c("conditions",ctrl,treated)
res <- results(dds, contrast, alpha = alpha_thrshold, lfcThreshold = fc_thrshold)
res <- results(dds, contrast)
# res <- lfcShrink(dds, coef=2, type="apeglm")
res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
res <- cbind(sampleA=ctrl, sampleB=treated, as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
# res <- subset(res, padj<alpha_thrshold)
res <- as.data.frame(res[order(res$pvalue),])
write.table(res, file=results_file, sep='\t', quote=FALSE, col.names=NA)

