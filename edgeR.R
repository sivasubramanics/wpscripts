#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 4) {
  stop("USAGE: Rscript deseq2.R sample_a sample_b counts.matrix dge_results.tsv", call.=FALSE)
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

data = read.table(counts_file, header=T, row.names=1, com='')
col_ordering = c(1,2,3,4,5,6)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = factor(c(rep(ctrl, 3), rep(treated, 3)))

exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateDisp(exp_study)
et = exactTest(exp_study, pair=c(ctrl, treated))
tTags = topTags(et,n=NULL)
result_table = tTags$table

# cpm_values <- cpm(exp_study)
# average_cpm_sampleA <- rowMeans(cpm_values[, exp_study$samples$group == ctrl])
# average_cpm_sampleB <- rowMeans(cpm_values[, exp_study$samples$group == treated])
# result_table$averageCPM_sampleA <- average_cpm_sampleA
# result_table$averageCPM_sampleB <- average_cpm_sampleB


result_table = data.frame(sampleA=ctrl, sampleB=treated, result_table)
result_table$logFC = -1 * result_table$logFC
write.table(result_table, file=results_file, sep='\t', quote=F, row.names=T, col.names=NA)
