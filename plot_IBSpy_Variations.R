#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

tsv_file <- NULL
chr <- NULL
name <- NULL
out_file <- NULL
cutoff <- 10

for (i in seq_along(args)) {
  if (args[i] == "--tsv" || args[i] == "-t") {
    tsv_file <- args[i + 1]
  } else if (args[i] == "--chr" || args[i] == "-c") {
    chr <- args[i + 1]
  } else if (args[i] == "--name" || args[i] == "-n") {
    name <- args[i + 1]
  } else if (args[i] == "--out" || args[i] == "-o") {
    out_file <- args[i + 1]
  } else if (args[i] == "--cutoff" || args[i] == "-p") {
    cutoff <- args[i + 1]
  } else if (args[i] == "--help" || args[i] == "-h") {
    cat("Usage: Rscript plot_IBSpy_variations.R --tsv variations.tsv --chr 1 --name LK002 --out LK002_chr2.jpeg --cutoff 100\n")
    q()
  }
}

if (is.null(tsv_file) || is.null(chr) || is.null(name) || is.null(out_file)) {
  stop("Usage: Rscript plot_IBSpy_variations.R --tsv variations.tsv --chr 1 --name LK002 --out LK002_chr2.jpeg --cutoff 100\n")
}

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(dplyr)) install.packages("dplyr")

# library(tidyverse)
df <- read.csv(tsv_file, sep = "\t", header = T)

print(chr)
# head(df)
if (chr == "all") {
  chr <- unique(df$seqname)
  is_facet <- TRUE
} else {
  # chr <- as.numeric(chr)
  is_facet <- FALSE
}
print(chr)
df$seqname <- as.factor(df$seqname)
df <- df[df$seqname %in% chr, ]
df$seqname <- as.factor(df$seqname)
df$variations <- as.numeric(df$variations)
# log2(max(df$variations))
print(head(df))

if (is_facet) {
  p <- ggplot(df, aes(x=start/1000000, y = variations, colour = variations < as.numeric(cutoff))) +
    geom_jitter(width = 0.3, size=0.2, alpha=0.9) +
    # scale_y_continuous(trans='log2') +
    labs(title=name,x ="Mbp", y = "variations") +
    theme(text=element_text(size=16,  family="Courier New")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none") +
    facet_grid(seqname ~ .)
    ggsave(out_file, plot = p, width = 15, height = 20, units = "in", dpi = 720)
} else {
  p <- ggplot(df, aes(x=start/1000000, y = variations, colour = variations < as.numeric(cutoff))) +
    geom_jitter(width = 0.3, size=0.2, alpha=0.9) +
    # scale_y_continuous(trans='log2') +
    labs(title=name,x ="Mbp", y = "variations") +
    theme(text=element_text(size=16,  family="Courier New")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")
    ggsave(out_file, plot = p, width = 15, height = 5, units = "in", dpi = 720)
}


