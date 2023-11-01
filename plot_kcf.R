#!/usr/bin/env Rscript

# required packages
requiredPackages = c('ggplot2')
suppressMessages(
  for (p in requiredPackages) {
    if (!require(p, character.only = TRUE)){
      install.packages(p)
    }
    library(p, character.only = TRUE)
  }
)

# extract the score from the SAMPLE column
extract_SC <- function(data, position) {
  as.numeric(sapply(strsplit(as.character(data$SAMPLE), ":"), "[", position))
}

# get the index of the tag from the FORMAT column
getIndex <- function(str, tag) {
  elements <- strsplit(as.character(str), ":")[[1]]
  index <- match(tag, elements)
  return(index)
}

# check if the sample is present in the dataframe
checkSample <- function(df, sample) {
  if (!(sample %in% colnames(df))) {
    cat("Sample is not present in the dataframe column names.")
    q()
  }
}

start <- NULL
end <- NULL
window <- NULL

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
# print(args)
for (i in seq_along(args)) {
  if (args[i] == "--kcf" || args[i] == "-k") {
    in_kcf <- args[i + 1]
  } else if (args[i] == "--chr" || args[i] == "-c") {
    chrom <- args[i + 1]
  } else if (args[i] == "--sample" || args[i] == "-s") {
    sample <- args[i + 1]
  } else if (args[i] == "--out" || args[i] == "-o") {
    out_file <- args[i + 1]
  } else if (args[i] == "--cutoff" || args[i] == "-p") {
    cutoff <- args[i + 1]
  } else if (args[i] == "--tag" || args[i] == "-t") {
    tag <- args[i + 1]
  } else if (args[i] == "--window" || args[i] == "-w") {
    window <- args[i + 1]
  } else if (args[i] == "--start" || args[i] == "-i") {
    start <- args[i + 1]
  } else if (args[i] == "--end" || args[i] == "-j") {
    end <- args[i + 1]
  } else if (args[i] == "--help" || args[i] == "-h") {
    cat("Usage: Rscript plot_kcf.R --kcf variations.kcf --chr 1 --sample LK002 --out LK002_chr2.jpeg --tag SC --cutoff 0.9\n")
    cat("Options:\n")
    cat("\t--kcf, -k\t\tKCF file\n")
    cat("\t--chr, -c\t\tChromosome\n")
    cat("\t--sample, -s\t\tSample name\n")
    cat("\t--out, -o\t\tOutput file name\n")
    cat("\t--tag, -t\t\tTag name\n")
    cat("\t--cutoff, -p\t\tCutoff value\n")
    cat("\t--window, -w\t\tWindow size\n")
    cat("\t--start, -i\t\tStart position\n")
    cat("\t--end, -j\t\tEnd position\n")
    cat("\t--help, -h\t\tPrint this help message\n")    
    q()
  }
}


# in_kcf <- "all.5k.ibs.kcf"
# chrom <- "Dovetail_reconstructed_14OCT2019_v4_lg_2"
# sample <- "LK421"
# tag <- "SC"
# cutoff <- 0.8 
# window <- 10000
# # start <- 50000000
# # end <- 52000000


# read KCF file
lines <- readLines(in_kcf)
# remove commented line from KCF file
lines <- lines[!startsWith(lines, "##")]
# read KCF data as dataframe
data <- read.delim(textConnection(lines), header = TRUE)
dim(data)
# check if the querid sample is there in the input KCF file
checkSample(data, sample)

# change the column name of #CHROM to CHROM
colnames(data)[colnames(data) == "X.CHROM"] <- "CHROM"
# change the column name
colnames(data)[colnames(data) == sample] <- "SAMPLE"

# subset the data for the chromosome
data <- data[data$CHROM == chrom, ]
chrom
dim(data)
if (is.null(start)) {
  start <- min(data$START)
}
if (is.null(end)) {
  end <- max(data$END)
}
# subset the data for the start and end position
data <- data[data$START >= as.numeric(start) & data$START <= as.numeric(end), ]
# get the index of SCORE tag from FORMAT column
sc_position <- getIndex(data$FORMAT[1], tag)
# extract data for the tag queried
data$SC <- extract_SC(data, sc_position)

# bin the data if windows have to be changed
if (is.null(window)) {
  plot_df <- data.frame(
    START = data$START,
    SC = data$SC,
    CHROM = data$CHROM
  )
} else {

  data$bin <- cut(data$START, breaks = seq(min(data$START), max(data$START), by = as.numeric(window)), include.lowest = TRUE)
  # get mean values for the tag
  if (tag == "SC"){
    SC <- tapply(data$SC, data$bin, mean, na.rm = TRUE)  
  }
  if (tag == "VA"){
    SC <- tapply(data$SC, data$bin, sum, na.rm = TRUE)  
  }
  # get mid points of the window for plotting
  midpoints <- sapply(strsplit(gsub("\\(|\\]|\\[", "", levels(data$bin)), ","), function(x) mean(as.numeric(x)))
  # constuct a new dataframe with START and SC
  plot_df <- data.frame(
    START = midpoints,
    SC = SC,
    CHROM = data$CHROM[match(names(SC), data$bin)]
  )
}

head(plot_df)
# remove the NA values
plot_df <- plot_df[!is.na(plot_df$CHROM), ]
# library("tidyverse")
# regions <- tibble(x1 = 20.4, x2 = 20.5, y1 = 0, y2 = max(plot_df$SC))

# plotting
p <- ggplot(plot_df, aes(x=START/1000000, y = SC, colour = SC < as.numeric(cutoff))) +
  geom_jitter(width = 0.3, size=0.5, alpha=0.9) +
  # scale_y_continuous(trans='log2') +
  labs(title=sample,x ="Mbp", y = "score [(obs_kmers - var)/tot_kmer]") +
  # theme(text=element_text(size=16,  family="Courier New")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
  # geom_rect(data=regions, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
  #           color = "transparent",
  #           fill = "blue",
  #           alpha=0.5,
  #           inherit.aes = FALSE)

# write the plot
ggsave(out_file, plot = p, device = "pdf", width = 20, height = 6, units = "in", dpi = 720)

