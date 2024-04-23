#!/usr/bin/env Rscript
rm(list = ls())

usage <- function(){
  cat("Usage      : vlookup.R <file1.tsv> <file2.tsv> <column_name> <out.tsv> [tobe_added]\n")
  cat("file1.tsv  : First file to be merged\n")
  cat("file2.tsv  : Second file to be merged\n")
  cat("column_name: Name of the column that is common in both files\n")
  cat("out.tsv    : Output file\n")
  cat("tobe_added : [Optional] Names of the columns to be added from the second file. If not provided, all columns from the second file will be added\n")
}

# if no argument is passed, print the usage
if (length(commandArgs(trailingOnly = TRUE)) == 0){
  usage()
  quit()
}

# read system arguments in this order: file1.tsv file2.tsv common_column_name out.tsv [tobe_added_col_names] and store in variables
sys_args <- commandArgs(trailingOnly = TRUE)
in_file_one <- sys_args[1]
in_file_two <- sys_args[2]
common_col_name <- sys_args[3]
out_file <- sys_args[4]
if (length(sys_args) >= 5){
  tobe_added_col_names <- sys_args[5:length(sys_args)]
} else {
  tobe_added_col_names <- NULL
}

in_data_one <- read.delim(in_file_one, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "#", quote = "")
in_data_two <- read.delim(in_file_two, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "#", quote = "")

# check if the common column name exists in both files
if (common_col_name %in% colnames(in_data_one) & common_col_name %in% colnames(in_data_two)){
  # check if the columns to be added exist in the second file
  if (length(tobe_added_col_names) > 0){
    # check if the columns to be added exist in the second file
    if (all(tobe_added_col_names %in% colnames(in_data_two))){
      in_data_two <- in_data_two[, c(common_col_name, tobe_added_col_names)]
    } else {
      stop("One or more columns to be added do not exist in the second file")
    }
  }
  # merge the two files
  merged_data <- merge(in_data_one, in_data_two, by = common_col_name)
} else {
  stop("The common column name does not exist in both files")
}

# write the merged data to the output file
write.table(merged_data, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
