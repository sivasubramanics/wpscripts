#!/bin/bash
# This script counts the number of reads in a directory of fastq.gz files
# Usage: bash untar.sh /path/to/directory

in_path=$(realpath $1)
for file in $(find $in_path -name "*.tar.");
do
    tar -xvf $file -C $(realpath $file)
done
