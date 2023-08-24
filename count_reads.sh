#!/bin/bash
# This script counts the number of reads in a directory of fastq.gz files
# Usage: bash count_reads.sh /path/to/directory

check_command() {
    if command -v $1 &> /dev/null
    then
        return 0
    else
        return 1
    fi
}

if [ $# -eq 0 ]
then
    echo "Usage: bash count_reads.sh </path/to/directory> <file_ext>"
    exit 1
fi

rapidgzip_flag=0

if check_command rapidgzip
then
    unzip_exe=$(which rapidgzip)
    rapidgzip_flag=1
elif check_command pigz
then
    unzip_exe=$(which pigz)
elif check_command gzip
then
    unzip_exe=$(which gzip)
else
    echo "No unzip command found"
    exit 1
fi

path=$1
file_ext=$2
count=0
n_files=0
for file in $(find $path -name "*.$file_ext.gz")
do
    n_files=$((n_files + 1))
    if [ $n_files -eq 1 ]
    then
        echo "file,count"
    fi
    if [ $rapidgzip_flag -eq 1 ]
    then
        file_count=$($unzip_exe --count-lines $file 2> /dev/null)
    else
        file_count=$($unzip_exe -kdc $file | wc -l)
    fi
    count=$((count + (file_count / 4)))
    echo -e "$file,$count"
done
# count=$((count / 2))
# echo -e "$path\t$count"
