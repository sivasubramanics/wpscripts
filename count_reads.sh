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
    echo "Usage: bash count_reads.sh /path/to/directory"
    exit 1
fi

if check_command rapidgzip
then
    unzip_exe=$(which rapidgzip)
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
count=0
for file in $(find $path -name "*.fastq.gz")
do
    file_count=$($unzip_exe -kdc $file | wc -l)
    count=$((count + (file_count / 4)))
done
count=$((count / 2))
echo -e "$path\t$count"
