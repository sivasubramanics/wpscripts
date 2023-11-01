#!/usr/bin/env bash
# Since the output of orthofinder contains directories with too many files, it makes it difficult to tranfer. So this script will tar.gz some folder within the othofinder output directory.
# usage: compress_orthofinder.sh <orthofinder output directory> 

# set -e
set -o pipefail
# set -u
# set -x

usage(){
    echo "usage: compress_orthofinder.sh <orthofinder output directory/Results_XYZ>"
    exit 1
}

dir_not_exist(){
    # if tar.gz already exists, then skip
    if [ -f "$1.tar.gz" ]; then
        echo "Warning: $1.tar.gz already exists. and $1 doesn't exists. So skipping"
        return
    fi
    echo "Error: $1 does not exist. Are you sure you are providing the correct orthofinder output directory? [eg. Orthofinder/Results_Oct28/]"
    exit 1
}

if [ $# -ne 1 ]; then
    usage
fi

# Check if the orthofinder output directory exists
if [ ! -d "$1" ]; then
    echo "Error: $1 does not exist"
    usage
fi

# Check if tar and gzip are installed
if ! command -v tar &> /dev/null
then
    echo "Error: tar could not be found"
    exit 1
fi
if ! command -v gzip &> /dev/null
then
    echo "Error: gzip could not be found"
    exit 1
fi

# compress the directory
tar_gz(){
    # resolve the path
    path=$(realpath $1)
    # remove all the / at the end of the directory name if it exists
    while [[ $path == */ ]]; do
        path=${path::-1}
    done
    echo "Compressing $path"
    tar -czvf $path.tar.gz $path --remove-files  > $path.tar.gz.log 2>&1
    echo "Removing $path"
    # rm -rf $path
}

# $1 is Orthofinder/Results_*/
# directories to compress
# Resolved_Gene_Trees
# MultipleSequenceAlignments
# Gene_Trees
# WorkingDirectory
# Orthogroup_Sequences
# Single_Copy_Orthologue_Sequences

if [ -d "$1/Orthogroup_Sequences" ]; then
    tar_gz $1/Orthogroup_Sequences
else
    dir_not_exist "$1/Orthogroup_Sequences"
fi

if [ -d "$1/Single_Copy_Orthologue_Sequences" ]; then
    tar_gz $1/Single_Copy_Orthologue_Sequences
else
    dir_not_exist "$1/Single_Copy_Orthologue_Sequences"
fi

if [ -d "$1/MultipleSequenceAlignments" ]; then
    tar_gz $1/MultipleSequenceAlignments
else
    dir_not_exist "$1/MultipleSequenceAlignments"
fi

if [ -d "$1/Gene_Trees" ]; then
    tar_gz $1/Gene_Trees
else
    dir_not_exist "$1/Gene_Trees"
fi

if [ -d "$1/Resolved_Gene_Trees" ]; then
    tar_gz $1/Resolved_Gene_Trees
else
    dir_not_exist "$1/Resolved_Gene_Trees"
fi

if [ -d "$1/WorkingDirectory" ]; then
    tar_gz $1/WorkingDirectory
else
    dir_not_exist "$1/WorkingDirectory"
fi

echo "Done"
