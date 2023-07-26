#!/bin/bash
# script: which_env.sh
# description: find the conda environment for a tool
# date: 2021-07-26
# version: 0.0.1

set -u
set -e

env_file="/lustre/BIF/nobackup/selva001/work/envs/env.tools"

echoerr() { echo "$@" 1>&2; }

if [ $# -eq 0 ]; then
    echo "Usage: $0 <tool1> <tool2> ... <toolN>"
    exit 1
fi

if [ ! -f $env_file ]; then
    echo "ERROR: $env not found"
    echo "ERROR: please run list_conda_tools.sh"
    exit 1
fi

random_filename=$(mktemp)

echo -e "ENVIRONMENT\tTOOL\tVERSION" > $random_filename
echo -e "-----------\t----\t-------" >> $random_filename
for tool in $@; do
    awk -v tool=$tool '{if($2==tool){print $1"\t"$2"\t"$3}}' $env_file
done | sort -k1,1 >> $random_filename

for tool in $@; do
    if [ $(cut -f2 $random_filename | grep -w $tool | wc -l) -eq 0 ]; then
        echoerr "ERROR: no environment found for $tool"
        echoerr "You may have to provide the conda package name instead. eg: ucsc-fasomerecords instead of faSomeRecords?"
        echoerr "--------------------------------------------------"
    fi
done


cat $random_filename | column -t -s $'\t'
rm $random_filename
