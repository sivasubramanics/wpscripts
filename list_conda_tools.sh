#!/bin/bash
# script: list_all_conda_tools.sh
# description: list all the tools installed in conda environments
# date: 2021-07-26
# version: 0.0.1

set -u
set -e

out_file="/lustre/BIF/nobackup/selva001/work/envs/env.tools"
out_file_bak="/lustre/BIF/nobackup/selva001/work/envs/env.tools.bak"

conda_envs=$(conda info --envs | grep -v "#" | awk '{print $1}')
for conda_env in $conda_envs; do
    if [[ $conda_env == "" ]]; then
        continue
    fi
    if [[ $conda_env =~ "/" ]]; then
        conda list -p $conda_env | grep -v "#" | awk -v env=$conda_env '{print env"\t"$1"\t"$2}'
        continue
    fi
    conda list -n $conda_env | grep -v "#" | awk -v env=$conda_env '{print env"\t"$1"\t"$2}'
done > $out_file_bak

mv $out_file_bak $out_file