#!/bin/bash
# script: list_all_conda_tools.sh
# description: list all the tools installed in conda environments
# date: 2021-07-26
# version: 0.0.2

set -u
set -e
set -o pipefail  # Fail script if any command in a pipeline fails

# check if /lustre/BIF/nobackup/selva001/ exists
if [[ ! -d "/lustre/BIF/nobackup/selva001/" ]]; then
    out_file="/Users/selva001/opt/envs/env.tools"
    out_file_bak="/Users/selva001/opt/envs/env.tools.bak"
else
    out_file="/lustre/BIF/nobackup/selva001/work/envs/env.tools"
    out_file_bak="/lustre/BIF/nobackup/selva001/work/envs/env.tools.bak"
fi

echo "Output file: $out_file"

# Get list of conda environments
conda_envs=$(conda info --envs | grep -v "#" | awk '{print $1}')

# Check if conda_envs is empty
if [[ -z "$conda_envs" ]]; then
    echo "No conda environments found."
    exit 1
fi

# Iterate over each conda environment
for conda_env in $conda_envs; do
    if [[ -z "$conda_env" ]]; then
        continue
    fi
    if [[ "$conda_env" =~ "/" ]]; then
        conda list -p "$conda_env" | grep -v "#" | awk -v env="$conda_env" '{print env"\t"$1"\t"$2}' || {
            echo "Error listing packages for environment: $conda_env"
            continue
        }
    else
        conda list -n "$conda_env" | grep -v "#" | awk -v env="$conda_env" '{print env"\t"$1"\t"$2}' || {
            echo "Error listing packages for environment: $conda_env"
            continue
        }
    fi
done > "$out_file_bak"

# Move the backup file to the final output file
mv "$out_file_bak" "$out_file"

echo "Tool list saved to: $out_file"
echo "Done."