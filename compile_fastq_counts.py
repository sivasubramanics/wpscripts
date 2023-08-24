#!/usr/bin/env python3

import sys
import os

def run_cmd(cmd):
    print(f"Running command: {cmd}")
    os.system(cmd)

def get_count(ifile):
    with open(ifile, 'r') as f:
        line_count = 0
        for line in f:
            line = line.strip().split()
            line_count += 1
            if line_count == 2:
                return(int(line[4].replace(',', '')))
    return 0

def get_base_count(file_list):
    base_count = 0
    with open(file_list, 'r') as f:
        for file in f:
            fq_name = file.strip()
            file = file.strip().replace('.fq.gz', '.seqkit.stats.txt')
            file = file.strip().replace('.fastq.gz', '.seqkit.stats.txt')
            if os.path.isfile(file):
                base_count += get_count(file)    
            else:
                # cmd = "seqkit stats " + fq_name + " > " + file
                # run_cmd(cmd)
                # base_count += get_count(file)
                print(f"File {file} does not exist for {fq_name}")
                continue
    return base_count

def main():
    if len(sys.argv) < 2:
        print("Usage: python compile_fastq_counts.py list_of_files.txt")
        sys.exit(1)
    
    for file in sys.argv[1:]:
        base_count = get_base_count(file)
        print(f"{file}\t{base_count}")

if __name__ == '__main__':
    main()