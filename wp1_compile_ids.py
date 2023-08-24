#!/usr/bin/env python3

import sys
import os

out_file = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/map_ids.tsv"
fo = open(out_file, "w")

file = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/all.trinity.renamed.fasta"
print(f"{file} processing")
with open(file, 'r') as f:
    for line in f:
        if line.startswith(">"):
            line = line.strip().split()
            from_id = line[0][1:]
            to_id = line[0][1:]
            fo.write(f"{from_id}\t{to_id}\t{file}\n")
        else:
            continue

file = "/lustre/BIF/nobackup/selva001/work/wp1/compile_data/all.trinity.renamed.fasta.pep"
print(f"{file} processing")
with open(file, 'r') as f:
    for line in f:
        if line.startswith(">"):
            line = line.strip().split()
            from_id = line[0][1:]
            to_id = from_id.split(".")[0]
            fo.write(f"{from_id}\t{to_id}\t{file}\n")
        else:
            continue


