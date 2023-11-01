#!/usr/bin/env python3
# script to rename fasta files from transcriptome assembly
# usage: rename_transcripts_fasta.py -f <fasta file> -p <prefix> -o <output file> -t <assembly type>
# example: rename_transcripts_fasta.py -f Trinity.fasta -p "Trinity_" -o Trinity_renamed.fasta -t Trinity

import argparse
import re
import sys

# create an instance of Argument Parser and add positional argument
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="input fasta file", required=True)
parser.add_argument("-p", "--prefix", help="prefix to add to fasta headers", required=True)
parser.add_argument("-o", "--output", help="output fasta file", required=True)
parser.add_argument("-t", "--type", help="assembly type", choices=["trinity", "rnaspades", "fasta"], required=True)
args = parser.parse_args()

def get_ids(line, asm_type):
    if asm_type == "trinity":
        # TRINITY_DN100000_c0_g1_i1
        g_id = re.search(r"(TRINITY_DN\d+_c\d+_g\d+)_i\d+", line).group(1)
        t_id = re.search(r"(TRINITY_DN\d+_c\d+_g\d+_i\d+)", line).group(1)
    elif asm_type == "rnaspades":
        # NODE_1_length_15401_cov_332.591075_g0_i0
        g_id = re.search(r"NODE_\d+_length_\d+_cov_\d+.\d+_(g\d+)_i\d+", line).group(1)
        t_id = re.search(r"(NODE_\d+_length_\d+_cov_\d+.\d+_g\d+_i\d+)", line).group(1)
    elif asm_type == "fasta":
        g_id = re.search(r"(.+)", line).group(1)
        t_id = re.search(r"(.+)", line).group(1)
    else:
        print("Assembly type not recognized")
        sys.exit(1)
    return g_id, t_id

# open fasta file and output file
fasta = open(args.fasta, "r")
output = open(args.output, "w")
out_map = open(args.output + ".map", "w")
in_gene_to_tmap = open(args.fasta + ".gene_to_transcript.map", "w")
out_gene_to_tmap = open(args.output + ".gene_to_transcript.map", "w")

fa_count = 1
with open(args.fasta, "r") as fasta:
    for line in fasta:
        if line.startswith(">"):
            line = line.strip()
            g_id, t_id = get_ids(line.split()[0][1:], args.type)
            output.write(">" + args.prefix + "_" + str(fa_count) + "\n")
            out_map.write(args.prefix + "_" + str(fa_count) + "\t" + line.split()[0][1:] + "\n")
            in_gene_to_tmap.write(g_id + "\t" + t_id + "\n")
            out_gene_to_tmap.write(g_id + "\t" + args.prefix + "_" + str(fa_count) + "\n")
            fa_count += 1
        else:
            output.write(line)
fasta.close()
output.close()
out_map.close()

