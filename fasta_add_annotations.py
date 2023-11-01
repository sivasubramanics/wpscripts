#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict

DBs = ["Coils", "Gene3D", "MobiDBLite", "PANTHER", "Pfam", "SUPERFAMILY", "TIGRFAM"]

parser = argparse.ArgumentParser(description='This script parses the interproscan tsv file and extracts the domain information for a given database and add those annotations to the fasta description line')
parser.add_argument('-i', '--input', help='interproscan tsv file', required=True)
parser.add_argument('-f', '--fasta', help='input fasta file', required=True)
parser.add_argument('-d', '--database', help='database name(s) comma seperated', required=True)
parser.add_argument('-o', '--output', help='output prefix', required=True)
parser.add_argument('-r', '--retain', help='retain original fasta description line', action='store_true')
args = parser.parse_args()

databases = []
if "," in args.database:
    databases = args.database.split(',')
    for database in databases:
        if database not in DBs:
            print("Database not recognized. Please choose from the following: " + str(DBs))
            sys.exit(1)
else:
    if args.databse not in DBs:
        print("Database not recognized. Please choose from the following: " + str(DBs))
        sys.exit(1)
    databases.append(args.databse)

ann_dict = defaultdict()
with open(args.input, 'r') as fh:
    print(f"Parsing interproscan tsv file {args.input}")
    for line in fh:
        data = line.rstrip().split("\t")
        if data[3] not in databases:
            continue
        if data[0] not in ann_dict:
            ann_dict[data[0]] = defaultdict()
        if data[3] not in ann_dict[data[0]]:
            ann_dict[data[0]][data[3]] = []
        if data[5] not in ann_dict[data[0]][data[3]] and data[5] != "-":
            ann_dict[data[0]][data[3]].append(data[5])

with open(args.fasta, 'r') as fh, open(f"{args.output}.fasta", 'w') as out_fh:
    print(f"Adding annotations to fasta file {args.fasta}")
    for line in fh:
        if line.startswith(">"):
            data = line.rstrip().split(" ")
            new_header = data[0]
            if data[0][1:] in ann_dict:
                for database in databases:
                    if database in ann_dict[data[0][1:]]:
                        new_header += f" {database}={';'.join(ann_dict[data[0][1:]][database])}"
                    else:
                        new_header += f" {database}=-"
            data[0] = new_header
            if args.retain:
                out_fh.write(" ".join(data) + "\n")
            else:
                out_fh.write(new_header + "\n")
        else:
            out_fh.write(line)

with open(f"{args.output}.ann", 'w') as out_fh:
    header = ["id", *databases]
    out_fh.write("# " + "\t".join(header) + "\n")
    for key, ann in ann_dict.items():
        out_fh.write(key + "\t")
        for database in databases:
            if database in ann:
                if ann[database] == []:
                    out_fh.write("-\t")
                else:
                    out_fh.write(";".join(ann[database]) + "\t")
            else:
                out_fh.write("-\t")
        out_fh.write("\n")
        




