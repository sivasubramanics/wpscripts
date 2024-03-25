#!/usr/bin/env python3


import argparse
from collections import defaultdict



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="input transcripts fasta file", required=True)
    parser.add_argument("-p", "--proteins", help="input proteins fasta file (evigene output: input.okay.aa)", required=True)
    parser.add_argument("-m", "--map", help="input map file (output of rename_isoforms. eg: input.map)", required=True)
    parser.add_argument("-o", "--output", help="output proteins fasta file", required=True)
    args = parser.parse_args()

    # read the map file and store the mapping in a dictionary
    map_dict = defaultdict()
    with open(args.map, "r") as map_file:
        for line in map_file:
            line = line.strip().split("\t")
            map_dict[line[1]] = line[0]

    # get the transcript ids from the fasta file and store them in a list
    trans_ids = []
    with open(args.fasta, "r") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                seq_id = line.split()[0][1:]
                if seq_id not in trans_ids:
                    trans_ids.append(seq_id)
                else:
                    print(f"Duplicate transcript id found: {seq_id}")
                    exit(1)

    # read the proteins fasta file, rename the proteins and write to the output file
    with open(args.proteins, "r") as prot_file, open(args.output, "w") as out_file:
        for line in prot_file:
            line = line.strip()
            if line.startswith(">"):
                seq_flag = False
                seq_id = line.split()[0][1:]
                if seq_id in map_dict:
                    new_id = map_dict[seq_id]
                if new_id in trans_ids:
                    out_file.write(f">{new_id}\n")
                    seq_flag = True
                    continue
            if seq_flag:
                out_file.write(f"{line}\n")
    print("Proteins renamed successfully")











if __name__ == "__main__":
    main()
