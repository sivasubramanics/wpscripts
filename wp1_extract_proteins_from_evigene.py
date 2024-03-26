#!/usr/bin/env python3


import argparse
from collections import defaultdict



def main():

    parser = argparse.ArgumentParser()
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

    # read the proteins fasta file, rename the proteins and write to the output file
    with open(args.proteins, "r") as prot_file, open(args.output, "w") as out_file:
        for line in prot_file:
            line = line.strip()
            if line.startswith(">"):
                seq_flag = False
                seq_id = line.split()[0][1:]
                if seq_id in map_dict:
                    new_id = map_dict[seq_id]
                    out_file.write(f">{new_id}\n")
                    seq_flag = True
                    continue
            if seq_flag:
                out_file.write(f"{line}\n")
    print("Proteins renamed successfully")











if __name__ == "__main__":
    main()
