#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict
def main():
    og_table = "/Users/selva001/projects/work/wp1/trinity_evigene/best_isoforms/on_expression/pantools_og.orthogroups.tsv"
    output_file = "/Users/selva001/projects/work/wp1/trinity_evigene/best_isoforms/on_expression/pantools_og.LFC.tsv"

    og_dict = {}
    deg_dict = defaultdict(dict)
    line_number = 0
    with open(og_table, 'r') as fh:
        for line in fh:
            line = line.strip()
            data = line.split()
            line_number += 1
            if line_number <= 2:
                continue
            for i in range(1, len(data)):
                genes = data[i].split(',')
                for gene in genes:
                    og_dict[gene] = data[0]
            # initialising deg_dict with list of 0s 9 elements long
            deg_dict[data[0]] = [0 for i in range(9)]

    parse_deg_table("/Users/selva001/projects/work/wp1/trinity_evigene/best_isoforms/on_expression/SAT_sat_DEG.tsv", deg_dict, og_dict, 0, 1, 2)
    parse_deg_table("/Users/selva001/projects/work/wp1/trinity_evigene/best_isoforms/on_expression/IL_il_DEG.tsv", deg_dict, og_dict, 3, 4, 5)
    parse_deg_table("/Users/selva001/projects/work/wp1/trinity_evigene/best_isoforms/on_expression/SAL_sal_DEG.tsv", deg_dict, og_dict, 6, 7, 8)

    with open(output_file, 'w') as fh:
        fh.write("OG")
        fh.write("\tSAT8\tSAT12\tSAT24")
        fh.write("\tIL8\tIL12\tIL24")
        fh.write("\tSAL8\tSAL12\tSAL24\n")
        for key in deg_dict:
            fh.write(key)
            for i in range(9):
                fh.write(f"\t{deg_dict[key][i]}")
            fh.write("\n")


def parse_deg_table(deg_file, deg_dict, og_dict, col_one, col_two, col_three):
    line_number = 0
    with open(deg_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            line_number += 1
            data = line.split()
            if line_number == 1:
                head_data = data
                continue
            if data[0] not in og_dict:
                continue
            if data[1] == 'NA':
                continue
            if float(data[1]) <= 0.05:
                if float(data[2]) <= 0.05:
                    deg_dict[og_dict[data[0]]][col_one] += float(data[3])
                if float(data[4]) <= 0.05:
                    deg_dict[og_dict[data[0]]][col_two] += float(data[5])
                if float(data[6]) <= 0.05:
                    deg_dict[og_dict[data[0]]][col_three] += float(data[7])


if __name__ == "__main__":
    main()