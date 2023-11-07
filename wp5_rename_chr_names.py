#!/usr/bin/env python
# script to rename extra contigs and scaffolds in a fai file
# usage: python rename_extra_contigs.py -i <input_fai> -e <extra_contigs_1> <extra_contigs_2> .. <extra_contigs_n> -o <output_fai>

import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Rename extra contigs and scaffolds in a fai file')
    parser.add_argument('-i', '--input', help='input fai file', required=True)
    parser.add_argument('-r', '--reference', help='reference chromosome names in list', required=True)
    parser.add_argument('-e', '--extra', nargs='+', help='extra contigs and scaffolds names in list', required=True)
    parser.add_argument('-o', '--output', help='output fai file', required=True)
    args = parser.parse_args()

    # chr_lengths = read_fai(args.input)
    new_labels = defaultdict(list)
    chr_no = 0
    with open(args.reference, 'r') as f:
        for line in f:
            line = line.strip()
            chr_no += 1
            new_labels[line] = chr_no
    for extra in args.extra:
        chr_no += 1
        with open(extra, 'r') as f:
            for line in f:
                line = line.strip()
                new_labels[line] = chr_no
    last_new_label = None
    with open(args.output, 'w') as f, open(args.input, 'r') as g:
        for line in g:
            line = line.strip().split('\t')
            if line[0] not in new_labels:
                continue
            if new_labels[line[0]] == last_new_label:
                new_length += int(line[1])
            else:
                new_length = int(line[1])
            f.write(f"{line[0]}\t{line[1]}\t{new_labels[line[0]]}\t{new_length-int(line[1])}\t{new_length}\n")
            last_new_label = new_labels[line[0]]














if __name__ == '__main__':
    main()