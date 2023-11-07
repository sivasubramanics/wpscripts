#!/usr/bin/env python
# script to rename extra contigs and scaffolds and adjust its length according to the dict file

import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Rename extra contigs and scaffolds in a fai file')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-d', '--dict', help='dict file (eg: lsalpg.chr_dict.tsv)', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-c', '--chrcol', help='column number of chromosome name (default: 0)', default=0, type=int)
    parser.add_argument('-p', '--poscol', help='column number(s) of position (default: 1)', default=1, type=int, nargs='+')

    args = parser.parse_args()

    chr_dict = defaultdict(list)
    with open(args.dict, 'r') as f:
        # old_name, old_length, new_name, new_start, new_end
        for line in f:
            line = line.strip().split('\t')
            chr_dict[line[0]] = line
    with open(args.output, 'w') as f, open(args.input, 'r') as g:
        for i, line in enumerate(g):
            if i == 0:
                f.write(line)
                continue
            if line.startswith('#'):
                f.write(line)
                continue
            line = line.strip().split('\t')
            if line[args.chrcol] not in chr_dict:
                print(f"Error: {line[args.chrcol]} not in dict file")
                exit(1)
            for poscol in args.poscol:
                line[poscol] = int(line[poscol]) + int(chr_dict[line[args.chrcol]][3])
            line[args.chrcol] = chr_dict[line[args.chrcol]][2]
            f.write('\t'.join(map(str, line)) + '\n')





if __name__ == '__main__':
    main()