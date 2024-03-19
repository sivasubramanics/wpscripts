#!/usr/bin/env python3

import sys
from collections import defaultdict
import argparse


def main():
    # for i in IBS.LK???.tsv;do sample=${i %.tsv};sample=${sample/IBS.}; awk -v s=$sample '{if($5>=400000 && NR>1){print s"\t"$0}}' ${i}; | awk '{if(NR==1){print "sample\tid\tseqname\tstart\tend\tlength\ttotal_blocks\tibs_blocks\tmean_score"};print}' > IBS_regions.tsv
    # awk '{if(NR==1){print $0"\tprop_blocks"}else{sc=($8/$7); if(sc>=0.75){print $0"\t"sc}}}' IBS_regions.tsv > IBS_regions.flt.tsv

    parser = argparse.ArgumentParser(description='This is a utility script to process IBS files, and perform simple operations')
    parser.add_argument('-i', '--input', help='input IBS file', required=True)
    parser.add_argument('-l', '--length', help='length file', required=True)
    parser.add_argument('-o', '--output', help='output table', required=True)
    parser.add_argument('-r', '--reverse', help='reverse the score cutoff', required=False, action='store_true')
    parser.add_argument('-s', '--score', help='score cutoff', required=True, type=float)
    args = parser.parse_args()

    len_dict = defaultdict()
    with open(args.length, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            len_dict[line[0]] = int(line[1])

    samples = []
    with open(args.input, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[1] == 'id':
                continue
            sample = line[0]
            if sample not in samples:
                samples.append(sample)

    out_matrix = defaultdict()
    for seqname in len_dict:
        for i in range(0, len_dict[seqname], 500000):
            win = int(i / 500000) * 0.5
            out_matrix[(seqname, win)] = [0] * len(samples)

    with open(args.input, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[1] == 'id':
                continue
            window_start = int(float(line[3]) / 500000)
            window_end = int(float(line[4]) / 500000)
            if args.reverse:
                if float(line[8]) <= args.score:
                    for i in range(window_start, window_end + 1):
                        out_matrix[(line[2], i * 0.5)][samples.index(line[0])] = line[8]
            else:
                if float(line[8]) >= args.score:
                    for i in range(window_start, window_end + 1):
                        out_matrix[(line[2], i * 0.5)][samples.index(line[0])] = line[8]

    with open(args.output, 'w') as f:
        f.write('seqname\tstart\t' + '\t'.join(samples) + '\n')
        for key in out_matrix:
            f.write(f'{key[0]}\t{key[1]}\t' + '\t'.join(map(str, out_matrix[key])) + '\n')

if __name__ == "__main__":
    main()