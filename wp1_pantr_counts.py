#!/usr/bin/env python3
# description: convert count matrix to pan matrix
# usage: wp1_pantr_counts.py

import sys
import os
import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='convert count matrix to pan matrix')
    parser.add_argument('-i', '--input', help='input count matrix', required=True)
    parser.add_argument('-m', '--map', help='PanID -> isoformID. Can be from mmseqs easy-clister tsv', required=True)
    parser.add_argument('-f', '--fai', help='fai file for the isoforms.fasta', required=True)
    parser.add_argument('-o', '--output', help='output pan matrix', required=True)
    args = parser.parse_args()
    
    # read fai
    isoform_len = defaultdict()
    with open(args.fai, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            isoform_len[line[0]] = int(line[1])
    f.close()

    # read map
    pan_map = defaultdict(list)
    pantr_count = 0
    prev = None
    with open(args.map, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[1] not in isoform_len:
                continue
            if prev != line[0] or prev == None:
                pantr_count += 1
                pantr = 'PanTR_' + str(pantr_count)       
            if pantr not in pan_map:
                pan_map[pantr] = []
            pan_map[pantr].append(line[1])
            prev = line[0] 
    f.close()

    # read count matrix
    count_matrix = defaultdict(list)
    line_num = 0
    with open(args.input, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            line_num += 1
            if line_num == 1:
                count_matrix['Geneid'] = line
                continue
            if line[0] not in isoform_len:
                continue
            line = normalize_count(line, isoform_len)
            count_matrix[line[0]] = list(map(float, line[1:]))
    f.close()

    # write pan matrix
    with open(args.output, 'w') as f:
        f.write('\t' + '\t'.join(count_matrix['Geneid']) + '\n')
        for pantr in pan_map:
            f.write(pantr + '\t' + join_elements(sum_coutns(pan_map[pantr], count_matrix), '\t') + '\n' )
    f.close()

    with open(args.output+'.map', 'w') as f:
        for pantr in pan_map:
            for isoform in pan_map[pantr]:
                if isoform in count_matrix and isoform in isoform_len:
                    f.write(pantr + '\t' + isoform + '\n')


def normalize_count(counts, isoform_len):
    iso_len = isoform_len[counts[0]]
    for idx,count in enumerate(counts[1:]):
        if float(count) == 0:
            counts[idx+1] = round(float(count), 2)
            continue
        counts[idx+1] = round((float(count) * 1000 / iso_len), 2)
    return counts

def sum_coutns(pantr, count_matrix):
    pantr_counts = []
    for isoform in pantr:
        if isoform in count_matrix:
            if len(pantr_counts) == 0:
                pantr_counts = count_matrix[isoform]
            else:
                for idx, count in enumerate(count_matrix[isoform]):
                    pantr_counts[idx] += count
    return list(map(round, pantr_counts))

def join_elements(ilist, sep):
    return sep.join(list(map(str, ilist)))

if __name__ == '__main__':
    main()

