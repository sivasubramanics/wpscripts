#!/usr/bin/env python

import argparse
from collections import defaultdict

TAB = '\t'
COMMA = ','
SEMICOLON = ';'
NEWLINE = '\n'


def get_node_to_tr(in_gfa):
    node_to_tr = defaultdict(list)
    with open(in_gfa, 'r') as f:
        for line in f:
            if line.startswith('P'):
                line = line.strip().split('\t')
                path_name = line[1]
                for node in line[2].split(','):
                    node_to_tr[node.rstrip('+-')].append(path_name)
    return node_to_tr


def parse_sam(in_sam, in_fasta, out_file):
    reads_to_tr = defaultdict(list)
    for sam_file in in_sam:
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                line = line.strip().split('\t')
                if line[2] == '*':
                    continue
                reads_to_tr[line[0]].append(line[2])

    with open(out_file, 'w') as f:
        for read in reads_to_tr:
            f.write(f"{read}\t{list_to_str(rm_dup(reads_to_tr[read]), COMMA)}\n")


def parse_gaf(in_gaf, in_gfa, out_file):
    node_to_tr = get_node_to_tr(in_gfa)
    reads_to_tr = defaultdict(list)
    for gaf_file in in_gaf:
        with open(gaf_file, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if ' ' in line[0]:
                    line[0] = line[0].split(' ')[0]
                if line[0] not in reads_to_tr:
                    reads_to_tr[line[0]] = []
                for node in get_nodes(line[5]):
                    reads_to_tr[line[0]] += node_to_tr[node]

    with open(out_file, 'w') as f:
        for read in reads_to_tr:
            f.write(f"{read}\t{list_to_str(rm_dup(reads_to_tr[read]), COMMA)}\n")


def rm_dup(ilist):
    return sorted(list(set(ilist)))


def get_nodes(aln_nodes):
    # replace < with > and split by >
    aln_nodes = aln_nodes.replace('<', '>')
    nodes = aln_nodes.split('>')
    return nodes[1:]


def list_to_str(ilist, sep):
    return sep.join(ilist)


def main():
    parser = argparse.ArgumentParser(
        description='parse sam/gaf files to get the number of reads mapped to each isofrom')
    parser.add_argument('-i', '--input', help='input sam/gaf file', required=True, type=str, nargs='+')
    parser.add_argument('-x', '--format', help='input file format', required=True, choices=['sam', 'gaf'])
    parser.add_argument('-r', '--reference', help='reference fasta/gfa', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    args = parser.parse_args()

    if args.format == 'sam':
        parse_sam(args.input, args.reference, args.output)
    else:
        parse_gaf(args.input, args.reference, args.output)


if __name__ == '__main__':
    main()
