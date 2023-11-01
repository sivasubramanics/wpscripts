#!/usr/bin/env python3

import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Get FC from DEG results')
    parser.add_argument('-i', '--in_file', type=str, help='path to input file. eg:isoforms_dge.tsv')
    parser.add_argument('-o', '--out_file', type=str, help='path to output file. eg:isoforms_dge_fc.tsv')
    parser.add_argument('-c', '--col', type=int, help='column number of FC', default=7)
    parser.add_argument('-g', '--genelist', type=str, help='file containing list of genes', default=None)
    parser.add_argument('-p', '--order', type=str, help='file containing order of samples', default=None)
    args = parser.parse_args()

    fc_dict = defaultdict(dict)
    line_num = 0
    sample_order = []
    gene_list = []
    with open(args.in_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            line_num += 1
            if line_num == 1:
                continue
            if line[0] not in fc_dict:
                fc_dict[line[0]] = defaultdict(dict)
            if line[2] not in fc_dict[line[0]]:
                fc_dict[line[0]][line[2]] = round(float(line[args.col]), 2)
            if line[2] not in sample_order:
                sample_order.append(line[2])
    f.close()

    gene_list = list(fc_dict.keys())

    if args.genelist != None:
        gene_list = []
        with open(args.genelist, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[0] == '' or line[0] == 'gene_id':
                    continue
                if line[0] not in gene_list:
                    gene_list.append(line[0])
        f.close()

    if args.order != None:
        sample_order = []
        with open(args.order, 'r') as f:
            for line in f:
                line = line.strip()
                if line == '' or line == 'gene_id':
                    continue
                if line not in sample_order:
                    sample_order.append(line)
        f.close()

    with open(args.out_file, 'w') as fo:
        fo.write('Geneid\t' + '\t'.join(sample_order) + '\n')
        for gene in gene_list:
            fo.write(gene)
            for sample in sample_order:
                if sample in fc_dict[gene]:
                    fo.write('\t' + str(fc_dict[gene][sample]))
                else:
                    fo.write('\tNA')
            fo.write('\n')


if __name__ == '__main__':
    main()