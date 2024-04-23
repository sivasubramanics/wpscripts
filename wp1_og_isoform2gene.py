#!/usr/bin/env python3

"""
We run orthofinder or pantools for the proteins from the isoform, hence we may expect multiple isoforms for a gene.
This script will try to compile the gene to isoform mapping.
"""

import argparse
from collections import defaultdict
import logging


def get_gene2isoform_map(gene2tr_map_file):
    """
    Get gene to isoform map
    """
    gene2isoform = defaultdict(list)
    with open(gene2tr_map_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            gene2isoform[line[0]].append(line[1])
    return gene2isoform


def get_isoform2gene_map(gene2tr_map_file):
    """
    Get isoform to gene map
    """
    isoform2gene = defaultdict(str)
    logging.info(f"Reading gene to isoform map from {gene2tr_map_file}")
    with open(gene2tr_map_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            isoform2gene[line[1]] = line[0]
    return isoform2gene


def convert_og_table(args):
    """
    Convert orthogroup table from isoform to gene
    """

    if len(args.map) != len(args.name):
        raise ValueError('Number of map files and names should be same')

    tr2gene = defaultdict(dict)
    for i in range(len(args.map)):
        tr2gene[args.name[i]] = get_isoform2gene_map(args.map[i])

    with open(args.input, 'r') as f, open(args.output, 'w') as o:
        line_num = 0
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                o.write(line + '\n')
                continue
            line = line.split('\t')
            line_num += 1
            if line[0] == 'Orthogroup':
                names = line[1:]
                for i in range(len(names)):
                    if names[i] not in tr2gene:
                        raise ValueError(f"Name {names[i]} not found in the map files")
                o.write('\t'.join(line) + '\n')
                continue
            if line_num == 1:
                raise ValueError('First line should be Orthogroup. check the input file is in orthotable format')
            og_id = line[0]
            o.write(f"{og_id}")
            for i in range(1, len(line)):
                if line[i] == '':
                    continue
                trs = line[i].split(',')
                gene_names = set()
                for tr in trs:
                    if tr in tr2gene[names[i-1]]:
                        gene_names.add(tr2gene[names[i-1]][tr])
                o.write(f"\t{','.join(gene_names)}")
            o.write('\n')


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    parser = argparse.ArgumentParser(description='convert isoform data to gene data')
    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    og_table_parser = subparsers.add_parser('og_table', help='convert orthogroup table from isoform to gene')
    og_table_parser.add_argument('-i', '--input', help='input orthogroup table. check wp1_comparative_transcriptomics.py for conversion', required=True)
    og_table_parser.add_argument('-m', '--map', help='map file with isoform to gene mapping', required=True, nargs='+')
    og_table_parser.add_argument('-n', '--name', help='name of the data (must be in the same order as map file atguments)', required=True, nargs='+')
    og_table_parser.add_argument('-o', '--output', help='output orthogroup table with gene names', required=True)

    args = parser.parse_args()

    if args.command == 'og_table':
        convert_og_table(args)




if __name__ == "__main__":
    main()