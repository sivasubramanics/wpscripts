#!/usr/bin/env python3

"""
This script reads through all the t_data.ctab files in the given directory and checks if cov is more than the cutoff and flags the transcript as expressed or not expressed
Then reads through the gtf file, and removed the transcripts that are not expressed
"""

import os
import sys
import argparse
from collections import defaultdict
import logging

class Feature:
    def __init__(self, line):
        self.line = line.strip().split('\t')
        self.seqname = self.line[0]
        self.source = self.line[1]
        self.feature = self.line[2]
        self.start = int(self.line[3])
        self.end = int(self.line[4])
        self.score = self.line[5]
        self.strand = self.line[6]
        self.frame = self.line[7]
        self.attribute = self.line[8]
        self.attribute_dict = {}
        for attr in self.attribute.split(';'):
            if attr:
                attr_list = attr.split()
                key = attr_list[0]
                value = ' '.join(attr_list[1:])
                self.attribute_dict[key] = value.strip('"')
        self.length = self.end - self.start + 1

    def __str__(self):
        return '\t'.join(self.line)

class TR:
    def __init__(self, tr_id, feature):
        self.tr_id = tr_id
        self.feature = feature
        self.parent = None
        self.children = []

    def __str__(self):
        output = [f'{self.feature}']
        output.extend([f'{child}' for child in self.children])
        return '\n'.join(output)


def main():
    parser = argparse.ArgumentParser(description='Filter Ballgown GTF file')
    parser.add_argument('-g', '--gtf', help='GTF file', required=True)
    parser.add_argument('-c', '--cutoff', type=float, help='Cutoff value for expression', default=1)
    parser.add_argument('-d', '--dir', help='Directory containing t_data.ctab files', required=True)
    parser.add_argument('-o', '--output', help='Output prefix', required=True)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    # get the expressed transcripts
    expressed_tr = set()
    tr_per_cov = defaultdict(dict)
    tr_per_fpkm = defaultdict(dict)
    # read through all the t_data.ctab files in the directory, recursively
    for root, dirs, files in os.walk(args.dir):
        for file in files:
            if file.endswith('t_data.ctab'):
                # get the sample name from the directory name of the file
                sample_name = os.path.basename(root)
                logging.info(f'Processing {sample_name} file {file}')
                with open(os.path.join(root, file), 'r') as fh:
                    for line in fh:
                        fields = line.strip().split('\t')
                        if fields[0] == 't_id':
                            continue
                        tr_id = fields[5]
                        if tr_id not in tr_per_cov:
                            tr_per_cov[tr_id] = defaultdict(float)
                            tr_per_fpkm[tr_id] = defaultdict(float)
                        if not sample_name in tr_per_cov[tr_id]:
                            tr_per_cov[tr_id][sample_name] = float(fields[10])
                            tr_per_fpkm[tr_id][sample_name] = float(fields[11])
                        else:
                            logging.warning(f'Duplicate entry for {tr_id} in {sample_name}')
    logging.info(f'Found {len(tr_per_cov)} transcripts')

    # check if the transcript is expressed in any of the samples
    logging.info(f'Filtering transcripts with coverage >= {args.cutoff}')
    for tr_id in tr_per_cov:
        count = 0
        for sample in tr_per_cov[tr_id]:
            if tr_per_cov[tr_id][sample] > args.cutoff:
                count += 1
        if count/len(tr_per_cov[tr_id]) >= 0.4:
            expressed_tr.add(tr_id)
    logging.info(f'Found {len(expressed_tr)} expressed transcripts')

    # read the GTF file and filter the transcripts
    gtf = defaultdict(dict)
    genes = defaultdict(dict)
    logging.info(f'Reading GTF file {args.gtf}')
    with open(args.gtf, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            feature = Feature(line)
            if feature.feature == 'gene':
                gene_id = feature.attribute_dict['gene_id']
                genes[gene_id] = feature
            if feature.feature == 'transcript':
                tr_id = feature.attribute_dict['transcript_id']
                gene_id = feature.attribute_dict['gene_id']
                gtf[tr_id] = TR(tr_id, feature)
            elif feature.feature == 'exon':
                tr_id = feature.attribute_dict['transcript_id']
                if not tr_id in gtf:
                    logging.error(f'Transcript {tr_id} not found in GTF file. File is not sorted? Please use agat to fix the gtf')
                    sys.exit(1)
                gtf[tr_id].children.append(feature)

    if len(genes) == 0:
        logging.warning('No genes found in GTF file')
    else:
        logging.info(f'{len(genes)} genes found in GTF file')
        for tr_id in gtf:
            gene_id = gtf[tr_id].feature.attribute_dict['gene_id']
            if gene_id not in genes:
                logging.error(f'Gene {gene_id} not found in GTF file')
                sys.exit(1)
            gtf[tr_id].parent = genes[gene_id]

    # write the filtered GTF file
    logging.info(f'Writing filtered GTF file {args.output}.gtf')
    with open(f'{args.output}.gtf', 'w') as fh:
        for tr_id in gtf:
            if tr_id in expressed_tr:
                fh.write(f'{gtf[tr_id]}\n')

    # write the cov table
    logging.info(f'Writing coverage table {args.output}.cov')
    with open(f'{args.output}.cov', 'w') as fh:
        fh.write('tr_id\t' + '\t'.join(sorted(tr_per_cov[list(tr_per_cov.keys())[0]])) + '\n')
        for tr_id in tr_per_cov:
            if tr_id in expressed_tr:
                fh.write(tr_id + '\t' + '\t'.join([str(tr_per_cov[tr_id][sample]) for sample in sorted(tr_per_cov[tr_id])]) + '\n')

    # write the fpkm table
    logging.info(f'Writing FPKM table {args.output}.fpkm')
    with open(f'{args.output}.fpkm', 'w') as fh:
        fh.write('tr_id\t' + '\t'.join(sorted(tr_per_fpkm[list(tr_per_fpkm.keys())[0]])) + '\n')
        for tr_id in tr_per_fpkm:
            if tr_id in expressed_tr:
                fh.write(tr_id + '\t' + '\t'.join([str(tr_per_fpkm[tr_id][sample]) for sample in sorted(tr_per_fpkm[tr_id])]) + '\n')



if __name__ == '__main__':
    main()
