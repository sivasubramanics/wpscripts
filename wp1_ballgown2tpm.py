#!/usr/bin/env python3

"""
This script will process the stringtie ballgown GTF files and output coverage, TPM and FPKM values in different tsv files
"""

import os
import sys
import argparse
from collections import defaultdict

def get_attributes(attr_str):
    """
    Get the attributes from the GTF line
    :param attr_str: attribute string
    :return: dictionary of attributes
    """
    attr_list = attr_str.split(';')
    attr_dict = {}
    for attr in attr_list:
        if attr:
            attr_list = attr.split()
            key = attr_list[0]
            value = ' '.join(attr_list[1:])
            attr_dict[key] = value.strip('"')
    return attr_dict

def main():
    parser = argparse.ArgumentParser(description='Convert Ballgown GTF to TPM')
    parser.add_argument('-g', '--gtf', help='Stringtie output GTF files (multiple files can be provided)', nargs='+', required=True)
    parser.add_argument('-o', '--output', help='Output prefix', required=True)
    args = parser.parse_args()

    fpkm = defaultdict(dict)
    tpm = defaultdict(dict)
    cov = defaultdict(dict)
    gene = defaultdict(dict)

    samples = []
    # read the GTF file
    for gtf_file in args.gtf:
        sample_name = os.path.basename(gtf_file).replace('.gtf', '')
        samples.append(sample_name)
        print(f'Processing {sample_name}')
        with open(gtf_file, 'r') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if fields[2] != 'transcript':
                    continue
                attributes = get_attributes(fields[8])
                tr_id = attributes['transcript_id']
                if tr_id not in fpkm:
                    fpkm[tr_id] = {}
                    tpm[tr_id] = {}
                    cov[tr_id] = {}
                if sample_name not in fpkm[tr_id]:
                    fpkm[tr_id][sample_name] = 0
                    tpm[tr_id][sample_name] = 0
                    cov[tr_id][sample_name] = 0
                if tr_id not in gene:
                    gene[tr_id] = attributes['gene_id']
                fpkm[tr_id][sample_name] = attributes['FPKM']
                tpm[tr_id][sample_name] = attributes['TPM']
                cov[tr_id][sample_name] = attributes['cov']


    # write the TPM, FPKM and coverage values to files
    with open(args.output + '.tpm.tsv', 'w') as tpm_fh, open(args.output + '.fpkm.tsv', 'w') as fpkm_fh, open(args.output + '.cov.tsv', 'w') as cov_fh:
        tpm_fh.write('gene_id\ttranscript_id\t' + '\t'.join(samples) + '\n')
        fpkm_fh.write('gene_id\ttranscript_id\t' + '\t'.join(samples) + '\n')
        cov_fh.write('gene_id\ttranscript_id\t' + '\t'.join(samples) + '\n')
        for tr_id in tpm:
            tpm_fh.write(gene[tr_id] + '\t' + tr_id + '\t' + '\t'.join([str(tpm[tr_id][sample]) for sample in samples]) + '\n')
            fpkm_fh.write(gene[tr_id] + '\t' + tr_id + '\t' + '\t'.join([str(fpkm[tr_id][sample]) for sample in samples]) + '\n')
            cov_fh.write(gene[tr_id] + '\t' + tr_id + '\t' + '\t'.join([str(cov[tr_id][sample]) for sample in samples]) + '\n')

    print(f'Output written to {args.output}.tpm, {args.output}.fpkm and {args.output}.cov')

if __name__ == '__main__':
    main()