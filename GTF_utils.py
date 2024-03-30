#!/usr/bin/env python3

import sys, os
import argparse
from collections import defaultdict
import logging

class GTFline:
    """
    Class object to hold GTF line
    """
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

def read_gtf(gtf_file):
    """
    Read a GTF file and yield GTFline objects
    :param gtf_file: GTF file
    :return: GTFline object
    """
    with open(gtf_file, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            yield GTFline(line)


def get_gene_to_tr_map(input, output):
    """
    Generate gene to transcript mapping
    :param input: input GTF file
    :param output: output table
    :return: None
    """
    # if output is not specified, create one
    if not output:
        output = input + '.gene_to_transcript_map'

    # create a dictionary to hold the gene to transcript mapping
    gene_to_tr = defaultdict(list)

    # read the GTF file and generate the gene to transcript mapping
    logging.info(f"reading GTF file: {input}")
    for record in read_gtf(input):
        if record.feature == 'transcript':
            gene_id = record.attribute_dict['gene_id']
            transcript_id = record.attribute_dict['transcript_id']
            gene_to_tr[gene_id].append(transcript_id)

    # write the gene to transcript mapping to the output table
    logging.info(f"writing gene to transcript mapping to: {output}")
    with open(output, 'w') as out_fh:
        for gene_id in gene_to_tr:
            for trans_id in gene_to_tr[gene_id]:
                out_fh.write(f"{gene_id}\t{trans_id}\n")




def main():
    # add command line arguments with subcommands
    parser = argparse.ArgumentParser(description='This is a utility script to process GTF files, and perform simple ooperations')
    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    # add subcommand for generating gene density per MB table
    density_parser = subparsers.add_parser('density', help='generate gene/transcript/exon density per MB table', description='generate gene/transcript/exon density per MB table')
    density_parser.add_argument('-i', '--input', help='input GTF file', required=True)
    density_parser.add_argument('-o', '--output', help='output table, if not specified, input.gtf.1mb.feature.density', required=False)
    density_parser.add_argument('-t', '--type', help='feature type (gene, transcript, exon)', required=True, type=str, choices=['gene', 'transcript', 'exon'])
    density_parser.add_argument('-m', '--mb', help='window size in MB', required=False, type=float, default=1)

    gene_to_tr_map_parser = subparsers.add_parser('gene_to_tr_map', help='generate gene to transcript mapping', description='generate gene to transcript mapping')
    gene_to_tr_map_parser.add_argument('-i', '--input', help='input GTF file', required=True)
    gene_to_tr_map_parser.add_argument('-o', '--output', help='output table, if not specified, input.gtf.gene_to_transcript_map', required=False)

    # parse the arguments and save to args
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    # set the logging level to INFO in format 'YYYY-MM-DD HH:MM:SS - LEVEL - MESSAGE'
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    # log the command line arguments provided, with the message "cmdline: ", and print on the stdout
    logging.info('cmdline: %s', ' '.join(sys.argv))




    # if subcommand is density, call the density function
    if args.command == 'density':
        get_density(args.input, args.output, args.type, args.mb)
    if args.command == 'gene_to_tr_map':
        get_gene_to_tr_map(args.input, args.output)

def get_density(in_gtf, out_table, feature, mb):
    """
    Generate gene/transcript/exon density per MB table
    :param in_gtf: input GTF file
    :param out_table: output table
    :param feature: feature type (gene, transcript, exon)
    :param mb: window size in MB
    :return: None
    """
    # if out_table is not specified, create one
    if not out_table:
        out_table = in_gtf + '.' + str(mb) + 'mb.' + feature + '.density'

    # create a dictionary to hold the density
    density = defaultdict(dict)

    # read the GTF file and calculate the density
    logging.info(f"reading GTF file: {in_gtf}")
    is_feature = False
    for record in read_gtf(in_gtf):
        if record.feature == feature:
            is_feature = True
            if record.seqname not in density:
                density[record.seqname] = defaultdict(dict)
            # get the window
            window = int(record.start / (mb * 1000000))
            if window not in density[record.seqname]:
                density[record.seqname][window] = 0
            density[record.seqname][window] += 1

    # if the feature is not found in the GTF file, exit with an error message
    if not is_feature:
        logging.error(f"feature: {feature} not found in the GTF file")
        sys.exit(1)

    # write the density to the output table
    logging.info(f"writing density to: {out_table}")
    with open(out_table, 'w') as out_fh:
        out_fh.write('seqname\twindow\t' + feature + '_density\n')
        for seqname in density:
            for window in density[seqname]:
                out_fh.write(seqname + '\t' + str(window * mb) + '\t' + str(density[seqname][window]) + '\n')





if __name__ == '__main__':
    main()





