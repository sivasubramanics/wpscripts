#!/usr/bin/env python3
# script to parse the GTF and GFF annotation files and return GTF or GFF object with its attributes

import sys
import os
import argparse
        

def strip_tag(tag):
    tag = tag.strip('"').replace('gene-', '').replace('transcript-', '').replace('protein-', '').replace('exon-', '').replace('CDS-', '').replace('";', '').replace('rna-', '')
    return tag

def get_attributes(attributes, format='gtf'):
    if format == 'gtf':
        attributes = attributes.split('; ')
        attributes = [attribute.split(' ') for attribute in attributes]
        attributes = {attribute[0]: strip_tag(attribute[1]) for attribute in attributes}
    elif format == 'gff':
        attributes = attributes.split(';')
        attributes = [attribute.split('=') for attribute in attributes]
        attributes = {attribute[0]: attribute[1] for attribute in attributes}
    else:
        print('ERROR: Unknown format.')
        sys.exit(1)
    return attributes


def parse_annotation_file(annotation_file, output_file, format='gtf'):
    out_attributes = {}
    with open(annotation_file, 'r') as fh, open(output_file, 'w') as out:
        for line in fh:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if format == 'gtf':
                if line[2] == 'CDS':
                    attributes = get_attributes(line[8], format='gtf')
                    if 'protein_id' not in attributes:
                        attributes['protein_id'] = attributes['transcript_id']
                    if (attributes['gene_id'], attributes['transcript_id'], attributes['protein_id']) in out_attributes:
                        continue
                    out_attributes[(attributes['gene_id'], attributes['transcript_id'], attributes['protein_id'])] = 1
                    out.write(f"{attributes['gene_id']}\t{attributes['transcript_id']}\t{attributes['protein_id']}\n")
            elif format == 'gff':
                attributes = get_attributes(line[8], format='gff')
            else:
                print('ERROR: Unknown format.')
                sys.exit(1)


def parse_gtf(annotation_file, output_file):
    tr2gene = {}
    with open(annotation_file, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            attributes = get_attributes(line[8], format='gtf')
            if 'transcript_id' in attributes:
                if line[2] == 'transcript':
                    tr2gene[attributes['transcript_id']] = {'chrom': line[0], 'start': line[3], 'end': line[4], 'strand': line[6]}
                    if 'protein_id' not in tr2gene[attributes['transcript_id']]:
                        tr2gene[attributes['transcript_id']]['protein_id'] = '-'
                    if 'gene_id' not in tr2gene[attributes['transcript_id']]:
                        tr2gene[attributes['transcript_id']]['gene_id'] = '-'
                if attributes['transcript_id'] not in tr2gene:
                    tr2gene[attributes['transcript_id']] = {'gene_id': '-', 'protein_id': '-'}
                if 'gene_id' in attributes:
                    tr2gene[attributes['transcript_id']]['gene_id'] = attributes['gene_id']
                if 'protein_id' in attributes:
                    tr2gene[attributes['transcript_id']]['protein_id'] = attributes['protein_id']
    return tr2gene


def main():
    parser = argparse.ArgumentParser(description='Parse GTF or GFF annotation files')
    # parser.add_argument('-t', '--task', help='Task to perform', choices=['map_ids', 'transcripts_bed', 'gene_bed'], required=True)
    parser.add_argument('-i', '--input', help='Input GTF or GFF file', required=True)
    parser.add_argument('-f', '--format', help='Input file format (GTF or GFF)', choices=['gtf', 'gff'], required=True)
    parser.add_argument('-o', '--output', help='Output file name', required=True)
    args = parser.parse_args()

    if args.format == 'gtf':
        tr2gene = parse_gtf(args.input, args.output)
        with open(args.output, 'w') as out:
            for tr in tr2gene:
                out.write(f"{tr2gene[tr]['gene_id']}\t{tr}\t{tr2gene[tr]['protein_id']}\t{tr2gene[tr]['chrom']}\t{tr2gene[tr]['start']}\t{tr2gene[tr]['end']}\t{tr2gene[tr]['strand']}\n")
    elif args.format == 'gff':
        parse_annotation_file(args.input, args.output, args.format)
    else:
        print('ERROR: Unknown format.')
        sys.exit(1)
    
if __name__ == '__main__':
    main()

