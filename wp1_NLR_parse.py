#!/usr/bin/env python3

"""
This script will parse the NLR genes from the WP1 NLR annotated output
"""

import argparse
from collections import defaultdict

RANKS = ['CC-TIR-NBARC-LRR',
         'TIR-NBARC-LRR', 'CC-NBARC-LRR',
         'TIR-NBARC', 'CC-NBARC',
         'TIR-LRR', 'CC-LRR',
         'TIR', 'CC',
         'NBARC-LRR',
         'NBARC']


def get_best_annotation(annotations):
    """
    Get the best annotation from the list of annotations
    """
    for rank in RANKS:
        if rank in annotations:
            return rank
    return 'NA'


def main():
    parser = argparse.ArgumentParser(description='Parse NLR annotater output and list NLR genes')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-m', '--map', help='gene to transcript map file', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)

    args = parser.parse_args()

    gene_transcript_map = defaultdict(list)
    with open(args.map, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0] not in gene_transcript_map:
                gene_transcript_map[line[0]] = []
            gene_transcript_map[line[0]].append(line[1])

    nlr_dict = defaultdict(dict)
    nlr_annotation = []
    last_transcript = ''
    with open(args.input, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            transcript = line[0]
            annotation = line[2]
            if transcript == last_transcript:
                nlr_annotation.append(annotation)
            else:
                if nlr_annotation:
                    nlr_dict[last_transcript] = get_best_annotation(nlr_annotation)
                    nlr_annotation = []
                nlr_annotation = [annotation]
            last_transcript = transcript
    if nlr_annotation:
        nlr_dict[last_transcript] = get_best_annotation(nlr_annotation)

    with open(args.output, 'w') as o:
        o.write("Geneid\tNLR\n")
        for gene in gene_transcript_map:
            annotations = []
            for transcript in gene_transcript_map[gene]:
                if transcript in nlr_dict:
                    annotations.append(nlr_dict[transcript])
            if annotations:
                o.write(f"{gene}\t{get_best_annotation(annotations)}\n")


if __name__ == '__main__':
    main()
