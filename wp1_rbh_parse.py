#!/usr/bin/env python3

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Parse mmseqs2 rbh output and convert to gene IDs based RBH')
    parser.add_argument('-i', '--input', help='Input file (mmseq2 easy-rbh with output format set as 2)', required=True)
    parser.add_argument('-a', '--annotation', help='Annotation file (genomic.tsv) columns: 1. gene_id 2. transcript_id, 3. protein_id', required=True, nargs='+')
    parser.add_argument('-o', '--output', help='Output file name', required=True)
    args = parser.parse_args()

    pr2gene = {}
    tr2gene = {}
    for annotation_file in args.annotation:
        with open(annotation_file, 'r') as fh:
            for line in fh:
                line = line.strip().split('\t')
                if line[0] == 'gene_id':
                    continue
                if line[2] == 'NA':
                    continue
                pr2gene[line[2]] = line[0]
                tr2gene[line[1]] = line[0]
    
    with open(args.input, 'r') as fh, open(args.output, 'w') as out:
        for line in fh:
            line = line.strip().split('\t')
            if line[0] == line[1]:
                continue
            if line[0] in pr2gene and line[1] in pr2gene:
                out.write(f"{pr2gene[line[0]]}\t{pr2gene[line[1]]}\n")
            elif line[0] in tr2gene and line[1] in tr2gene:
                out.write(f"{tr2gene[line[0]]}\t{tr2gene[line[1]]}\n")
            else:
                print(f"ERROR: {line[0]} or {line[1]} not found in annotation file.")
                # sys.exit(1)



if __name__ == '__main__':
    main()
