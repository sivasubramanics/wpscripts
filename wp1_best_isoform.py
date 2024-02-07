#!/usr/bin/env python3
# script to pick the best isoforms from the assembled transcripts fasta file, based on either length or expression
# usage: python wp1_best_isoform.py -i <input_fasta> -o <output> -m <gene_to_trans_map> -t <len/exp>

import argparse
from collections import defaultdict
from typing import NoReturn
import sys


class FASTA:
    """
    Class to store the fasta sequence
    """
    def __init__(self, name: str, sequence: str, description: str = None, readcounts: float = None):
        self.name = name
        self.sequence = sequence
        self.description = description
        self.readcounts = None

    def __str__(self):
        fold_seq = self.fold(60)
        return fold_seq

    def fold(self, width: int) -> str:
        out_str = f">{self.name}"
        if self.description:
            out_str += f" {self.description}"
        if self.readcounts:
            out_str += f" readcounts={self.readcounts}"
        out_str += "\n"
        out_str += '\n'.join([self.sequence[i:i+width] for i in range(0, len(self.sequence), width)])
        return out_str

    def set_readcounts(self, readcounts: float):
        self.readcounts = readcounts

    def __len__(self):
        return len(self.sequence)


def parse_fasta(fasta_file):
    """
    Parse the fasta file
    Parameters
    ----------
    fasta_file

    Returns
    -------
    Yields the FASTA object
    """
    with open(fasta_file, 'r') as f:
        name = ""
        sequence = ""
        begun = False
        for line in f:
            if line.startswith(">"):
                line = line.strip()
                line = line.split(" ")
                if begun:
                    yield FASTA(name, sequence, description)
                name = line[0].replace('>', '')
                if len(line) > 1:
                    description = " ".join(line[1:])
                else:
                    description = None
                sequence = ""
                begun = True
            else:
                sequence += line.strip()
        yield FASTA(name, sequence, description)


def read_gene_to_tr(in_file: str) -> defaultdict:
    """
    Read the gene to transcript map file
    Parameters
    ----------
    in_file

    Returns
    -------
    gene_to_tr: defaultdict
    """
    print(f"Reading gene to transcript map file {in_file}", file=sys.stderr)
    gene_to_tr = defaultdict(list)
    with open(in_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0] not in gene_to_tr:
                gene_to_tr[line[0]] = []
            if line[1] in gene_to_tr[line[0]]:
                print(f"Duplicate transcript {line[1]} for gene {line[0]} in the map file.", file=sys.stderr)
                continue
            gene_to_tr[line[0]].append(line[1])
    return gene_to_tr


def main() -> NoReturn:
    parser = argparse.ArgumentParser(
        description='Pick the best isoforms from the assembled transcripts fasta file, based on either length or '
                    'expression')
    parser.add_argument('-i', '--input', help='input fasta file', required=True)
    parser.add_argument('-o', '--output', help='output prefix file', required=True)
    parser.add_argument('-m', '--map', help='gene to transcript map file', required=True)
    parser.add_argument('-t', '--type', help='type of selection: len or exp', required=True)
    parser.add_argument('-e', '--exp', help='expression count matrix file', required=False)
    args = parser.parse_args()

    # read the gene to transcript map file
    gene_to_tr = read_gene_to_tr(args.map)

    # read the fasta file and read the sequence length
    transcripts = defaultdict(int)
    for fasta in parse_fasta(args.input):
        transcripts[fasta.name] = fasta

    # read expression count matrix file if provided
    if args.exp:
        line_num = 0
        with open(args.exp, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                line_num += 1
                if line_num == 1:
                    continue
                if line[0] not in transcripts:
                    print(f"Transcript {line[0]} not found in the fasta file", file=sys.stderr)
                    continue
                transcripts[line[0]].set_readcounts(sum([float(i) for i in line[1:]]))

    if args.type == 'len':
        # select the longest transcript for each gene
        for gene in gene_to_tr:
            lengths = [len(transcripts[tr]) for tr in gene_to_tr[gene]]
            best_tr = gene_to_tr[gene][lengths.index(max(lengths))]
            out_tr = FASTA(gene, transcripts[best_tr].sequence, transcripts[best_tr].description,
                           transcripts[best_tr].readcounts)
            with open(f"{args.output}.len.fa", 'a') as f:
                f.write(f"{str(out_tr)}\n")
            with open(f"{args.output}.len.map", 'a') as g:
                g.write(f"{gene}\t{best_tr}\n")
        f.close()
        g.close()

    elif args.type == 'exp':
        # select the transcript with the highest expression for each gene
        for gene in gene_to_tr:
            exps = [transcripts[tr].readcounts for tr in gene_to_tr[gene]]
            best_tr = gene_to_tr[gene][exps.index(max(exps))]
            out_tr = FASTA(gene, transcripts[best_tr].sequence, transcripts[best_tr].description,
                           transcripts[best_tr].readcounts)
            with open(f"{args.output}.exp.fa", 'a') as f:
                f.write(f"{str(out_tr)}\n")
            with open(f"{args.output}.exp.map", 'a') as g:
                g.write(f"{gene}\t{best_tr}\n")
        f.close()
        g.close()

    else:
        print(f"Invalid type {args.type}. Please choose len or exp", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
