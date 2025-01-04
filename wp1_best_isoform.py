#!/usr/bin/env python3
# script to pick the best isoforms from the assembled transcripts fasta file, based on either length or expression
# usage: python wp1_best_isoform.py -i <input_fasta> -o <output> -m <gene_to_trans_map> -t <len/exp>

import argparse
import os
from collections import defaultdict
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


def read_gene_to_tr(in_file: str):
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


def main():
    parser = argparse.ArgumentParser(
        description='Pick the best isoforms from the assembled transcripts fasta file, based on either length or '
                    'expression')
    parser.add_argument('-i', '--input', help='input fasta file', required=True)
    parser.add_argument('-o', '--output', help='output prefix file', required=True)
    parser.add_argument('-m', '--map', help='gene to transcript map file', required=True)
    parser.add_argument('-t', '--type', help='type of selection: len or exp or len_pep', required=True)
    parser.add_argument('-e', '--exp', help='expression count matrix file', required=False)
    parser.add_argument('-p', '--pep', help='extract proteins file', required=False)
    args = parser.parse_args()

    # read the gene to transcript map file
    gene_to_tr = read_gene_to_tr(args.map)

    # read the fasta file and read the sequence length
    transcripts = defaultdict()
    for fasta in parse_fasta(args.input):
        transcripts[fasta.name] = fasta

    for gene in gene_to_tr:
        for tr in gene_to_tr[gene]:
            if tr not in transcripts:
                # remove the transcript from the gene to transcript map
                gene_to_tr[gene].remove(tr)
                print(f"Transcript {tr} not found in the fasta file. Removing from the gene to transcript map.",
                      file=sys.stderr)

    if args.type == 'len':
        # remove the output files if they exist
        os.remove(f"{args.output}.len.fa") if os.path.exists(f"{args.output}.len.fa") else None
        os.remove(f"{args.output}.len.map") if os.path.exists(f"{args.output}.len.map") else None

        # output filenames
        tr_fasta_out = f"{args.output}.len.fa"
        tr_map_out = f"{args.output}.len.map"

        # open the output files
        try:
            f = open(tr_fasta_out, 'w')
            g = open(tr_map_out, 'w')
        except IOError:
            print(f"Error opening the output files {tr_fasta_out} or {tr_map_out}", file=sys.stderr)
            sys.exit(1)

        # select the longest transcript for each gene
        for gene in gene_to_tr:
            next_best_tr = None
            lengths = [len(transcripts[tr]) for tr in gene_to_tr[gene]]
            write_isoform(f, g, gene, gene_to_tr, lengths, next_best_tr, transcripts)
        f.close()
        g.close()

    elif args.type == 'exp':
        # read the expression count matrix file
        transcripts = parse_exp_file(args.exp, transcripts)

        # remove the output files if they exist
        os.remove(f"{args.output}.exp.fa") if os.path.exists(f"{args.output}.exp.fa") else None
        os.remove(f"{args.output}.exp.map") if os.path.exists(f"{args.output}.exp.map") else None

        # output filenames
        tr_fasta_out = f"{args.output}.exp.fa"
        tr_map_out = f"{args.output}.exp.map"

        # open the output files
        try:
            f = open(tr_fasta_out, 'w')
            g = open(tr_map_out, 'w')
        except IOError:
            print(f"Error opening the output files {tr_fasta_out} or {tr_map_out}", file=sys.stderr)
            sys.exit(1)

        # select the transcript with the highest expression for each gene
        for gene in gene_to_tr:
            next_best_tr = None
            exps = [transcripts[tr].readcounts for tr in gene_to_tr[gene]]
            write_isoform(f, g, gene, gene_to_tr, exps, next_best_tr, transcripts)
        f.close()
        g.close()

        count_matrix_out = f"{args.output}.exp.count.matrix"
        map_dict = defaultdict(list)
        with open(tr_map_out, 'r') as mf:
            for line in mf:
                line = line.strip().split('\t')
                if line[1] not in map_dict:
                    map_dict[line[1]] = line[0]
                else:
                    print(f"Duplicate transcript {line[1]} in the map file", file=sys.stderr)
                    continue
        with open(count_matrix_out, 'w') as g, open(args.exp, 'r') as f:
            line_num = 0
            for line in f:
                line = line.rstrip()
                line_num += 1
                if line_num == 1:
                    g.write(f"{line}\n")
                    continue
                line = line.split('\t')
                if line[0] in map_dict:
                    line[0] = map_dict[line[0]]
                    out_line = "\t".join(line)
                    g.write(f"{out_line}\n")

    elif args.type == 'len_pep':
        if args.pep is None:
            print("Protein fasta file is required for type len_pep", file=sys.stderr)
            sys.exit(1)

        tr_fasta_out = f"{args.output}.len_orf.fasta"
        tr_map_out = f"{args.output}.len_orf.map"
        pep_fasta_out = f"{args.output}.len_orf.pep.fasta"

        proteins = defaultdict()
        for protein in parse_fasta(args.pep):
            proteins[protein.name] = protein

        # open the output files
        try:
            f = open(tr_fasta_out, 'w')
            g = open(tr_map_out, 'w')
            h = open(pep_fasta_out, 'w')
        except IOError:
            print(f"Error opening the output files {tr_fasta_out} or {tr_map_out} or {pep_fasta_out}", file=sys.stderr)
            sys.exit(1)

        # select the longest protein for each gene
        for gene in gene_to_tr:
            next_best_tr = None
            lengths = [len(proteins[tr]) for tr in gene_to_tr[gene]]
            write_isoform(f, g, gene, gene_to_tr, lengths, next_best_tr, transcripts, h, proteins)

    else:
        print(f"Invalid type {args.type}. Please choose len or exp", file=sys.stderr)
        sys.exit(1)

    if args.pep and args.type != 'len_pep':
        if args.type == 'len':
            map_file = f"{args.output}.len.map"
            pep_file = f"{args.output}.len.pep"
        if args.type == 'exp':
            map_file = f"{args.output}.exp.map"
            pep_file = f"{args.output}.exp.pep"
        map_dict = defaultdict(list)
        with open(map_file, 'r') as mf:
            for line in mf:
                line = line.strip().split('\t')
                if line[1] not in map_dict:
                    map_dict[line[1]] = line[0]
                else:
                    print(f"Duplicate transcript {line[1]} in the map file", file=sys.stderr)
                    continue
        with open(pep_file, 'w') as g:
            for faa in parse_fasta(args.pep):
                if faa.name in map_dict:
                    faa.name = map_dict[faa.name]
                    g.write(f"{str(faa)}\n")


def parse_exp_file(count_matrix_file, transcripts):
    """
    Parse the expression count matrix file
    Parameters
    ----------
    args
    transcripts: defaultdict

    Returns
    -------
    transcripts: defaultdict
    """
    # read the expression count matrix file
    line_num = 0
    with open(count_matrix_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            line_num += 1
            if line_num == 1:
                continue
            if line[0] not in transcripts:
                print(f"Transcript {line[0]} not found in the fasta file", file=sys.stderr)
                continue
            transcripts[line[0]].set_readcounts(sum([float(i) for i in line[1:]]))
    return transcripts


def write_isoform(f, g, gene, gene_to_tr, values, next_best_tr, transcripts, h=None, proteins=None):
    """
    Write the best isoform to the output file
    Parameters
    ----------
    f: file handle for the output fasta file
    g: file handle for the gene to transcript map file
    gene: str
    gene_to_tr: defaultdict
    values: list of lengths or expression values
    next_best_tr: str
    transcripts: dictionary of isoforms FASTA objects

    Returns
    -------
    None
    """
    best_tr = gene_to_tr[gene][values.index(max(values))]
    if len(gene_to_tr[gene]) > 1:
        next_best_tr = gene_to_tr[gene][values.index(sorted(values)[-2])]
    out_tr = FASTA(gene, transcripts[best_tr].sequence,
                   f"old_name={best_tr} next_best={next_best_tr} {transcripts[best_tr].description}",
                   transcripts[best_tr].readcounts)

    f.write(f"{str(out_tr)}\n")
    g.write(f"{gene}\t{best_tr}\n")

    if h is not None:
        out_prot = FASTA(gene, proteins[best_tr].sequence, f"old_name={best_tr} next_best={next_best_tr}")
        h.write(f"{str(out_prot)}\n")


if __name__ == "__main__":
    main()
