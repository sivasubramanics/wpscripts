#!/usr/bin/env python3
# script to parse one or more transcriptome fasta files and rename the sequences


import argparse
import os
from collections import defaultdict

class PSL:
    def __init__(self, psl_line):
        psl = psl_line.strip().split("\t")
        self.matches = psl[0]
        self.misMatches = psl[1]
        self.repMatches = psl[2]
        self.nCount = psl[3]
        self.qNumInsert = psl[4]
        self.qBaseInsert = psl[5]
        self.tNumInsert = psl[6]
        self.tBaseInsert = psl[7]
        self.strand = psl[8]
        self.qName = psl[9]
        self.qSize = psl[10]
        self.qStart = psl[11]
        self.qEnd = psl[12]
        self.tName = psl[13]
        self.tSize = psl[14]
        self.tStart = psl[15]
        self.tEnd = psl[16]
        self.blockCount = psl[17]
        self.blockSizes = psl[18]
        self.qStarts = psl[19]
        self.tStarts = psl[20]
        self.qcov = (int(self.qEnd) - int(self.qStart)) / int(self.qSize)
        self.tcov = (int(self.tEnd) - int(self.tStart)) / int(self.tSize)

    def __str__(self):
        return f"{self.qName}\t{self.tName}\t{self.strand}\t{self.qStart}\t{self.qEnd}\t{self.tStart}\t{self.tEnd}\t{self.blockCount}\t{self.blockSizes}\t{self.qStarts}\t{self.tStarts}"

    def __repr__(self):
        return self.__str__()

class Annotation:
    def __init__(self, eggnog_line):
        eggnog = eggnog_line.strip().split("\t")
        self.seqname, self.seed_ortholog, self.evalue, self.score, self.eggNOG_OGs, 
        self.max_annot_lvl, self.cog_category, self.description, self.preferred_name, 
        self.GOs, self.ec, self.kegg_ko, self.kegg_Pathway, self.kegg_Module, self.kegg_Reaction, 
        self.kegg_rclass, self.brite, self.kegg_TC, self.cazy, self.bigg_Reaction, self.PFAMs = eggnog_line.strip().split("\t")

    def __str__(self):
        return f"{self.seqname}\t{self.cog_category}\t{self.description}\t{self.preferred_name}\t{self.GOs}\t{self.PFAMs}"

    def __repr__(self):
        return self.__str__()
    
    def add_psl(self, psl, gname):
        if not hasattr(self, "psl"):
            self.psl = {}
        if not gname in self.psl:
            self.psl[gname] = []
        self.psl[gname].append(psl)


def parse_psl(pslfile, gname, ann):
    with open(pslfile, 'r') as f:
        for line in f:
            psl = PSL(line)
            if psl.tcov < 0.95:
                continue
            if psl.qName in ann:
                ann[psl.qName].add_psl(psl, gname)
            else:
                print(f"Could not find {psl.qName} in annotation file")
                exit(1)
    return ann

def parse_blat(args, ann):
    if not args.genomes is None and not args.psl is None:
        print("Please provide either genomes or psl files. Not both")
        exit(1)
    elif not args.genomes is None:
        if args.names is None:
            print("Please provide names for the genomes")
            exit(1)
        if len(args.genomes) != len(args.names):
            print("Number of genomes and names are not equal. please check")
            exit(1)
        for i,genome in enumerate(args.genomes):
            blat_out = f"{args.ouput}_blat_{genome}.psl"
            blat_cmd = f"pblat -minIdentity=95 -noHead -threads={args.threads} {genome} {args.fasta} {blat_out}"
            print(blat_cmd)
            os.system(blat_cmd)
            ann = parse_psl(blat_out, args.names[i], ann)
    elif not args.psl is None:
        if args.names is None:
            print("Please provide names for the genomes")
            exit(1)
        if len(args.psl) != len(args.names):
            print("Number of genomes and names are not equal. please check")
            exit(1)
        for i,psl in enumerate(args.psl):
            ann = parse_psl(psl, args.names[i], ann)
    else:
        print("Please provide either genomes or psl files")
        exit(1)
    return ann
        

def parse_annotation(ann):
    ann_dict = {}
    with open(ann, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                a = Annotation(line)
                ann_dict[a.seqname] = a
    return ann_dict


def parse_transtogene(trans_to_gene):
    trans_to_gene_map = {}
    with open(trans_to_gene, 'r') as f:
        for line in f:
            line = line.strip().split("\t")
            if line[0] not in trans_to_gene_map:
                trans_to_gene_map[line[0]] = []
            trans_to_gene_map[line[0]].append(line[1])
    return trans_to_gene_map

def parse_count_matrix(matrix, conditions):
    """
    function to parse count matrix and conditions file, and return matrix with mean values per sample
    """




def main():
    parser = argparse.ArgumentParser(description="This script will process the RNAseq data and many more")
    parser.add_argument('-f', '--fasta', help='Input transcriptome fasta file', required=True)
    parser.add_argument('-m', '--trans_to_gene', help='Input transcript to gene map file', required=True)
    parser.add_argument('-a', '--annotation', help='Input annotation file (from EGGNOG output)', required=True)
    parser.add_argument('-x', '--matrix', help='Input gene.count.matrix file ', required=True)
    parser.add_argument('-c', '--conditions', help='Input conditions file (columns should contain rep and name)', required=True)
    parser.add_argument('-g', '--genomes', nargs='+', help='Input genome files (or multiple files)')
    parser.add_argument('-p', '--psl', nargs='+', help='Input psl files (or multiple files)')
    parser.add_argument('-n', '--names', nargs='+', help='Input genome names (or multiple files)'
    parser.add_argument('-t', '--threads', help='Number of threads to use', default=2)
    parser.add_argument('-o', '--output', help='Output prefix', required=True)
    args = parser.parse_args()

    trans_to_gene_map = parse_transtogene(args.trans_to_gene)
    
    ann = parse_annotation(args.annotation)

    ann = parse_blat(args, ann)

    counts = parse_count_matrix(args.matrix, args.conditions)


if __name__ == "__main__":
    main()
