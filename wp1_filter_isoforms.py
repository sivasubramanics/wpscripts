#!/usr/bin/env python3

"""
This script will help in filtering the isoforms based on the uniref annotation
"""

import logging
import os
import sys
from collections import defaultdict
import argparse


class Fasta:
    def __init__(self, id, seq, desc: str = None):
        self.id = id
        self.seq = seq
        self.desc = desc


class Aln:
    def __init__(self, qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qcovhsp, scovhsp, qlen, slen):
        self.qseqid = qseqid
        self.sseqid = sseqid
        self.pident = pident
        self.length = length
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.evalue = evalue
        self.bitscore = bitscore
        self.qcovhsp = qcovhsp
        self.scovhsp = scovhsp
        self.qlen = qlen
        self.slen = slen

    def is_complete(self, min_pident, min_scovhsp):
        if self.pident >= min_pident and self.scovhsp >= min_scovhsp:
            return True
        return False


def parse_fasta(in_fa_file):
    """
    Parse the fasta file
    """
    seq = ""
    id = ""
    desc = ""
    with open(in_fa_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if seq:
                    yield Fasta(id, seq, desc)
                    seq = ""
                id = line.strip().split()[0][1:]
                desc = " ".join(line.strip().split()[1:])
            else:
                seq += line.strip()
        if seq:
            yield Fasta(id, seq, desc)


def parse_uniref_desc(desc):
    """
    Parse the uniref description "titin isoform X1 n=7 Tax=Crotalus tigris TaxID=88082 RepID=UPI00192FAB0"
    return
    desc = titin isoform X1
    n= 7
    Tax = Crotalus tigris
    TaxID = 88082
    RepID = UPI00192FAB0
    """
    parts = desc.split()
    parsed_data = {}

    desc_parts = []
    for part in parts:
        if part.startswith("n="):
            break
        desc_parts.append(part)
    parsed_data['desc'] = ' '.join(desc_parts)
    # Iterate over the remaining parts and parse them
    for part in parts[len(desc_parts):]:
        key, value = part.split('=')
        parsed_data[key] = value
    return parsed_data


def main():
    parser = argparse.ArgumentParser(description='Filter isoforms based on the uniref annotation')
    parser.add_argument('-i', '--input', help='Input isoforms fasta file', required=True)
    parser.add_argument('-o', '--output', help='Output prefix', required=True)
    parser.add_argument('-m', '--gene_to_tr_map', help='Gene to transcript map file', required=True)
    parser.add_argument('-d', '--diamond', help='Diamond output file', required=True)
    parser.add_argument('-u', '--uniref', help='Uniref fasta file', required=True)
    parser.add_argument('-p', '--percent', help='Minimum percent identity', required=False, default=60)
    parser.add_argument('-c', '--coverage', help='Minimum coverage', required=False, default=80)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

    logging.info(f'Reading the uniref90 fasta file {args.uniref}')
    uniref = defaultdict()
    for fasta in parse_fasta(args.uniref):
        uniref[fasta.id] = parse_uniref_desc(fasta.desc)

    logging.info(f'Reading the gene to transcript map file {args.gene_to_tr_map}')
    gene_to_tr = defaultdict()
    tr_to_gene = defaultdict()
    with open(args.gene_to_tr_map, 'r') as f:
        for line in f:
            gene, tr = line.strip().split()
            tr_to_gene[tr] = gene
            if gene not in gene_to_tr:
                gene_to_tr[gene] = []
            gene_to_tr[gene].append(tr)

    logging.info(f'Reading the diamond output file {args.diamond}')
    aln_dict = defaultdict()
    with open(args.diamond, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            aln = Aln(*line.strip().split())
            if aln.qseqid not in aln_dict:
                aln_dict[aln.qseqid] = []
            else:
                logging.error(f"Duplicate alignment found for {aln.qseqid}. Make sure you provide only the best hits in the diamond output file")
            aln_dict[aln.qseqid].append(aln)
            aln_dict[aln.qseqid].append(uniref[aln.sseqid])

    logging.info(f'Filtering the alignments')
    with open(args.output + ".classified.tsv", 'w') as f:
        for gene_id in gene_to_tr:
            for tr_id in gene_to_tr[gene_id]:
                if tr_id not in aln_dict:
                    f.write(f"{gene_id}\t{tr_id}\tUnmapped\n")
                    continue
                aln = aln_dict[tr_id][0]
                ann = aln_dict[tr_id][1]
                f.write(
                    f"{gene_id}\t{tr_id}\t{aln.sseqid}\t{aln.pident}\t{aln.qcovhsp}\t{aln.scovhsp}\t{ann['Tax']}\t{ann['desc']}")
                if aln.is_complete(args.percent, args.coverage):
                    f.write("\tComplete\n")
                else:
                    f.write("\tIncomplete\n")


if __name__ == "__main__":
    main()
