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
    def __init__(self, qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore,
                 qcovhsp, scovhsp, qlen, slen):
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
        if float(self.pident) >= min_pident and float(self.scovhsp) >= min_scovhsp:
            return True
        return False


def fold_seq(seq, n=60):
    """
    Fold the sequence into n characters
    """
    return '\n'.join([seq[i:i + n] for i in range(0, len(seq), n)])


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
                seq += line.strip().replace("*", "")
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
        if '=' not in part:
            # merge the last part with the current part
            parsed_data[list(parsed_data.keys())[-1]] += ' ' + part
            continue
        key, value = part.split('=')
        parsed_data[key] = value
    return parsed_data


def main():
    parser = argparse.ArgumentParser(description='Filter isoforms based on the uniref annotation')
    parser.add_argument('-i', '--input', help='Input isoforms fasta file', required=True)
    parser.add_argument('-a', '--proteins', help='Proteins fasta file', required=True)
    parser.add_argument('-o', '--output', help='Output prefix', required=True)
    parser.add_argument('-m', '--gene_to_tr_map', help='Gene to transcript map file', required=True)
    parser.add_argument('-d', '--diamond', help='Diamond output file', required=True)
    parser.add_argument('-u', '--uniref', help='Uniref fasta file', required=True)
    parser.add_argument('-p', '--percent', help='Minimum percent identity', required=False, default=60)
    parser.add_argument('-c', '--coverage', help='Minimum coverage', required=False, default=80)
    args = parser.parse_args()

    # Set up the logging in the format to print the time: [2024-02-02 12:12:12] message
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

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

    logging.info(f'Reading the proteins fasta file {args.proteins}')
    proteins = defaultdict()
    for fasta in parse_fasta(args.proteins):
        proteins[fasta.id] = fasta

    logging.info(f'Reading the isoforms fasta file {args.input}')
    isoforms = defaultdict()
    for fasta in parse_fasta(args.input):
        isoforms[fasta.id] = fasta

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
                logging.error(
                    f"Duplicate alignment found for {aln.qseqid}. Make sure you provide only the best hits in the diamond output file")
            aln_dict[aln.qseqid].append(aln)
            aln_dict[aln.qseqid].append(uniref[aln.sseqid])

    logging.info(f'Filtering the alignments')
    with open(args.output + ".classified.tsv", 'w') as raw, open(args.output + ".classified.filtered.tsv", 'w') as flt:
        for gene_id in gene_to_tr:
            scores = defaultdict()
            for tr_id in gene_to_tr[gene_id]:
                if tr_id not in proteins:
                    raw.write(f"{gene_id}\t{tr_id}\tNoProtein\n")
                    continue
                if tr_id not in aln_dict:
                    raw.write(f"{gene_id}\t{tr_id}\tUnmapped\n")
                    continue
                aln = aln_dict[tr_id][0]
                ann = aln_dict[tr_id][1]
                scores[tr_id] = aln.bitscore
                raw.write(
                    f"{gene_id}\t{tr_id}\t{aln.sseqid}\t{aln.pident}\t{aln.qcovhsp}\t{aln.scovhsp}\t{aln.bitscore}\t{ann['Tax']}\t{ann['desc']}\t{len(isoforms[tr_id].seq)}")
                if aln.is_complete(args.percent, args.coverage):
                    raw.write("\tComplete\n")
                else:
                    raw.write("\tIncomplete\n")
            if scores:
                tr_id = max(scores, key=scores.get)
                # check if other transcripts have the same score
                if list(scores.values()).count(scores[tr_id]) > 1:
                    # get the list of transcripts with the same score
                    tr_ids = [tr for tr in scores if scores[tr] == scores[tr_id]]
                    # get the longest transcript
                    tr_id = max(tr_ids, key=lambda x: len(isoforms[x].seq))
                aln = aln_dict[tr_id][0]
                ann = aln_dict[tr_id][1]
                flt.write(
                    f"{gene_id}\t{tr_id}\t{aln.sseqid}\t{aln.pident}\t{aln.qcovhsp}\t{aln.scovhsp}\t{aln.bitscore}\t{ann['Tax']}\t{ann['desc']}\t{len(isoforms[tr_id].seq)}")
                if aln.is_complete(args.percent, args.coverage):
                    flt.write("\tComplete\n")
                else:
                    flt.write("\tIncomplete\n")

    logging.info(f'Writing the filtered isoforms')
    with (open(f"{args.output}.fasta", 'w') as tr_fa,
          open(f"{args.output}.gene_to_tr_map", 'w') as map,
          open(f"{args.output}.faa", 'w') as faa,
            open(f"{args.output}.classified.filtered.tsv", 'r') as f):
        for line in f:
            gene_id, tr_id, uniref_id, pident, qcovhsp, scovhsp, bitscore, tax, desc, length, status = line.strip().split('\t')
            tr_fa.write(f">{gene_id}\n{fold_seq(isoforms[tr_id].seq)}\n")
            faa.write(f">{gene_id}\n{fold_seq(proteins[tr_id].seq)}\n")
            map.write(f"{gene_id}\t{gene_id}\n")


if __name__ == "__main__":
    main()
