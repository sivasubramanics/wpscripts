#!/usr/bin/env python3

"""
This script will process the genome files and prepare the genome for further analysis
"""

import argparse
import sys
from collections import defaultdict
import logging
import os
import subprocess


def run_cmd(cmd) -> list:
    logging.info(f"Running command: {cmd}")
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        logging.error(f"Command failed: {cmd}")
        logging.error(f"Error: {stderr.decode('utf-8')}")
        sys.exit(1)
    return [stdout.decode('utf-8'), stderr.decode('utf-8')]


def main():
    parser = argparse.ArgumentParser(
        description='This is a utility script to process genome files and prepare the genome for further analysis')
    parser.add_argument('-f', '--fasta', help='input genome fasta file', required=True)
    parser.add_argument('-d', '--dict', help='input dictionary tsv file with <old_name>\t<new_name>', required=True)
    parser.add_argument('-g', '--gff', help='input gff file', required=True)
    parser.add_argument('-o', '--outdir', help='output directory', required=True)
    parser.add_argument('-r', '--rna', help='input rna fasta file', required=False)
    parser.add_argument('-p', '--protein', help='input protein fasta file', required=False)
    parser.add_argument('--force', help='force overwrite existing files', action='store_true')
    # if no argument is provided, print help
    args = parser.parse_args(sys.argv[1:] if len(sys.argv) > 1 else ['-h'])

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

    TAB = '\t'


    # read dictionary file
    name_dict = defaultdict()
    with open(args.dict, 'r') as f:
        for line in f:
            line = line.strip().split(TAB)
            if len(line) != 2:
                logging.error(f"Dictionary file must have two columns: {args.dict}")
                sys.exit(1)
            name_dict[line[0]] = line[1]

    # remove output directory if force is enabled
    if args.force:
        logging.info("Force overwrite enabled")
        os.remove(args.outdir)

    # create output directory
    if not os.path.exists(args.outdir):
        logging.info(f"Creating output directory: {args.outdir}")
        os.makedirs(args.outdir)

    # read genome fasta file
    genome_fasta = args.outdir + '/genome.fa'
    genome_gff = args.outdir + '/genomic.gff'
    genome_gtf = args.outdir + '/genomic.gtf'
    genome_rna = args.outdir + '/rna.fa'
    genome_protein = args.outdir + '/protein.fa'
    genome_tsv = args.outdir + '/genomic.tsv'

    if not os.path.exists(f"{genome_fasta}.ok") or not os.path.exists(f"{genome_fasta}"):
        logging.info(f"Creating genome fasta file: {genome_fasta}")
        with open(args.fasta, 'r') as f, open(genome_fasta, 'w') as out:
            for line in f:
                if line.startswith('>'):
                    line = line.strip().split(' ')
                    name = line[0][1:]
                    if name in name_dict:
                        name = name_dict[name]
                    out.write(f">{name}\n")
                else:
                    out.write(line)
        open(f"{genome_fasta}.ok", 'w').close()

    if not os.path.exists(f"{genome_gff}.ok") or not os.path.exists(f"{genome_gff}"):
        logging.info(f"Creating genome gff file: {genome_gff}")
        with open(args.gff, 'r') as f, open(genome_gff, 'w') as out:
            for line in f:
                if line.startswith('#'):
                    out.write(line)
                else:
                    line = line.strip().split(TAB)
                    if line[0] in name_dict:
                        line[0] = name_dict[line[0]]
                    out.write(f"{TAB.join(line)}\n")
        open(f"{genome_gff}.ok", 'w').close()

    if not os.path.exists(f"{genome_gtf}.ok") or not os.path.exists(f"{genome_gtf}"):
        logging.info(f"Converting gff to gtf using agat : {genome_gtf}")
        cmd = f"agat_convert_sp_gff2gtf.pl --gff {genome_gff} -o {genome_gtf} > {genome_gtf}.log 2>&1"
        run_cmd(cmd)
        open(f"{genome_gtf}.ok", 'w').close()

    if not os.path.exists(f"{genome_tsv}.ok") or not os.path.exists(f"{genome_tsv}"):
        logging.info(f"Converting gff to tsv using agat : {genome_tsv}")
        cmd = f"agat_convert_sp_gff2tsv.pl --gff {genome_gff} -o {genome_tsv} > {genome_tsv}.log 2>&1"
        run_cmd(cmd)
        open(f"{genome_tsv}.ok", 'w').close()

    logging.info(f"Reading genome tsv file: {genome_tsv}")
    gene_ids = defaultdict()
    rna_ids = defaultdict()
    protein_ids = defaultdict()
    with open(genome_tsv, 'r') as f:
        for line in f:
            line = line.strip().split(TAB)
            if line[0] == 'seq_id':
                id_idx = line.index('ID')
                continue
            if line[2] == 'gene':
                line[id_idx] = line[id_idx].replace('gene-', '')
                if not line[id_idx] in gene_ids:
                    gene_ids[line[id_idx]] = 1
                continue
            if line[2] == 'mRNA' or line[2] == 'RNA':
                line[id_idx] = line[id_idx].replace('rna-', '')
                if not line[id_idx] in rna_ids:
                    rna_ids[line[id_idx]] = 1
                continue
            if line[2] == 'CDS':
                line[id_idx] = line[id_idx].replace('cds-', '')
                if not line[13] in protein_ids:
                    protein_ids[line[id_idx]] = 1
                continue

    if not os.path.exists(f"{genome_rna}.ok") or not os.path.exists(f"{genome_rna}"):
        if not args.rna:
            if not os.path.exists(f"{genome_rna}.bak") or not os.path.exists(f"{genome_rna}.bak.ok"):
                logging.info("RNA fasta file not provided. Extracting RNA sequences from genome fasta and gff file")
                cmd = f"agat_sp_extract_sequences.pl --gff {genome_gff} --fasta {genome_fasta} --mrna --output {genome_rna}.bak > {genome_rna}.log 2>&1"
                run_cmd(cmd)
                open(f"{genome_rna}.bak.ok", 'w').close()
            args.rna = f"{genome_rna}.bak"
        if args.rna:
            logging.info(f"Creating rna fasta file: {genome_rna}")
            with open(args.rna, 'r') as f, open(genome_rna, 'w') as out:
                for line in f:
                    if line.startswith('>'):
                        line = line.strip().split(' ')
                        name = line[0][1:].replace('rna-', '')
                        if not name in rna_ids:
                            logging.warning(f"RNA id {name} not found in genome tsv file. Maybe the ID is different in the gff file.")
                        out.write(f">{name}\n")
                    else:
                        out.write(line)
            open(f"{genome_rna}.ok", 'w').close()
            os.remove(f"{genome_rna}.bak")
            os.remove(f"{genome_rna}.bak.ok")

    if not os.path.exists(f"{genome_protein}.ok") or not os.path.exists(f"{genome_protein}"):
        if not args.protein:
            if not os.path.exists(f"{genome_protein}.bak") or not os.path.exists(f"{genome_protein}.bak.ok"):
                logging.info("Protein fasta file not provided. Extracting protein sequences from genome fasta and gff file")
                cmd = f"agat_sp_extract_sequences.pl --gff {genome_gff} --fasta {genome_fasta} --protein --output {genome_protein}.bak > {genome_protein}.log 2>&1"
                run_cmd(cmd)
                open(f"{genome_protein}.bak.ok", 'w').close()
            args.protein = f"{genome_protein}.bak"
        if args.protein:
            logging.info(f"Creating protein fasta file: {genome_protein}")
            with open(args.protein, 'r') as f, open(genome_protein, 'w') as out:
                for line in f:
                    if line.startswith('>'):
                        line = line.strip().split(' ')
                        name = line[0][1:].replace('rna-', '')
                        if not name in protein_ids:
                            logging.warning(f"Protein id {name} not found in genome tsv file. Maybe the ID is different in the gff file.")
                        out.write(f">{name}\n")
                    else:
                        if line.endswith('*'):
                            line = line[:-1]
                        out.write(line)
            open(f"{genome_protein}.ok", 'w').close()
            os.remove(f"{genome_protein}.bak")
            os.remove(f"{genome_protein}.bak.ok")


if __name__ == '__main__':
    main()
