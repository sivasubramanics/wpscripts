#!/usr/bin/env python3
# script to parse blat output and create clusters
# usage: blat_to_clusters.py -b <blat output> -o <output file> -c <coverage cut-off>
# example: blat_to_clusters.py -b Trinity_renamed.fasta_vs_Trinity_renamed.fasta.blat -o Trinity_renamed.fasta_vs_Trinity_renamed.fasta.blat.clusters -c 0.8
# version: 1.0
# contact: siva.selvanayagam[at]wur.nl
# date: 2023-07-27

import argparse
import os
import sys
import subprocess
import pandas as pd
import networkx as nx
from datetime import datetime
import pysam
from collections import defaultdict

blat_head="matches\tmismatches\trep.matches\tNs\tQgapcount\tQgapbases\tTgapcount\tTgapbases\tstrand\tQname\tQsize\tQstart\tQend\tTname\tTsize\tTstart\tTend\tblockcount\tblocksizes\tqStarts\ttStarts"

class CheckPoint:
    def __init__(self):
        self.out_dir = False
        self.groups_dir = False
        self.clusters_done = False
        self.groups_done = False
        self.sigletons_dir = False
        self.blat_done = False
        self.edges_file = False
        self.aln_done = False
        self.alnsum_done = False

    def set_alnsum_done(self):
        self.alnsum_done = True

    def set_aln_done(self):
        self.aln_done = True
    
    def set_out_dir(self):
        self.out_dir = True
    
    def set_groups_dir(self):
        self.groups_dir = True
    
    def set_clusters_done(self):
        self.clusters_done = True
    
    def set_sigletons_dir(self):
        self.sigletons_dir = True
    
    def set_blat_done(self):
        self.blat_done = True
    
    def set_edges_file(self):
        self.edges_file = True
    
    def __str__(self):
        return f"out_dir: {self.out_dir}\n" \
               f"groups_dir: {self.groups_dir}\n" \
               f"clusters_done: {self.clusters_done}\n" \
               f"sigletons_dir: {self.sigletons_dir}\n" \
               f"blat_done: {self.blat_done}\n" \
               f"edges_file: {self.edges_file}\n" \
               f"aln_done: {self.aln_done}\n" \
               f"alnsum_done: {self.alnsum_done}\n"
    
def get_check_points(out_dir):
    """
    Get the check points
    """
    check_points = CheckPoint()
    if os.path.exists(out_dir):
        check_points.set_out_dir()
        if os.path.exists(out_dir + "/groups"):
            check_points.set_groups_dir()
        if os.path.exists(out_dir + "/sigletons"):
            check_points.set_sigletons_dir()
        if os.path.exists(out_dir + "/blat.psl.ok"):
            check_points.set_blat_done()
        if os.path.exists(out_dir + "/edges.tsv.ok"):
            check_points.set_edges_file()
        if os.path.exists(out_dir + "/clusters.tsv.ok"):
            check_points.set_clusters_done()
        if os.path.exists(out_dir + "/aln.ok"):
            check_points.set_aln_done()
        if os.path.exists(out_dir + "/aln.summary.ok"):
            check_points.set_alnsum_done()
    return check_points

def reverse_complement(seq):
    """
    Reverse complement a sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def fold_seq(seq, length):
    """
    Fold a sequence
    """
    return "\n".join(seq[i:i+length] for i in range(0, len(seq), length))


class FASTA(object):
    def __init__(self, name, sequence, description=None):
        self.name = name
        self.sequence = sequence
        self.description = description
    
    def __str__(self):
        return f">{self.name}\n{self.sequence}"
    
    def __len__(self):
        return len(self.sequence)
    
    def rev_complement(self):
        return FASTA(self.name, reverse_complement(self.sequence))

    def write_seq(self, handle):
        handle.write(f">{self.name}\n{self.sequence}\n")
    
    def write_aln(self, handle, strand='+'):
        if len(self.name) > 20:
            name = self.name[:20]
        else:
            name = self.name + " "*(20-len(self.name)) 
        handle.write(f"{name}\t{strand}\t{self.sequence}\n")


def parse_fasta(fasta_file):
    """
    Parse a fasta file
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


def correct_seq_direction(transcript_direction, in_fasta, out_fasta):
    """
    Correct the sequence direction
    """
    sequences = parse_fasta(in_fasta)
    with open(out_fasta, 'w') as f:
        for seq in sequences:
            if transcript_direction[seq.name] == '-':
                f.write(f">{seq.name} converted_direction=true\n{fold_seq(reverse_complement(seq.sequence), 80)}\n")
            else:
                f.write(f">{seq.name} converted_direction=false\n{fold_seq(seq.sequence, 80)}\n")
        f.close()
    

def extract_seq(in_fasta, out_fasta, seq_ids):
    """
    Extract sequences from a fasta file
    """
    fasta_file = pysam.FastaFile(in_fasta)
    with open(out_fasta, 'w') as f:
        for seq_id in seq_ids:
            f.write(f">{seq_id}\n{fasta_file.fetch(seq_id)}\n")
    fasta_file.close()


def print_log(msg):
    """
    Print log messages
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", file=sys.stderr)


def run_cmd(cmd_line):
    """
    Run a command line
    """
    p = subprocess.Popen(cmd_line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        raise Exception(err)
    return out    


def check_tools(tools):
    """
    Check if the required tools are installed
    """
    for tool in tools:
        cmd_line = f"which {tool}"
        try:
            run_cmd(cmd_line)
        except:
            sys.exit(f"{tool} is not installed. Please install blat and add it to PATH")


def parse_arguments():
    # create an instance of Argument Parser and add positional argument
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="input fasta file", required=True)
    parser.add_argument("-o", "--outdir", help="output directory", required=True)
    parser.add_argument("-x", "--rewrite", help="rewrite the output directory", action="store_true")
    parser.add_argument("-c", "--cutoff", type=float, help="coverage cut-off", default=0.8)
    parser.add_argument("-t", "--threads", type=int, help="number of threads", default=2)
    args = parser.parse_args()
    return args


def check_fasta_index(args):
    # check if the fasta file is indexed
    if not os.path.exists(args.fasta + ".fai"):
        print_log(f"Indexing {args.fasta}")
        cmd_line = f"samtools faidx {args.fasta}"
        run_cmd(cmd_line)


def main():
    # start the clock
    start_time = datetime.now()

    # parse arguments
    args = parse_arguments()

    # check if the tools are installed
    check_tools(["samtools", "faSomeRecords", "clustalo"])

    # check if the fasta file is indexed
    check_fasta_index(args)
    
    print(f"{get_check_points(args.outdir)}")

    # check if the output directory exists
    if args.rewrite:
        # remove the old files
        print_log(f"Removing old files if any")
        cmd_line = f"rm -rf {args.outdir}"
        run_cmd(cmd_line)
        check_points = get_check_points(args.outdir)
    elif os.path.exists(args.outdir + "/.checkpoints"):
        check_points = get_check_points(args.outdir)
    else:
        check_points = get_check_points(args.outdir)
    
    # create output directory
    if not check_points.out_dir:
        # create output directory
        print_log(f"Creating output directory {args.outdir}")
        cmd_line = f"mkdir -p {args.outdir}"
        run_cmd(cmd_line)
        cmd_line = f"mkdir -p {args.outdir}/logs"
        run_cmd(cmd_line)
        cmd_line = f"mkdir -p {args.outdir}/sh/"
        run_cmd(cmd_line)    

    in_seq_desc = {}
    in_seqs = parse_fasta(args.fasta)
    for seq in in_seqs:
        in_seq_desc[seq.name] = seq.description

    in_seq_len = {}
    in_seqs = parse_fasta(args.fasta)
    for seq in in_seqs:
        in_seq_len[seq.name] = len(seq.sequence)

    if not check_points.blat_done:
        # run blat
        run_blat(args)
    else:
        print_log(f"blat is already done")
    
    if not check_points.edges_file:
        # open blat output file and output file
        print_log(f"Creating edges file")
        combinatons = defaultdict(dict)
        transcript_direction = {}
        blat = open(args.outdir + '/blat.psl', "r")
        output = open(args.outdir + "/edges.tsv", "w")
        output.write(f"Query\tQuery_length\tQuery_coverage\tTarget\tTarget_length\tTarget_coverage\tNo_blocks\t{blat_head}\n")
        for line in blat:
            if not line.startswith("#"):
                line = line.strip()
                in_line = line
                line = line.split("\t")
                if line[9] == line[13]:
                    continue
                if line[9] in combinatons and line[13] in combinatons[line[9]]:
                    continue
                combinatons[line[9]][line[13]] = 1
                combinatons[line[13]][line[9]] = 1
                q_cov = round(float(line[0]) / float(line[10]), 2)
                t_cov = round(float(line[0]) / float(line[14]), 2)
                if q_cov >= args.cutoff or t_cov >= args.cutoff:
                    if float(line[17]) <= 4:
                        output.write(f"{line[9]}\t{line[10]}\t{q_cov}\t{line[13]}\t{line[14]}\t{t_cov}\t{line[17]}\t{in_line}\n")
                        transcript_direction = update_strand_info(transcript_direction, line[9], line[13], line[8])
        blat.close()
        output.close()
        strand_info = open(args.outdir + "/strand_info.tsv", "w")
        for k, v in transcript_direction.items():
            strand_info.write(f"{k}\t{v}\n")
        strand_info.close()
        cmd_line = f"touch {args.outdir + '/edges.tsv.ok'}"
        run_cmd(cmd_line)
    else:
        print_log(f"edges file is already created")

    transcript_direction = {}
    strand_info = open(args.outdir + "/strand_info.tsv", "r")
    for line in strand_info:
        line = line.strip()
        line = line.split("\t")
        transcript_direction[line[0]] = line[1]
    strand_info.close()
    
    if not check_points.clusters_done:
        # create clusters file
        print_log(f"Creating clusters file")
        edges = pd.read_csv(args.outdir + "/edges.tsv", sep='\t')
        G = nx.from_pandas_edgelist(edges, source='Query', target='Target')
        components = nx.connected_components(G)
        output = open(args.outdir + "/clusters.tsv", 'w')
        # Output nodes of each component into separate files
        for i, component in enumerate(components):
            g_id = f'group{i}'
            for node in component:
                output.write(f'{g_id}\t{node}\t{transcript_direction[node]}\t{in_seq_len[node]}\n')
        output.close()
        cmd_line = f"touch {args.outdir + '/clusters.tsv.ok'}"
        run_cmd(cmd_line)
    else:
        print_log(f"clusters file is already created")

    print_log(f"parse the clusters file")
    clust_dict = {}
    clusters = open(args.outdir + "/clusters.tsv", 'r')
    for line in clusters:
        line = line.strip()
        line = line.split("\t")
        if line[0] not in clust_dict:
            clust_dict[line[0]] = []
        clust_dict[line[0]].append(line[1])
    clusters.close()
    print_log(f"Total number of clusters: {len(clust_dict)}")

    # split the fasta file into groups
    if not check_points.groups_dir:
        # create groups directory
        print_log(f"Creating groups directory")
        cmd_line = f"mkdir -p {args.outdir}/groups"
        run_cmd(cmd_line)
        sh = open(args.outdir + "/sh/faSomeRecords.sh", 'w')
        for k, v in clust_dict.items():
            # Its commented out because it takes a lot of time to create the fasta file, so we decided to use the faSomeRecords tool in parallel
            # extract_seq(args.fasta, args.outdir + f"/groups/{k}.fasta", v)
            with open(args.outdir + f"/groups/{k}.ids", 'w') as f:
                f.write("\n".join(v) + "\n")
            cmd_line = f"faSomeRecords {args.fasta} {args.outdir}/groups/{k}.ids {args.outdir}/groups/{k}.fasta"
            sh.write(cmd_line + "\n")
        sh.close()
        cmd_line = f"parallel -j {args.threads} < {args.outdir}/sh/faSomeRecords.sh"
        run_cmd(cmd_line)

        for k, v in clust_dict.items():
            correct_seq_direction(transcript_direction, args.outdir + f"/groups/{k}.fasta", args.outdir + f"/groups/{k}.corrected.fasta")
        cmd_line = f"touch {args.outdir + '/groups.ok'}"
        run_cmd(cmd_line)
    else:
        print_log(f"groups directory is already created")
    
    if not check_points.aln_done:
        print_log(f"Creating alignment for each group")
        sh = open(args.outdir + "/sh/clustalo.sh", 'w')
        # align the fasta files
        for k, v in clust_dict.items():
            cmd_line = f"clustalo -i {args.outdir}/groups/{k}.corrected.fasta -o {args.outdir}/groups/{k}.aln.fasta --threads=2 --force > {args.outdir}/logs/{k}.aln.log 2>&1"
            sh.write(cmd_line + "\n")
        sh.close()
        cmd_line = f"parallel -j {int(args.threads/2)} < {args.outdir}/sh/clustalo.sh"
        run_cmd(cmd_line)
        cmd_line = f"touch {args.outdir + '/aln.ok'}"
        run_cmd(cmd_line)
    else:
        print_log(f"alignment is already done")
    
    if not check_points.alnsum_done:
        print_log(f"Creating summary file")
        summary = open(args.outdir + "/aln.summary.tsv", 'w')
        summary.write("Group\tSeq_id\tGroup_len\tSeq_len\tStrand\n")
        for k, v in clust_dict.items():
            sequences = parse_fasta(args.outdir + f"/groups/{k}.aln.fasta")
            aln_seq = open(args.outdir + f"/groups/{k}.aln", 'w')
            summary.write(f"--------------------------------------------\n")
            for seq in sequences:
                summary.write(f"{k}\t{seq.name}\t{len(seq.sequence)}\t{len(seq.sequence.replace('-', ''))}\t{transcript_direction[seq.name]}\t{in_seq_desc[seq.name]}\n")
                seq.write_aln(aln_seq, transcript_direction[seq.name])
            aln_seq.close()
        summary.close()
        cmd_line = f"touch {args.outdir + '/aln.summary.ok'}"
        run_cmd(cmd_line)
    else:
        print_log(f"alignment summary is already done")
        
    # print the run time
    run_time = str(datetime.now() - start_time).split(":")
    print_log(f"Total time taken: {run_time[0]}:{run_time[1]}:{round(float(run_time[2]))} ")

def run_blat(args):
    """
    function to run blat
    """
    print_log(f"running blat in parallel")
    cmd_line = f"pblat -t=dna -q=dna -threads={args.threads} -noHead {args.fasta} {args.fasta} {args.outdir + '/blat.psl'} > {args.outdir + '/logs/blat.log'} 2>&1"    
    run_cmd(cmd_line)
    cmd_line = f"touch {args.outdir + '/blat.psl.ok'}"
    run_cmd(cmd_line)

def update_strand_info(transcript_direction, qName, tName, strand):
    """
    Update the strand information
    """
    # Directionality of both transcripts already defined
    if((tName in transcript_direction) and (qName in transcript_direction)):
        return transcript_direction
    # Transcript directiionality already defined
    elif(tName in transcript_direction): 
        tdir = transcript_direction[tName]
        if(tdir == strand): #If both the transcript and strand same
            transcript_direction[qName] = '+'
        else:
            transcript_direction[qName] = '-'
    # Query transcript already defined
    elif(qName in transcript_direction): 
        qdir = transcript_direction[qName]
        if(qdir == strand): # If both the query and strand same
            transcript_direction[tName] = '+'
        else:
            transcript_direction[tName] = '-'
    # Neither of the sequences defined
    else:
        if(strand=='+'): # Both transcripts in the same direction
            transcript_direction[qName] = '+'
            transcript_direction[tName] = '+'

        else: # arbitrarily define transcript as postive and query as negative
            transcript_direction[tName] = '+'
            transcript_direction[qName] = '-'
    return transcript_direction

if __name__ == "__main__":
    main()