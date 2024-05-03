#!/usr/bin/env python3

import os
import sys
import logging
import argparse
import subprocess
import gzip
import concurrent.futures
from collections import defaultdict

TR_FQ = "in.fastq"
TR_FA = "in.fasta"
TR_RENAMED_FA = "renamed.fasta"
ALN_SAM = "aln.sam"
CLT_DIR = "clusters"
LOGDIR = "logs"
TAB = "\t"


def fold_seq(seq, n):
    """
    Fold the sequence to n characters
    """
    return '\n'.join([seq[i:i + n] for i in range(0, len(seq), n)])


class FASTA:
    def __init__(self, name, seq, description):
        self.name = name
        self.seq = seq
        self.description = description

    def __str__(self):
        return f">{self.name} {self.description}\n{fold_seq(self.seq, 80)}"

    def to_upper(self):
        self.seq = self.seq.upper()
        return self.seq

    def to_fastq(self):
        return f"@{self.name}\n{self.to_upper()}\n+\n{'I' * len(self.seq)}"


def parse_fasta(fasta):
    """
    Parse the fasta file and yeild the FASTA object
    """
    with open(fasta) as f:
        name = None
        seq = ""
        description = ""
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield FASTA(name, seq, description)
                name = line.split(" ")[0][1:]
                seq = ""
                # get the description as string excluding the name
                description = " ".join(line.split(" ")[1:])
            else:
                seq += line
        yield FASTA(name, seq, description)


def fasta_to_fastq(fasta, fastq):
    """
    Convert fasta to fastq
    """
    with open(fastq, "w") as f:
        for seq in parse_fasta(fasta):
            f.write(f"{seq.to_fastq()}\n")


def sam_to_fastq_clustrs(sam, cluster_dir):
    """
    Convert the sam file to fastq files for each cluster based on the reference name
    """


def run_cmd(cmd, log_file=None, prtlog=True):
    """
    Run the command and return the stdout, and stderr
    """
    if prtlog:
        if log_file:
            logging.info(f"Running command: {cmd} and logging to {log_file}")
        else:
            logging.info(f"Running command: {cmd}")
    try:
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        # get system signal for the subprocess
        if p.returncode != 0:
            logging.error(f"Error: {stderr.decode()}")
            sys.exit(1)
        if log_file:
            with open(log_file, "w") as f:
                f.write(stdout.decode())
        return stdout.decode()
    except Exception as e:
        logging.error(e)
        sys.exit(1)


class FASTQ:
    def __init__(self, name):
        self.name = name
        self.fseq = ""
        self.rseq = ""
        self.fqual = ""
        self.rqual = ""

    def __str__(self):
        return f"@{self.name} 1\n{self.fseq}\n+\n{self.fqual}\n@{self.name} 2\n{self.rseq}\n+\n{self.rqual}"


def explain_sam_flag(sam_flag):
    """
    This function process the SAM flag and returns the information about the flag
    """
    flag_dict = {
        0x1: "read paired",
        0x2: "read mapped in proper pair",
        0x4: "read unmapped",
        0x8: "mate unmapped",
        0x10: "read reverse strand",
        0x20: "mate reverse strand",
        0x40: "first in pair",
        0x80: "second in pair",
        0x100: "not primary alignment",
        0x200: "read fails platform/vendor quality checks",
        0x400: "read is PCR or optical duplicate",
        0x800: "supplementary alignment"
    }
    return [flag_dict[flag] for flag in flag_dict if sam_flag & flag]


class SAM:
    def __init__(self, line):
        self.qname = line[0]
        self.flag = int(line[1])
        self.rname = line[2]
        self.pos = int(line[3])
        self.mapq = int(line[4])
        self.cigar = line[5]
        self.rnext = line[6]
        self.pnext = int(line[7])
        self.tlen = int(line[8])
        self.seq = line[9]
        self.qual = line[10]
        self.opt = line[11:]
        self.flag_info = explain_sam_flag(self.flag)

    def is_fread(self):
        return "first in pair" in self.flag_info

    def is_rread(self):
        return "second in pair" in self.flag_info

    def is_primary(self):
        if "not primary alignment" in self.flag_info or "supplementary alignment" in self.flag_info:
            return False
        return True

    def is_unmapped(self):
        return "read unmapped" in self.flag_info

    def __str__(self):
        return f"{self.qname}\t{self.flag}\t{self.rname}\t{self.pos}\t{self.mapq}\t{self.cigar}\t{self.rnext}\t{self.pnext}\t{self.tlen}\t{self.seq}\t{self.qual}\t{TAB.join(self.opt)}"


def process_sams(sams):
    """
    Process the SAM objects and return the FASTQ objects
    """
    for sam in sams:
        if sam.is_primary():
            print(sam.qname, sam)


def sam_to_fq(in_sam, nthreads=1):
    """
    Add sequences to the sam file 10th column for secondary alignements
    """
    out_sam = in_sam.replace(".sam", ".seq.tsv")
    base_dir = os.path.dirname(in_sam)
    with open(in_sam) as f, open(out_sam, "w") as fo:
        read_set = set()
        last_ref = None
        for line in f:
            if line.startswith("@"):
                continue
            line = line.rstrip().split("\t")
            line[2] = line[2].split("_")[0]
            if line[2] == "*":
                line[2] = "unmapped"
            if line[2] != last_ref and last_ref is not None:
                readid_list_file = os.path.join(base_dir, 'clusters', f"{last_ref}.readids")
                with open(readid_list_file, "w") as f:
                    for readid in read_set:
                        f.write(f"{readid}\n")
                read_set = set()
            last_ref = line[2]
            if line[0] in read_set:
                continue
            read_set.add(line[0])
        readid_list_file = os.path.join(base_dir, 'clusters', f"{last_ref}.readids")
        with open(readid_list_file, "w") as f:
            for readid in read_set:
                f.write(f"{readid}\n")
        logging.info(f"Extracting fastq per cluster")
        with concurrent.futures.ThreadPoolExecutor(max_workers=nthreads) as executor:
            for readid_list_file in os.listdir(os.path.join(base_dir, 'clusters')):
                if readid_list_file.endswith(".readids"):
                    executor.submit(extract_fasta, readid_list_file, base_dir)


def extract_fasta(readid_list_file, base_dir):
    clust_id = readid_list_file.split(".")[0]
    # base_dir = os.path.dirname(readid_list_file).rsrip("clusters")
    run_cmd(f"faSomeRecords {os.path.join(base_dir, 'fwd.fasta')} {os.path.join(base_dir, 'clusters', readid_list_file)} {os.path.join(base_dir, 'clusters', f'{clust_id}.1.fa')}", None, False)
    run_cmd(f"faSomeRecords {os.path.join(base_dir, 'rev.fasta')} {os.path.join(base_dir, 'clusters', readid_list_file)} {os.path.join(base_dir, 'clusters', f'{clust_id}.2.fa')}", None, False)


def fastq_to_fasta(in_fq, out_fa, read="fwd", nthreads=1):
    """
    Convert the fastq to fasta and rename the reads with incrementing numbers
    """
    if read == "fwd":
        pair = 1
    elif read == "rev":
        pair = 2
    else:
        raise Exception("Read should be either fwd or rev")
    if in_fq.endswith(".gz"):
        fq = gzip.open(in_fq, 'rb')
    else:
        fq = open(in_fq, "r")
    fo = open(out_fa, "w")
    line_no = 0
    fq_num = 0
    for line in fq:
        line = line.strip()
        if in_fq.endswith(".gz"):
            line = line.decode()
        if line_no % 4 == 0:
            fq_num += 1
        elif line_no % 4 == 1:
            seq = line
        elif line_no % 4 == 3:
            fo.write(f">{fq_num} {pair}\n{seq}\n")
        line_no += 1


def assemble_clusters(base_dir, nthreads, spades_threads=4, spades_mem=12):
    """
    read all the fasta file in the clusters directory and assemble them using rnaspades.py
    """
    for readids_file in os.listdir(os.path.join(base_dir, "clusters")):
        if not readids_file.endswith(".readids"):
            continue
        asm_out = os.path.join(base_dir, "assemblies", readids_file.split(".")[0], "transcripts.fasta")
        if os.path.exists(asm_out):
            continue
        clust_id = readids_file.split(".")[0]
        mem = int(((os.path.getsize(f"{os.path.join(base_dir, 'clusters', clust_id)}.1.fa") / 1024 / 1024) + (os.path.getsize(f"{os.path.join(base_dir, 'clusters', clust_id)}.2.fa") / 1024 / 1024)) * 24) + 1
        num_runs = int(nthreads / spades_threads) + 1
        cmd = (f"rnaspades.py -o {os.path.join(base_dir, 'assemblies', clust_id)} "
               f"-1 {os.path.join(base_dir, 'clusters', clust_id)}.1.fa "
               f"-2 {os.path.join(base_dir, 'clusters', clust_id)}.2.fa "
               f"--threads {spades_threads} --memory {mem}")
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_runs) as executor:
            executor.submit(run_cmd, cmd, os.path.join(base_dir, LOGDIR, f"asm_{clust_id}.log"), False)


def main():
    parser = argparse.ArgumentParser(description="Assembly based on the clustering")
    parser.add_argument("-t", "--transcript", help="Reference transcriptome", required=True)
    parser.add_argument("-o", "--output", help="Output directory", required=True)
    parser.add_argument("-f", "--fwd", help="Forward reads", required=True)
    parser.add_argument("-r", "--rev", help="Reverse reads", required=True)
    parser.add_argument("-p", "--threads", help="Number of threads", default=2, type=int)
    parser.add_argument("--force", help="Force the clustering step", action="store_true")

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # if output directory does not exist, create it
    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)
    if not os.path.exists(os.path.join(args.output, LOGDIR)):
        os.makedirs(os.path.join(args.output, LOGDIR), exist_ok=True)
    if not os.path.exists(os.path.join(args.output, CLT_DIR)):
        os.makedirs(os.path.join(args.output, CLT_DIR), exist_ok=True)

    if os.path.exists(os.path.join(args.output, TR_FQ)) and args.force is False:
        logging.info(f"Fastq file {os.path.join(args.output, TR_FQ)} already exists. Skipping the conversion step")
    else:
        # convert fasta to fastq
        logging.info("Converting fasta to fastq")
        fasta_to_fastq(args.transcript, os.path.join(args.output, TR_FQ))

    # check if args.output/isONclust/final_clusters.tsv exists
    if os.path.exists(os.path.join(args.output, "isONclust/final_clusters.tsv")) and args.force is False:
        logging.info("Final clusters already exists. Skipping the clustering step")
    else:
        # construct the clusters
        logging.info("Clustering the transcripts")
        cmd = f"isONclust --isoseq --t {args.threads} --fastq {os.path.join(args.output, TR_FQ)} --outfolder {os.path.join(args.output, 'isONclust')}"
        run_cmd(cmd, os.path.join(args.output, LOGDIR, "isONclust.log"))

    # rename the transcripts in final_clusters.tsv
    logging.info("Add new names to the transcripts")
    with open(os.path.join(args.output, "isONclust/final_clusters.tsv")) as f, open(
            os.path.join(args.output, "clusters.tsv"), "w") as fo:
        last_clust = None
        seq_count = 0
        for line in f:
            line = line.rstrip().split("\t")
            if line[0] != last_clust and last_clust:
                seq_count = 0
            fo.write(f"{line[0]}\t{line[1]}\t{line[0]}_{seq_count}\n")
            seq_count += 1
            last_clust = line[0]

    # rename fasta file with the third column in clusters.tsv
    logging.info("Renaming the fasta file")
    with open(os.path.join(args.output, "clusters.tsv")) as f:
        clust_dict = {}
        for line in f:
            line = line.rstrip().split("\t")
            clust_dict[line[1]] = line[2]
    with open(os.path.join(args.output, TR_RENAMED_FA), "w") as fo:
        for seq in parse_fasta(args.transcript):
            seq.description += f" old_name={seq.name}"
            seq.name = clust_dict[seq.name]
            fo.write(f"{seq}\n")

    if os.path.exists(os.path.join(args.output, "fwd.fasta")) and os.path.exists(os.path.join(args.output, "rev.fasta")) and args.force is False:
        logging.info("Fasta files already exists. Skipping the conversion step")
    else:
        # convert the fastq files to fasta
        logging.info(f"Converting the fastq files to fasta: {args.fwd}")
        fastq_to_fasta(args.fwd, os.path.join(args.output, "fwd.fasta"), "fwd", args.threads)
        logging.info(f"Converting the fastq files to fasta: {args.rev}")
        fastq_to_fasta(args.rev, os.path.join(args.output, "rev.fasta"), "rev", args.threads)

    if os.path.exists(os.path.join(args.output, ALN_SAM)) and args.force is False:
        logging.info(f"Sam file {os.path.join(args.output, ALN_SAM)} already exists. Skipping the alignment step")
    else:
        # align the reads to the transcriptome
        logging.info("Aligning the reads to the transcriptome")
        cmd = (f"minimap2 -2 -ax sr --secondary=yes -N1000 -t {args.threads} "
               f"{os.path.join(args.output, TR_RENAMED_FA)} "
               f"{os.path.join(args.output, 'fwd.fasta')} "
               f"{os.path.join(args.output, 'rev.fasta')} | "
               f"samtools sort -@ {args.threads} --output-fmt SAM -o {os.path.join(args.output, ALN_SAM)}")
        run_cmd(cmd, os.path.join(args.output, LOGDIR, "minimap2.log"))

    if os.path.exists(os.path.join(args.output, "clusters", "done")) and args.force is False:
        logging.info("SAM to FASTQ conversion already done. Skipping the conversion step")
    else:
        # convert the sam file to fastq files
        logging.info("SAM to FASTQ conversion")
        sam_to_fq(os.path.join(args.output, ALN_SAM), args.threads)
        run_cmd(f"touch {os.path.join(args.output, 'clusters', 'done')}", None, False)

    if os.path.exists(os.path.join(args.output, "assemblies", 'done')) and args.force is False:
        logging.info("Assemblies already exists. Skipping the assembly step")
    else:
        assemble_clusters(args.output, args.threads)
        run_cmd(f"touch {os.path.join(args.output, 'assemblies', 'done')}", None, False)




if __name__ == "__main__":
    main()
