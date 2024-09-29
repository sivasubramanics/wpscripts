#!/usr/bin/env python3

"""
This script is a try-out script to convert the psl file to gfa file
"""
import sys

class FASTA:
    def __init__(self, id, seq, desc: str = None):
        self.id = id
        self.seq = seq
        self.desc = desc

    def __repr__(self):
        return f">{self.id}\n{self.seq}\n"


class PSL:
    def __init__(self, psl_line):
        psl_entries = psl_line.strip().split("\t")
        self.matches = int(psl_entries[0])
        self.mismatches = int(psl_entries[1])
        self.rep_matches = int(psl_entries[2])
        self.n_count = int(psl_entries[3])
        self.q_num_insert = int(psl_entries[4])
        self.q_base_insert = int(psl_entries[5])
        self.t_num_insert = int(psl_entries[6])
        self.t_base_insert = int(psl_entries[7])
        self.strand = psl_entries[8]
        self.q_name = psl_entries[9]
        self.q_size = int(psl_entries[10])
        self.q_start = int(psl_entries[11])
        self.q_end = int(psl_entries[12])
        self.t_name = psl_entries[13]
        self.t_size = int(psl_entries[14])
        self.t_start = int(psl_entries[15])
        self.t_end = int(psl_entries[16])
        self.block_count = int(psl_entries[17])
        self.block_sizes = psl_entries[18]
        self.q_starts = psl_entries[19]
        self.t_starts = psl_entries[20]
        self.psl_score = self.matches + (self.rep_matches / 2) - self.mismatches - self.q_num_insert - self.t_num_insert
        self.q_score = self.psl_score / self.q_size * 1000
        self.qcov = (self.q_end - self.q_start) / self.q_size
        self.tcov = (self.t_end - self.t_start) / self.t_size

    def __repr__(self):
        return (f"{self.matches}\t"
                f"{self.mismatches}\t"
                f"{self.rep_matches}\t"
                f"{self.n_count}\t"
                f"{self.q_num_insert}\t"
                f"{self.q_base_insert}\t"
                f"{self.t_num_insert}\t"
                f"{self.t_base_insert}\t"
                f"{self.strand}\t"
                f"{self.q_name}\t"
                f"{self.q_size}\t"
                f"{self.q_start}\t"
                f"{self.q_end}\t"
                f"{self.t_name}\t"
                f"{self.t_size}\t"
                f"{self.t_start}\t"
                f"{self.t_end}\t"
                f"{self.block_count}\t"
                f"{self.block_sizes}\t"
                f"{self.q_starts}\t"
                f"{self.t_starts}\n")

    def __str__(self):
        return self.__repr__()

    def get_blocks(self):
        q_starts = self.q_starts.split(",")
        t_starts = self.t_starts.split(",")
        block_sizes = self.block_sizes.split(",")
        blocks = []
        for i in range(len(block_sizes)):
            blocks.append((int(q_starts[i]), int(t_starts[i]), int(block_sizes[i])))
        return blocks


def parse_fasta(in_fasta):
    id = None
    seq = None
    desc = None
    with open(in_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                if seq:
                    yield FASTA(id, seq, desc)
                id = line.strip().split()[0][1:]
                desc = " ".join(line.strip().split()[1:])
                seq = ""
            else:
                seq += line.strip()
        yield FASTA(id, seq, desc)


def main():
    in_fasta = sys.argv[1]
    in_psl = sys.argv[2]
    out_gfa = sys.argv[3]

    transcripts = {}
    for seq in parse_fasta(in_fasta):
        transcripts[seq.id] = seq

    with open(in_psl, "r") as f, open(out_gfa, "w") as out:
        for line in f:
            psl = PSL(line)
            if psl.q_name == psl.t_name:
                continue
            if psl.qcov < 0.8 or psl.tcov < 0.8:
                continue
            print(psl)

if __name__ == "__main__":
    main()
