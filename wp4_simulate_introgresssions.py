#!/usr/bin/env python3

"""
This script takes two fasta files and introduces introgressions from one fasta file to another fasta file, also returns the locations where it has introduced the introgressions
"""

import argparse
import random


class FASTA:
    def __init__(self, name: str, sequence: str, description: str = None):
        self.name = name
        self.sequence = sequence
        self.description = description

    def __str__(self):
        fold_seq = self.fold(60)
        return fold_seq

    def fold(self, width: int) -> str:
        out_str = f">{self.name}"
        if self.description:
            out_str += f" {self.description}"
        out_str += "\n"
        out_str += '\n'.join([self.sequence[i:i + width] for i in range(0, len(self.sequence), width)])
        return out_str

    def __len__(self):
        return len(self.sequence)


def parse_fasta(file):
    name = None
    description = None
    sequence = ""
    with open(file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield FASTA(name, sequence, description)
                name = line[1:].split()[0]
                description = line[len(name) + 2:]
                sequence = ""
            else:
                sequence += line
        yield FASTA(name, sequence, description)


def introduce_introgression(donor_seq, recurrent_seq, introgression_size):
    donor_start = random.randint(0, len(donor_seq) - introgression_size)
    donor_end = donor_start + introgression_size

    recurrent_start = random.randint(0, len(recurrent_seq) - introgression_size)
    recurrent_end = recurrent_start + introgression_size

    introgressed_seq = (
            recurrent_seq[:recurrent_start] +
            donor_seq[donor_start:donor_end] +
            recurrent_seq[recurrent_end:]
    )

    return introgressed_seq, (donor_start, donor_end, recurrent_start, recurrent_end)


def main():
    parser = argparse.ArgumentParser(description='Introduce introgressions from one fasta file to another fasta file')
    parser.add_argument('-d', '--donor', help='Donor fasta file', required=True)
    parser.add_argument('-e', '--donor_gtf', help='Donor gtf file', required=True)
    parser.add_argument('-r', '--recurrent', help='Recurrent fasta file', required=True)
    parser.add_argument('-s', '--recurrent_gtf', help='Recurrent gtf file', required=True)
    parser.add_argument('-o', '--output', help='Output prefix', required=True)
    args = parser.parse_args()

    donor_sequences = list(parse_fasta(args.donor))
    recurrent_sequences = list(parse_fasta(args.recurrent))

    introgressed_sequences = []
    introgression_info = []

    # if input file has more thant 1 sequence then quite the program
    if len(donor_sequences) > 1:
        print("Input donor file has more than 1 sequence, please provide only 1 sequence in donor file")
        quit()
    if len(recurrent_sequences) > 1:
        print("Input recurrent file has more than 1 sequence, please provide only 1 sequence in recurrent file")
        quit()

    donor_seq = donor_sequences[0]
    recurrent_seq = recurrent_sequences[0]

    if recurrent_seq:
        introgression_size = random.randint(10000000, 60000000)  # 1Mb to 6Mb
        new_seq, (donor_start, donor_end, recurrent_start, recurrent_end) = introduce_introgression(
            donor_seq.sequence, recurrent_seq.sequence, introgression_size
        )

        introgressed_sequences.append(FASTA(recurrent_seq.name, new_seq, recurrent_seq.description))
        introgression_info.append((
            donor_seq.name, donor_start, donor_end,
            recurrent_seq.name, recurrent_start, recurrent_end
        ))
    else:
        print(f"Warning: No matching sequence found in recurrent for {donor_seq.name}")

    # Write introgressed sequences to FASTA file
    with open(f"{args.output}.fasta", "w") as f:
        for seq in introgressed_sequences:
            seq.name = "intr"
            f.write(str(seq) + "\n")

    # Write introgression information to TSV file
    with open(f"{args.output}.tsv", "w") as f:
        f.write("Donor_Chr\tDonor_Start\tDonor_End\tRecurrent_Chr\tRecurrent_Start\tRecurrent_End\tlength\n")
        for info in introgression_info:
            f.write("\t".join(map(str, info)) + f"\t{info[2] - info[1]}\n")

    # Write genes from donor gtf file for the regions between donor_start and donor_end
    with open(f"{args.output}.gtf", "w") as f:
        flag = False
        with open(args.donor_gtf) as gtf:
            for line in gtf:
                line = line.strip()
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                start = int(fields[3])
                end = int(fields[4])
                if fields[0] == donor_seq.name and fields[2] == "gene":
                    if start >= donor_start and end <= donor_end:
                        flag = True
                    else:
                        flag = False
                if flag:
                    # convert the start of the gene relative to the actual position in the introgresed sequence
                    seq_name = "intr"
                    start = int(fields[3]) - donor_start + recurrent_start
                    end = int(fields[4]) - donor_start + recurrent_start
                    f.write(f"{seq_name}\t{fields[1]}\t{fields[2]}\t{start}\t{end}\t{fields[5]}\t{fields[6]}\t{fields[7]}\t{fields[8]}\n")
        flag = False
        with open(args.recurrent_gtf) as gtf:
            for line in gtf:
                line = line.strip()
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                start = int(fields[3])
                end = int(fields[4])
                if fields[0] == recurrent_seq.name and fields[2] == "gene":
                    if end >= recurrent_start or start >= recurrent_end:
                        flag = True
                    else:
                        flag = False
                if flag:
                    seq_name = "intr"
                    f.write(f"{seq_name}\t{fields[1]}\t{fields[2]}\t{start}\t{end}\t{fields[5]}\t{fields[6]}\t{fields[7]}\t{fields[8]}\n")

    print(f"Introgressed sequences written to {args.output}_introgressed.fasta")
    print(f"Introgression information written to {args.output}_introgressions.tsv")


if __name__ == "__main__":
    main()