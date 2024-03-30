#!/usr/bin/env python3

import argparse

class FASTA:
    """
    Class to store the fasta sequence
    """
    def __init__(self, name: str, sequence: str, description: str = None, readcounts: float = None):
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
        out_str += '\n'.join([self.sequence[i:i+width] for i in range(0, len(self.sequence), width)])
        return out_str

    def __len__(self):
        return len(self.sequence)

    def get_seq(self, start, end):
        return self.sequence[start:end]

    def get_gc_content(self):
        return (self.sequence.count("G") + self.sequence.count("C")) / len(self.sequence)

    def get_at_content(self):
        return (self.sequence.count("A") + self.sequence.count("T")) / len(self.sequence)

    def get_n_content(self):
        return self.sequence.count("N")

    def get_repeat_content(self, is_hardmasked: bool = False):
        if is_hardmasked:
            return self.sequence.count("N")
        else:
            return self.sequence.count("a") + self.sequence.count("t") + self.sequence.count("g") + self.sequence.count("c")


def parse_fasta(fasta_file) -> None:
    """
    Parse the fasta file
    """
    with open(fasta_file, 'r') as f:
        name = ""
        sequence = ""
        description = ""
        begun = False
        for line in f:
            if line.startswith(">"):
                if begun:
                    yield FASTA(name, sequence, description)
                line = line.strip()
                name = line[1:].split()[0]
                description = line[1:]
                sequence = ""
                begun = True
            else:
                sequence += line.strip()
        yield FASTA(name, sequence, description)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', help='input fasta file', required=True)
    parser.add_argument('-b', '--bed', help='input bed file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-H', '--hardmask', help='input fasta soft masked for repeats', action='store_true', default=False)
    parser.add_argument('-B', '--bedheader', help='bed file has header', action='store_true', default=False)

    args = parser.parse_args()

    fasta_dict = {}
    for fasta in parse_fasta(args.fasta):
        fasta_dict[fasta.name] = fasta

    with open(args.bed, 'r') as bed, open(args.output, 'w') as out:
        out.write("seq_id\tstart\tend\tlength\trepeat_content\n")
        if args.bedheader:
            next(bed)
        for line in bed:
            line = line.strip().split('\t')
            seq_id = line[0]
            start = int(line[1])
            end = int(line[2])
            if seq_id in fasta_dict:
                seq = FASTA(f"{seq_id}:{start}-{end}", fasta_dict[seq_id].get_seq(start, end))
                out.write(f"{seq_id}\t{start}\t{end}\t{len(seq)}\t{seq.get_repeat_content(args.hardmask)/len(seq):.2f}\n")
            else:
                print(f"Sequence {seq_id} not found in fasta file")


if __name__ == "__main__":
    main()
