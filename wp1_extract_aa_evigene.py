#!/usr/bin/env python3

"""After evigene we may run kraken2 decontamination. Then we do rename the contigs to the preferred format. This
script will help to extract the protein sequences from the evigene output.
"""

import argparse

class Fasta:
    def __init__(self, header, sequence, description=None):
        self.header = header
        self.sequence = sequence
        self.description = description

    def __str__(self):
        return f">{self.header} {self.description}\n{self.fold()}"

    def fold(self, width=60):
        return '\n'.join([self.sequence[i:i+width] for i in range(0, len(self.sequence), width)])


def parse_fasta(fasta_file):
    with open(fasta_file) as f:
        header = None
        sequence = ''
        for line in f:
            if line.startswith('>'):
                if header:
                    yield Fasta(header, sequence, description)
                header = line.strip().split()[0][1:]
                if len(line.strip().split()) > 1:
                    description = ' '.join(line.strip().split()[1:])
                else:
                    description = ''
                sequence = ''
            else:
                sequence += line.strip()
        yield Fasta(header, sequence, description)


def main():
    parser = argparse.ArgumentParser(description='Extract protein sequences from evigene output')
    parser.add_argument('-i', '--evigene_output', help='Evigene output aa file', required=True)
    parser.add_argument('-m', '--map_file', help='Map file to rename the protein', required=True)
    parser.add_argument('-o', '--output', help='Output file to store protein sequences', required=True)
    args = parser.parse_args()

    map_dict = {}
    with open(args.map_file) as f:
        for line in f:
            new, old = line.strip().split()
            if old in map_dict:
                print(f"Duplicate entry found for {old}")
                exit(1)
            map_dict[old] = new

    with open(args.output, 'w') as out:
        for fasta in parse_fasta(args.evigene_output):
            if fasta.header in map_dict:
                out.write(f">{map_dict[fasta.header]} {fasta.description}\n{fasta.fold()}\n")

if __name__ == '__main__':
    main()


