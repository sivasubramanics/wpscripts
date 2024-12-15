#!/usr/bin/env python

import argparse
import re
import sys
from collections    import defaultdict

class Fasta:
    """
    Class to represent a fasta sequence
    """
    def __init__(self, name, seq, desc=None):
        self.name = name
        self.seq = seq.upper()
        self.desc = desc

    def __str__(self):
        return ">{} {}\n{}".format(self.name, self.desc, fold_fasta(self))
    
    def __repr__(self):
        return str(self)
    
    def __len__(self):
        return len(self.seq)
    

def fold_fasta(fasta, width=80):
    """
    Fold a fasta sequence to a given width
    """
    return "\n".join(fasta.seq[i:i+width] for i in range(0, len(fasta.seq), width))

def parse_fasta(fasta_file):
    """
    parse fasta and yield Fasta objects
    """
    begin = False
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                if begin:
                    yield Fasta(name, seq, desc)
                name = line.strip().split()[0][1:]
                desc = " ".join(line.strip().split()[1:])
                seq = ""
                begin = True
            else:
                seq += line.strip()
        yield Fasta(name, seq, desc)


def get_ids(fa_name, asm_type):
    """
    Get gene and transcript ids from fasta name
    """
    if asm_type == "trinity":
        # TRINITY_DN100000_c0_g1_i1
        g_id = re.search(r"(TRINITY_DN\d+_c\d+_g\d+)_i\d+", fa_name).group(1)
        t_num = re.search(r"TRINITY_DN\d+_c\d+_g\d+_i(\d+)", fa_name).group(1)
        t_id = re.search(r"(TRINITY_DN\d+_c\d+_g\d+_i\d+)", fa_name).group(1)
    elif asm_type == "rnaspades":
        # NODE_1_length_15401_cov_332.591075_g0_i0
        g_id = re.search(r"NODE_\d+_length_\d+_cov_\d+.\d+_(g\d+)_i\d+", fa_name).group(1)
        t_num = re.search(r"NODE_\d+_length_\d+_cov_\d+.\d+_g\d+_i(\d+)", fa_name).group(1)
        t_id = re.search(r"(NODE_\d+_length_\d+_cov_\d+.\d+_g\d+_i\d+)", fa_name).group(1)
    elif asm_type == "fasta":
        g_id = re.search(r"(.+)", fa_name).group(1)
        t_id = re.search(r"(.+)", fa_name).group(1)
        t_num = -1
    elif asm_type == "evigene":
        g_id = re.search(r"(NonamEVm\d+)t\d+", fa_name).group(1)
        t_num = re.search(r"NonamEVm\d+t(\d+)", fa_name).group(1)
        t_id = re.search(r"(NonamEVm\d+t\d+)", fa_name).group(1)
    else:
        print("Assembly type not recognized")
        sys.exit(1)
    return g_id, t_id, t_num
    

def main():
    parser = argparse.ArgumentParser(description='Rename isoforms in a fasta file')
    parser.add_argument('-f', '--fasta', help='input fasta file', required=True)
    parser.add_argument('-p', '--prefix', help='prefix to add to fasta headers', required=True)
    parser.add_argument('-o', '--output', help='output file prefix', required=True)
    parser.add_argument('-t', '--type', help='assembly type', choices=['trinity', 'rnaspades', 'fasta', 'evigene'], required=True)
    parser.add_argument('-a', '--protein', help='proteins fasta file', required=False)
    args = parser.parse_args()

    if args.protein:
        prot_dict = defaultdict()
        for prot in parse_fasta(args.protein):
            prot_dict[prot.name] = prot

    fo = open(f"{args.output}.fasta", 'w')
    fm = open(f"{args.output}.map", 'w')
    fgtr = open(f"{args.output}.gene_to_transcript.map", 'w')
    if args.protein:
        fp = open(f"{args.output}.faa", 'w')
    genes = defaultdict(list)
    n_genes = 0
    for fasta in parse_fasta(args.fasta):
        g_id, t_id, t_num = get_ids(fasta.name, args.type)
        if g_id not in genes:
            n_genes += 1
            genes[g_id] = []
            genes[g_id].append(n_genes)
        genes[g_id].append(t_id)
        # get index of g_id in genes and add 1 to get the gene number
        if t_num != -1:
            new_gene_name = f"{args.prefix}_g{genes[g_id][0]}"
            new_isoforms_name = f"{new_gene_name}.i{t_num}"
        else:
            new_gene_name = f"{args.prefix}_g{genes[g_id][0]}"
            new_isoforms_name = f"{new_gene_name}.i{genes[g_id].index(t_id)}"
        fm.write(f"{new_isoforms_name}\t{fasta.name}\n")
        fasta.name = new_isoforms_name
        fo.write(f"{str(fasta)}\n")
        fgtr.write(f"{new_gene_name}\t{new_isoforms_name}\n")
        if args.protein:
            if t_id in prot_dict:
                prot_dict[t_id].name = new_isoforms_name
                fp.write(f"{str(prot_dict[t_id])}\n")
            else:
                print(f"Protein {t_id} not found in protein fasta file", file=sys.stderr)


if __name__ == '__main__':
    main()
        

