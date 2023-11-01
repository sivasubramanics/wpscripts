#!/usr/bin/env python3
# script to compare different orthogroups files and output a table with the number of shared orthogroups
# usage: wp1_compare_orthogroups.py -i orthogroups_1.tsv orthogroups_2.tsv -f orthofinder pantools -o orthogroup_comparison.tsv

import argparse
import time
import sys
import os
from collections    import defaultdict
from datetime       import datetime

class DE(object):
    def __init__(self, gene_id, basemeanA, basemeanB, log2fc, p_value, q_value):
        self.gene_id = gene_id
        self.basemeanA = basemeanA
        self.basemeanB = basemeanB
        self.log2fc = log2fc
        self.p_value = p_value
        self.q_value = q_value
    
        

class Annotation(object):
    """
    class to store annotation information
    """
    def __init__(self, attr_id):
        self.id = attr_id
        self.ann = defaultdict(list)

    def add_ann(self, ann_type, ann):
        if ann_type not in self.ann:
            self.ann[ann_type] = []
        if ',' in ann:
            ann = ann.split(',')
        else:
            ann = [ann]
        for a in ann:
            if a not in self.ann[ann_type]:
                self.ann[ann_type].append(a)

    def __str__(self):
        return f"{self.id}\t{'; '.join(self.ann['desc'])}\t{'; '.join(self.ann['go'])}\t{'; '.join(self.ann['ko'])}\t{'; '.join(self.ann['pfam'])}"

class Fasta(object):
    def __init__(self, name, sequence, description=None):
        self.name = name
        self.sequence = sequence
        self.description = description
    
    def __str__(self):
        return f">{self.name} {self.description}\n{fold_sequence(self.sequence)}"

class OG(object):
    """
    class to store orthogroups information: contains a dictionary of orthogroups with gene ids as values, and a dictionary of annotations with orthogroup ids as keys
    """
    def __init__(self, filename=None, iformat=None, annotations=None, species_fasta=None):
        self.filename = filename
        self.format = iformat
        self.orthogroups = defaultdict(list)
        self.annotations = defaultdict(Annotation)
        self.og_ann = defaultdict(Annotation)
        # self.species_fasta = defaultdict(lambda: defaultdict(Fasta))
        # self.og_dict_genes = defaultdict(lambda: defaultdict(list))
        # self.og_dict_numbers = defaultdict(lambda: defaultdict(int))
        # self.classification = defaultdict(list)
        # self.core_ogs = []
        # self.acc_ogs = []
        # self.scc_ogs = []
        # self.uniq_ogs = []
        if not iformat is None:
            self.orthogroups = self.read_orthogroups()
            self.summary()
        if not annotations is None:    
            self.read_annotations(annotations)
            self.update_og_ann()
        if not species_fasta is None:
            self.read_species_fasta(species_fasta)
            self.og_dict_genes()
            self.og_dict_numbers()
            self.classify_ogs()            

    def read_orthogroups(self):
        """
        switch case to read orthogroups file in different formats
        """
        if self.format == "orthofinder":
            return self.read_orthofinder()
        elif self.format == "pantools":
            return self.read_pantools()
        elif self.format == "blatclusters":
            return self.read_blatclusters()
        elif self.format == "ogtable":
            return self.read_ogtable()
        else:
            print("Error: unknown orthogroups format")
            exit(1)
    
    def read_ogtable(self):
        """"
        read ogtable ortho group file, where the first column is the orthogroup id, and the rest of the columns are the gene ids per species.
        """
        orthogroups = defaultdict(list)
        print_log(f"Reading {self.filename} in ogtable format")
        with open(self.filename) as f:
            for line in f:
                line = line.rstrip().split()
                og = line[0]
                if og not in orthogroups:
                    orthogroups[og] = []
                orthogroups[og] += line[1:]
        return orthogroups


    def read_pantools(self):
        """
        read pantools orthogroups file and return a dictionary of orthogroups with gene ids as values
        file name example: pantools_homology_groups.txt
        """
        orthogroups = defaultdict(list)
        print_log(f"Reading {self.filename} in pantools format")
        warning = False
        with open(self.filename) as f:
            for line in f:
                line = line.rstrip().split()
                if line[0].startswith('#'):
                    continue
                og = line[0].replace(':', '')
                if og.startswith('OG') and warning == False:
                    print_log(f"Warning: Are you sure the file is in pantools format? check the -f option, before continuing with the output")
                    warning = True
                if og not in orthogroups:
                    orthogroups[og] = []
                for i in range(1, len(line)):
                    line[i] = line[i].split('#')[0]
                    orthogroups[og].append(line[i])
        return orthogroups

    def read_orthofinder(self):
        """
        read orthofinder orthogroups file and return a dictionary of orthogroups with gene ids as values
        file name example: Orthogroups.txt (not the Orthogroups.tsv file, because Orthogroups.tsv does not contain singletons)
        """
        orthogroups = defaultdict(list)
        warning = False
        print_log(f"Reading {self.filename} in orthofinder format")
        with open(self.filename) as f:
            for line in f:
                line = line.rstrip().split()
                if line[0].startswith('#'):
                    continue
                og = line[0].replace(':', '')
                if not og.startswith('OG') and warning == False:
                    print_log(f"Warning: Are you sure the file is in orthofinder format? check the -f option, before continuing with the output")
                    warning = True
                if og not in orthogroups:
                    orthogroups[og] = []
                orthogroups[og] += line[1:]
        return orthogroups
    
    def read_blatclusters(self):
        """
        read blatclusters orthogroups file and return a dictionary of orthogroups with gene ids as values
        way to get this file: wp2_test_gfa.py psl_to_gfa -i lettuce.psl -f lettuce.fna -o lettuce.gfa -b 6 > psl2gfa.log 2>&1
        you may use the output cluster file as the input for this script: lettuce.clusters.tsv
        """
        orthogroups = defaultdict(list)
        print_log(f"Reading {self.filename} in blatclusters format")
        warning = False
        with open(self.filename) as f:
            for line in f:
                line = line.rstrip().split()
                if line[0][1] == "#":
                    continue
                og = line[0]
                if og.startswith('OG') and warning == False:
                    print_log(f"Warning: Are you sure the file is in blatclusters format? check the -f option, before continuing with the output")
                    warning = True
                if og not in orthogroups:
                    orthogroups[og] = []
                orthogroups[og].append(line[1])
        return orthogroups

    def orthogroups_to_pairs(self):
        """
        convert orthogroups dictionary to gene pairs dictionary, and singletons list
        """
        orthogroups = self.orthogroups
        gene_pairs = defaultdict(int)
        singletons = []
        for og in orthogroups:
            if len(orthogroups[og]) == 1:
                singletons.append(orthogroups[og][0])
            else:
                for i in range(0, len(orthogroups[og])):
                    for j in range(i+1, len(orthogroups[og])):
                        if (orthogroups[og][i], orthogroups[og][j]) not in gene_pairs:
                            gene_pairs[(orthogroups[og][i], orthogroups[og][j])] = 1
                        else:
                            gene_pairs[(orthogroups[og][i], orthogroups[og][j])] += 1
                        if (orthogroups[og][j], orthogroups[og][i]) not in gene_pairs:
                            gene_pairs[(orthogroups[og][j], orthogroups[og][i])] = 1
                        else:
                            gene_pairs[(orthogroups[og][j], orthogroups[og][i])] += 1
        return gene_pairs, singletons
    
    def read_annotations(self, annotations):
        """
        read annotation files and store them in a dictionary, this function is called only if the annotation file is provided
        this function is written specifically for eggnog.mapper.annotations files
        cmd to get this file: eggnog-mapper -i input.faa -m diamond --itype protein -d /path/to/eggnog/db -o output --cpu 8 --target_taxa 33090 --temp_dir tmp
        file name example: eggnog.mapper.annotations
        """
        self.annotations = defaultdict(Annotation)
        for ann_file in annotations:
            print_log(f"Reading annotation file {ann_file}")
            with open(ann_file) as f:
                for line in f:
                    line = line.rstrip().split("\t")
                    if line[0].startswith('##'):
                        continue
                    if line[0].startswith('#'):
                        ann_header = line
                    attr_id = line[ann_header.index('#query')]
                    if attr_id not in self.annotations:
                        self.annotations[attr_id] = Annotation(attr_id)
                    if line[ann_header.index('Description')] != '-':
                        self.annotations[attr_id].add_ann('desc', line[ann_header.index('Description')])
                    if line[ann_header.index('GOs')] != '-':
                        self.annotations[attr_id].add_ann('go', line[ann_header.index('GOs')])
                    if line[ann_header.index('KEGG_Pathway')] != '-':
                        self.annotations[attr_id].add_ann('kegg', line[ann_header.index('KEGG_Pathway')])
                    if line[ann_header.index('PFAMs')] != '-':
                        self.annotations[attr_id].add_ann('pfam', line[ann_header.index('PFAMs')])

    def update_og_ann(self):
        """
        update the og_ann dictionary with combined annotations of all the genes in an orthogroup
        """
        for og in self.orthogroups:
            anns = []
            for gene in self.orthogroups[og]:
                if gene in self.annotations:
                    anns.append(self.annotations[gene])
            combined_ann = combine_annotations(og, anns)
            self.og_ann[og] = combined_ann

    def summary(self):
        """
        print summary of the orthogroups
        """
        n_genes = len(set([gene for og in self.orthogroups for gene in self.orthogroups[og]]))
        n_orthogroups = len(self.orthogroups)
        n_singletons = len([og for og in self.orthogroups if len(self.orthogroups[og]) == 1])
        outtable = [
            ["File name", self.filename],
            ["Format", self.format],
            ["# genes", n_genes],
            ["# orthogroups", n_orthogroups],
            ["# singletons", f"{n_singletons} ({round(n_singletons/n_orthogroups * 100, 2)}%)" ]
        ]
        tabulate(outtable, header=True)
        return

    def read_species_fasta(self, species_fasta):
        """
        read species fasta file and store them in a dictionary
        """
        self.species_fasta = defaultdict(lambda: defaultdict(Fasta))
        for fasta_file in species_fasta:
            print_log(f"Reading species fasta file {fasta_file}")
            basename = strip_extension(os.path.basename(fasta_file))
            print_log(f"Species name: {basename}")
            if basename not in self.species_fasta:
                self.species_fasta[basename] = defaultdict(Fasta)
            else:
                print_log(f"Error: {basename} already exists in species fasta dictionary, please rename the filename, since the script takes the species name from the filename")
                sys.exit(1)
            
            for fasta in parse_fasta(fasta_file):
                if fasta.name not in self.species_fasta[basename]:
                    self.species_fasta[basename][fasta.name] = fasta
                else:
                    print_log(f"Error: {fasta.name} already exists in {basename} fasta dictionary, please rename the fasta header, since the script takes the gene name from the fasta header")
                    sys.exit(1)

    def og_dict_genes(self):
        """
        returns the dictionary of [og_id][species][list of genes]
        """
        self.og_dict_genes = defaultdict(lambda: defaultdict(list))
        for og in self.orthogroups:
            for gene in self.orthogroups[og]:
                for species in self.species_fasta:
                    if gene in self.species_fasta[species]:
                        self.og_dict_genes[og][species].append(gene)

    def og_dict_numbers(self):
        """
        returns the dictionary of [og_id][species][number of genes]
        """
        self.og_dict_numbers = defaultdict(lambda: defaultdict(int))
        for og in self.orthogroups:
            for gene in self.orthogroups[og]:
                for species in self.species_fasta:
                    if gene in self.species_fasta[species]:
                        self.og_dict_numbers[og][species] += 1
    
    def core_ogs(self):
        """
        returns the list of core genes, meaning the og should have at least one gene from each species
        """
        self.core_ogs = []
        for og in self.og_dict_numbers:
            if len(self.og_dict_numbers[og]) == len(self.species_fasta):
                self.core_ogs.append(og)
        return self.core_ogs
    
    def acc_ogs(self):
        """
        returns the list of accessory genes, meaning the og should have at least one gene from one species
        """
        self.acc_ogs = []
        for og in self.og_dict_numbers:
            if len(self.og_dict_numbers[og]) > 1:
                self.acc_ogs.append(og)
        return self.acc_ogs
    
    def scc_ogs(self):
        """
        returns the list of single copy genes, meaning the og should have only one gene from each species
        """
        self.scc_ogs = []
        for og in self.og_dict_numbers:
            if len(self.og_dict_numbers[og]) == len(self.species_fasta) and sum(self.og_dict_numbers[og].values()) == len(self.species_fasta):
                self.scc_ogs.append(og)
        return self.scc_ogs
    
    def uniq_ogs(self):
        """
        returns the list of unique genes, meaning the og should have only one gene from one species
        """
        self.uniq_ogs = []
        for og in self.og_dict_numbers:
            if len(self.og_dict_numbers[og]) == 1:
                self.uniq_ogs.append(og)
        return self.uniq_ogs
    
    def classify_ogs(self):
        self.classification = defaultdict(list)
        self.classification['core'] = self.core_ogs()
        self.classification['acc'] = self.acc_ogs()
        self.classification['scc'] = self.scc_ogs()
        self.classification['uniq'] = self.uniq_ogs()
        total_ogs = list(set(self.core_ogs + self.acc_ogs + self.scc_ogs + self.uniq_ogs))
        tabulate([[f"# {key}", f"{len(self.classification[key])} ({round(len(self.classification[key])/len(total_ogs) * 100, 2)}%)"] for key in self.classification] + [["# total OGs (species)", len(total_ogs)]] + [["# total OGs (overall)", len(self.orthogroups)]], header=False)

    def write_ogtable(self, out_prefix):
        """
        write orthogroups table in the format og_id\tspecies1\tspecies2\t...
        """
        print_log(f"Writing orthogroups table to {out_prefix}.ogtable")
        if not self.species_fasta:
            print_log("Error: species fasta files not provided, or empty")
            sys.exit(1)
        with open(out_prefix + ".ogtable", 'w') as f:
            f.write("OG\t")
            for species in self.species_fasta:
                f.write(f"{species}\t")
            f.write("\n")
            for og in self.og_dict_genes:
                f.write(f"{og}\t")
                for species in self.species_fasta:
                    if species in self.og_dict_genes[og]:
                        f.write(f"{','.join(self.og_dict_genes[og][species])}\t")
                    else:
                        f.write("\t")
                f.write("\n")

        # write upset count table   
        with open(out_prefix + ".ogtable.counts", 'w') as f:
            f.write("OG\t")
            for species in self.species_fasta:
                f.write(f"{species}\t")
            f.write("total\n")
            for og in self.og_dict_numbers:
                total = 0
                f.write(f"{og}\t")
                for species in self.species_fasta:
                    if species in self.og_dict_numbers[og]:
                        total += self.og_dict_numbers[og][species]
                        f.write(f"{self.og_dict_numbers[og][species]}\t")
                    else:
                        f.write("0\t")
                f.write(f"{total}\n")

class DEG(object):
    def __init__(self, reference, sampleA, sampleB):
        self.reference = reference
        self.sampleA = sampleA
        self.sampleB = sampleB
        self.genes = defaultdict(DE)
    
    def read_dge(self, dge_file, method):
        """
        read DEG output tsv file and update self.genes dictionary
        """
        with open(dge_file) as f:
            for line in f:
                line = line.rstrip().split("\t")
                if line[0].startswith('#'):
                    continue
                if line[0].startswith(''):
                    continue
                if method == "deseq2":
                    # self.genes[gene_id] = DE(gene_id, basemeanA, basemeanB, log2fc, p_value, q_value)
                    if line[0] in self.genes:
                        print_log(f"Error: {line[0]} already exists in {self.genes}")
                        sys.exit(1)
                    self.genes[line[0]] = DE(line[0], line[3], line[4], line[6], line[9], line[10])
                else:
                    print_log(f"Error: unknown method {method}")
                    sys.exit(1)


def write_upset_tsv(og_obj, output):
    """
    function to write upset data in format og_id\tspecies1\tspecies2\t...
    """
    print_log(f"Writing upset data to {output}")
    with open(output, 'w') as f:
        f.write("OG\t")
        for species in og_obj.species_fasta:
            f.write(f"{species}\t")
        f.write("\n")
        for og in og_obj.og_dict_genes:
            f.write(f"{og}\t")
            for species in og_obj.species_fasta:
                if species in og_obj.og_dict_genes[og]:
                    f.write(f"{','.join(og_obj.og_dict_genes[og][species])}\t")
                else:
                    f.write("\t")
            f.write("\n")
    
    print_log(f"Writing upset data to {output}.counts")
    with open(output + ".counts", 'w') as f:
        f.write("OG\t")
        for species in og_obj.species_fasta:
            f.write(f"{species}\t")
        f.write("total\n")
        for og in og_obj.og_dict_numbers:
            total = 0
            f.write(f"{og}\t")
            for species in og_obj.species_fasta:
                if species in og_obj.og_dict_numbers[og]:
                    total += og_obj.og_dict_numbers[og][species]
                    f.write(f"{og_obj.og_dict_numbers[og][species]}\t")
                else:
                    f.write("0\t")
            f.write(f"{total}\n")
        
def parse_fasta(fasta_file):
    """
    Parse a fasta file
    """
    with open(fasta_file, 'r') as f:
        begun = False
        for line in f:
            if line.startswith(">"):
                if begun:
                    yield Fasta(name[1:], sequence, " ".join(description) if description else None)
                name, *description = line.strip().split()
                sequence = ""
                begun = True
            else:
                sequence += line.strip()
                if not sequence:
                    continue
        yield Fasta(name[1:], sequence, " ".join(description) if description else None)

def fold_sequence(seq, width=80):
    """
    fold sequence to a given width
    """
    return '\n'.join(seq[i:i+width] for i in range(0, len(seq), width))

def combine_annotations(og_id, anns):
    """
    combine annotations of all the genes in an orthogroup
    """
    combined_ann = Annotation(og_id)
    for ann in anns:
        for ann_type in ann.ann:
            combined_ann.ann[ann_type] = list(set(combined_ann.ann[ann_type] + ann.ann[ann_type]))
    return combined_ann

def print_log(msg):
    """
    print log messages
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", file=sys.stderr)

def compare_orthogroups(args):
    """
    compare orthogroups files and output a table with the number of shared gene-pairs
    """
    # read in orthogroups files
    og_one = OG(args.input[0], args.format[0])
    og_two = OG(args.input[1], args.format[1])
    genes_pairs_1, singletons_1 = og_one.orthogroups_to_pairs()
    genes_pairs_2, singletons_2 = og_two.orthogroups_to_pairs()
    total_pairs_1 = int(len(genes_pairs_1)/2)
    total_pairs_2 = int(len(genes_pairs_2)/2)
    
    # compare common pairs
    common_pairs = set(genes_pairs_1.keys()) & set(genes_pairs_2.keys())
    true_pairs = int(len(common_pairs)/2)
    # remove common_pairs from both the dictionaries
    for key in common_pairs:
        try:
            del genes_pairs_1[key]
            del genes_pairs_2[key]
            del genes_pairs_1[(key[1], key[0])]
            del genes_pairs_2[(key[1], key[0])]
        except KeyError:
            pass
    
    unique_pairs_1 = genes_pairs_1.keys()
    unique_pairs_2 = genes_pairs_2.keys()

    genes_in_pairs_1 = list(set([key[0] for key in unique_pairs_1] + [key[1] for key in unique_pairs_1]))
    genes_in_pairs_2 = list(set([key[0] for key in unique_pairs_2] + [key[1] for key in unique_pairs_2]))

    total_singletons_1 = len(singletons_1)
    total_singletons_2 = len(singletons_2)

    # compare singletons
    true_singletons = len(list(set(singletons_1) & set(singletons_2)))

    pair_single = list(set(genes_in_pairs_1) & set(singletons_2))
    single_pair = list(set(genes_in_pairs_2) & set(singletons_1))

    out_table = [
        ["", args.input[0], args.input[1]],
        ["Total Pairs", total_pairs_1, total_pairs_2],
        ["True Pairs", f"{true_pairs} ({round(true_pairs/total_pairs_1 * 100, 2)}%)", f"{true_pairs} ({round(true_pairs/total_pairs_2 * 100, 2)}%)"],
        ["Singletons", total_singletons_1, total_singletons_2],
        ["True Singletons", f"{true_singletons} ({round(true_singletons/total_singletons_1 * 100, 2)}%)", f"{true_singletons} ({round(true_singletons/total_singletons_2 * 100, 2)}%)"],
        ["Pair-Single", len(set(pair_single)), len(set(single_pair))]
    ]
    tabulate(out_table, header=True)

def tabulate(itable, header=False):
    """
    print a table with the given list of list of strings
    """
    # get the largest string in the list of list of strings    
    max_len = 0
    for row in itable:
        for item in row:
            if len(str(item)) > max_len:
                max_len = len(str(item))
    
    # add padding to each item in the list of list of strings
    for row in itable:
        for i in range(0, len(row)):
            row[i] = str(row[i]).ljust(max_len)

    # print the table with first row as header
    i = 0
    for row in itable:
        if i == 0:
            print_log('-+-'.join(['-'*max_len]*len(row)))
            print_log(' | '.join(row))
            if header:
                print_log('-+-'.join(['-'*max_len]*len(row)))
        else:
            print_log(' | '.join(row))
        i += 1
    print_log('-+-'.join(['-'*max_len]*len(row)))

def convert_orthogroups(args):
    """
    convert orthofinder/pantools/blatclusters to gene2trans_map format
    """
    og_data = OG(args.input, args.format, args.annotation)
    with open(args.output, 'w') as f:
        f.write(f"#OG\tGene\tDescription\n")
        for og in og_data.orthogroups:
            for gene in sorted(og_data.orthogroups[og]):
                if og_data.annotations:
                    if gene in og_data.annotations:
                        desc = '; '.join(og_data.annotations[gene].ann['desc'])
                        f.write(f"{og}\t{gene}\t{desc}\n")
                    else:
                        f.write(f"{og}\t{gene}\t-\n")
                else:
                    f.write(f"{og}\t{gene}\n")
    
    with open(args.output + ".og.ann", 'w') as f:
        f.write(f"#OG\tDescription\tGO\tKEGG\tPFAM\n")
        for og in og_data.og_ann:
            f.write(f"{og_data.og_ann[og]}\n")

def orthogroups_to_upset(args):
    """
    plot an upset plot of the orthogroups
    """
    og_data = OG(args.input, args.format, species_fasta=args.fasta)
    og_data.write_ogtable(args.output)
    # write_upset_tsv(og_data, args.output)
    print_log(f"Orthogroups table written to {args.output}")

def strip_extension(filename):
    """
    strip the extension of a file
    """
    return os.path.splitext(filename)[0]


def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description='Compare orthogroups files and output a table with the number of shared orthogroups')
    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    compare = subparsers.add_parser('compare_og', help='compare orthogroups files and output a table with the number of shared gene-pairs')
    compare.add_argument('-i', '--input', nargs=2, help='input orthogroups files', required=True)
    compare.add_argument('-f', '--format', nargs=2, help='format of the input orthogroups files', required=True, choices=['orthofinder', 'pantools', 'blatclusters'])

    convert = subparsers.add_parser('convert_og', help='convert orthogroups file to a gene-trans format')
    convert.add_argument('-i', '--input', help='input orthogroups file', required=True)
    convert.add_argument('-f', '--format', help='format of the input orthogroups file', required=True, choices=['orthofinder', 'pantools', 'blatclusters'])
    convert.add_argument('-o', '--output', help='output prefix', required=True)
    convert.add_argument('-a', '--annotation', help='annotation file (eggnog.mapper.annotations)', required=False, nargs='*')


    upset = subparsers.add_parser('upset_og', help='convert input orthogroup file to species specifc genelist and count file')
    upset.add_argument('-i', '--input', help='input orthogroups file', required=True)
    upset.add_argument('-f', '--format', help='format of the input orthogroups file', required=True, choices=['orthofinder', 'pantools', 'blatclusters'])
    upset.add_argument('-o', '--output', help='output prefix', required=True)
    upset.add_argument('-s', '--fasta', help='species fasta file, output will be generated only for the gene_ids present in the provided fasta file(s)', required=True, nargs='*')


    args = parser.parse_args(args=(sys.argv[1:] or ['--help'] or ['-h']))

    if args.command == 'compare_og':
        compare_orthogroups(args)
    if args.command == 'convert_og':
        convert_orthogroups(args)
    if args.command == 'upset_og':
        orthogroups_to_upset(args)

    

    end_time = time.time()
    print_log(f"Total time taken: {round(end_time - start_time, 2)} seconds")

if __name__ == '__main__':
    main()

# EOF