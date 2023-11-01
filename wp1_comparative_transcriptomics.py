#!/usr/bin/env python3


import argparse
import sys
import os
from datetime       import datetime
from collections    import defaultdict


class Annotation(object):
    """
    class to store annotation information
    """
    def __init__(self, gene_id):
        self.gene_id = gene_id
        self.ann = defaultdict(list)

    def add_ann(self, ann_type, ann):
        if ann_type not in self.ann:
            self.ann[ann_type] = []
        if ann_type == 'desc':
            self.ann[ann_type].append(ann)
            return
        if ',' in ann:
            ann = ann.split(',')
        else:
            ann = [ann]
        for a in ann:
            if a not in self.ann[ann_type]:
                self.ann[ann_type].append(a)

    def __str__(self):
        if len(self.ann['desc']) == 0:
            self.ann['desc'].append('-')
        if len(self.ann['go']) == 0:
            self.ann['go'].append('-')
        if len(self.ann['ko']) == 0:
            self.ann['ko'].append('-')
        if len(self.ann['pfam']) == 0:
            self.ann['pfam'].append('-')
        return f"{'; '.join(self.ann['desc'])}\t{'; '.join(self.ann['go'])}\t{'; '.join(self.ann['ko'])}\t{'; '.join(self.ann['pfam'])}"
    
class AnnotationParser(object):
    def __init__(self, reference, annotation_file):
        self.reference = reference
        self.annotation_file = annotation_file
        self.parse_annotation()

    def parse_annotation(self):
        """
        read annotation files and store them in a dictionary, this function is called only if the annotation file is provided
        this function is written specifically for eggnog.mapper.annotations files
        cmd to get this file: eggnog-mapper -i input.faa -m diamond --itype protein -d /path/to/eggnog/db -o output --cpu 8 --target_taxa 33090 --temp_dir tmp
        file name example: eggnog.mapper.annotations
        """
        print_log(f"Reading annotation file: {self.annotation_file}")
        self.annotations = defaultdict(Annotation)
        with open(self.annotation_file) as f:
            for line in f:
                line = line.rstrip().split("\t")
                if line[0].startswith('##'):
                    continue
                if line[0].startswith('#'):
                    ann_header = line
                gene_id = line[ann_header.index('#query')]
                if gene_id not in self.annotations:
                    self.annotations[gene_id] = Annotation(gene_id)
                if line[ann_header.index('Description')] != '-':
                    self.annotations[gene_id].add_ann('desc', line[ann_header.index('Description')])
                if line[ann_header.index('GOs')] != '-':
                    self.annotations[gene_id].add_ann('go', line[ann_header.index('GOs')])
                if line[ann_header.index('KEGG_Pathway')] != '-':
                    self.annotations[gene_id].add_ann('kegg', line[ann_header.index('KEGG_Pathway')])
                if line[ann_header.index('PFAMs')] != '-':
                    self.annotations[gene_id].add_ann('pfam', line[ann_header.index('PFAMs')])

class Fasta(object):
    def __init__(self, name, sequence, description=None):
        self.name = name
        self.sequence = sequence.upper()
        self.description = description
    
    def __str__(self):
        return f">{self.name} {self.description}\n{fold_sequence(self.sequence)}"

class FastaParser(object):
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.parse_fasta()

    def parse_fasta(self):
        self.fasta = defaultdict(Fasta)
        begun = False
        print_log(f"Reading fasta file: {self.fasta_file}")
        with open(self.fasta_file) as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('>'):
                    if begun:
                        if name in self.fasta:
                            print_log(f"ERROR: Duplicate sequence name found: {name}")
                            sys.exit(1)
                        self.fasta[name] = Fasta(name, sequence, description)
                    begun = True
                    name = line.split()[0][1:]
                    description = ' '.join(line.split()[1:])
                    sequence = ''
                else:
                    sequence += line
            if name in self.fasta:
                print_log(f"ERROR: Duplicate sequence name found: {name}")
                sys.exit(1)
            self.fasta[name] = Fasta(name, sequence, description)
        print_log(f"Total sequences read ({self.fasta_file}): {len(self.fasta)}")

class OrthogroupsParser(object):
    def __init__(self, orthogroups_file, iformat, fasta_dict=None, annotation_dict=None):
        self.orthogroups_file = orthogroups_file
        self.iformat = iformat
        self.fasta_dict = fasta_dict
        self.parse_orthogroups()
        self.summary()
        self.annotation_dict = annotation_dict
        if self.annotation_dict is not None:
            self.annotate_orthogroups()
        

    def parse_orthogroups(self):
        if self.iformat == 'orthofinder':
            self.parse_orthofinder()
        elif self.iformat == 'pantools':
            self.parse_pantools()
        elif self.iformat == 'orthotable':
            self.parse_orthotable()
        else:
            print_log(f"ERROR: Unknown orthogroups file format: {self.iformat}")
            sys.exit(1)
        
    def parse_orthofinder(self):
        """
        read orthogroups file and store them in a dictionary
        this function is written specifically for OrthoFinder results
        """
        self.orthogroups = defaultdict(lambda: defaultdict(list))
        if self.fasta_dict is None:
            print_log("ERROR: Fasta dictionary is required for orthofinder format. please provide fasta files")
            sys.exit(1)
        self.references = list(self.fasta_dict.keys())
        print_log(f"Reading orthogroups file: {self.orthogroups_file}")
        with open(self.orthogroups_file) as f:
            for line in f:
                line = line.rstrip().split(": ")
                gene_ids = line[1].split()
                for gene_id in gene_ids:
                    for reference in self.fasta_dict:
                        if gene_id in self.fasta_dict[reference].fasta:
                            self.orthogroups[line[0]][reference].append(gene_id)
                            break

    def parse_pantools(self):
        """
        read orthogroups file and store them in a dictionary
        this function is written specifically for PanTools results
        """
        self.orthogroups = defaultdict(lambda: defaultdict(list))
        self.references = []
        if self.fasta_dict is None:
            print_log("ERROR: Fasta dictionary is required for pantools format. please provide fasta files")
            sys.exit(1)
        with open(self.orthogroups_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip().split()
                line[0] = line[0].replace(':', '')
                for i in range(1, len(line)):
                    line[i] = line[i].split('#')[0]
                    for reference in self.fasta_dict:
                        if reference not in self.references:
                            self.references.append(reference)
                        if line[i] in self.fasta_dict[reference].fasta:
                            self.orthogroups[line[0]][reference].append(line[i])
                            break
    
    def parse_orthotable(self):
        """
        read orthogroups file and store them in a dictionary
        this function is written specifically for OrthoTable results
        the table contains, first column is the orthogroup name, and the rest of the columns are the species names
        """
        self.orthogroups = defaultdict(lambda: defaultdict(list))
        self.references = []
        print_log(f"Reading orthogroups file: {self.orthogroups_file}")
        first_line = True
        n_genes = 0
        og_format = None
        with open(self.orthogroups_file) as f:
            for line in f:
                line = line.rstrip().split("\t")
                if line[0].startswith('#'):
                    og_format = line[0].split('=')[1]
                    continue
                if og_format != 'orthotable':
                    print(f"ERROR: Orthogroups file format is not orthotable\n"
                                    f"ERROR: Please provide orthogroups file in orthotable format\n"
                                    f"ERROR: See {sys.argv[0]} convert --help for more details")
                    sys.exit(1)
                if first_line:
                    first_line = False
                    references = line[1:]
                    self.references.extend([reference for reference in references if reference not in self.references])
                    continue
                for i in range(1, len(line)):
                    genes = []
                    if ',' in line[i]:
                        genes = line[i].split(',')
                    elif line[i] != "":
                        genes = [line[i]]
                    else:
                        continue
                    n_genes += len(genes)
                    self.orthogroups[line[0]][references[i-1]] = genes

    def write_orthogroup_table(self, output_prefix, references=None):
        """
        write orthogroups table
        """
        if references is None:
            references = self.references
        print_log(f"Writing orthogroups table: {output_prefix}.orthogroups.tsv")
        with open(f"{output_prefix}.orthogroups.tsv", 'w') as f:
            f.write(f"# format=orthotable\n")
            f.write(f"Orthogroup")
            for reference in references:
                f.write(f"\t{reference}")
            f.write("\n")
            for orthogroup in self.orthogroups:
                f.write(f"{orthogroup}")
                for reference in references:
                    f.write(f"\t{','.join(self.orthogroups[orthogroup][reference])}")
                f.write("\n")

        print_log(f"Writing orthogroups counts: {output_prefix}.orthogroups.counts")
        with open(f"{output_prefix}.orthogroups.counts", 'w') as f:
            f.write(f"Orthogroup")
            for reference in self.fasta_dict:
                f.write(f"\t{reference}")
            f.write("\ttotal\n")
            for orthogroup in self.orthogroups:
                total = 0
                f.write(f"{orthogroup}")
                for reference in references:
                    total += len(self.orthogroups[orthogroup][reference])
                    f.write(f"\t{len(self.orthogroups[orthogroup][reference])}")
                f.write(f"\t{total}\n")
    
    def annotate_orthogroups(self):
        """
        annotate orthogroups with annotations
        """
        if self.annotation_dict is None:
            raise Exception("Annotation dictionary is required for annotating orthogroups. please provide annotation files")
        print_log(f"Annotating orthogroups")
        self.orthogroups_annotations = defaultdict(Annotation)
        for orthogroup in self.orthogroups:
            og_anns = []
            for reference in self.orthogroups[orthogroup]:
                for gene_id in self.orthogroups[orthogroup][reference]:
                    if reference not in self.annotation_dict:
                        continue
                    if gene_id not in self.annotation_dict[reference].annotations:
                        continue
                    og_anns.append(self.annotation_dict[reference].annotations[gene_id])
            self.orthogroups_annotations[orthogroup] = combine_annotations(orthogroup, og_anns)
    
    def write_og2tr_table(self, output_prefix):
        """
        write orthogroups to transcript table {og}\t{tr}\t{annotation}
        """
        if self.annotation_dict is None:
            raise Exception("Annotation dictionary is required for writing orthogroups to transcript table. please provide annotation files")
        print_log(f"Writing orthogroups to transcript table: {output_prefix}.orthogroups.map.tsv")
        with open(f"{output_prefix}.orthogroups.map.tsv", 'w') as f:
            f.write(f"Orthogroup\tTranscript\tDescription\tGO\tKEGG\tPFAM\n")
            for orthogroup in sorted(self.orthogroups):
                for reference in self.orthogroups[orthogroup]:
                    if reference not in self.annotation_dict:
                        continue
                    for gene_id in self.orthogroups[orthogroup][reference]:
                        f.write(f"{orthogroup}\t{gene_id}")
                        if gene_id in self.annotation_dict[reference].annotations:
                            f.write(f"\t{self.annotation_dict[reference].annotations[gene_id]}\n")
                        else:
                            f.write(f"\t-\t-\t-\t-\n")

        
    def write_orthogroups_annotations(self, output_prefix):
        """
        write orthogroups annotations {og}\t{annotation}
        """
        if self.annotation_dict is None:
            raise Exception("Annotation dictionary is required for writing orthogroups annotations. please provide annotation files")
        print_log(f"Writing orthogroups annotations: {output_prefix}.orthogroups.annotations.tsv")
        with open(f"{output_prefix}.orthogroups.annotations.tsv", 'w') as f:
            f.write(f"Orthogroup\tDescription\tGO\tKEGG\tPFAM\n")
            for orthogroup in sorted(self.orthogroups_annotations):
                f.write(f"{self.orthogroups_annotations[orthogroup]}\n")
    
    def summary(self):
        """
        print summary of the orthogroups
        """
        print_log(f"Summary of orthogroups file: {self.orthogroups_file}")
        num_genes = 0
        num_orthogroups = 0
        num_species = len(self.references)
        core_orthogroups = 0
        accessory_orthogroups = 0
        unique_orthogroups = 0
        single_copy_orthogroups = 0
        for orthogroup in self.orthogroups:
            num_orthogroups += 1
            num_genes_in_orthogroup = 0
            for reference in self.orthogroups[orthogroup]:
                num_genes_in_orthogroup += len(self.orthogroups[orthogroup][reference])
                num_genes += len(self.orthogroups[orthogroup][reference])
            if len(self.orthogroups[orthogroup]) == num_species:
                core_orthogroups += 1
            elif len(self.orthogroups[orthogroup]) == 1:
                unique_orthogroups += 1
            else:
                accessory_orthogroups += 1
            if num_genes_in_orthogroup == num_species and len(self.orthogroups[orthogroup]) == num_species:
                single_copy_orthogroups += 1
        
        outtable = []
        outtable.append(['Orthogroup file', os.path.basename(self.orthogroups_file)])
        outtable.append(['# references', num_species])
        outtable.append(['# genes', num_genes])
        outtable.append(['# orthogroups', num_orthogroups])
        outtable.append(['# core orthogroups', core_orthogroups])
        outtable.append(['# accessory orthogroups', accessory_orthogroups])
        outtable.append(['# unique orthogroups', unique_orthogroups])
        outtable.append(['# single copy orthogroups', single_copy_orthogroups])
        tabulate(outtable, header=True)

    def pairs(self):
        """
        returns a dictionary of gene pairs from all orthogroups where gene_pairs[gene_a][gene_b] = 1, and list of singletons
        """
        gene_pairs = defaultdict(lambda: defaultdict(int))
        singletons = []
        self.og_genes = defaultdict(list)
        for orthogroup in self.orthogroups:
            for reference in self.references:
                self.og_genes[orthogroup] += self.orthogroups[orthogroup][reference]
        
        for orthogroup in self.og_genes:
            if len(self.og_genes[orthogroup]) == 1:
                singletons.append(self.og_genes[orthogroup][0])
                continue
            for i in range(0, len(self.og_genes[orthogroup])):
                for j in range(i+1, len(self.og_genes[orthogroup])):
                    gene_pairs[self.og_genes[orthogroup][i]][self.og_genes[orthogroup][j]] += 1
                    gene_pairs[self.og_genes[orthogroup][j]][self.og_genes[orthogroup][i]] += 1
        return gene_pairs, singletons

class DEG(object):
    def __init__(self, gene_id, log2fc, pvalue, padj):
        self.gene_id = gene_id
        self.log2fc = log2fc
        self.pvalue = pvalue
        self.padj = padj

    def is_true(self, lfc_cutoff, pvalue_cutoff):
        if self.log2fc == 'NA' or self.pvalue == 'NA' or self.padj == 'NA':
            return False
        if self.log2fc == 'Inf' or self.log2fc == '-Inf':
            if float(self.padj) <= pvalue_cutoff:
                return True
            else:
                return False
        if abs(float(self.log2fc)) >= lfc_cutoff and float(self.padj) <= pvalue_cutoff:
            return True
        else:
            return False
    
    def __str__(self):
        return f"{self.log2fc}\t{self.pvalue}\t{self.padj}"

class DEGParser(object):
    def __init__(self, deg_file, iformat, annotation_dict=None):
        self.deg_file = deg_file
        self.iformat = iformat
        self.annotation_dict = annotation_dict
        self.parse_deg()
    
    def parse_deg(self):
        if self.iformat == 'edgeR':
            self.parse_edgeR()
        elif self.iformat == 'deseq2':
            self.parse_deseq2()
        else:
            print_log(f"ERROR: Unknown deg file format: {self.iformat}")
            sys.exit(1)
    
    def parse_deseq2(self):
        """
        read deg file and store them in a dictionary
        this function is written specifically for DESeq2 results
        columns are :""\tsampleA\tsampleB\tbaseMeanA\tbaseMeanB\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj
        """
        self.deg = defaultdict(DEG)
        print_log(f"Reading deg file: {self.deg_file}")
        with open(self.deg_file) as f:
            first_line = True
            for line in f:
                line = line.rstrip().split("\t")
                if line[0].startswith('#'):
                    continue
                if first_line:
                    first_line = False
                    columns = line
                    continue
                gene_id = line[0]
                log2fc = line[columns.index('log2FoldChange')]
                pvalue = line[columns.index('pvalue')]
                padj = line[columns.index('padj')]
                self.deg[gene_id] = DEG(gene_id, log2fc, pvalue, padj)
        
    def write_deg(self, output_prefix, lfc_cutoff=1, pvalue_cutoff=0.05):
        """
        write deg file
        """
        print_log(f"Writing deg file: {output_prefix}.deg.ann.tsv")
        with open(f"{output_prefix}.deg.ann.tsv", 'w') as f:
            f.write(f"Gene\tlog2FoldChange\tpvalue\tpadj")
            if self.annotation_dict is not None:
                f.write(f"\tDescription\tGO\tKEGG\tPFAM\n")
            else:
                f.write("\n")
            for gene_id in self.deg:
                if not self.deg[gene_id].is_true(lfc_cutoff, pvalue_cutoff):
                    continue
                f.write(f"{gene_id}\t{self.deg[gene_id]}")
                if self.annotation_dict is not None:
                    for reference in self.annotation_dict:
                        if gene_id in self.annotation_dict[reference].annotations:
                            f.write(f"\t{self.annotation_dict[reference].annotations[gene_id]}\n")
                            break
                        else:
                            f.write(f"\t-\t-\t-\t-\n")
                else:
                    f.write("\n")
                    

class GeneToTransMap(object):
    """
    class object that holds, isoforms to gene id, or protein to gene id
    """
    def __init__(self, iso_id, gene_id):
        self.gene = gene_id
        self.isoform = iso_id

class GeneToTransMapParser(object):
    def __init__(self,)


def fold_sequence(sequence, width=60):
    return "\n".join(sequence[i:i+width] for i in range(0, len(sequence), width))

def num2human(num):
    """
    convert integer to human readable format. Every 3 digits add a comma
    """
    return "{:,}".format(num)

def tabulate(itable, header=False):
    """
    print a table with the given list of list of strings
    """
    # get the largest string in the list of list of strings    
    max_len = 0
    for row in itable:
        for item in row:
            if type(item) is int or type(item) is float:
                # convert int to human readable format. Every 3 digits add a comma
                itable[itable.index(row)][row.index(item)] = num2human(item)
                len_item = len(num2human(item))
            else:
                len_item = len(str(item))
            if len_item > max_len:
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

def strip_extension(filename):
    """
    strip the extension of a file
    """
    return os.path.splitext(filename)[0]

def compare_orthogroups(og_one, og_two):
    """
    compare orthogroups files and output a table with the number of shared gene-pairs
    """
    genes_pairs_1, singletons_1 = og_one.pairs()
    genes_pairs_2, singletons_2 = og_two.pairs()
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
        ["Orthogroups files", os.path.basename(og_one.orthogroups_file), os.path.basename(og_two.orthogroups_file)],
        ["Total Pairs", total_pairs_1, total_pairs_2],
        ["True Pairs", f"{true_pairs} ({round(true_pairs/total_pairs_1 * 100, 2)}%)", f"{true_pairs} ({round(true_pairs/total_pairs_2 * 100, 2)}%)"],
        ["Singletons", total_singletons_1, total_singletons_2],
        ["True Singletons", f"{true_singletons} ({round(true_singletons/total_singletons_1 * 100, 2)}%)", f"{true_singletons} ({round(true_singletons/total_singletons_2 * 100, 2)}%)"],
        ["Pair-Single", len(set(pair_single)), len(set(single_pair))]
    ]
    tabulate(out_table, header=True)



def main():
    start_time = datetime.now()
    parser = argparse.ArgumentParser(description='Compare orthogroups files and output a table with the number of shared orthogroups')
    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    convert = subparsers.add_parser('convert', help='convert orthogroups file to a table', description='convert orthogroups file to a table')
    convert.add_argument('-i', '--input', help='input orthogroups file', required=True)
    convert.add_argument('-x', '--format', help='input orthogroups file format', choices=['orthofinder', 'pantools', 'orthotable'], required=True)
    convert.add_argument('-f', '--fasta', help='input fasta file(s). Make sure to provide fasta for all the species involved', required=True, nargs='+')
    convert.add_argument('-o', '--output', help='output prefix', required=True)

    og2map = subparsers.add_parser('og2map', help='convert orthogroups file to a gene-transcript map format', description='convert orthogroups file to a gene-transcript map format')
    og2map.add_argument('-i', '--input', help='input orthogroups table file', required=True)
    og2map.add_argument('-a', '--annotation', help='input annotation file(s). Make sure to provide annotation for all the species involved', required=True, nargs='+')
    og2map.add_argument('-o', '--output', help='output prefix', required=True)

    summary = subparsers.add_parser('summary', help='print summary of orthogroups file', description='print summary of orthogroups file')
    summary.add_argument('-i', '--input', help='input orthogroups file', required=True)

    compare = subparsers.add_parser('compare', help='compare orthogroups files and output a table with the number of shared gene pairs', description='compare orthogroups files and output a table with the number of shared gene pairs')
    compare.add_argument('-i', '--input', help='input orthogroups file (in orthotable format)', required=True, nargs=2)

    annotate_deg = subparsers.add_parser('annotate_deg', help='annotate differentially expressed genes', description='annotate differentially expressed genes')
    annotate_deg.add_argument('-i', '--input', help='input deg file', required=True)
    annotate_deg.add_argument('-x', '--format', help='input deg file format', choices=['edgeR', 'deseq2'], required=True)
    annotate_deg.add_argument('-a', '--annotation', help='input annotation file. Make sure to provide annotation for all the species involved', required=True)
    annotate_deg.add_argument('-o', '--output', help='output prefix', required=True)
    annotate_deg.add_argument('-l', '--lfc', help='log2 fold change cutoff', type=float, default=1)
    annotate_deg.add_argument('-p', '--pvalue', help='pvalue cutoff', type=float, default=0.05)
    annotate_deg.add_argument('-m', '--map', help='transcript to gene map file')
    
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if args.command == 'convert':
        fasta_dict = defaultdict(FastaParser)
        for fasta_file in args.fasta:
            reference = strip_extension(os.path.basename(fasta_file))
            if reference in fasta_dict:
                raise Exception(f"Duplicate reference name found: {reference}")
            fasta_dict[reference] = FastaParser(os.path.realpath(fasta_file))
        orthogroups = OrthogroupsParser(os.path.realpath(args.input), args.format, fasta_dict)
        orthogroups.write_orthogroup_table(args.output)
    elif args.command == 'og2map':
        annotation_dict = defaultdict(AnnotationParser)
        for annotation_file in args.annotation:
            reference = strip_extension(strip_extension(os.path.basename(annotation_file)))
            if reference in annotation_dict:
                raise Exception(f"Duplicate reference name found: {reference}")
            annotation_dict[reference] = AnnotationParser(reference, os.path.realpath(annotation_file))
        orthogroups = OrthogroupsParser(os.path.realpath(args.input), 'orthotable', None, annotation_dict)
        orthogroups.write_og2tr_table(args.output)
        orthogroups.write_orthogroups_annotations(args.output)
    elif args.command == 'summary':
        orthogroups = OrthogroupsParser(args.input, 'orthotable')
    elif args.command == 'compare':
        og_one = OrthogroupsParser(args.input[0], 'orthotable')
        og_two = OrthogroupsParser(args.input[1], 'orthotable')
        compare_orthogroups(og_one, og_two)
    elif args.command == 'annotate_deg':
        annotation_dict = defaultdict(AnnotationParser)
        reference = strip_extension(strip_extension(os.path.basename(args.annotation)))
        annotation_dict[reference] = AnnotationParser(reference, os.path.realpath(args.annotation))
        map_dict = defaultdict(str)
        if args.map is not None:
            with open(args.map) as f:
                for line in f:
                    line = line.rstrip().split()
                    map_dict[line[1]] = line[0]
        # add transcript to gene map to annotation dictionary
        deg = DEGParser(args.input, args.format, annotation_dict)
        
        deg.write_deg(args.output, args.lfc, args.pvalue)
    else:
        parser.print_help()
        sys.exit(1)


    print_log(f"Total time taken: {datetime.now() - start_time}")
        
        
        

    

if __name__ == '__main__':
    main()
