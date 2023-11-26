#!/usr/bin/env python3
# script to compare orthogroups files and output a table with the number of shared orthogroups
# contact: siva.selvanayagam[at]wur.nl
# data: 2023-10-20
# version: 0.1


import argparse
import sys
import os
import re
from datetime import datetime
from collections import defaultdict

# set gloabal constants
PVALUE_CUTOFF = 0.05
LFC_CUTOFF = 1
EMPTY_ANNOTATION = f'-\t-\t-\t-\t-'
ANNOTATION_HEAD = f"Description\tGO\tKEGG_ko\tKEGG_path\tPFAM"


class Annotation(object):
    """
    class to store annotation information
    """

    def __init__(self, attr_id):
        self.attr_id = attr_id
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
        if len(self.ann['kegg_ko']) == 0:
            self.ann['kegg_ko'].append('-')
        if len(self.ann['kegg_path']) == 0:
            self.ann['kegg_path'].append('-')
        if len(self.ann['pfam']) == 0:
            self.ann['pfam'].append('-')
        return f"{'; '.join(self.ann['desc'])}" \
               f"\t{'; '.join(self.ann['go'])}" \
               f"\t{'; '.join(self.ann['kegg_ko'])}" \
               f"\t{'; '.join(self.ann['kegg_path'])}" \
               f"\t{'; '.join(self.ann['pfam'])}"


class AnnotationParser(object):
    def __init__(self, reference, annotation_file, map_dict=None):
        self.annotations = None
        self.reference = reference
        self.annotation_file = annotation_file
        self.parse_annotation()
        if map_dict is not None:
            self.map_dict = map_dict
            self.gene_annotations()

    def parse_annotation(self):
        """
        read annotation files and store them in a dictionary, this function is called only if the annotation file is provided
        this function is written specifically for eggnog.mapper.annotations files
        cmd to get this file: eggnog-mapper -i input.faa -m diamond --itype protein -d /path/to/eggnog/db -o output --cpu 8 --target_taxa 33090 --temp_dir tmp
        file name example: eggnog.mapper.annotations
        """
        print_log(f"Reading annotation file [{self.reference}]: {self.annotation_file}")
        self.annotations = defaultdict(Annotation)
        with open(self.annotation_file) as f:
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
                    self.annotations[attr_id].add_ann('kegg_path', line[ann_header.index('KEGG_Pathway')])
                if line[ann_header.index('KEGG_ko')] != '-':
                    self.annotations[attr_id].add_ann('kegg_ko', line[ann_header.index('KEGG_ko')])
                if line[ann_header.index('PFAMs')] != '-':
                    self.annotations[attr_id].add_ann('pfam', line[ann_header.index('PFAMs')])


class Fasta(object):
    def __init__(self, name, sequence, description=None):
        self.name = name
        self.sequence = sequence.upper()
        self.description = description

    def __str__(self):
        return f">{self.name} {self.description}\n{fold_sequence(self.sequence)}"


class FastaParser(object):
    def __init__(self, fasta_file):
        self.fasta = None
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
        self.og_attributes = None
        self.orthogroups_annotations = None
        self.references = None
        self.orthogroups = None
        self.orthogroups_file = orthogroups_file
        self.iformat = iformat
        self.fasta_dict = fasta_dict
        self.parse_orthogroups()
        self.summary()
        self.annotation_dict = annotation_dict
        if self.annotation_dict is not None:
            self.annotate_orthogroups()
        if self.fasta_dict is not None:
            self.intersect_orthogroups()

    def intersect_orthogroups(self):
        """
        removed the gene ids from the orthogroups that are not present in the fasta files
        """
        for orthogroup in self.orthogroups:
            for reference in self.orthogroups[orthogroup]:
                for attr_id in self.orthogroups[orthogroup][reference]:
                    if attr_id not in self.fasta_dict[reference].fasta:
                        self.orthogroups[orthogroup][reference].remove(attr_id)

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
                attr_ids = line[1].split()
                for attr_id in attr_ids:
                    for reference in self.fasta_dict:
                        if attr_id in self.fasta_dict[reference].fasta:
                            self.orthogroups[line[0]][reference].append(attr_id)
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
        n_attributes = 0
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
                    attributes = []
                    if ',' in line[i]:
                        attributes = line[i].split(',')
                    elif line[i] != "":
                        attributes = [line[i]]
                    else:
                        continue
                    n_attributes += len(attributes)
                    self.orthogroups[line[0]][references[i - 1]] = attributes

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
            raise Exception(
                "Annotation dictionary is required for annotating orthogroups. please provide annotation files")
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
            raise Exception(
                "Annotation dictionary is required for writing orthogroups to transcript table. please provide "
                "annotation files")
        print_log(f"Writing orthogroups to transcript table: {output_prefix}.orthogroups.map.tsv")
        with open(f"{output_prefix}.orthogroups.map.tsv", 'w') as f:
            f.write(f"Orthogroup\tTranscript\t{ANNOTATION_HEAD}\n")
            for orthogroup in sorted(self.orthogroups):
                for reference in self.orthogroups[orthogroup]:
                    if reference not in self.annotation_dict:
                        continue
                    for attr_id in self.orthogroups[orthogroup][reference]:
                        f.write(f"{orthogroup}\t{attr_id}")
                        if attr_id in self.annotation_dict[reference].annotations:
                            f.write(f"\t{self.annotation_dict[reference].annotations[attr_id]}\n")
                        else:
                            f.write(f"\t{EMPTY_ANNOTATION}\n")

    def write_orthogroups_annotations(self, output_prefix):
        """
        write orthogroups annotations {og}\t{annotation}
        """
        if self.annotation_dict is None:
            raise Exception(
                "Annotation dictionary is required for writing orthogroups annotations. please provide annotation files")
        print_log(f"Writing orthogroups annotations: {output_prefix}.orthogroups.annotations.tsv")
        with open(f"{output_prefix}.orthogroups.annotations.tsv", 'w') as f:
            f.write(f"Orthogroup\t{ANNOTATION_HEAD}\n")
            for orthogroup in sorted(self.orthogroups_annotations):
                f.write(f"{orthogroup}\t{self.orthogroups_annotations[orthogroup]}\n")

    def summary(self):
        """
        print summary of the orthogroups
        """
        print_log(f"Summary of orthogroups file: {self.orthogroups_file}")
        num_attributes = 0
        num_orthogroups = 0
        num_species = len(self.references)
        core_orthogroups = 0
        accessory_orthogroups = 0
        unique_orthogroups = 0
        single_copy_orthogroups = 0
        for orthogroup in self.orthogroups:
            num_orthogroups += 1
            num_attributes_in_orthogroup = 0
            for reference in self.orthogroups[orthogroup]:
                num_attributes_in_orthogroup += len(self.orthogroups[orthogroup][reference])
                num_attributes += len(self.orthogroups[orthogroup][reference])
            if len(self.orthogroups[orthogroup]) == num_species:
                core_orthogroups += 1
            elif len(self.orthogroups[orthogroup]) == 1:
                unique_orthogroups += 1
            else:
                accessory_orthogroups += 1
            if num_attributes_in_orthogroup == num_species and len(self.orthogroups[orthogroup]) == num_species:
                single_copy_orthogroups += 1

        outtable = [['Orthogroup file', os.path.basename(self.orthogroups_file)], ['# references', num_species],
                    ['# genes', num_attributes], ['# orthogroups', num_orthogroups],
                    ['# core orthogroups', core_orthogroups], ['# accessory orthogroups', accessory_orthogroups],
                    ['# unique orthogroups', unique_orthogroups],
                    ['# single copy orthogroups', single_copy_orthogroups]]
        tabulate(outtable, header=True)

    def pairs(self):
        """
        returns a dictionary of gene pairs from all orthogroups where attribute_pairs[gene_a][gene_b] = 1, and list of singletons
        """
        attribute_pairs = defaultdict(lambda: defaultdict(int))
        singletons = []
        self.og_attributes = defaultdict(list)
        for orthogroup in self.orthogroups:
            for reference in self.references:
                self.og_attributes[orthogroup] += self.orthogroups[orthogroup][reference]

        for orthogroup in self.og_attributes:
            if len(self.og_attributes[orthogroup]) == 1:
                singletons.append(self.og_attributes[orthogroup][0])
                continue
            for i in range(0, len(self.og_attributes[orthogroup])):
                for j in range(i + 1, len(self.og_attributes[orthogroup])):
                    attribute_pairs[self.og_attributes[orthogroup][i]][self.og_attributes[orthogroup][j]] += 1
                    attribute_pairs[self.og_attributes[orthogroup][j]][self.og_attributes[orthogroup][i]] += 1
        return attribute_pairs, singletons

    def summarize_DEGs(self, out_prefix, deg_dict, pvlaue_cutoff=PVALUE_CUTOFF, lfc_cutoff=LFC_CUTOFF):
        """
        compare the ogs and degs and write summary for each og in a tsv file for each sample
        """
        for reference in self.references:
            if reference not in deg_dict:
                continue
            diff_dict = defaultdict(lambda: defaultdict(float))
            for sample in deg_dict[reference]:
                fo = open(f"{out_prefix}.{reference}_{sample}.orthogroups.summary.tsv", 'w')
                fo.write(f"Orthogroup\tnum_attributes\tnum_deg\tnum_up\tnum_down\tratio_deg\tratio_up\tattributes"
                         f"\tratio_down\tdiff")
                if self.annotation_dict is not None:
                    fo.write(f"\t{ANNOTATION_HEAD}\n")
                else:
                    fo.write("\n")
                for orthogroup in self.orthogroups:
                    num_attributes = len(self.orthogroups[orthogroup][reference])
                    num_deg = 0
                    num_up = 0
                    num_down = 0
                    for attr_id in self.orthogroups[orthogroup][reference]:
                        if attr_id in deg_dict[reference][sample].deg:
                            if deg_dict[reference][sample].deg[attr_id].is_true(pvlaue_cutoff):
                                num_deg += 1
                                if deg_dict[reference][sample].deg[attr_id].is_up(lfc_cutoff):
                                    num_up += 1
                                elif deg_dict[reference][sample].deg[attr_id].is_down(lfc_cutoff):
                                    num_down += 1
                    if num_deg == 0:
                        continue
                    ratio_deg = round(num_deg / num_attributes, 2)
                    ratio_up = round(num_up / num_attributes, 2)
                    ratio_down = round(num_down / num_attributes, 2)
                    diff = round(ratio_up - ratio_down, 2)
                    fo.write(f"{orthogroup}\t"
                             f"{num_attributes}\t"
                             f"{num_deg}\t"
                             f"{num_up}\t"
                             f"{num_down}\t"
                             f"{ratio_deg}\t"
                             f"{ratio_up}\t"
                             f"{ratio_down}\t"
                             f"{diff}"
                             f"\t{','.join(self.orthogroups[orthogroup][reference])}")
                    if self.annotation_dict is not None:
                        if orthogroup in self.orthogroups_annotations:
                            fo.write(f"\t{self.orthogroups_annotations[orthogroup]}\n")
                        else:
                            fo.write(f"\t{EMPTY_ANNOTATION}\n")
                    else:
                        fo.write("\n")
                    diff_dict[orthogroup][sample] = diff
                fo.close()
            fr = open(f"{out_prefix}.{reference}.orthogroups.up-down.tsv", 'w')
            fr.write(f"Orthogroup")
            for sample in deg_dict[reference]:
                fr.write(f"\t{sample}")
            fr.write("\n")
            if self.orthogroups_annotations is not None:
                fr.write(f"{ANNOTATION_HEAD}\n")
            for orthogroup in diff_dict:
                fr.write(f"{orthogroup}")
                for sample in deg_dict[reference]:
                    fr.write(f"\t{diff_dict[orthogroup][sample]}")
                if orthogroup in self.orthogroups_annotations:
                    fr.write(f"\t{self.orthogroups_annotations[orthogroup]}\n")
                else:
                    fr.write(f"\t{EMPTY_ANNOTATION}\n")
            fr.close()

    def write_orthogroups_tpm(self, output, tpms, metadata=None):
        """
        write orthogroups tpm table
        """
        for reference in self.references:
            print_log(f"Writing orthogroups tpm table: {output}.{reference}.orthogroups.tpm.tsv")
            fo = open(f"{output}.{reference}.orthogroups.tpm.tsv", 'w')
            fo.write(f"Orthogroup")
            if metadata is not None:
                for sample in metadata.data:
                    is_sample = False
                    for replicate in metadata.data[sample].replicates:
                        if replicate not in tpms.samples:
                            continue
                        is_sample = True
                        break
                    if is_sample:
                        fo.write(f"\t{sample}")
            else:
                for sample in tpms.data:
                    fo.write(f"\t{sample}")
            fo.write("\n")
            for orthogroup in self.orthogroups:
                fo.write(f"{orthogroup}")
                if metadata is not None:
                    for sample in metadata.data:
                        is_sample = False
                        for replicate in metadata.data[sample].replicates:
                            if replicate not in tpms.samples:
                                continue
                            is_sample = True
                            break
                        if is_sample:
                            tpm = 0.00
                            for attr_id in self.orthogroups[orthogroup][reference]:
                                tpm += tpms.get_tpm(attr_id, metadata.data[sample].replicates)
                            fo.write(f"\t{round(tpm, 2)}")
                else:
                    for sample in tpms.data:
                        tpm = 0.00
                        for attr_id in self.orthogroups[orthogroup][reference]:
                            tpm += tpms.get_tpm(attr_id, sample)
                        fo.write(f"\t{tpm}")
                fo.write("\n")
            fo.close()


class DEG(object):
    def __init__(self, attr_id, log2fc, pvalue, padj):
        self.attr_id = attr_id
        self.log2fc = log2fc
        self.pvalue = pvalue
        self.padj = padj

    def is_true(self, pvalue_cutoff=PVALUE_CUTOFF):
        if self.padj == 'NA':
            return False
        if float(self.padj) <= pvalue_cutoff:
            return True
        else:
            return False

    def is_up(self, lfc_cutoff=LFC_CUTOFF):
        if self.log2fc == 'NA' or self.pvalue == 'NA' or self.padj == 'NA':
            return False
        if self.log2fc == 'Inf':
            return True
        if self.log2fc == '-Inf':
            return False
        if float(self.log2fc) >= lfc_cutoff:
            return True
        else:
            return False

    def is_down(self, lfc_cutoff=LFC_CUTOFF):
        if self.log2fc == 'NA' or self.pvalue == 'NA' or self.padj == 'NA':
            return False
        if self.log2fc == 'Inf':
            return False
        if self.log2fc == '-Inf':
            return True
        if float(self.log2fc) <= -lfc_cutoff:
            return True
        else:
            return False

    def __str__(self):
        return f"{self.log2fc}\t{self.pvalue}\t{self.padj}"


class DEGParser(object):
    def __init__(self, deg_file, iformat, annotation_dict=None, lfc_cutoff=LFC_CUTOFF, pvalue_cutoff=PVALUE_CUTOFF):
        self.deg = None
        self.deg_file = deg_file
        self.iformat = iformat
        self.annotation_dict = annotation_dict
        self.lfc_cutoff = lfc_cutoff
        self.pvalue_cutoff = pvalue_cutoff
        self.parse_deg()
        self.summary()

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
                attr_id = line[0]
                log2fc = line[columns.index('log2FoldChange')]
                pvalue = line[columns.index('pvalue')]
                padj = line[columns.index('padj')]
                self.deg[attr_id] = DEG(attr_id, log2fc, pvalue, padj)

    def summary(self):
        """
        print summary of the deg file,
        number of genes/isoforms
        number of degs (percentage)
        number of up-regulated genes (percentage)
        number of down-regulated genes (percentage)
        """
        print_log(f"Summary of deg file: {self.deg_file}")
        num_attributes = len(self.deg)
        num_deg = 0
        num_up = 0
        num_down = 0
        for attr_id in self.deg:
            if self.deg[attr_id].is_true(self.pvalue_cutoff):
                num_deg += 1
                if self.deg[attr_id].is_up(self.lfc_cutoff):
                    num_up += 1
                elif self.deg[attr_id].is_down(self.lfc_cutoff):
                    num_down += 1
        percentage_deg = int(num_deg / num_attributes * 100)
        percentage_up = int(num_up / num_attributes * 100)
        percentage_down = int(num_down / num_attributes * 100)
        outtable = [['DEG file', os.path.basename(self.deg_file)],
                    ['# attributes', f"{num2human(num_attributes)}"],
                    ['# DEGs', f"{num2human(num_deg)} ({percentage_deg} %))"],
                    ['# up-regulated', f"{num2human(num_up)} ({percentage_up} %)"],
                    ['# down-regulated', f"{num2human(num_down)} ({percentage_down} %)"]]
        tabulate(outtable, header=True)

    def write_deg(self, output_prefix, lfc_cutoff=1, pvalue_cutoff=0.05):
        """
        write deg file
        """
        print_log(f"Writing deg file: {output_prefix}.deg.ann.tsv")
        with open(f"{output_prefix}.deg.ann.tsv", 'w') as f:
            f.write(f"Gene\tlog2FoldChange\tpvalue\tpadj")
            if self.annotation_dict is not None:
                f.write(f"\t{ANNOTATION_HEAD}\n")
            else:
                f.write("\n")
            for attr_id in self.deg:
                if not self.deg[attr_id].is_true(pvalue_cutoff):
                    continue
                if not self.deg[attr_id].is_up(lfc_cutoff) and not self.deg[attr_id].is_down(lfc_cutoff):
                    continue
                f.write(f"{attr_id}\t{self.deg[attr_id]}")
                if self.annotation_dict is not None:
                    for reference in self.annotation_dict:
                        if attr_id in self.annotation_dict[reference].annotations:
                            f.write(f"\t{self.annotation_dict[reference].annotations[attr_id]}\n")
                            break
                        else:
                            f.write(f"\t{EMPTY_ANNOTATION}\n")
                else:
                    f.write("\n")

    def parse_edgeR(self):
        pass


class GeneToTransMap(object):
    """
    class object that holds, isoforms to gene id, or protein to gene id
    """

    def __init__(self, iso_id, gene_id):
        self.gene = gene_id
        self.isoform = iso_id

    def __str__(self):
        return f"{self.isoform}\t{self.gene}"

    def gene_id(self):
        return self.gene

    def isoform_id(self):
        return self.isoform


class GeneToTransMapParser(object):
    def __init__(self, mapfile):
        self.map = None
        self.mapfile = mapfile
        self.parse_mapfile()

    def parse_mapfile(self):
        """
        read map file and store them in a dictionary
        map file format: isoform_id\tgene_id
        """
        self.map = defaultdict(list)
        print_log(f"Reading map file: {self.mapfile}")
        with open(self.mapfile) as f:
            for line in f:
                line = line.rstrip().split()
                if line[0] not in self.map:
                    # print(f"Map: {line[0]}")
                    self.map[line[0]] = []
                self.map[line[0]].append(GeneToTransMap(line[1], line[0]))

    def get_attrs(self, gene_id):
        if gene_id not in self.map:
            return []
        for gene_map in self.map[gene_id]:
            yield gene_map.isoform_id()


def fold_sequence(sequence, width=60):
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def num2human(num):
    """
    convert integer to human-readable format. Every 3 digits add a comma
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
                # convert int to human-readable format. Every 3 digits add a comma
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
            print_log('-+-'.join(['-' * max_len] * len(row)))
            print_log(' | '.join(row))
            if header:
                print_log('-+-'.join(['-' * max_len] * len(row)))
        else:
            print_log(' | '.join(row))
        i += 1
    print_log('-+-'.join(['-' * max_len] * len(row)))


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
    attribute_pairs_1, singletons_1 = og_one.pairs()
    attribute_pairs_2, singletons_2 = og_two.pairs()
    total_pairs_1 = int(len(attribute_pairs_1) / 2)
    total_pairs_2 = int(len(attribute_pairs_2) / 2)

    # compare common pairs
    common_pairs = set(attribute_pairs_1.keys()) & set(attribute_pairs_2.keys())
    true_pairs = int(len(common_pairs) / 2)
    # remove common_pairs from both the dictionaries
    for key in common_pairs:
        try:
            del attribute_pairs_1[key]
            del attribute_pairs_2[key]
            del attribute_pairs_1[(key[1], key[0])]
            del attribute_pairs_2[(key[1], key[0])]
        except KeyError:
            pass

    unique_pairs_1 = attribute_pairs_1.keys()
    unique_pairs_2 = attribute_pairs_2.keys()

    attributes_in_pairs_1 = list(set([key[0] for key in unique_pairs_1] + [key[1] for key in unique_pairs_1]))
    attributes_in_pairs_2 = list(set([key[0] for key in unique_pairs_2] + [key[1] for key in unique_pairs_2]))

    total_singletons_1 = len(singletons_1)
    total_singletons_2 = len(singletons_2)

    # compare singletons
    true_singletons = len(list(set(singletons_1) & set(singletons_2)))

    pair_single = list(set(attributes_in_pairs_1) & set(singletons_2))
    single_pair = list(set(attributes_in_pairs_2) & set(singletons_1))

    out_table = [
        ["Orthogroups files", os.path.basename(og_one.orthogroups_file), os.path.basename(og_two.orthogroups_file)],
        ["Total Pairs", total_pairs_1, total_pairs_2],
        ["True Pairs", f"{true_pairs} ({round(true_pairs / total_pairs_1 * 100, 2)}%)",
         f"{true_pairs} ({round(true_pairs / total_pairs_2 * 100, 2)}%)"],
        ["Singletons", total_singletons_1, total_singletons_2],
        ["True Singletons", f"{true_singletons} ({round(true_singletons / total_singletons_1 * 100, 2)}%)",
         f"{true_singletons} ({round(true_singletons / total_singletons_2 * 100, 2)}%)"],
        ["Pair-Single", len(set(pair_single)), len(set(single_pair))]
    ]
    tabulate(out_table, header=True)


class Sample(object):
    def __init__(self, sample_name):
        self.name = sample_name
        self.replicates = []
        self.traits = {}

    def add_replicate(self, replicate):
        self.replicates.append(replicate)

    def add_trait(self, trait, value):
        self.traits[trait] = value


class MetadataParser:
    def __init__(self, metadata_file):
        self.metadata_file = metadata_file
        self.data = {}
        self.parse_metadata()

    def parse_metadata(self):
        print_log(f"Reading metadata file: {self.metadata_file}")
        with open(self.metadata_file) as f:
            for i, line in enumerate(f):
                line = line.rstrip().split("\t")
                if i == 0:
                    header = line
                    continue
                sample = line[header.index('name')]
                rep = line[header.index('run')]
                if sample not in self.data:
                    self.data[sample] = Sample(sample)
                self.data[sample].add_replicate(rep)
                for trait in header:
                    if trait in ['name', 'run']:
                        continue
                    self.data[sample].add_trait(trait, line[header.index(trait)])


class TPMParser:
    def __init__(self, tpm_file):
        self.tpm_file = tpm_file
        self.data = defaultdict(lambda: defaultdict(float))
        self.samples = []
        self.parse_tpm()

    def parse_tpm(self):
        print_log(f"Reading TPM file: {self.tpm_file}")
        with open(self.tpm_file) as f:
            for i, line in enumerate(f):
                line = line.rstrip().split("\t")
                if i == 0:
                    header = line
                    self.samples = header[1:]
                    continue
                attr_id = line[0]
                if attr_id not in self.data:
                    self.data[attr_id] = defaultdict(float)
                for sample in header[1:]:
                    self.data[attr_id][sample] = float(line[header.index(sample)])

    def get_tpm(self, attr_id, samples):
        tpm = []
        for sample in samples:
            tpm.append(self.data[attr_id][sample])
        return sum(tpm) / len(tpm)




class CountMatrixParser:
    # TODO: add metadata parser and make changes to the constructor
    def __init__(self, count_matrix_file, metadata_file):
        self.count_matrix_file = count_matrix_file
        self.metadata_file = metadata_file
        self.data = defaultdict(lambda: defaultdict(float))
        self.samples = []
        self.parse_count_matrix()

    def parse_count_matrix(self):
        print_log(f"Reading count matrix file: {self.count_matrix_file}")
        with open(self.count_matrix_file) as f:
            for i, line in enumerate(f):
                line = line.rstrip().split("\t")
                if i == 0:
                    header = line
                    self.samples = header[1:]
                    continue
                attr_id = line[0]
                if attr_id not in self.data:
                    self.data[attr_id] = defaultdict(float)
                for sample in header[1:]:
                    self.data[attr_id][sample] = float(line[header.index(sample)])





def main():
    start_time = datetime.now()
    parser = argparse.ArgumentParser(
        description='Compare orthogroups files and output a table with the number of shared orthogroups')
    subparsers = parser.add_subparsers(help='sub-command help',
                                       dest='command')

    convert = subparsers.add_parser('convert',
                                    help='convert orthogroups file to a table',
                                    description='convert orthogroups file to a table')
    convert.add_argument('-i', '--input',
                         help='input orthogroups file',
                         required=True)
    convert.add_argument('-x', '--format',
                         help='input orthogroups file format',
                         choices=['orthofinder', 'pantools', 'orthotable'],
                         required=True)
    convert.add_argument('-f', '--fasta',
                         help='input fasta file(s). Make sure to provide fasta for all the species involved',
                         required=True,
                         nargs='+')
    convert.add_argument('-o', '--output',
                         help='output prefix',
                         required=True)

    og2map = subparsers.add_parser('og2map',
                                   help='convert orthogroups file to a gene-transcript map format',
                                   description='convert orthogroups file to a gene-transcript map format')
    og2map.add_argument('-i', '--input',
                        help='input orthogroups table file',
                        required=True)
    og2map.add_argument('-a', '--annotation',
                        help='input annotation file(s). Make sure to provide annotation for all the species involved',
                        required=True,
                        nargs='+')
    og2map.add_argument('-o', '--output',
                        help='output prefix',
                        required=True)

    summary = subparsers.add_parser('summary',
                                    help='print summary of orthogroups file',
                                    description='print summary of orthogroups file')
    summary.add_argument('-i', '--input',
                         help='input orthogroups file',
                         required=True)

    compare = subparsers.add_parser('compare',
                                    help='compare orthogroups files and output a table with the number of shared gene '
                                         'pairs',
                                    description='compare orthogroups files and output a table with the number of '
                                                'shared gene pairs')
    compare.add_argument('-i', '--input',
                         help='input orthogroups file (in orthotable format)',
                         required=True,
                         nargs=2)

    annotate_deg = subparsers.add_parser('annotate_deg',
                                         help='annotate differentially expressed genes',
                                         description='annotate differentially expressed genes')
    annotate_deg.add_argument('-i', '--input',
                              help='input deg file',
                              required=True)
    annotate_deg.add_argument('-x', '--format',
                              help='input deg file format',
                              choices=['edgeR', 'deseq2'],
                              required=True)
    annotate_deg.add_argument('-a', '--annotation',
                              help='input annotation file. Make sure to provide annotation for all the species involved',
                              required=True)
    annotate_deg.add_argument('-o', '--output',
                              help='output prefix',
                              required=True)
    annotate_deg.add_argument('-l', '--lfc',
                              help='log2 fold change cutoff',
                              type=float,
                              default=LFC_CUTOFF)
    annotate_deg.add_argument('-p', '--pvalue',
                              help='pvalue cutoff',
                              type=float,
                              default=PVALUE_CUTOFF)
    annotate_deg.add_argument('-m', '--map',
                              help='transcript to gene map file in {gene_id}\\t{transcript_id} format')

    og_deg = subparsers.add_parser('og_deg',
                                   help='compare orthogrouping and DEGs and return summary',
                                   description='compare orthogrouping and DEGs and return summary')
    og_deg.add_argument('-i', '--input',
                        help='input orthogroups file (in orthotable format)',
                        required=True)
    og_deg.add_argument('-d', '--deg',
                        help='input deg file(s). Filename should in {SPECIES}_{SAMPLE}.deg.tsv format',
                        required=True,
                        nargs='+')
    og_deg.add_argument('-o', '--output',
                        help='output prefix',
                        required=True)
    og_deg.add_argument('-a', '--annotation',
                        help='input annotation file(s). Make sure to provide annotation for all the species involved',
                        required=True,
                        nargs='+')
    og_deg.add_argument('-p', '--pvalue',
                        help='pvalue cutoff',
                        type=float,
                        default=PVALUE_CUTOFF)
    og_deg.add_argument('-l', '--lfc',
                        help='log2 fold change cutoff',
                        type=float,
                        default=LFC_CUTOFF)

    og_ann = subparsers.add_parser('og_ann',
                                   help='annotate orthogroups',
                                   description='annotate orthogroups')
    og_ann.add_argument('-i', '--input',
                        help='input orthogroups file (in orthotable format)',
                        required=True)
    og_ann.add_argument('-a', '--annotation',
                        help='input annotation file(s). Make sure to provide annotation for all the species involved',
                        required=True,
                        nargs='+')
    og_ann.add_argument('-o', '--output',
                        help='output prefix',
                        required=True)

    og_tpm = subparsers.add_parser('og_tpm',
                                   help='export TPM for orthogroups from count matrix',
                                   description='export TPM for orthogroups from count matrix')
    og_tpm.add_argument('-i', '--input',
                        help='input orthogroups file (in orthotable format)',
                        required=True)
    og_tpm.add_argument('-c', '--count',
                        help='input count matrix file',
                        required=True)
    og_tpm.add_argument('-m', '--metadata',
                        help='input metadata file')
    og_tpm.add_argument('-o', '--output',
                        help='output prefix',
                        required=True)

    longest_isofrom = subparsers.add_parser('longest_isoform',
                                            help='extract longest isoform from each gene',
                                            description='extract longest isoform from each gene')
    longest_isofrom.add_argument('-i', '--input',
                                 help='input fasta file',
                                 required=True)
    longest_isofrom.add_argument('-o', '--output',
                                 help='output fasta file',
                                 required=True)
    longest_isofrom.add_argument('-m', '--map',
                                 help='transcript to gene map file in {gene_id}\\t{transcript_id} format',
                                 required=True)

    summarize_isoforms = subparsers.add_parser('summarize_isoforms',
                                               help='summarize isoforms',
                                               description='summarize isoforms')
    summarize_isoforms.add_argument('-m', '--map',
                                    help='transcript to gene map file in {gene_id}\\t{transcript_id} format',
                                    required=True)
    summarize_isoforms.add_argument('-c', '--count',
                                    help='input count matrix file',
                                    required=True)
    summarize_isoforms.add_argument('-f', '--fasta',
                                    help='input isoforms fasta file',
                                    required=True)
    summarize_isoforms.add_argument('-p', '--phenotype',
                                    help='input phenotype file',
                                    required=True)
    summarize_isoforms.add_argument('-o', '--output',
                                    help='output prefix',
                                    required=True)

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
        OrthogroupsParser(args.input, 'orthotable')
    elif args.command == 'compare':
        og_one = OrthogroupsParser(args.input[0], 'orthotable')
        og_two = OrthogroupsParser(args.input[1], 'orthotable')
        compare_orthogroups(og_one, og_two)
    elif args.command == 'annotate_deg':
        annotation_dict = defaultdict(AnnotationParser)
        map_dict = defaultdict(GeneToTransMapParser)
        reference = strip_extension(strip_extension(os.path.basename(args.annotation)))
        annotation_dict[reference] = AnnotationParser(reference, os.path.realpath(args.annotation))
        if args.map is not None:
            map_dict[reference] = GeneToTransMapParser(os.path.realpath(args.map))
            annotation_dict[reference] = isoforms_to_gene_ann(annotation_dict[reference], map_dict[reference])
        deg = DEGParser(args.input, args.format, annotation_dict)
        deg.write_deg(args.output, args.lfc, args.pvalue)
    elif args.command == 'og_deg':
        if args.annotation is None:
            orthogroups = OrthogroupsParser(args.input, 'orthotable')
        else:
            annotation_dict = defaultdict(AnnotationParser)
            for annotation_file in args.annotation:
                reference, ann_format = get_names(annotation_file, 'annotation')
                if reference in annotation_dict:
                    raise Exception(f"Duplicate reference name found: {reference}")
                annotation_dict[reference] = AnnotationParser(reference, os.path.realpath(annotation_file))
            orthogroups = OrthogroupsParser(args.input, 'orthotable', None, annotation_dict)
        deg_dict = defaultdict(lambda: defaultdict(DEGParser))
        for deg_file in args.deg:
            reference, sample, deg_method = get_names(deg_file, 'deg')
            deg_dict[reference][sample] = DEGParser(deg_file, deg_method, annotation_dict, args.lfc, args.pvalue)
        orthogroups.summarize_DEGs(args.output, deg_dict, args.pvalue, args.lfc)
    elif args.command == 'og_ann':
        annotation_dict = defaultdict(AnnotationParser)
        for annotation_file in args.annotation:
            reference, ann_format = get_names(annotation_file, 'annotation')
            if reference in annotation_dict:
                raise Exception(f"Duplicate reference name found: {reference}")
            annotation_dict[reference] = AnnotationParser(reference, os.path.realpath(annotation_file))
        orthogroups = OrthogroupsParser(args.input, 'orthotable', None, annotation_dict)
        orthogroups.write_orthogroups_annotations(args.output)
    elif args.command == 'og_tpm':
        orthogroups = OrthogroupsParser(args.input, 'orthotable')
        if args.metadata is not None:
            metadata = MetadataParser(args.metadata)
        else:
            metadata = None
        tpms = TPMParser(args.count)
        orthogroups.write_orthogroups_tpm(args.output, tpms, metadata)
    elif args.command == 'longest_isoform':
        # check if the input file exists
        if os.path.isfile(f'{args.input}.fai'):
            # read the fasta index file
            print_log(f"Reading fasta index file: {args.input}.fai")
            with open(f'{args.input}.fai') as f:
                fasta_index = {}
                for line in f:
                    line = line.rstrip().split("\t")
                    fasta_index[line[0]] = int(line[1])
        else:
            # read fasta file and get length of each sequence
            print_log(f"Reading fasta file: {args.input}")
            fasta_index = {}
            with open(args.input) as f:
                for line in f:
                    line = line.rstrip()
                    if line.startswith('>'):
                        attr_id = line[1:].split()[0]
                        fasta_index[attr_id] = 0
                    else:
                        fasta_index[attr_id] += len(line)

        # read the map file
        map_dict = GeneToTransMapParser(args.map)
        # read the fasta file
        fasta = FastaParser(args.input)
        # write the longest isoform fasta file
        print_log(f"Writing longest isoform fasta file: {args.output}")
        with open(args.output, 'w') as f:
            for gene_id in map_dict.map:
                longest_isoform = None
                longest_isoform_len = 0
                for attr_id in map_dict.get_attrs(gene_id):
                    if attr_id not in fasta_index:
                        continue
                    if fasta_index[attr_id] > longest_isoform_len:
                        longest_isoform_len = fasta_index[attr_id]
                        longest_isoform = attr_id
                if longest_isoform is None:
                    print_log(f"WARNING: No isoform found for gene: {gene_id}")
                    continue
                f.write(f">{longest_isoform}\n{fold_sequence(fasta.fasta[longest_isoform].sequence)}\n")
    elif args.command == 'summarize_isoforms':
        # read the map file
        map_dict = GeneToTransMapParser(args.map)
        # read the fasta file
        fasta = FastaParser(args.fasta)
        # read the count matrix file
        count_matrix = CountMatrixParser(args.count, args.phenotype)

    else:
        parser.print_help()
        sys.exit(1)

    print_log(f"Total time taken: {datetime.now() - start_time}")


def get_names(ifile, itype):
    """
    get reference and sample names from the input file name
    """
    if itype == 'fasta':
        # check if the fasta name is in {reference}.fasta format or {reference}.fa format or {reference}.fna format or {reference}.faa format
        if re.match(r'^.+\.fasta$', os.path.basename(ifile)):
            return os.path.basename(ifile).split('.')[0]
        elif re.match(r'^.+\.fa$', os.path.basename(ifile)):
            return os.path.basename(ifile).split('.')[0]
        elif re.match(r'^.+\.fna$', os.path.basename(ifile)):
            return os.path.basename(ifile).split('.')[0]
        elif re.match(r'^.+\.faa$', os.path.basename(ifile)):
            return os.path.basename(ifile).split('.')[0]
        else:
            print_log(f"ERROR: Input fasta file name should be in {{reference}}.fasta format")
            sys.exit(1)
    elif itype == 'annotation':
        # check if the annotation name is in {reference}.emapper.annotations format
        if re.match(r'^.+\.emapper\.annotations$', os.path.basename(ifile)):
            return [os.path.basename(ifile).split('.')[0], 'emapper']
        else:
            print_log(f"ERROR: Input annotation file name should be in {{reference}}.emapper.annotations format")
            sys.exit(1)
    elif itype == 'deg':
        # check if the deg name is in {reference}_{sample}.{deg_method}.deg.tsv format
        if re.match(r'^.+_.+\.(edgeR|deseq2)\.deg\.tsv$', os.path.basename(ifile)):
            reference, sample = os.path.basename(ifile).split('.')[0].split('_')
            deg_method = os.path.basename(ifile).split('.')[1]
            return [reference, sample, deg_method]
        else:
            print_log(f"ERROR: Input deg file name should be in {{reference}}_{{sample}}.{{dge_method}}.deg.tsv format")
            sys.exit(1)
    else:
        print_log(f"ERROR: Unknown input type: {itype}")
        sys.exit(1)


def isoforms_to_gene_ann(annotation_dict, map_dict):
    """
    change annotations from isofrom to gene, combine the isoforms annotation to gene annotation
    """
    # change annotations from isofrom to gene, combine the isoforms annotation to gene annotation
    for gene_id in map_dict.map:
        print(gene_id)
        for attr_id in map_dict.get_attrs(gene_id):
            ann_list = []
            if attr_id not in annotation_dict.annotations:
                continue
            ann_list.append(annotation_dict.annotations[attr_id])
        annotation_dict.annotations[gene_id] = combine_annotations(gene_id, ann_list)
    return annotation_dict


if __name__ == '__main__':
    main()
