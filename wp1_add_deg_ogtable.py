#!/usr/bin/env python3

"""
The script can be used to add DEG results to the Orthogroups table
"""

import argparse
from collections import defaultdict
import logging

class DGE:
    def __init__(self, gene, lrtpval, lfc, padj):
        self.gene = gene
        if lrtpval == 'NA':
            self.lrtpval = 'NA'
        else:
            self.lrtpval = float(lrtpval)
        self.lfc = float(lfc)
        if padj == 'NA':
            self.padj = 'NA'
        else:
            self.padj = float(padj)

    def is_significant(self):
        if self.lrtpval == 'NA':
            return False
        # use global LFC and padj cutoffs
        if self.padj <= PADJ_CUTOFF and abs(self.lfc) >= LFC_CUTOFF and self.lrtpval <= PADJ_CUTOFF:
            return True
        return False

    def __str__(self):
        return f"{self.gene}\t{self.lrtpval}\t{self.lfc}\t{self.padj}"


class OG_stat:
    def __init__(self, og_id):
        self.og_id = og_id
        self.genes = defaultdict(dict)

    def add_gene(self, dge):
        if dge.gene in self.genes:
            raise ValueError(f"Gene {dge.gene} already present in the orthogroup {self.og_id}")
        self.genes[dge.gene] = dge

    def total_genes(self):
        return len(self.genes)

    def significant_genes(self):
        count = 0
        for gene in self.genes:
            if self.genes[gene].is_significant():
                count += 1
        return count

    def upregulated_genes(self):
        count = 0
        for gene in self.genes:
            if self.genes[gene].is_significant() and self.genes[gene].lfc >= LFC_CUTOFF:
                count += 1
        return count

    def downregulated_genes(self):
        count = 0
        for gene in self.genes:
            if self.genes[gene].is_significant() and self.genes[gene].lfc <= -LFC_CUTOFF:
                count += 1
        return count

    def lfc(self):
        lfc = 0
        for gene in self.genes:
            # if self.genes[gene].is_significant():
            #     lfc += self.genes[gene].lfc
            lfc += self.genes[gene].lfc
        return lfc


def parse_deg(deg_file):
    """
    Parse DEG results
    """
    deg = defaultdict(dict)
    with open(deg_file, 'r') as f:
        line_num = 0
        timepoints = []
        for line in f:
            line = line.strip().split('\t')
            line_num += 1
            if line_num == 1:
                header = line
                for i in range(1, len(header)):
                    if header[i].endswith('_Wald_padj'):
                        timepoints.append(header[i].replace('_Wald_padj', ''))
                continue
            else:
                gene = line[0]
                lrtpval = line[1]
                for i in range(1, len(line)):
                    if header[i].endswith('_Wald_padj'):
                        tp = header[i].replace('_Wald_padj', '')
                        deg[gene][tp] = DGE(gene, lrtpval, line[i+1], line[i])
    print(f"Parsed {line_num} lines from DEG file")
    return deg, timepoints


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(description='Add DEG results to Orthogroups table')
    parser.add_argument('-i', '--input', help='Orthogroups table', required=True)
    parser.add_argument('-d', '--deg', help='DEG results', required=True)
    parser.add_argument('-n', '--name', help='Species name present in orthotable to add the DEG', required=True)
    parser.add_argument('-o', '--output', help='Output file prefix', required=True)
    parser.add_argument('-a', '--annotation', help='Annotation file', required=False)
    parser.add_argument('-r', '--rlk', help='RLK/RLP annotation out file', required=False)
    parser.add_argument('-l', '--lfc', help='Log2 fold change cutoff', required=False, default=1, type=float)
    parser.add_argument('-p', '--padj', help='padj cutoff', required=False, default=0.05, type=float)

    args = parser.parse_args()

    # update global LFC and padj cutoffs
    global LFC_CUTOFF
    global PADJ_CUTOFF
    LFC_CUTOFF = args.lfc
    PADJ_CUTOFF = args.padj

    deg, timepoints = parse_deg(args.deg)

    print(f"Timepoints: {timepoints}")
    print(f"no of genes in DEG: {len(deg)}")

    if args.annotation:
        print("Annotation file is provided")
        og_annotation = defaultdict(dict)
        with open(args.annotation, 'r') as f:
            line_num = 0
            for line in f:
                line = line.strip().split('\t')
                line_num += 1
                if line_num == 1:
                    ann_header = line
                    continue
                og_id = line[0]
                for i in range(1, len(line)):
                    og_annotation[og_id][ann_header[i]] = line[i]
        ann_header.remove('GO')
        ann_header.remove('KEGG_ko')
        # ann_header.remove('KEGG_path')


    if args.rlk:
        rlk_dict = defaultdict(dict)
        rlk_tags = []
        with open(args.rlk, 'r') as f:
            line_num = 0
            for line in f:
                line = line.strip().split('\t')
                if line[0].startswith('#'):
                    continue
                line_num += 1
                if line_num == 1:
                    header = line
                    continue
                tag = f"{line[header.index('Receptor')]}_{line[header.index('Type')]}"
                if tag == '-_-':
                    continue
                if tag not in rlk_tags:
                    rlk_tags.append(tag)
                rlk_dict[line[header.index('Protein')]] = tag
        rlk_tags.sort()

    with open(args.input, 'r') as f, open(f"{args.output}.summary", 'w') as o:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#'):
                continue
            line = line.split('\t')
            if line[0] == 'Orthogroup':
                o.write(f"{line[0]}")
                for tp in timepoints:
                    o.write(f"\t{tp}_n_genes\t{tp}_n_sign\t{tp}_n_up\t{tp}_n_down\t{tp}_lfc")
                if args.annotation:
                    for i in range(1, len(ann_header)):
                        o.write(f"\t{ann_header[i]}")

                if args.rlk:
                    for tag in rlk_tags:
                        o.write(f"\t{tag}")
                o.write('\n')
                if args.name not in line:
                    raise ValueError(f"Name {args.name} not found in the orthogroup table")
                name_idx = line.index(args.name)
                header = line
                continue
            og_id = line[0]
            og_data = defaultdict(dict)
            for tp in timepoints:
                og_data[tp] = OG_stat(og_id)

            geneids = line[name_idx].split(',')
            for geneid in geneids:
                if geneid in deg:
                    for tp in timepoints:
                        if tp in deg[geneid]:
                            og_data[tp].add_gene(deg[geneid][tp])

            o.write(f"{og_id}")
            for tp in timepoints:
                o.write(f"\t{og_data[tp].total_genes()}"
                        f"\t{og_data[tp].significant_genes()}"
                        f"\t{og_data[tp].upregulated_genes()}"
                        f"\t{og_data[tp].downregulated_genes()}"
                        f"\t{og_data[tp].lfc():.2f}")
            if args.annotation:
                for i in range(1, len(ann_header)):
                    if ann_header[i] in og_annotation[og_id]:
                        o.write(f"\t{og_annotation[og_id][ann_header[i]]}")
                    else:
                        o.write('\t')
            if args.rlk:
                # initialize rlks with len(rlk_tags) elements with 0
                rlks = [0] * len(rlk_tags)
                for geneid in geneids:
                    if geneid in rlk_dict:
                        tag = rlk_dict[geneid]
                        rlks[rlk_tags.index(tag)] += 1
                for rlk in rlks:
                    o.write(f"\t{rlk}")
            o.write('\n')


if __name__ == "__main__":
    main()



