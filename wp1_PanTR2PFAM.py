#!/usr/bin/env python3

import argparse
from collections import defaultdict

def parse_pfam_hmm(in_pfam_hmm):
    print(f"{in_pfam_hmm} processing")
    pfam_hmm_dict = {}
    pfam_acc = ""
    pfam_name = ""
    pfam_desc = ""
    with open(in_pfam_hmm, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            pfam_acc = line[0]
            pfam_name = line[1]
            pfam_desc = line[2]
            if pfam_acc not in pfam_hmm_dict:
                pfam_hmm_dict[pfam_acc] = {}
                pfam_hmm_dict[pfam_acc]['name'] = pfam_name
                pfam_hmm_dict[pfam_acc]['desc'] = pfam_desc
    return pfam_hmm_dict

def main():
    parser = argparse.ArgumentParser(description='Add Pfam annotation to PanTR map file')
    parser.add_argument('-i', '--in_map', type=str, help='path to input map file eg: pantr.counts.matrix.map')
    parser.add_argument('-a', '--ann', type=str, help='path to InterPro annotation file eg: isoforms.ann.tsv')
    parser.add_argument('-o', '--out', type=str, help='path to output file')
    args = parser.parse_args()

    pfam_hmm = "/lustre/BIF/nobackup/selva001/work/dbs/Pfam/Pfam-A.hmm.tsv"
    pfam_hmm_dict = parse_pfam_hmm(pfam_hmm)

    ann_map = defaultdict(list)
    with open(args.ann, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[3] != 'Pfam':
                continue
            if line[0] not in ann_map:
                ann_map[line[0]] = []
            ann_map[line[0]].append(line[4])

    gene_map = defaultdict(list)
    with open(args.in_map, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0] not in gene_map:
                gene_map[line[0]] = []
            gene_map[line[0]].append(line[1])

    print('number of genes: ' + str(len(gene_map)))

    gene_ann_map = defaultdict(list)
    with open(args.out, 'w') as fo:
        fo.write('Geneid\tPfam\n')
        for gene in gene_map:
            annotation = []
            for isoform in gene_map[gene]:
                if isoform in ann_map:
                    for pfam_id in ann_map[isoform]:
                        if pfam_id in pfam_hmm_dict:
                            if pfam_hmm_dict[pfam_id]['name'] not in annotation:
                                annotation.append(pfam_hmm_dict[pfam_id]['name'])
            if len(annotation) > 0:
                fo.write(gene + '\t' + ','.join(annotation) + '\n')
            else:
                fo.write(gene + '\t-\n')
    fo.close()

                            
                    
if __name__ == "__main__":
    main()                
    


