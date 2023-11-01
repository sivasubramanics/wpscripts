#!/usr/bin/env python3

import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='filter DEG results')
    parser.add_argument('-i', '--in_file', type=str, help='path to input file. eg:isoforms_dge.tsv')
    parser.add_argument('-o', '--out_prefix', type=str, help='out prefix', default='out')
    parser.add_argument('-l', '--log2fc', type=float, help='log2 fold change cutoff', default=1)
    parser.add_argument('-p', '--pval', type=float, help='p-value cutoff', default=0.05)
    parser.add_argument('-c', '--contrast', type=str, help='contrast type', default='treatment')
    parser.add_argument('-s', '--order', type=str, help='sample order file', default=None)
    parser.add_argument('-g', '--genelist', type=str, help='gene list file', default=None)
    args = parser.parse_args()

    deg_dict = defaultdict(dict)
    line_num = 0
    sample_order = []
    with open(args.in_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            line_num += 1
            if line_num == 1:
                fc_col = line.index('log2FoldChange')
                pval_col = line.index('padj')
                continue
            if line[1] != args.contrast:
                continue
            if line[0] not in deg_dict:
                deg_dict[line[0]] = defaultdict(dict)
            if line[2] not in deg_dict[line[0]]:
                deg_dict[line[0]][line[2]] = [round(float(line[fc_col]), 2), round(float(line[pval_col]), 2)]
            if line[2] not in sample_order:
                sample_order.append(line[2])
    f.close()

    gene_list = list(deg_dict.keys())

    if args.genelist != None:
        gene_list = []
        with open(args.genelist, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[0] == '' or line[0] == 'gene_id':
                    continue
                if line[0] not in gene_list:
                    gene_list.append(line[0])
        f.close()

    if args.order != None:
        sample_order = []
        with open(args.order, 'r') as f:
            for line in f:
                line = line.strip()
                if line == '' or line == 'gene_id':
                    continue
                if line not in sample_order:
                    sample_order.append(line)
        f.close()

    out_sig = open(args.out_prefix + '.sig', 'w')
    out_sig.write('gene_id\t' + '\t'.join(sample_order) + '\n')
    out_down = open(args.out_prefix + '.down', 'w')
    out_down.write('gene_id\t' + '\t'.join(sample_order) + '\n')
    out_up = open(args.out_prefix + '.up', 'w')
    out_up.write('gene_id\t' + '\t'.join(sample_order) + '\n')
    out_sig_updown = open(args.out_prefix + '.sig.updown', 'w')
    out_sig_updown.write('gene_id\t' + '\t'.join(sample_order) + '\n')
    out_stats = open(args.out_prefix + '.stats', 'w')
    out_stats.write('sample\tup\tdown\n')
    
    # filter deg_dict
    stats = {}
    for gene in gene_list:
        sig = False
        updown_flag = False
        updown_list = [0]*len(sample_order)
        if gene not in deg_dict:
                continue
        for sample in sample_order:
            if sample not in deg_dict[gene]:
                continue
            deg = deg_dict[gene][sample]
            if deg[1] <= args.pval:
                sig = True
                if deg[0] >= args.log2fc:
                    updown_list[sample_order.index(sample)] = '1'
                    if sample not in stats:
                        stats[sample] = [0, 0]
                    stats[sample][0] += 1
                elif deg[0] <= (args.log2fc * -1):
                    updown_list[sample_order.index(sample)] = '-1'
                    if sample not in stats:
                        stats[sample] = [0, 0]
                    stats[sample][1] += 1
                else:
                    updown_list[sample_order.index(sample)] = '0'
        if sig:
            out_sig.write(gene)
            for sample in sample_order:
                if sample in deg_dict[gene]:
                    out_sig.write('\t' + str(deg_dict[gene][sample][0]))
                else:
                    out_sig.write('\tNA')
            out_sig.write('\n')
        if '1' in updown_list:
            out_up.write(gene)
            for fc in updown_list:
                if fc == '1':
                    out_up.write('\t1')
                else:
                    out_up.write('\t0')
            out_up.write('\n')
            updown_flag = True
        if '-1' in updown_list:
            out_down.write(gene)
            for fc in updown_list:
                if fc == '-1':
                    out_down.write('\t1')
                else:
                    out_down.write('\t0')
            out_down.write('\n')
            updown_flag = True
        if updown_flag:
            out_sig_updown.write(gene)
            for sample in sample_order:
                if sample in deg_dict[gene]:
                    out_sig_updown.write('\t' + str(deg_dict[gene][sample][0]))
                else:
                    out_sig_updown.write('\tNA')
            out_sig_updown.write('\n')
            
            
    
    for sample in sample_order:
        if sample in stats:
            out_stats.write(sample + '\t' + str(stats[sample][0]) + '\t' + str(stats[sample][1]) + '\n')
        else:
            out_stats.write(sample + '\t0\t0\n')



if __name__ == '__main__':
    main()
