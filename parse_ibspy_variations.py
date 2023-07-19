#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from datetime import datetime
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# python combine_variations.py -i data/metadata.csv -w 50000 -o test -c variations -a LK087 -b LK453 -q LK079 -s 0.2

def plot_histogram(dataframe, column_name, filename):
    """
    Plot histogram of a column in a dataframe
    :param dataframe:
    :param column_name:
    :param filename:
    :return:
    """
    plt.hist(dataframe[column_name])
    plt.xlabel(column_name)
    plt.ylabel('Frequency')
    plt.title('Histogram of ' + column_name)
    plt.savefig(filename + ".png")


def add_variations(in_var_file, in_sample_name, dict_var, column):
    """
    Add variations to a dictionary
    :param in_var_file:
    :param in_sample_name:
    :param dict_var:
    :param column:
    :return:
    """
    with open(in_var_file) as f:
        for line in f:
            line = line.strip()
            data = line.split('\t')
            if data[0] == 'seqname':
                col_idx = data.index(column)
                continue
            else:
                seqname = data[0]
                start = data[1]
                end = data[2]
                variations = data[col_idx]
                # key = seqname + '_' + start + '_' + end
                if (seqname, start, end) not in dict_var:
                    dict_var[(seqname, start, end)] = defaultdict()
                dict_var[(seqname, start, end)][in_sample_name] = variations
    return dict_var


def main():
    parser = argparse.ArgumentParser("Combine multiple IBSpy output files by window & reference")
    parser.add_argument('-i', '--input', help='File containing the list of variations tsv to combine', required=True)
    parser.add_argument('-w', '--window_size', type=int, help='Window size to combine the variations', required=True)
    parser.add_argument('-o', '--output', help='Output file prefix', required=True)
    parser.add_argument('-c', '--column', help='Columns to combine the variations', required=True)
    parser.add_argument('-a', '--ref_a', help='Recipient genome accession', required=False)
    parser.add_argument('-b', '--ref_b', help='Donor genome accession', required=False)
    parser.add_argument('-q', '--query', help='Query accession', required=False)
    parser.add_argument('-s', '--sd', help='Std Dev', type=float, required=False)
    args = parser.parse_args()

    start_time = datetime.now()
    dict_var = defaultdict()
    # add vaitaions from each files to a dictionary
    with open(args.input) as f:
        for line in f:
            line = line.strip().split(',')
            in_var_file = line[1]
            in_sample_name = line[0]
            print(f"Processing {in_var_file}")
            dict_var = add_variations(in_var_file, in_sample_name, dict_var, args.column)

    print(f"Writing output file {args.output}.combined.tsv")
    with open(args.output + '.combined.tsv', 'w') as f:
        f.write('seqname\tstart\tend\t' + '\t'.join(dict_var[list(dict_var.keys())[0]].keys()) + '\n')
        for k, v in dict_var.items():
            f.write('\t'.join(k) + '\t' + '\t'.join(v.values()) + '\n')


    print(f"Writing output file for the heatmap {args.output}.combined.heatmap.tsv")
    with open(args.output + '.combined.heatmap.tsv', 'w') as f:
        f.write('seqname\twindow\t' + '\t'.join(dict_var[list(dict_var.keys())[0]].keys()) + '\n')
        i = 0
        p_chr = ""
        for k, v in dict_var.items():
            i += 1
            if k[0] != p_chr:
                i = 1
            f.write(k[0] + '\t' + str(i) + '\t' + '\t'.join(v.values()) + '\n')
            p_chr = k[0]


    df = pd.read_csv(args.output + '.combined.tsv', delimiter='\t')
    print(df.head())
    out_df = df
    out_df.insert(3, 'min', df.iloc[:, 3:].min(axis=1))
    out_df.to_csv(args.output + '.combined.with_min.tsv', sep='\t', index=False)

    plot_histogram(out_df, 'min', args.output + '.combined.with_min.hist')

    if args.ref_a and args.ref_b and args.query:
        df = pd.read_csv(args.output + '.combined.tsv', delimiter='\t')
        out_df = df[['seqname', 'start', 'end', args.ref_a, args.query, args.ref_b]]
        out_df['aq'] = round(out_df[args.ref_a] / out_df[args.query], 2)
        out_df['bq'] = round(out_df[args.ref_b] / out_df[args.query], 2)
        out_df['introgression'] = np.where((out_df['aq'] <= 1 + args.sd) & (out_df['bq'] <= 1 + args.sd), 'b',
                                           np.where((out_df['aq'] >= 1 - args.sd) & (out_df['bq'] >= 1 - args.sd), 'a',
                                                    'n'))
        out_df['ibd'] = np.where((out_df['bq'] <= 1 + args.sd) & (out_df['bq'] >= 1 - args.sd), 'yes', 'no')
        print(f"Writing output file {args.output}.{args.ref_a}_{args.query}_{args.ref_b}.tsv")
        out_df.to_csv(args.output + '.' + args.ref_a + '_' + args.query + '_' + args.ref_b + '.tsv', sep='\t',
                      index=False)

    print(f"Done! Time elapsed: {datetime.now() - start_time}")


if __name__ == '__main__':
    main()
