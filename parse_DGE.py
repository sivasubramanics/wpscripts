#!/usr/bin/env python3
# script: parse_DGE.py
# contact: siva.selvanayagam@wur.nl
# description: parse DGE file from runDGE.py
# input: DGE file from runDGE.py
# output: DGE file for a specific contrast/condition
# version: 1.0
# date: 11-05-2023

import pandas as pd
import numpy as np
from collections import defaultdict
import argparse


def extract_contrast_df(in_df, contrast):
    """
    Extract dataframe for a specific condition from DGE file
    Parameters
    ----------
    in_df
    contrast

    Returns
    -------

    """
    df = in_df.loc[in_df['contrast'] == contrast]
    return df


def flatten_dict(nested_dict):
    """
    Flatten a nested dictionary
    Parameters
    ----------
    nested_dict

    Returns
    -------
    res : dict
    """
    res = {}
    if isinstance(nested_dict, dict):
        for k in nested_dict:
            flattened_dict = flatten_dict(nested_dict[k])
            for key, val in flattened_dict.items():
                key = list(key)
                key.insert(0, k)
                res[tuple(key)] = val
    else:
        res[()] = nested_dict
    return res


def dict_to_df(values_dict):
    """
    Convert a nested dictionary to a dataframe
    Parameters
    ----------
    values_dict : dict

    Returns
    -------
    out_df : pd.DataFrame
    """
    flat_dict = flatten_dict(values_dict)
    out_df = pd.DataFrame.from_dict(flat_dict, orient="index")
    out_df.index = pd.MultiIndex.from_tuples(out_df.index)
    out_df = out_df.unstack(level=-1)
    out_df.columns = out_df.columns.map("{0[1]}".format)
    return out_df


def replace_empty_cells(data_frame):
    """
    Replace empty cells with NA
    Parameters
    ----------
    data_frame: pd.DataFrame

    Returns
    -------
    data_frame: pd.DataFrame

    """
    data_frame.replace('', np.nan, inplace=True)
    data_frame.fillna('NA', inplace=True)
    return data_frame


parser = argparse.ArgumentParser(description='Parse DGE file')
parser.add_argument('-i', '--input', dest='in_deg', help='DGE file result of runDGE.py', required=True)
parser.add_argument('-c', '--contrast', dest='contrast', help='condition/contrast e.g. "treatment" or "time"',
                    required=True)
parser.add_argument('-o', '--output', dest='out_prefix', help='output prefix', required=True)
parser.add_argument('-l', '--log2fc', dest='log2fc', help='log2 fold change cutoff', required=False, default=1, type=float)
parser.add_argument('-p', '--pval', dest='pval', help='p-value cutoff', required=False, default=0.05, type=float)
args = parser.parse_args()

in_deg = args.in_deg # eg: /Users/selva001/projects/bremia_rnaseq/DEG/counts/rb_genome/lsal_featurecounts_dge.tsv
contrast = args.contrast
out_prefix = args.out_prefix
log2fc_cutoff = args.log2fc
pval_cutoff = float(args.pval)

# read input file into dataframe
df = pd.read_csv(in_deg, sep='\t')
# extract dataframe for a specific condition/contast
df = extract_contrast_df(df, contrast)

dge_dict = defaultdict()
gene_list = []
sample_list = []
for index, row in df.iterrows():
    gene_id = row['trans_id']
    sample = row['sampleA']
    lfc = row['log2FoldChange']
    pval = row['padj']
    if pval > pval_cutoff:
        continue
    if gene_id not in gene_list:
        gene_list.append(gene_id)
    if sample not in sample_list:
        sample_list.append(sample)
    if gene_id not in dge_dict:
        dge_dict[gene_id] = defaultdict()
    if sample not in dge_dict[gene_id]:
        dge_dict[gene_id][sample] = round(lfc, 2)

# convert nested dictionary to dataframe
dge_df = dict_to_df(dge_dict)
# add gene_id as column name for index column
dge_df.index.name = 'gene_id'
# replace empty cells with NA
dge_df = replace_empty_cells(dge_df)
# write dataframe to file
dge_df.to_csv(out_prefix + '_' + contrast + '.summary', sep='\t')

# open output stats file handle
out_fh = open(out_prefix + '_' + contrast + '.stats', 'w')
# write header
out_fh.write(f"sample,up,down\n")
for (sample_name, sample_values) in dge_df.items():
    up = 0
    down = 0
    for (gene_id, lfc) in sample_values.items():
        if lfc == 'NA':
            continue
        elif lfc >= log2fc_cutoff:
            up = up + 1
        elif lfc <= -log2fc_cutoff:
            down = down + 1
        else:
            continue
    out_fh.write(f"{sample_name},{up},{down}\n")
out_fh.close()
