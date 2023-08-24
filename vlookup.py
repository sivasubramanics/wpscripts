#!/usr/bin/env python3
# script: vlookup.py
# contact: siva.selvanayagam@wur.nl
# description: vlookup between two files
# input: Input file, Lookup file
# output: Output file
# version: 1.0
# date: 16-05-2023
import sys

import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='vlookup between two files')
parser.add_argument('-i', '--input', help='Input file', required=True)
parser.add_argument('-l', '--lookup', help='Lookup file', required=True)
parser.add_argument('-o', '--output', help='Output file', required=True)
parser.add_argument('-c', '--column', help='Column to lookup (in the Input file)', required=True)
parser.add_argument('-k', '--key', help='Key column (in the Lookup file)', required=True)
parser.add_argument('-v', '--value', help='Value column (in the Lookup file)', required=True)
args = parser.parse_args()


def get_dict_lookup(lookup_df, key, value):
    """
    Convert a lookup dataframe to a dictionary
    Parameters
    ----------
    lookup_df
    key
    value

    Returns
    -------
    lookup_dict : dict
    """
    lookup_dict = defaultdict()
    for index, row in lookup_df.iterrows():
        if row[key] in lookup_dict:
            print('Duplicate key in the lookup file: ' + row[key])
            sys.exit(1)
        lookup_dict[row[key]] = row[value]
    return lookup_dict


def vlookup(in_df, lookup_dict, column):
    """
    Perform a vlookup between two dataframes
    Parameters
    ----------
    in_df
    lookup_df
    key
    value

    Returns
    -------
    out_df : dataframe
    """
    in_df[column] = in_df[column].map(lookup_dict)
    return in_df


in_df = pd.read_csv(args.input, sep='\t', keep_default_na=False)
lookup_df = pd.read_csv(args.lookup, sep='\t')
lookup_dict = get_dict_lookup(lookup_df, args.key, args.value)
out_df = vlookup(in_df, lookup_dict, args.column)
out_df.to_csv(args.output, sep='\t', index=False)
