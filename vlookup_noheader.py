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
parser.add_argument('-c', '--column', type=int, help='Column to lookup (in the Input file)', required=True)
parser.add_argument('-k', '--key', type=int, help='Key column (in the Lookup file)', required=True)
parser.add_argument('-v', '--value', type=int, help='Value column (in the Lookup file)', required=True)
args = parser.parse_args()


def get_dict_lookup(lookup_file, key, value):
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
    with open(lookup_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[key] in lookup_dict:
                print('Duplicate key in the lookup file: ' + line[key])
                sys.exit(1)
            lookup_dict[line[key]] = line[value]
    return lookup_dict



lookup_dict = get_dict_lookup(args.lookup, args.key, args.value)
out_fo = open(args.output, "w")
with open(args.input, 'r') as f:
    for line in f:
        line = line.strip()
        data = line.split('\t')
        if data[args.column] in lookup_dict:
            out_fo.write(f"{line}\t{lookup_dict[data[args.column]]}\n")
        else:
            out_fo.write(f"{line}\tNA\n")

