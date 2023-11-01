#!/usr/bin/env python3

import argparse
from collections import defaultdict

sep = { 'comma': ',', 'tab': '\t', 'colon': ':', 'semi-colon': ';', 'pipe': '|' }

def join_elements(elements, sep):
    return sep.join(map(str, elements))

def main():
    parser = argparse.ArgumentParser(description='Rearrange TSV file based on order from another list file')
    parser.add_argument('-i', '--in_table', type=str, help='path to input table file')
    parser.add_argument('-l', '--in_list', type=str, help='path to list file')
    parser.add_argument('-s', '--sep', type=str, help='delimiter', choices = {'comma', 'tab', 'colon', 'semi-colon', 'pipe'}, default='tab', required=False)
    parser.add_argument('-o', '--out_table', type=str, help='path to output table file')
    args = parser.parse_args()

    # read list file
    list_order = []
    with open(args.in_list, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue
            if line in list_order:
                print('Error: duplicate entry in list file: ' + line)
                exit(1)
            list_order.append(line.strip())
    f.close()

    print('List order: ' + str(list_order))

    # read table file
    table = defaultdict(list)
    lin_num = 0
    with open(args.in_table, 'r') as f, open(args.out_table, 'w') as fo:
        for line in f:
            line = line.strip().split(sep[args.sep])
            if line[0] == '':
                continue
            lin_num += 1
            if lin_num == 1:
                col_order = line
            olist = []
            for item in list_order:
                
                olist.append(line[col_order.index(item)])
            if lin_num == 1:
                print(join_elements(olist, sep[args.sep]))
            fo.write(join_elements(olist, sep[args.sep]) + '\n')


    

if __name__ == '__main__':
    main()
