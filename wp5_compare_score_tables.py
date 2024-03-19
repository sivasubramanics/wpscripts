#!/usr/bin/env python3

import argparse
from collections import defaultdict
import sys
import pandas as pd

def main():

    parser = argparse.ArgumentParser(description='This is a utility script to process IBS scores tables and compare them between the files')
    parser.add_argument('-i', '--input', help='input TSV file', required=True)
    parser.add_argument('-c', '--compare', help='compare TSV file(s)', required=True, nargs='+')
    parser.add_argument('-m', '--map', help='map file', required=True)
    parser.add_argument('-o', '--output', help='output prefix', required=True)
    parser.add_argument('-r', '--reference', help='reference sample name', required=False)
    parser.add_argument('-n', '--names', help='falg to print sample names', action='store_true')
    parser.add_argument('-b', '--binary', help='0=no sample from compare file, 1=sample from compare file (it will overide -n option)', action='store_true')

    args = parser.parse_args()

    tab = '\t'
    comma = ','

    # read map file
    map_dict = defaultdict()
    chr_dict = defaultdict()
    with open(args.map, 'r') as f:
        print(f"Reading {args.map}")
        for line in f:
            line = line.strip().split('\t')
            if line[0] == 'WIN':
                continue
            if float(line[0]) not in map_dict:
                map_dict[float(line[0])] = (line[1], line[2], line[3])
            if line[1] not in chr_dict:
                chr_dict[line[1]] = []
            chr_dict[line[1]].append(line[0])

    n_windows = len(map_dict)

    # read compare input files
    target_dict = defaultdict()
    for file in args.compare:
        print(f"Reading {file}")
        ref_name = file.split('/')[-1].split('.')[0]
        target_dict[ref_name] = defaultdict()
        target_win_count = 0
        with open(file, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[0] == 'WIN':
                    samples = line[1:]
                    continue
                target_win_count += 1
                if line[0] not in target_dict[ref_name]:
                    target_dict[ref_name][line[0]] = defaultdict()
                for i in range(1, len(line)):
                    target_dict[ref_name][line[0]][samples[i-1]] = round(float(line[i]), 3)
        if target_win_count != n_windows:
            print(f"ERROR: {file} has different number of windows than the previous file")
            sys.exit(1)

    # read input file
    input_dict = defaultdict()
    with open(args.input, 'r') as f:
        print(f"Reading {args.input}")
        input_win_count = 0
        for line in f:
            line = line.strip().split('\t')
            if line[0] == 'WIN':
                samples = line[1:]
                for sample in samples:
                    if sample not in input_dict:
                        input_dict[sample] = defaultdict()
                continue
            input_win_count += 1
            for i in range(1, len(line)):
                input_dict[samples[i-1]][line[0]] = round(float(line[i]), 3)

    if target_win_count != input_win_count:
        print(f"ERROR: {args.input} has different number of windows than the compare file(s)")
        sys.exit(1)

    references = list(target_dict.keys())

    # compare
    if args.binary:
        for chr_name in chr_dict:
            for ref_name in references:
                is_first = True
                out_file = args.output + '.' + ref_name + '_' + chr_name + '.IBS.tsv'
                for sample in input_dict:
                    if is_first:
                        # create dataframe with win as first column and index
                        df = pd.DataFrame(chr_dict[chr_name], columns=['WIN'])
                        is_first = False
                    data = []
                    for win in chr_dict[chr_name]:
                        if win in target_dict[ref_name]:
                            tar_samples = []
                            for tar_sample in target_dict[ref_name][win]:
                                if input_dict[sample][win] == target_dict[ref_name][win][tar_sample]:
                                    tar_samples.append(tar_sample)
                            if len(tar_samples) > 0:
                                data.append(1)
                            else:
                                data.append(0)
                    df[sample] = data
                df.to_csv(out_file, sep='\t', index=False)

            for ref_name in references:
                df = pd.read_csv(args.output + '.' + ref_name + '_' + chr_name + '.IBS.tsv', sep='\t')
                # outdf with win as first column and index
                out_df = pd.DataFrame(chr_dict[chr_name], columns=['WIN'])
                for sample in input_dict:
                    data = df[sample].tolist()
                    out_data = []
                    block_no = 0
                    na_count = 0
                    for i in data:
                        if i == 0:
                            na_count += 1
                            out_data.append(0)
                        else:
                            if na_count < 6:
                                # replace last na_count elements to 1 in out_data
                                for j in range(na_count):
                                    idx = len(out_data) - j - 1
                                    out_data[idx] = 1
                            na_count = 0
                            out_data.append(1)
                    out_df[sample] = out_data
                out_df.to_csv(args.output + '.' + ref_name + '_' + chr_name + '.IBS.blocks.tsv', sep='\t', index=False)

            for ref_name in references:
                df = pd.read_csv(args.output + '.' + ref_name + '_' + chr_name + '.IBS.blocks.tsv', sep='\t')
                out_ibs_tsv = args.output + '.' + ref_name + '_' + chr_name + '.IBS.regions.tsv'
                f = open(out_ibs_tsv, 'w')
                f.write(f"sample\tseqname\tstart\tend\tlength\n")
                wins = df['WIN'].tolist()
                for sample in input_dict:
                    blocks = df[sample].tolist()
                    block_num = 0
                    started = False
                    for win in wins:
                        if blocks[wins.index(win)] == 1:
                            block_num += 1
                            if not started:
                                started = True
                                start = map_dict[win][1]
                            end = map_dict[win][2]
                        if blocks[wins.index(win)] == 0:
                            if started:
                                started = False
                                f.write(f"{sample}\t{chr_name}\t{start}\t{end}\t{int(end)-int(start)+1}\n")
                            block_num = 0
                    if started:
                        f.write(f"{sample}\t{chr_name}\t{start}\t{end}\t{int(end)-int(start)+1}\n")
                f.close()
    else:
        for sample in input_dict:
            out_file = args.output + '.' + sample + '.comp.tsv'
            with open(out_file, 'w') as f:
                f.write(f"WIN\t{join_list(references, tab)}\n")
                for win in input_dict[sample]:
                    f.write(f"{win}")
                    for reference in references:
                        tar_samples = []
                        if win in target_dict[reference]:
                            for tar_sample in target_dict[reference][win]:
                                if input_dict[sample][win] == target_dict[reference][win][tar_sample]:
                                    tar_samples.append(tar_sample)
                        f.write(f"\t{len(tar_samples)}")
                        if args.names:
                            f.write(f":{join_list(tar_samples, comma)}")
                    f.write('\n')

def join_list(ilist, sep):
    return sep.join(ilist)





if __name__ == '__main__':
    main()