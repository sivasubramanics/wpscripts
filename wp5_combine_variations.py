#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from datetime import datetime
from collections import defaultdict
import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
# from plotnine import *


# python combine_variations.py -i data/metadata.csv -w 50000 -o test -c variations -a LK087 -b LK453 -q LK079 -s 0.2
class Window:
    def __init__(self, chr_name, start, end):
        self.chr_name = chr_name
        self.start = start
        self.end = end
        self.total_kmers = 0
        self.variations = defaultdict()
        self.observed_kmers = defaultdict()
        self.blocks = defaultdict()

    def get_sample_names(self):
        return list(self.variations.keys())

    def set_sample_names(self, sample_names):
        for sample_name in sample_names:
            if sample_name not in self.variations:
                self.variations[sample_name] = 0
            if sample_name not in self.observed_kmers:
                self.observed_kmers[sample_name] = 0

    def set_total_kmers(self, total_kmers):
        self.total_kmers = int(total_kmers)

    def get_total_kmers(self):
        return self.total_kmers

    def add_total_kmers(self, total_kmers):
        self.total_kmers += int(total_kmers)

    def set_variation(self, sample_name, variations):
        if sample_name not in self.variations:
            self.variations[sample_name] = 0
        self.variations[sample_name] = variations

    def get_variation(self, sample_name):
        return self.variations[sample_name]

    def add_variation(self, sample_name, variations):
        if sample_name not in self.variations:
            self.set_variation(sample_name, 0)
        self.variations[sample_name] += int(variations)

    def set_observed_kmer(self, sample_name, observed_kmers):
        if sample_name not in self.observed_kmers:
            self.observed_kmers[sample_name] = 0
        self.observed_kmers[sample_name] = int(observed_kmers)

    def get_observed_kmer(self, sample_name):
        return self.observed_kmers[sample_name]

    def add_observed_kmer(self, sample_name, observed_kmers):
        if sample_name not in self.observed_kmers:
            self.set_observed_kmer(sample_name, 0)
        self.observed_kmers[sample_name] += int(observed_kmers)

    def get_variations(self, sample_list):
        variations = []
        for sample in sample_list:
            variations.append(self.variations[sample])
        return variations

    def set_variations(self, sample_list, variations_list):
        for sample, variation in zip(sample_list, variations_list):
            if sample not in self.variations:
                self.variations[sample] = 0
            self.variations[sample] = int(variation)

    def add_variations(self, sample_list, variations_list):
        for sample, variation in zip(sample_list, variations_list):
            self.variations[sample] += int(variation)

    def get_observed_kmers(self, sample_list):
        observed_kmers = []
        for sample in sample_list:
            observed_kmers.append(self.observed_kmers[sample])
        return observed_kmers

    def set_observed_kmers(self, sample_list, observed_kmers_list):
        for sample, observed_kmers in zip(sample_list, observed_kmers_list):
            if sample not in self.observed_kmers:
                self.observed_kmers[sample] = 0
            self.observed_kmers[sample] = int(observed_kmers)

    def add_observed_kmers(self, sample_list, observed_kmers_list):
        for sample, observed_kmers in zip(sample_list, observed_kmers_list):
            if sample not in self.observed_kmers:
                self.observed_kmers[sample] = 0
            self.observed_kmers[sample] += int(observed_kmers)

    def get_all_variations(self):
        return self.variations.keys(), self.variations.values()

    def get_all_observed_kmers(self):
        return self.observed_kmers.keys(), self.observed_kmers.values()

    def print(self):
        print(self.chr_name, self.start, self.end, self.total_kmers, self.variations)

    def min_variations(self):
        return min(self.variations.values())

    def max_variations(self):
        return max(self.variations.values())

    def min_observed_kmers(self):
        return min(self.observed_kmers.values())

    def max_observed_kmers(self):
        return max(self.observed_kmers.values())

    def __str__(self):
        return self.chr_name + '\t' + str(self.start) + '\t' + str(self.end)

    def __repr__(self):
        return self.chr_name + '\t' + str(self.start) + '\t' + str(self.end)

    def get_block(self, sample_name):
        return self.blocks[sample_name]

    def set_block(self, sample_name, block):
        if sample_name not in self.blocks:
            self.blocks[sample_name] = 'NA'
        self.blocks[sample_name] = block

    def get_blocks(self, sample_list):
        blocks = []
        for sample in sample_list:
            blocks.append(self.blocks[sample])
        return blocks

    def set_blocks(self, sample_list, blocks_list):
        for sample, block in zip(sample_list, blocks_list):
            if sample not in self.blocks:
                self.blocks[sample] = 'NA'
            self.blocks[sample] = block

    def get_all_blocks(self):
        return self.blocks.keys(), self.blocks.values()
    
    def get_score(self, sample):
        return round((self.observed_kmers[sample]-self.variations[sample])/self.total_kmers, 2)

    def write_window(self, sample):
        return (f"{str(self.variations[sample])}"
                f":{str(self.observed_kmers[sample])}"
                f":{str(self.get_score(sample))}")


class Block:
    def __init__(self, chr_name, start, end):
        self.chr_name = chr_name
        self.start = int(start)
        self.end = int(end)
        self.block_nums = []
        self.windows = []
        self.variations = []
        self.observed_kmers = []

    def add_block_num(self, block_num):
        self.block_nums.append(str(block_num))

    def get_block_nums(self):
        self.trim_blocks_for_na()
        return self.block_nums

    def get_start(self):
        return self.start

    def set_start(self, start):
        if self.start > int(start):
            self.start = int(start)

    def get_end(self):
        return self.end

    def set_end(self, end):
        if self.end < int(end):
            self.end = int(end)

    def get_length(self):
        return self.end - self.start

    def get_na_count(self):
        self.trim_blocks_for_na()
        na_count = 0
        for block_num in self.block_nums:
            if block_num == 'F':
                na_count += 1
        return na_count

    def get_window_count(self):
        self.trim_blocks_for_na()
        return len(self.block_nums)

    def trim_blocks_for_na(self):
        for i in reversed(range(len(self.block_nums))):
            if self.block_nums[i] == 'F':
                self.block_nums.pop(i)
            else:
                break

    def get_variations(self):
        variations = []
        for i in range(len(self.get_block_nums())):
            variations.append(self.variations[i])
        self.variations = variations
        return self.variations

    def add_variation(self, variation):
        self.variations.append(str(variation))

    def get_windows(self):
        windows = []
        for i in range(len(self.get_block_nums())):
            windows.append(self.windows[i])
        self.windows = windows
        return self.windows

    def add_window(self, window):
        self.windows.append(str(window))

    def add_observed_kmer(self, observed_kmer):
        self.observed_kmers.append(str(observed_kmer))

    def get_observed_kmers(self):
        observed_kmers = []
        for i in range(len(self.get_block_nums())):
            observed_kmers.append(self.observed_kmers[i])
        self.observed_kmers = observed_kmers
        return self.observed_kmers

    def __repr__(self):
        return (f"{str(self.chr_name)}"
                f"\t{str(self.start)}"
                f"\t{str(self.end)}"
                f"\t{str(self.get_length())}"
                f"\t{str(self.get_window_count())}"
                f"\t{str(self.get_na_count())}"
                f"\t{str(','.join(self.get_block_nums()))}"
                f"\t{str(','.join(self.get_windows()))}"
                f"\t{str(','.join(self.get_variations()))}"
                f"\t{str(','.join(self.get_observed_kmers()))}")


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


def convert_windows(windows, chr_lenths, o_win_size, n_win_size, kmer_size):
    """

    :param output:
    :param windows:
    :param chr_lenths:
    :param window_size:
    :return:
    """
    new_windows = defaultdict(dict)
    for chr in chr_lenths:
        n_win_start = 0
        win_starts, win_ends = get_window(0, chr_lenths[chr], n_win_size, kmer_size)
        for win_start, win_end in zip(win_starts, win_ends):
            new_window = Window(chr, win_start, win_end)
            new_key = (str(chr), str(win_start), str(win_end))
            o_win_starts, o_win_ends = get_window(n_win_start, win_end, o_win_size, kmer_size)
            for o_win_start, o_win_end in zip(o_win_starts, o_win_ends):
                # print(win_start, win_end, o_win_start, o_win_end)
                key = (str(chr), str(o_win_start), str(o_win_end))
                if key not in windows:
                    n_win_start = o_win_start
                    continue
                old_window = windows[key]
                samples = old_window.get_sample_names()
                new_window.set_sample_names(samples)
                new_window.add_observed_kmers(samples, old_window.get_observed_kmers(samples))
                new_window.add_total_kmers(old_window.get_total_kmers())
                new_window.add_variations(samples, old_window.get_variations(samples))
            new_windows[new_key] = new_window
    return new_windows


def plot_variations(windows, chr_name, window_size, sample, prefix):
    f = open(prefix + "_" + sample + "_" + chr_name + ".tsv", "w")
    f.write(f"window\tmin\t{sample}\n")
    for key in windows:
        window = windows[key]
        if window.chr_name == chr_name:
            f.write(f"{str(window.end)}\t{str(window.min_variations())}\t{str(window.variations[sample])}\n")
    f.close()
    df = pd.read_csv(prefix + "_" + sample + "_" + chr_name + ".tsv", sep='\t')
    df_long = df.melt(id_vars='window', value_vars=[sample, 'min'], var_name='variable', value_name='variations')
    df_long['color'] = np.where(df_long['variations'] <= 60, 'red', 'blue')
    df_long['variations'] = np.log2(df_long['variations'])

    # Create the plot
    plt = (ggplot(df_long)
           + aes(x='window', y='variations', color='color')
           + geom_point()
           + facet_wrap('~ variable', ncol=1)
           + ggtitle('window vs LK001 and min')
           + theme(figure_size=(20, 5)))
    plt.save(prefix + "_" + sample + "_" + chr_name + ".svg")


def find_IBS_ref(windows, chr_name, output):
    f = open(output + "_" + chr_name + ".tsv", "w")
    f.write("window\tmin\tmax\n")
    itr = 0
    for key in windows:
        window = windows[key]
        if window.chr_name == chr_name:
            f.write(f"{str(window.end)}\t{str(window.min_variations())}\t{str(window.max_variations())}\n")
            itr += 1
    f.close()


def get_chr_windows(windows):
    chr_windows = defaultdict()
    for key in windows:
        window = windows[key]
        if not window.chr_name in chr_windows:
            chr_windows[window.chr_name] = defaultdict()
        chr_windows[window.chr_name][key] = window
    return chr_windows


def main():
    parser = argparse.ArgumentParser("Combine multiple IBSpy output files by window & reference")
    parser.add_argument('-i', '--input', help='File containing the list of variations tsv to combine', required=True)
    parser.add_argument('-w', '--window_size', type=int, help='Window size to combine the variations', required=False)
    parser.add_argument('-k', '--kmer_size', type=int, help='Kmer size', required=False)
    parser.add_argument('-o', '--output', help='Output file prefix', required=True)
    parser.add_argument('-t', '--task', help='task', choices=['combine_variations', 'find_IBS_cross', 'find_IBS_ref',
                                                              'plot_heatmap', 'plot_hist', 'convert_windows',
                                                              'plot_variations'],
                        required=False)
    parser.add_argument('-l', '--length', help='Length of the chromosome lengths', required=False)
    parser.add_argument('-a', '--ref_a', help='Recipient genome accession', required=False)
    parser.add_argument('-b', '--ref_b', help='Donor genome accession', required=False)
    parser.add_argument('-q', '--query', help='Query accession', required=False)
    parser.add_argument('-c', '--cut_off', type=int, help='Cut off for the number of variations', required=False,
                        default=10)
    parser.add_argument('-s', '--step_size', type=int, help='Step size to include windows within the given concordant windows', required=False, default=10)
    args = parser.parse_args()

    start_time = datetime.now()

    print(f"Reading chromosome lengths from {args.length}")
    chr_lenths = get_chr_length(args.length)

    if args.task == 'combine_variations':
        task_combine_variations(args, chr_lenths)

    if args.task == 'convert_windows':
        task_convert_windows(args, chr_lenths)

    if args.task == 'plot_variations':
        task_plot_variations(args, chr_lenths)

    if args.task == 'find_IBS_ref':
        if not os.path.exists(args.output + ".variations.tsv"):
            task_combine_variations(args, chr_lenths)
        windows, current_window_size, samples = read_variations(args.output, chr_lenths)
        if args.window_size != current_window_size:
            prefix = args.output + ".win" + str(args.window_size)
            if not os.path.exists(prefix + ".variations.tsv"):
                print(f"Converting windows from {current_window_size} to {args.window_size}")
                task_convert_windows(args, chr_lenths)
            windows, current_window_size, samples = read_variations(prefix, chr_lenths)
        windows = get_chr_windows(windows)
        for_plot = defaultdict()
        for sample in samples:
            print(f"Finding IBS for {sample}")
            block_no = 0
            win_count = 0
            for chr_name in chr_lenths:
                itr = 0
                na = 0
                first_in_chromosome = True
                for key in windows[chr_name]:
                    win_count += 1
                    if str(win_count) not in for_plot:
                        for_plot[str(win_count)] = defaultdict()
                    if sample not in for_plot[str(win_count)]:
                        for_plot[str(win_count)][sample] = [0, 'F']
                    if 'chr' not in for_plot[str(win_count)]:
                        for_plot[str(win_count)]['chr'] = chr_name
                    if 'start' not in for_plot[str(win_count)]:
                        for_plot[str(win_count)]['start'] = windows[chr_name][key].start
                    if 'end' not in for_plot[str(win_count)]:
                        for_plot[str(win_count)]['end'] = windows[chr_name][key].end
                    window = windows[chr_name][key]
                    var_count = window.get_variation(sample)
                    if var_count > args.cut_off or window.get_observed_kmer(sample) == 0:
                        window.set_block(sample, 'NA')
                        na += 1
                    else:
                        for_plot[str(win_count)][sample] = [1, 'F']
                        if na >= args.step_size or first_in_chromosome:
                            block_no += 1
                            first_in_chromosome = False
                        na = 0
                        window.set_block(sample, block_no)
                    itr += 1
        for sample in samples:
            out_prefix = args.output
            if current_window_size != args.window_size:
                out_prefix = args.output + ".win" + str(args.window_size)
            out_file = out_prefix + "_" + sample + ".blocks"
            out_sum_file = out_prefix + "_" + sample + ".blocks.summary"
            out_sum_flt_file = out_prefix + "_" + sample + ".blocks.summary.filtered"
            f = open(out_file, "w")
            fs = open(out_sum_file, "w")
            fsf = open(out_sum_flt_file, "w")
            print(f"Writing blocks to {out_file}")
            f.write("win_num\tseqname\twin_start\twin_end\tvariations\tblock\n")
            fs.write(f"block_num\tseqname\tstart\tend\tlength\twin_count\t>{str(args.cut_off)}\tconditions\twin_nums\tvariations\tobserved_kmers\n")
            fsf.write(f"block_num\tseqname\tstart\tend\tlength\twin_count\t>{str(args.cut_off)}\tconditions\twin_nums\tvariations\tobserved_kmers\n")
            win_count = 0
            blocks_dict = defaultdict()
            for chr_name in chr_lenths:
                last_block_no = 'NA'
                for key in windows[chr_name]:
                    win_count += 1
                    window = windows[chr_name][key]
                    f.write(
                        f"{str(win_count)}"
                        f"\t{str(window.chr_name)}"
                        f"\t{str(window.start)}"
                        f"\t{str(window.end)}"
                        f"\t{str(window.get_variation(sample))}"
                        f"\t{str(window.get_block(sample))}\n")
                    block_no = window.get_block(sample)
                    if block_no != 'NA':
                        if block_no not in blocks_dict:
                            blocks_dict[block_no] = Block(chr_name, window.start, window.end)
                        blocks_dict[block_no].set_start(window.start)
                        blocks_dict[block_no].set_end(window.end)
                        blocks_dict[block_no].add_block_num('T')
                        blocks_dict[block_no].add_variation(window.get_variation(sample))
                        blocks_dict[block_no].add_window(win_count)
                        blocks_dict[block_no].add_observed_kmer(window.get_observed_kmer(sample))
                        last_block_no = block_no
                    else:
                        if last_block_no != 'NA':
                            blocks_dict[last_block_no].add_block_num('F')
                            blocks_dict[last_block_no].add_variation(window.get_variation(sample))
                            blocks_dict[last_block_no].add_window(win_count)
                            blocks_dict[last_block_no].add_observed_kmer(window.get_observed_kmer(sample))
            for block_no in blocks_dict:
                block = blocks_dict[block_no]
                fs.write(f"{block_no}\t{block}\n")
                true_count = block.get_window_count() - block.get_na_count()
                na_count = block.get_na_count()
                if true_count > na_count and true_count > 2:
                    fsf.write(f"{block_no}\t{block}\n")
                    for win_num in block.get_windows():
                        for_plot[win_num][sample] = [1, 'T']
            f.close()
            fs.close()
            fsf.close()

        for_plot_file = args.output + ".for_plot"
        print(f"Writing for plot to {for_plot_file}")
        f = open(for_plot_file, "w")
        f.write("win_num\tseqname\tstart\tend\t" + "\t".join(samples) + "\n")
        for win_num in for_plot:
            f.write(f"{str(win_num)}\t{str(for_plot[win_num]['chr'])}\t{str(for_plot[win_num]['start'])}\t{str(for_plot[win_num]['end'])}")
            for sample in samples:
                f.write("\t" + str(for_plot[win_num][sample][0]))
            f.write("\n")
        f.close()
        for_plot_file = args.output + ".for_plot.filtered"
        print(f"Writing for plot to {for_plot_file}")
        f = open(for_plot_file, "w")
        f.write("win_num\tseqname\tstart\tend\t" + "\t".join(samples) + "\n")
        for win_num in for_plot:
            f.write(f"{str(win_num)}\t{str(for_plot[win_num]['chr'])}\t{str(for_plot[win_num]['start'])}\t{str(for_plot[win_num]['end'])}")
            for sample in samples:
                if for_plot[win_num][sample][1] == 'T':
                    f.write("\t" + str(for_plot[win_num][sample][0]))
                else:
                    f.write("\t0")
            f.write("\n")

    print(f"Done! Time elapsed: {datetime.now() - start_time}")


def task_plot_variations(args, chr_lenths):
    windows, current_window_size, samples = read_variations(args.output, chr_lenths)
    for chr_name in chr_lenths:
        print(f"Plotting variations for {chr_name}")
        for sample in samples:
            plot_variations(windows, chr_name, current_window_size, sample, args.output)


def task_convert_windows(args, chr_lenths):
    windows, current_window_size, samples = read_variations(args.output, chr_lenths)
    print(f"Converting windows from {current_window_size} to {args.window_size}")
    new_windows = convert_windows(windows, chr_lenths, current_window_size, args.window_size, args.kmer_size)
    write_variations(args.output + ".win" + str(args.window_size), samples, new_windows)


def task_combine_variations(args, chr_lenths):
    windows = defaultdict()
    samples = []
    current_window_size = 0
    with open(args.input) as f:
        for line in f:
            line = line.strip().split(',')
            in_var_file = line[1]
            in_sample_name = line[0]
            samples.append(in_sample_name)
            print(f"Processing {in_var_file}")
            windows, current_window_size = parse_variations(chr_lenths, current_window_size, in_sample_name,
                                                            in_var_file, windows)
    print(f"Writing combined variations to {args.output}")
    write_variations(args.output, samples, windows)


def parse_variations(chr_lenths, current_window_size, in_sample_name, in_var_file, windows):
    """
    parse input variations file from IBSpy_cpp
    :param chr_lenths:
    :param current_window_size:
    :param in_sample_name:
    :param in_var_file:
    :param windows:
    :return:
    """
    with open(in_var_file) as fh:
        for line_fh in fh:
            line_fh = line_fh.strip()
            data = line_fh.split('\t')
            if data[0] == 'seqname':
                continue
            else:
                seqname = data[0]
                start = data[1]
                end = data[2]
                total_kmers = data[3]
                observed_kmers = data[4]
                variations = data[5]
                window_size = int(end) - int(start)
                if current_window_size == 0:
                    current_window_size = window_size
                if current_window_size != window_size:
                    if int(end) != int(chr_lenths[seqname]):
                        print(
                            f"Window size is not same for all the files. Current window size: {current_window_size}, \n "
                            f"window size for {in_var_file}: {window_size}")
                        sys.exit(1)
                if (seqname, start, end) not in windows:
                    windows[(seqname, start, end)] = Window(seqname, start, end)
                windows[(seqname, start, end)].add_variation(in_sample_name, variations)
                windows[(seqname, start, end)].add_observed_kmer(in_sample_name, observed_kmers)
                windows[(seqname, start, end)].set_total_kmers(total_kmers)
    return windows, current_window_size


def read_variations(prefix, chr_lenths):
    """
    read variations from the files and return the windows dictionary
    :param prefix:
    :return:
    """
    windows = defaultdict()
    samples = []
    current_window_size = 0
    print(f"Reading {prefix}.variations.tsv")
    with open(prefix + '.variations.tsv') as f:
        for line in f:
            line = line.strip()
            data = line.split('\t')
            if data[0] == 'seqname':
                header = data
                samples = header[3:]
                continue
            else:
                seqname = data[0]
                start = data[1]
                end = data[2]
                window_size = int(end) - int(start)
                key = (str(seqname), str(start), str(end))
                if current_window_size == 0:
                    current_window_size = window_size
                if current_window_size != window_size:
                    if int(end) != int(chr_lenths[seqname]):
                        print(f"Window size is not same for all the files. Current window size: {current_window_size}, "
                              f"window size for {prefix}.variations.tsv: {window_size}")
                        sys.exit(1)
                for i in range(3, len(data)):
                    variations = data[i]
                    # key = seqname + '_' + start + '_' + end
                    if key not in windows:
                        windows[key] = Window(seqname, start, end)
                    windows[key].add_variation(header[i], variations)
                    # windows[(seqname, start, end)][header[i]] = variations

    print(f"Reading {prefix}.observed_kmers.tsv")
    with open(prefix + '.observed_kmers.tsv') as f:
        for line in f:
            line = line.strip()
            data = line.split('\t')
            if data[0] == 'seqname':
                header = data
                continue
            else:
                seqname = data[0]
                start = data[1]
                end = data[2]
                window_size = int(end) - int(start)
                if current_window_size != window_size:
                    if int(end) != int(chr_lenths[seqname]):
                        print(f"Window size is not same for all the files. Current window size: {current_window_size}, "
                              f"window size for {prefix}.observed_kmers.tsv: {window_size}")
                        sys.exit(1)
                for i in range(3, len(data)):
                    observed_kmers = data[i]
                    windows[(seqname, start, end)].add_observed_kmer(header[i], observed_kmers)

    print(f"Reading {prefix}.total_kmers.tsv")
    with open(prefix + '.total_kmers.tsv') as f:
        for line in f:
            line = line.strip()
            data = line.split('\t')
            if data[0] == 'seqname':
                header = data
                continue
            else:
                seqname = data[0]
                start = data[1]
                end = data[2]
                window_size = int(end) - int(start)
                if current_window_size != window_size:
                    if int(end) != int(chr_lenths[seqname]):
                        print(f"Window size is not same for all the files. Current window size: {current_window_size}, "
                              f"window size for {prefix}.total_kmers.tsv: {window_size}")
                        sys.exit(1)
                for i in range(3, len(data)):
                    total_kmers = data[i]
                    # key = seqname + '_' + start + '_' + end
                    windows[(seqname, start, end)].set_total_kmers(total_kmers)
    return windows, current_window_size, samples


def write_variations(prefix, samples, windows):
    """
    Write variations to 3 different files. variations.tsv, observed_kmers.tsv, total_kmers.tsv
    :param prefix:
    :param samples:
    :param windows:
    :return:
    """
    print(f"Writing output file {prefix}.variations.kvf")
    with open(prefix + '.variations.kvf', 'w') as f:
        f.write('seqname\tstart\tend\ttotal_kmers\tFORMAT' + '\t'.join(samples) + '\n')
        for window in windows:
            f.write(f"{str(windows[window])}\t{str(windows[window].get_total_kmers())}\t"
                    f"V:O:S\t")
            for sample in samples:
                f.write(f"{str(windows[window].write_window(sample))}\t")
        f.close()

    print(f"Writing output file {prefix}.variations.tsv")
    with open(prefix + '.variations.tsv', 'w') as f:
        f.write('seqname\tstart\tend\t' + '\t'.join(samples) + '\n')
        for window in windows:
            f.write(str(windows[window]) + '\t')
            f.write('\t'.join(map(str, windows[window].get_variations(samples))) + '\n')
        f.close()
    print(f"Writing output file {prefix}.observed_kmers.tsv")
    with open(prefix + '.observed_kmers.tsv', 'w') as f:
        f.write('seqname\tstart\tend\t' + '\t'.join(samples) + '\n')
        for window in windows:
            f.write(str(windows[window]) + '\t')
            f.write('\t'.join(map(str, windows[window].get_observed_kmers(samples))) + '\n')
        f.close()
    print(f"Writing output file {prefix}.total_kmers.tsv")
    with open(prefix + '.total_kmers.tsv', 'w') as f:
        f.write('seqname\tstart\tend\ttotal_kmers\n')
        for window in windows:
            f.write(str(windows[window]) + '\t')
            f.write(str(windows[window].get_total_kmers()) + '\n')
        f.close()


def get_window(start, end, window_size, kmer_size):
    """
    return a list of windows
    :param start:
    :param end:
    :param window_size:
    :param kmer_size:
    :return:
    """
    win_start = []
    win_end = []
    t_win_start = start
    t_win_end = start + window_size
    while t_win_end < end:
        win_start.append(t_win_start)
        win_end.append(t_win_end)
        t_win_start = t_win_end - kmer_size
        t_win_end = t_win_start + window_size
    if t_win_end >= end:
        win_start.append(t_win_start)
        win_end.append(end)
    return win_start, win_end


def get_chr_length(chr_length_file):
    """
    return a dictionary with key as chromosome name and value as chromosome length
    :param chr_length_file:
    :return:
    """
    chr_len = defaultdict()
    with open(chr_length_file) as f:
        for line in f:
            line = line.strip().split('\t')
            if line:
                # print(line)
                chr_len[line[0]] = int(line[1])
    f.close()
    return chr_len


if __name__ == '__main__':
    main()
