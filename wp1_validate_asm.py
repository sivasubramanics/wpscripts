#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd

class PSL:
    def __init__(self, line):
        self.line = line.strip().split('\t')
        self.matches = int(self.line[0])
        self.mismatches = int(self.line[1])
        self.rep_matches = int(self.line[2])
        self.n_count = int(self.line[3])
        self.q_num_insert = int(self.line[4])
        self.q_base_insert = int(self.line[5])
        self.t_num_insert = int(self.line[6])
        self.t_base_insert = int(self.line[7])
        self.strand = self.line[8]
        self.q_name = self.line[9]
        self.q_size = int(self.line[10])
        self.q_start = int(self.line[11])
        self.q_end = int(self.line[12])
        self.t_name = self.line[13]
        self.t_size = int(self.line[14])
        self.t_start = int(self.line[15])
        self.t_end = int(self.line[16])
        self.block_count = int(self.line[17])
        self.block_sizes = [int(x) for x in self.line[18].split(',') if x]
        self.q_starts = [int(x) for x in self.line[19].split(',') if x]
        self.t_starts = [int(x) for x in self.line[20].split(',') if x]

    def get_query_coverage(self):
        return (self.matches + self.mismatches + self.rep_matches) / self.q_size

    def get_target_coverage(self):
        return (self.matches + self.mismatches + self.rep_matches) / self.t_size

    def get_identity(self):
        return (self.matches + self.rep_matches) / (self.matches + self.mismatches + self.rep_matches)

    def get_max_blocks(self):
        return max(self.block_sizes)

    def __str__(self):
        return '\t'.join(self.line)

    def get_blocks_count(self):
        return self.block_count

    def psl_score(self):
        # my $pslScore = $sizeMul * ($matches + ( $repMatches >> 1) ) - $sizeMul * $misMatches - $qNumInsert - $tNumInsert;
        return self.matches + (self.rep_matches / 2) - self.mismatches - self.q_num_insert - self.t_num_insert


def parse_psl(psl_file):
    psl_list = []
    q_name = None
    with open(psl_file, 'r') as fh:
        for line in fh:
            psl = PSL(line)
            if q_name is None:
                q_name = psl.q_name
            if psl.q_name != q_name:
                sorted_psl_list = sorted(psl_list, key=lambda x: x.psl_score(), reverse=True)
                yield sorted_psl_list
                psl_list = []
                q_name = psl.q_name
            psl_list.append(psl)
    sorted_psl_list = sorted(psl_list, key=lambda x: x.psl_score(), reverse=True)
    yield sorted_psl_list


def main():
    parser = argparse.ArgumentParser(description='Validate the assembly based on the mapping with reference transcripts')
    parser.add_argument('-i', '--input', help='input psl file', required=True)
    parser.add_argument('-o', '--outdir', help='output directory', required=True)
    parser.add_argument('-t', '--t_cov', help='target coverage threshold', type=int, default=80)
    parser.add_argument('-q', '--q_cov', help='query coverage threshold', type=int, default=80)
    args = parser.parse_args()

    # create output directory if not exists
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    best_hit_psl = os.path.join(args.outdir, 'best_hit_withscores.psl')
    best_hit_psl_fh = open(best_hit_psl, 'w')

    best_hit_tsv = os.path.join(args.outdir, 'best_hit_withscores.tsv')
    best_hit_tsv_fh = open(best_hit_tsv, 'w')
    best_hit_tsv_fh.write("query\tq_len\ttarget\tt_len\tscore\tq_overage\tt_coverage\tidentity\tn_blocks\n")
    for psl_list in parse_psl(args.input):
        best_hit_psl_fh.write(f"{str(psl_list[0])}\t{psl_list[0].psl_score()}\n")
        best_hit_tsv_fh.write(f"{psl_list[0].q_name}"
                              f"\t{psl_list[0].q_size}"
                              f"\t{psl_list[0].t_name}"
                              f"\t{psl_list[0].t_size}"
                              f"\t{psl_list[0].psl_score()}"
                              f"\t{round(psl_list[0].get_query_coverage() * 100, 2)}"
                              f"\t{round(psl_list[0].get_target_coverage() * 100, 2)}"
                              f"\t{round(psl_list[0].get_identity() * 100, 2)}"
                              f"\t{psl_list[0].get_blocks_count()}"
                              f"\n")
    best_hit_psl_fh.close()
    best_hit_tsv_fh.close()

    get_coverage_summary(best_hit_tsv, 'target')
    get_coverage_summary(best_hit_tsv, 'query')

    complete_asm(best_hit_tsv, args.t_cov, args.q_cov)


def get_coverage_summary(best_hit_tsv, coverage_type):
    if coverage_type == 'query':
        coverage_type = 'q_overage'
        # get directory of the best_hit_tsv and create a new file with the coverage summary
        outfile = os.path.join(os.path.dirname(best_hit_tsv), 'query_coverage.tsv')
    if coverage_type == 'target':
        coverage_type = 't_coverage'
        # get directory of the best_hit_tsv and create a new file with the coverage summary
        outfile = os.path.join(os.path.dirname(best_hit_tsv), 'target_coverage.tsv')

    # create a dataframe from the tsv file
    df = pd.read_csv(best_hit_tsv, sep='\t')
    # get percentage of querys for each 10% of the target coverage
    coverage = df[coverage_type]
    coverage = coverage.apply(lambda x: int(x / 10) * 10)
    coverage = coverage.value_counts()
    coverage = coverage.sort_index()
    coverage = coverage / coverage.sum() * 100
    coverage = coverage.round(2)
    coverage.to_csv(outfile, sep='\t', header=False)


def complete_asm(best_hit_tsv, t_cov_threshold=80, q_cov_threshold=80):
    """
    Extract the t_name and q_name from the best_hit_tsv where the q_coverage and t_coverage are more than cut_off
    """
    out_file = os.path.join(os.path.dirname(best_hit_tsv), 'complete_asm.tsv')
    df = pd.read_csv(best_hit_tsv, sep='\t')
    df = df[(df['q_overage'] >= q_cov_threshold) & (df['t_coverage'] >= t_cov_threshold)]
    df.to_csv(out_file, sep='\t', index=False)



if __name__ == '__main__':
    main()


