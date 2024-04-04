#!/usr/bin/env python3
# script to prepare data for transcriptome assembly pipeline

import os
import argparse


def get_fastqs(in_dir):
    """
    Get fastq filenames from a sample directory. search recursively for fastq files, yeilding pairs of fastq files.
    yield sample_name, run_name, rep_name, fwd, rev
    Ideally the path for each fastq will be {sample_name}/{run_name}/{rep_name}_R1.fastq.gz
    """
    for root, dirs, files in os.walk(in_dir):
        for file in files:
            if file.endswith('_R1_001.fastq.gz'):
                fwd = os.path.join(root, file)
                rev = fwd.replace('_R1_001.fastq.gz', '_R2_001.fastq.gz')
                if os.path.exists(rev):
                    # ideally the path for each fastq will be {sample_name}/{run_name}/{rep_name}_R1.fastq.gz
                    acc_name = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(fwd))))
                    sample_name = os.path.basename(os.path.dirname(os.path.dirname(fwd)))
                    run_name = os.path.basename(os.path.dirname(fwd))
                    rep_name = os.path.basename(fwd).replace('_R1_001.fastq.gz', '')
                    yield [acc_name, sample_name, run_name, rep_name, fwd, rev]


def main():
    parser = argparse.ArgumentParser(description='Prepare data for transcriptome assembly pipeline')
    parser.add_argument('-i', '--input', help='Input directory containing raw data', required=True)
    parser.add_argument('-o', '--output', help='Output directory to store processed data', required=True)
    args = parser.parse_args()

    # create output directory
    os.makedirs(args.output, exist_ok=True)

    for meta_info in get_fastqs(args.input):
        print(f"{meta_info[0]}\t{meta_info[1]}\t{meta_info[2]}\t{meta_info[3]}\t{meta_info[4]}\t{meta_info[5]}")


if __name__ == '__main__':
    main()
