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

    with open(os.path.join(args.output, 'metadata.tsv'), 'w') as out:
        out.write('#accession\tsample\trep\tfq1\tfq2\n')
        for acc,sample,run,rep,fwd,rev in get_fastqs(args.input):
            os.makedirs(os.path.join(args.output, acc, sample, run), exist_ok=True)
            os.symlink(fwd, os.path.join(args.output, acc, sample, run, f'{rep}_R1.fastq.gz'))
            os.symlink(rev, os.path.join(args.output, acc, sample, run, f'{rep}_R2.fastq.gz'))
            out.write(f'{acc}\t{sample}\t{rep}\t{os.path.join(acc, sample, run, f"{rep}_R1.fastq.gz")}\t{os.path.join(acc, sample, run, f"{rep}_R2.fastq.gz")}\n')



if __name__ == '__main__':
    main()
